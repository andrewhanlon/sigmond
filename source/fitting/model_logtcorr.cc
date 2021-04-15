#include "model_logtcorr.h"
#include "model_tcorr.h"
#include <cmath>
#include <string>
#include "task_utils.h"
using namespace std;


// ******************************************************************************


void create_logtcorr_model(const string& modeltype, uint in_Tperiod,
                        LogTemporalCorrelatorModel* &mptr)
{
 if (modeltype=="LogTimeForwardSingleExponential"){        
    mptr=new LogTimeForwardSingleExponential(in_Tperiod);}
 else if (modeltype=="LogTimeForwardTwoExponential"){
    mptr=new LogTimeForwardTwoExponential(in_Tperiod);}
 else{
    mptr=0;
    throw(std::invalid_argument(string("Invalid Model in LogRealTemporalCorrelatorFit: ")+modeltype));}
}

// ******************************************************************************

void LogTemporalCorrelatorModel::setupInfos(XMLHandler& xmlm, vector<MCObsInfo>& fitparam_info, int taskcount) const
{
 try{
 fitparam_info.resize(m_nparams);
 int param_count=0;
 for (vector<string>::const_iterator param_name_it=param_names.begin(); param_name_it!=param_names.end(); ++param_name_it, ++param_count) {
    XMLHandler xmlparam(xmlm,*param_name_it);
    string name; int index;
    xmlreadchild(xmlparam,"Name",name);
    if (name.empty()) throw(std::invalid_argument("Must provide name for parameter " + *param_name_it));
    index=taskcount;
    xmlreadifchild(xmlparam,"IDIndex",index);
    fitparam_info[param_count]=MCObsInfo(name,index);
 }

 for (uint k=0;k<m_nparams;k++)
 for (uint l=k+1;l<m_nparams;l++)
    if (fitparam_info[k]==fitparam_info[l])
        throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string(model_name)+" -- "+string(errmsg.what())));}

}

void LogTemporalCorrelatorModel::simpleSetFitInfo(
                          const std::vector<MCObsInfo>& fitparams_info,  
                          const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                          uint fit_tmax, double chisq_dof, double qual,
                          TCorrFitInfo& fitinfo) const
{
 fitinfo.tmin=fit_tmin;
 fitinfo.tmax=fit_tmax;
 fitinfo.meff_approach.clear();
 fitinfo.energy_mean=fitparams[0].getFullEstimate();
 fitinfo.energy_err=fitparams[0].getSymmetricError();
 fitinfo.chisq_dof=chisq_dof;
 fitinfo.quality=qual;
 fitinfo.energy_key=fitparams_info[0];
 fitinfo.amplitude_key=fitparams_info[1];
}



void LogTemporalCorrelatorModel::approachSetFitInfo(
                          const std::vector<MCObsInfo>& fitparams_info,  
                          const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                          uint fit_tmax, uint meff_step, double chisq_dof, double qual,
                          TCorrFitInfo& fitinfo, bool added_constant) const
{
 fitinfo.tmin=fit_tmin;
 fitinfo.tmax=fit_tmax;
 fitinfo.energy_mean=fitparams[0].getFullEstimate();
 fitinfo.energy_err=fitparams[0].getSymmetricError();
 fitinfo.chisq_dof=chisq_dof;
 fitinfo.quality=qual;
 fitinfo.energy_key=fitparams_info[0];
 fitinfo.amplitude_key=fitparams_info[1];

 vector<double> fitmeans(fitparams.size());
 for (int k=0;k<int(fitparams.size());k++)
    fitmeans[k]=fitparams[k].getFullEstimate();
 double subt_const=0.0;
 if (added_constant){
    subt_const=fitmeans[fitparams.size()-1];}
 int npoints=400;
 fitinfo.meff_approach.resize(npoints);
 double curvestep=(double(fit_tmax)-double(fit_tmin))/(double(npoints)-1.0);
 double meffstep=double(meff_step);
 double tval=double(fit_tmin);
 double corr,corrstep,yval;
 double corrback=0.0;
 double tmin=fit_tmax;
 int efftype=m_effmasstype;
// if (efftype>1){     // subtract fit constant
//    efftype-=2;}     // efftypes 2 and 3 remove constant, but noisy
 EffectiveEnergyCalculator Feff(meff_step,T_period,efftype);
 double shift=(Feff.needsBackStep())? 0.0 : 0.5*double(meff_step);
 for (int k=0;k<npoints;k++){
    evaluate(fitmeans,tval,corr);
    evaluate(fitmeans,tval+meffstep,corrstep);
    if (Feff.needsBackStep())
       evaluate(fitmeans,tval-meffstep,corrback);
    if (added_constant){
       corr-=subt_const;
       corrstep-=subt_const;
       corrback-=subt_const;}
    bool flag=Feff.calculate(yval,tval,corr,corrstep,corrback);
    if (flag){
       fitinfo.meff_approach[k]=XYPoint(tval+shift,yval);
       if ((fabs(yval-fitinfo.energy_mean)<=2.0*fitinfo.energy_err)
            &&(tval<tmin)) tmin=tval;}
    tval+=curvestep;}
 fitinfo.tmin=tmin;
}

void LogTemporalCorrelatorModel::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("LogModel",model_name);
}

void LogTemporalCorrelatorModel::setFitInfo(
                   const std::vector<MCObsInfo>& fitparams_info,
                   const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                   uint fit_tmax, bool show_approach,
                   uint meff_step, double chisq_dof, double qual,
                   TCorrFitInfo& fitinfo) const
{
 if (show_approach)
    approachSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,meff_step,chisq_dof,qual,fitinfo);
 else
    simpleSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,chisq_dof,qual,fitinfo);
}

uint LogTemporalCorrelatorModel::getParameterIndex(const string& param_name) const
{
 vector<string>::const_iterator param_it = find(param_names.begin(), param_names.end(), param_name);
 if (param_it == param_names.end())
    throw(std::invalid_argument("Param name does not exist"));
 return (param_it-param_names.begin());
}


// ******************************************************************************


      // Fitting function is single exponential time-forward only:
      //
      //       f(t) = log(A) - m*t
      //
      // where 
      //           m = fitparams[0]
      //      log(A) = fitparams[1].


void LogTimeForwardSingleExponential::evaluate(const vector<double>& fitparams, 
                                               double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],tval,value);
}


void LogTimeForwardSingleExponential::evalGradient(const vector<double>& fitparams, 
                         double tval, vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],tval,grad[1],grad[0]);
}


void LogTimeForwardSingleExponential::guessInitialParamValues(
                   const vector<double>& data, const vector<uint>& tvals,
                   vector<double>& fitparams) const
{
 if (data.size()<2)
    throw(std::invalid_argument("SingleExponential -- Error: at least two data points needed! in exponential guess"));
 vector<double> exp_data(data.size());
 for (uint i=0;i<data.size();i++)
   exp_data[i] = exp(data[i]);
 TimeForwardSingleExponential::get_exp_guess(tvals,exp_data,fitparams[0],fitparams[1]);
 fitparams[1] = log(fitparams[1]);
}



      //       f(t) = log(A) - m*t

 
void LogTimeForwardSingleExponential::eval_func(double logA, double m, double tf, double& funcval) const
{
 funcval=logA-m*tf;
}


void LogTimeForwardSingleExponential::eval_grad(double logA, double m, double tf, double& dlogAval, 
                                               double& dmval) const
{
 dlogAval=1.0;
 dmval=-tf;
}





// ******************************************************************************

      // The fitting function is the log of a sum of two exponentials,
      // time-forward:
      //
      //    f(t) = log(A) - m*t + log[ 1 + B*exp(-Delta^2*t) ]
      //
      //  where 
      //          m = fitparams[0]
      //     log(A) = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //


void LogTimeForwardTwoExponential::evaluate(
            const vector<double>& fitparams, double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,value);
}


void LogTimeForwardTwoExponential::evalGradient(
                const vector<double>& fitparams, double tval, 
                vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,
           grad[1],grad[0],grad[3],grad[2]);
}




void LogTimeForwardTwoExponential::guessInitialParamValues(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams) const
{
 vector<double> exp_data(data.size());
 for (uint i=0;i<data.size();i++)
   exp_data[i] = exp(data[i]);
 double tasymfrac=0.33;
 TimeForwardTwoExponential::get_two_exp_guess(tvals,exp_data,fitparams[0],fitparams[1],fitparams[2],fitparams[3],tasymfrac);
 fitparams[1]=log(fitparams[1]);
}


      //    f(t) = log(A) - m*t + log[ 1 + B*exp(-DD^2*t) ] 


void LogTimeForwardTwoExponential::eval_func(
              double logA, double m, double B, double DD,
              double tf, double& funcval) const
{
 funcval=logA-m*tf+(1.0+B*exp(-DD*DD*tf));
}


void LogTimeForwardTwoExponential::eval_grad(
              double logA, double m, double B, double DD,
              double tf, double& dlogAval, double& dmval,
              double& dBval, double& dDDval) const
{
 double gap=DD*DD;
 double r2=exp(-gap*tf);
 dlogAval=1.0;
 dmval=-tf;
 dBval=r2/(1.0+B*r2);
 dDDval=-2.0*tf*B*DD*dBval;
}



 // ***********************************************************************************

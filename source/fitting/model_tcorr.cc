#include "model_tcorr.h"
#include <cmath>
#include <string>
#include "task_utils.h"
using namespace std;


// ******************************************************************************


void create_tcorr_model(const string& modeltype, uint in_Tperiod,
                        TemporalCorrelatorModel* &mptr)
{
 if (modeltype=="TimeForwardSingleExponential"){        
    mptr=new TimeForwardSingleExponential(in_Tperiod);}
 else if (modeltype=="TimeSymSingleExponential"){
    mptr=new TimeSymSingleExponential(in_Tperiod);}
 else if (modeltype=="TimeForwardSingleExponentialPlusConstant"){
    mptr=new TimeForwardSingleExponentialPlusConstant(in_Tperiod);}
 else if (modeltype=="TimeSymSingleExponentialPlusConstant"){
    mptr=new TimeSymSingleExponentialPlusConstant(in_Tperiod);}
 else if (modeltype=="TimeForwardTwoExponential"){
    mptr=new TimeForwardTwoExponential(in_Tperiod);}
 else if (modeltype=="TimeSymTwoExponential"){
    mptr=new TimeSymTwoExponential(in_Tperiod);}
 else if (modeltype=="TimeForwardTwoExponentialPlusConstant"){
    mptr=new TimeForwardTwoExponentialPlusConstant(in_Tperiod);}
 else if (modeltype=="TimeSymTwoExponentialPlusConstant"){
    mptr=new TimeSymTwoExponentialPlusConstant(in_Tperiod);}
 else if (modeltype=="TimeForwardGeomSeriesExponential"){
    mptr=new TimeForwardGeomSeriesExponential(in_Tperiod);}
 else if (modeltype=="TimeSymGeomSeriesExponential"){
    mptr=new TimeSymGeomSeriesExponential(in_Tperiod);}
 else{
    mptr=0;
    throw(std::invalid_argument(string("Invalid Model in RealTemporalCorrelatorFit: ")+modeltype));}
}

// ******************************************************************************


void TemporalCorrelatorModel::simpleSetFitInfo(
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



void TemporalCorrelatorModel::approachSetFitInfo(
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
 if (efftype>1){     // subtract fit constant
    efftype-=2;}     // efftypes 2 and 3 remove constant, but noisy
 EffectiveEnergyCalculator Feff(meff_step,T_period,efftype);
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
       fitinfo.meff_approach[k]=XYPoint(tval,yval);
       if ((fabs(yval-fitinfo.energy_mean)<=2.0*fitinfo.energy_err)
            &&(tval<tmin)) tmin=tval;}
    tval+=curvestep;}
 fitinfo.tmin=tmin;
}

// ******************************************************************************


      // Fitting function is single exponential time-forward only:
      //
      //       f(t) = A * exp( -m*t ) 
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1].


void TimeForwardSingleExponential::setupInfos(XMLHandler& xmlm, 
                          vector<MCObsInfo>& fitparam_info, int taskcount) const
{
 setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeForwardSingleExponential::setup(XMLHandler& xmlm, 
                 vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount)
{
 try{
 fitparam_info.resize(nparam);
 XMLHandler xmlen(xmlm,"Energy");
 string name; int index;
 xmlreadchild(xmlen,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlen,"IDIndex",index);
 fitparam_info[0]=MCObsInfo(name,index);

 XMLHandler xmla(xmlm,"Amplitude");
 name.clear();
 xmlreadchild(xmla,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for amplitude parameter"));
 index=taskcount;
 xmlreadifchild(xmla,"IDIndex",index);
 fitparam_info[1]=MCObsInfo(name,index);

 if (fitparam_info[0]==fitparam_info[1])
     throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("SingleExponential -- ")+string(errmsg.what())));}
}


void TimeForwardSingleExponential::evaluate(const vector<double>& fitparams, 
                                            double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],tval,value);
}


void TimeForwardSingleExponential::evalGradient(const vector<double>& fitparams, 
                         double tval, vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],tval,grad[1],grad[0]);
}


void TimeForwardSingleExponential::guessInitialParamValues(
                   const vector<double>& data, const vector<uint>& tvals,
                   vector<double>& fitparams) const
{
 if (data.size()<2)
    throw(std::invalid_argument("SingleExponential -- Error: at least two data points needed! in exponential guess"));
 eval_guess(tvals[0],data[0],tvals[1],data[1],fitparams[1],fitparams[0]);
}


void TimeForwardSingleExponential::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardSingleExponential");
}


void TimeForwardSingleExponential::setFitInfo(
                   const std::vector<MCObsInfo>& fitparams_info,
                   const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                   uint fit_tmax, bool show_approach,
                   uint meff_step, double chisq_dof, double qual,
                   TCorrFitInfo& fitinfo) const
{ 
 simpleSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,chisq_dof,qual,fitinfo);
}

      //       f(t) = A * exp( -m*t ) 

 
void TimeForwardSingleExponential::eval_func(double A, double m, double tf, double& funcval) const
{
 funcval=A*exp(-m*tf);
}


void TimeForwardSingleExponential::eval_grad(double A, double m, double tf, double& dAval, 
                                             double& dmval) const
{
 dAval=exp(-m*tf); 
 dmval=-tf*A*dAval;
}


void TimeForwardSingleExponential::eval_guess(int tval, double corrt, int tnext, double corrtnext, 
                                              double& A, double& m)
{
 double s=corrt/corrtnext;
 if ((s<=0.0)||(tval==tnext))
    throw(std::invalid_argument("SingleExponential -- could not compute a guess for exponential"));
 m=log(s)/(double(tnext)-double(tval));       // guess for m
 A=exp(m*double(tval))*corrt;                 // guess for A
}



// ******************************************************************************


      // Fitting function is single exponential time-symmetric
      // (forwards and backwards):
      //
      //       f(t) = A * {exp( -m*t ) +  exp( -m*(T_period-t) )}
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1].

void TimeSymSingleExponential::setupInfos(XMLHandler& xmlm, 
                        vector<MCObsInfo>& fitparam_info, int taskcount) const
{
 TimeForwardSingleExponential::setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeSymSingleExponential::evaluate(const vector<double>& fitparams, 
                                        double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],tval,T_period,value);
}


void TimeSymSingleExponential::evalGradient(const vector<double>& fitparams, 
                         double tval, vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],tval,T_period,grad[1],grad[0]);
}


void TimeSymSingleExponential::guessInitialParamValues(
                   const vector<double>& data, const vector<uint>& tvals,
                   vector<double>& fitparams) const
{
 if (data.size()<2)
    throw(std::invalid_argument("SingleExponential -- Error: at least two data points needed! in exponential guess"));
 TimeForwardSingleExponential::eval_guess(tvals[0],data[0],tvals[1],data[1],
                                          fitparams[1],fitparams[0]);
}


void TimeSymSingleExponential::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeSymSingleExponential");
}


void TimeSymSingleExponential::setFitInfo(
                   const std::vector<MCObsInfo>& fitparams_info,
                   const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                   uint fit_tmax, bool show_approach,
                   uint meff_step, double chisq_dof, double qual,
                   TCorrFitInfo& fitinfo) const
{
 simpleSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,chisq_dof,qual,fitinfo);
}


      //    f(t) = A * { exp(-m*t)  + exp(-m*(Nt-t))  } 

void TimeSymSingleExponential::eval_func(double A, double m,
                                         double tf, int Nt, double& funcval) const
{
 double tb=double(Nt)-tf;
 funcval=A*(exp(-m*tf)+exp(-m*tb));
}


void TimeSymSingleExponential::eval_grad(double A, double m, 
                       double tf, int Nt, double& dAval, double& dmval) const
{
 dAval=exp(-m*tf); 
 dmval=-tf*A*dAval;
 double tb=double(Nt)-tf;
 double r1=exp(-m*tb);
 dAval+=r1;
 dmval-=tb*A*r1;
}


// ******************************************************************************


      // The fitting function is a single exponential time-forward
      // with an added constant:
      //
      //       f(t) = A * exp( -m*t ) + c0
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1]
      //          c0 = fitparams[2]
      //


void TimeForwardSingleExponentialPlusConstant::setupInfos(XMLHandler& xmlm, 
                                      vector<MCObsInfo>& fitparam_info, int taskcount) const
{
 setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeForwardSingleExponentialPlusConstant::setup(XMLHandler& xmlm, 
                        vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount)
{
 try{
 fitparam_info.resize(nparam);
 XMLHandler xmlen(xmlm,"Energy");
 string name; int index;
 xmlreadchild(xmlen,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlen,"IDIndex",index);
 fitparam_info[0]=MCObsInfo(name,index);

 XMLHandler xmla(xmlm,"Amplitude");
 name.clear();
 xmlreadchild(xmla,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for amplitude parameter"));
 index=taskcount;
 xmlreadifchild(xmla,"IDIndex",index);
 fitparam_info[1]=MCObsInfo(name,index);

 XMLHandler xmlc(xmlm,"AddedConstant");
 name.clear();
 xmlreadchild(xmlc,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for amplitude parameter"));
 index=taskcount;
 xmlreadifchild(xmlc,"IDIndex",index);
 fitparam_info[2]=MCObsInfo(name,index);

 for (uint k=0;k<nparam;k++)
 for (uint l=k+1;l<nparam;l++)
    if (fitparam_info[k]==fitparam_info[l])
        throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("SingleExponentialPlusConst -- ")+string(errmsg.what())));}

}


void TimeForwardSingleExponentialPlusConstant::evaluate(
                             const vector<double>& fitparams, 
                             double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[2],tval,value);
}


void TimeForwardSingleExponentialPlusConstant::evalGradient(
                         const vector<double>& fitparams, 
                         double tval, vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],tval,grad[1],grad[0],grad[2]);
}


void TimeForwardSingleExponentialPlusConstant::guessInitialParamValues(
                   const vector<double>& data, const vector<uint>& tvals,
                   vector<double>& fitparams) const
{
 if (data.size()<3)
    throw(std::invalid_argument("SingleExponentialPlusConst -- Error: at least three data points needed! in exponential+const guess"));
 eval_guess(tvals[0],data[0],tvals[1],data[1],tvals[2],data[2],fitparams[1],fitparams[0],fitparams[2]);
}


void TimeForwardSingleExponentialPlusConstant::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardSingleExponentialPlusConstant");
}


void TimeForwardSingleExponentialPlusConstant::setFitInfo(
                   const std::vector<MCObsInfo>& fitparams_info,
                   const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                   uint fit_tmax, bool show_approach,
                   uint meff_step, double chisq_dof, double qual,
                   TCorrFitInfo& fitinfo) const
{
 simpleSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,chisq_dof,qual,fitinfo);
}



      //    f(t) = A * exp(-m*t)  + c0

void TimeForwardSingleExponentialPlusConstant::eval_func(
                double A, double m, double c0, double tf, double& funcval) const
{
 funcval=A*exp(-m*tf)+c0;
}


void TimeForwardSingleExponentialPlusConstant::eval_grad(
                double A, double m, double tf, double& dAval, double& dmval,
                double& dc0val) const
{
 dAval=exp(-m*tf); 
 dmval=-tf*A*dAval;
 dc0val=1.0;
}


void TimeForwardSingleExponentialPlusConstant::eval_guess(
                int tval, double corrt, int tp1, double corrtp1, int tp2, double corrtp2,
                double& A, double& m, double& c0)
{
 if (((tp2-tp1)!=(tp1-tval))||(tval==tp1))
    throw(std::invalid_argument("SingleExponentialPlusConst -- could not compute a guess for exponential"));
 double cor0=corrtp1-corrt;
 double cor1=corrtp2-corrtp1;
 double s=cor0/cor1;
 double d=double(tp1)-double(tval);
 if (s<=0.0)
    throw(std::invalid_argument("SingleExponentialPlusConst -- could not compute a guess for exponential"));
 m=log(s)/d;
 double rr=exp(-m*double(tval));
 A=cor0/(rr*(exp(-m*d)-1.0));
 c0=corrt-A*rr;  
}


// ******************************************************************************


      // The fitting function is a single exponential time-symmetric
      // (forwards and backwards) with an added constant:
      //
      //       f(t) = A * {exp( -m*t ) +  exp( -m*(T_period-t) )} + c0
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1]
      //          c0 = fitparams[2]
      //


void TimeSymSingleExponentialPlusConstant::setupInfos(XMLHandler& xmlm, 
                                      vector<MCObsInfo>& fitparam_info, int taskcount) const
{
 TimeForwardSingleExponentialPlusConstant::setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeSymSingleExponentialPlusConstant::evaluate(
                             const vector<double>& fitparams, 
                             double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[2],tval,T_period,value);
}

    
void TimeSymSingleExponentialPlusConstant::evalGradient(
                         const vector<double>& fitparams, 
                         double tval, vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],tval,T_period,grad[1],grad[0],grad[2]);
}


void TimeSymSingleExponentialPlusConstant::guessInitialParamValues(
                   const vector<double>& data, const vector<uint>& tvals,
                   vector<double>& fitparams) const
{
 if (data.size()<3)
    throw(std::invalid_argument("SingleExponentialPlusConst -- Error: at least three data points needed! in exponential+const guess"));
 TimeForwardSingleExponentialPlusConstant::eval_guess(
      tvals[0],data[0],tvals[1],data[1],tvals[2],data[2],fitparams[1],fitparams[0],fitparams[2]);
}


void TimeSymSingleExponentialPlusConstant::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeSymSingleExponentialPlusConstant");
}


void TimeSymSingleExponentialPlusConstant::setFitInfo(
                   const std::vector<MCObsInfo>& fitparams_info,
                   const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                   uint fit_tmax, bool show_approach,
                   uint meff_step, double chisq_dof, double qual,
                   TCorrFitInfo& fitinfo) const
{
 simpleSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,chisq_dof,qual,fitinfo);
}



      //    f(t) = A * { exp(-m*t)  + exp(-m*(Nt-t))  }   + c0

void TimeSymSingleExponentialPlusConstant::eval_func(
               double A, double m, double c0, double tf, int Nt, double& funcval) const
{
 double tb=double(Nt)-tf;
 funcval=A*(exp(-m*tf)+exp(-m*tb))+c0;
}


void TimeSymSingleExponentialPlusConstant::eval_grad(
               double A, double m, double tf, int Nt, double& dAval, double& dmval,
               double& dc0val) const
{
 dAval=exp(-m*tf); 
 dmval=-tf*A*dAval;
 double tb=double(Nt)-tf;
 double r1=exp(-m*tb);
 dAval+=r1;
 dmval-=tb*A*r1;
 dc0val=1.0;
}


// ******************************************************************************

      // The fitting function is a sum of two exponentials, time-forward:
      //
      //    f(t) = A * exp(-m*t) * [ 1 + B*exp(-Delta^2*t) ]
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //


void TimeForwardTwoExponential::setupInfos(XMLHandler& xmlm, 
                                vector<MCObsInfo>& fitparam_info, int taskcount) const
{
 setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeForwardTwoExponential::setup(XMLHandler& xmlm, 
                   vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount)
{
 try{
 fitparam_info.resize(nparam);
 XMLHandler xmlen(xmlm,"FirstEnergy");
 string name; int index;
 xmlreadchild(xmlen,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for first energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlen,"IDIndex",index);
 fitparam_info[0]=MCObsInfo(name,index);

 XMLHandler xmla(xmlm,"FirstAmplitude");
 name.clear();
 xmlreadchild(xmla,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for first amplitude parameter"));
 index=taskcount;
 xmlreadifchild(xmla,"IDIndex",index);
 fitparam_info[1]=MCObsInfo(name,index);

 XMLHandler xmlg(xmlm,"SqrtGapToSecondEnergy");
 name.clear();
 xmlreadchild(xmlg,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for sqrt gap to second energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlg,"IDIndex",index);
 fitparam_info[2]=MCObsInfo(name,index);

 XMLHandler xmlb(xmlm,"SecondAmplitudeRatio");
 name.clear();
 xmlreadchild(xmlb,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for second amplitude ratio parameter"));
 index=taskcount;
 xmlreadifchild(xmlb,"IDIndex",index);
 fitparam_info[3]=MCObsInfo(name,index);

 for (uint k=0;k<nparam;k++)
 for (uint l=k+1;l<nparam;l++)
    if (fitparam_info[k]==fitparam_info[l])
        throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("TwoExponential -- ")+string(errmsg.what())));}

}


void TimeForwardTwoExponential::evaluate(
            const vector<double>& fitparams, double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,value);
}


void TimeForwardTwoExponential::evalGradient(
                const vector<double>& fitparams, double tval, 
                vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,
           grad[1],grad[0],grad[3],grad[2]);
}

/*
void TimeForwardTwoExponential::guessInitialParamValues(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams) const
{
 if (data.size()<4)
    throw(std::invalid_argument("TwoExponential -- Error: at least four data points needed! in two exponential guess"));
 for (int nd=3;nd<int(data.size()-3);nd++)
    if (eval_guess(tmin,data[0],data[1],data[nd-1],data[nd],
                   fitparams[1],fitparams[0],fitparams[3],fitparams[2])) return;
 fitparams[2]=0.0;
 fitparams[3]=0.0;
 TimeForwardSingleExponential::eval_guess(tmin,data[2],data[3],
                                          fitparams[1],fitparams[0]);
} */


void TimeForwardTwoExponential::guessInitialParamValues(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams) const
{
 eval_guess(data,tvals,fitparams);
}


void TimeForwardTwoExponential::eval_guess(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams)
{
 if (data.size()<4)
    throw(std::invalid_argument("TwoExponential -- Error: at least four data points needed! in two exponential guess"));
 uint kfar=2*data.size()/3;
 if (kfar<2) kfar=2;
 if (kfar==data.size()-1) kfar--;
 eval_guess(tvals[kfar],data[kfar],tvals[kfar+1],data[kfar+1], 
            tvals[0],data[0],tvals[1],data[1], 
            fitparams[1],fitparams[0],fitparams[3],fitparams[2]);
}


void TimeForwardTwoExponential::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardTwoExponential");
}


void TimeForwardTwoExponential::setFitInfo(
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


      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-DD^2*t) ]  }


void TimeForwardTwoExponential::eval_func(
              double A, double m, double B, double DD,
              double tf, double& funcval) const
{
 funcval=A*(exp(-m*tf)*(1.0+B*exp(-DD*DD*tf)));
}


void TimeForwardTwoExponential::eval_grad(
              double A, double m, double B, double DD,
              double tf, double& dAval, double& dmval,
              double& dBval, double& dDDval) const
{
 double gap=DD*DD;
 double r1=exp(-m*tf); 
 double r2=exp(-gap*tf);
 dAval=r1*(1.0+B*r2);
 dmval=-tf*A*dAval;
 dBval=A*r1*r2;
 dDDval=-2.0*tf*B*DD*dBval;
}

/*
bool TimeForwardTwoExponential::eval_guess(
              int tval, double f0, double f1, double f2, double f3,
              double& A, double& m, double& B, double& DD)
{
 double a=f1*f1-f0*f2;
 double b=f3*f0-f1*f2;
 double c=f2*f2-f1*f3;
 double x1=b*b-4.0*a*c;
 if (x1>0.0){
    x1=sqrt(x1);
    double x2=0.5*(-b-x1)/a;
    x1=0.5*(-b+x1)/a;
    if ((x1>0.0)&&(x1<=1.0)&&(x2>0.0)&&(x2<=1.0)){
       double x;
       if (x1>x2) x=x1;
       else x=x2;
       double y=(f2-x*f1)/(x*(f1-x*f0));
       if ((x>0.0)&&(y>0.0)){
          m=-log(x);
          double gap=-log(y);
          a=x*y*f0-f1;
          A=a/(x*(y-1.0));
          B=(f1-x*f0)/a;
          A*=exp(double(tval)*m);
          B*=exp(double(tval)*gap);
          DD=sqrt(gap);
          return true;}
       }}
 return false;
}
*/
     //  initial guess for a two-exponential fit 
     //     A * exp( -m*t ) * (1 + B * exp(-DD^2*t) )
     //
     //  choose tfar in large time region where second
     //  exponential is negligible, then choose tnear
     //  in small time region where second exponential
     //  can be exposed;  ffar = corr(tfar), ffarnext=corr(tfarnext)
     //  and fnear=corr(tnear), fnearnext=corr(tnearnext)

void TimeForwardTwoExponential::eval_guess(
              int tfar, double ffar, int tfarnext, double ffarnext, 
              int tnear, double fnear, int tnearnext, double fnearnext,
              double& A, double& m, double& B, double& DD)
{
 double s=ffar/ffarnext;
 if ((s<=0.0)||(tfarnext==tfar))
    throw(std::invalid_argument("Two exponential -- could not compute an initial guess"));
 m=log(s)/(double(tfarnext)-double(tfar)); // guess for m
 A=exp(m*double(tfar))*ffar;               // guess for A
 double r1=exp(m*double(tnear))*fnear/A-1.0;
 double r2=exp(m*double(tnearnext))*fnearnext/A-1.0;
 s=r1/r2;
 if ((s<=1.0)||(tnearnext==tnear))
    throw(std::invalid_argument("Two exponential -- could not compute an initial guess"));
 double DDsq=log(s)/(double(tnearnext)-double(tnear)); // guess for DD
 DD=sqrt(DDsq);
 B=exp(DDsq*double(tnear))*r1;       // guess for B
}


 // ***********************************************************************************


      // The fitting function is a sum of two exponentials, time-symmetric:
      //
      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-Delta^2*t) ]
      //          + exp(-m*(T_period-t)) * [ 1 + B*exp(-Delta^2*(T_period-t)) ] }
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //


void TimeSymTwoExponential::setupInfos(XMLHandler& xmlm, 
                                vector<MCObsInfo>& fitparam_info, int taskcount) const
{
 TimeForwardTwoExponential::setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeSymTwoExponential::evaluate(const vector<double>& fitparams, double tval, 
                                     double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,T_period,value);
}


void TimeSymTwoExponential::evalGradient(const vector<double>& fitparams, double tval, 
                                         vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,T_period,
           grad[1],grad[0],grad[3],grad[2]);
}

/*
void TimeSymTwoExponential::guessInitialParamValues(
                              const vector<double>& data, const vector<uint>& tvals,
                              vector<double>& fitparams) const
{
 if (data.size()<4)
    throw(std::invalid_argument("TwoExponential -- Error: at least four data points needed! in two exponential guess"));
 for (int nd=3;nd<int(data.size()-3);nd++)
    if (TimeForwardTwoExponential::eval_guess(tmin,data[0],data[1],data[nd-1],data[nd],
                   fitparams[1],fitparams[0],fitparams[3],fitparams[2])) return;
 fitparams[2]=0.0;
 fitparams[3]=0.0;
 TimeForwardSingleExponential::eval_guess(tmin,data[2],data[3],
                                          fitparams[1],fitparams[0]);
} */

void TimeSymTwoExponential::guessInitialParamValues(
                              const vector<double>& data, const vector<uint>& tvals,
                              vector<double>& fitparams) const
{
 TimeForwardTwoExponential::eval_guess(data,tvals,fitparams);
}



void TimeSymTwoExponential::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeSymTwoExponential");
}


void TimeSymTwoExponential::setFitInfo(
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



      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-DD^2*t) ]
      //          + exp(-m*(Nt-t)) * [ 1 + B*exp(-DD^2*(Nt-t)) ] }

void TimeSymTwoExponential::eval_func(
                double A, double m, double B, double DD,
                double tf, int Nt, double& funcval) const
{
 double tb=double(Nt)-tf;
 funcval=A*(exp(-m*tf)*(1.0+B*exp(-DD*DD*tf))+exp(-m*tb)*(1.0+B*exp(-DD*DD*tb)));
}


void TimeSymTwoExponential::eval_grad(
                double A, double m, double B, double DD,
                double tf, int Nt, double& dAval, double& dmval,
                double& dBval, double& dDDval) const
{
 double gap=DD*DD;
 double r1=exp(-m*tf); 
 double r2=exp(-gap*tf);
 dAval=r1*(1.0+B*r2);
 dmval=-tf*A*dAval;
 dBval=A*r1*r2;
 dDDval=-2.0*tf*B*DD*dBval;
 double tb=double(Nt)-tf;
 r1=exp(-m*tb); r2=exp(-gap*tb);
 double s=r1*(1.0+B*r2);
 dAval+=s;
 dmval-=tb*A*s;
 s=A*r1*r2;
 dBval+=s;
 dDDval-=2.0*tb*B*DD*s;
}


 // ***********************************************************************************


      // The fitting function is a sum of two exponentials, time-forward + constant
      //
      //    f(t) = A * exp(-m*t) * [ 1 + B*exp(-Delta^2*t) ] + c0
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //         c0 = fitparams[4]
      //


void TimeForwardTwoExponentialPlusConstant::setupInfos(XMLHandler& xmlm, 
                                vector<MCObsInfo>& fitparam_info, int taskcount) const
{
 setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeForwardTwoExponentialPlusConstant::setup(XMLHandler& xmlm,
                       vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount)
{
 try{
 fitparam_info.resize(nparam);
 XMLHandler xmlen(xmlm,"FirstEnergy");
 string name; int index;
 xmlreadchild(xmlen,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for first energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlen,"IDIndex",index);
 fitparam_info[0]=MCObsInfo(name,index);

 XMLHandler xmla(xmlm,"FirstAmplitude");
 name.clear();
 xmlreadchild(xmla,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for first amplitude parameter"));
 index=taskcount;
 xmlreadifchild(xmla,"IDIndex",index);
 fitparam_info[1]=MCObsInfo(name,index);

 XMLHandler xmlg(xmlm,"SqrtGapToSecondEnergy");
 name.clear();
 xmlreadchild(xmlg,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for sqrt gap to second energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlg,"IDIndex",index);
 fitparam_info[2]=MCObsInfo(name,index);

 XMLHandler xmlb(xmlm,"SecondAmplitudeRatio");
 name.clear();
 xmlreadchild(xmlb,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for second amplitude ratio parameter"));
 index=taskcount;
 xmlreadifchild(xmlb,"IDIndex",index);
 fitparam_info[3]=MCObsInfo(name,index);

 XMLHandler xmlc(xmlm,"AddedConstant");
 name.clear();
 xmlreadchild(xmlc,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for amplitude parameter"));
 index=taskcount;
 xmlreadifchild(xmlc,"IDIndex",index);
 fitparam_info[4]=MCObsInfo(name,index);

 for (uint k=0;k<nparam;k++)
 for (uint l=k+1;l<nparam;l++)
    if (fitparam_info[k]==fitparam_info[l])
        throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("TwoExponentialPlusConst -- ")+string(errmsg.what())));}
}



void TimeForwardTwoExponentialPlusConstant::evaluate(
                        const vector<double>& fitparams, double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[3],fitparams[2],fitparams[4],tval,value);
}


void TimeForwardTwoExponentialPlusConstant::evalGradient(
                        const vector<double>& fitparams, double tval, 
                        vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,
           grad[1],grad[0],grad[3],grad[2],grad[4]);
}

/*
void TimeForwardTwoExponentialPlusConstant::guessInitialParamValues(
                        const vector<double>& data, const vector<uint>& tvals,
                        vector<double>& fitparams) const
{
 if (data.size()<5)
    throw(std::invalid_argument("TwoExponentialPlusConst -- Error: at least five data points needed! in two exponential guess"));
 for (int nd=3;nd<int(data.size()-4);nd++)
    if (eval_guess(tmin,data[0],data[1],data[nd-1],data[nd],data[nd+1],
                   fitparams[1],fitparams[0],fitparams[3],fitparams[2],fitparams[4])) return;
 fitparams[2]=0.0;
 fitparams[3]=0.0;
 TimeForwardSingleExponentialPlusConstant::eval_guess(tmin,data[2],data[3],data[4],
                                          fitparams[1],fitparams[0],fitparams[4]);
}
*/

void TimeForwardTwoExponentialPlusConstant::guessInitialParamValues(
                        const vector<double>& data, const vector<uint>& tvals,
                        vector<double>& fitparams) const
{
 eval_guess(data,tvals,fitparams);
}


void TimeForwardTwoExponentialPlusConstant::eval_guess(
                        const vector<double>& data, const vector<uint>& tvals,
                        vector<double>& fitparams)
{
 if (data.size()<5)
    throw(std::invalid_argument("TwoExponentialPlusConst -- Error: at least five data points needed! in two exponential guess"));
 uint kfar=2*data.size()/3;
 if (kfar<2) kfar=2;
 if (kfar>=data.size()-2) kfar=data.size()-3;
 eval_guess(tvals[kfar],data[kfar],tvals[kfar+1],data[kfar+1],tvals[kfar+2],data[kfar+2],
            tvals[0],data[0],tvals[1],data[1], 
            fitparams[1],fitparams[0],fitparams[3],fitparams[2],fitparams[4]);
}


void TimeForwardTwoExponentialPlusConstant::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardTwoExponentialPlusConstant");
}



void TimeForwardTwoExponentialPlusConstant::setFitInfo(
                   const std::vector<MCObsInfo>& fitparams_info,
                   const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                   uint fit_tmax, bool show_approach,
                   uint meff_step, double chisq_dof, double qual,
                   TCorrFitInfo& fitinfo) const
{
 if (show_approach)
    approachSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,meff_step,chisq_dof,qual,fitinfo,true);
 else
    simpleSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,chisq_dof,qual,fitinfo);
}



      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-DD^2*t) ] }   + c0

void TimeForwardTwoExponentialPlusConstant::eval_func(
                    double A, double m, double B, double DD, double c0,
                    double tf, double& funcval) const
{
 funcval=A*(exp(-m*tf)*(1.0+B*exp(-DD*DD*tf)))+c0;
}


void TimeForwardTwoExponentialPlusConstant::eval_grad(
                    double A, double m, double B, double DD, 
                    double tf, double& dAval, double& dmval,
                    double& dBval, double& dDDval, double& dc0val) const
{
 double gap=DD*DD;
 double r1=exp(-m*tf); 
 double r2=exp(-gap*tf);
 dAval=r1*(1.0+B*r2);
 dmval=-tf*A*dAval;
 dBval=A*r1*r2;
 dDDval=-2.0*tf*B*DD*dBval;
 dc0val=1.0;
}

/*
bool TimeForwardTwoExponentialPlusConstant::eval_guess(
                   int tval, double f0, double f1, double f2, double f3,
                   double f4, double& A, double& m, double& B, 
                   double& DD, double& c0)
{
 double a=f0*f2-f0*f3-f1*f1+f1*f2+f1*f3-f2*f2;
 double b=-f0*f3+f0*f4+f1*f2-f1*f4-f2*f2+f2*f3;
 double c=f1*f3-f1*f4-f2*f2+f2*f3+f2*f4-f3*f3;
 double x1=b*b-4.0*a*c;
 if (x1>0.0){
    x1=sqrt(x1);
    double x2=0.5*(-b-x1)/a;
    x1=0.5*(-b+x1)/a;
    if ((x1>0.0)&&(x1<1.0)&&(x2>0.0)&&(x2<1.0)){
       double x=(x1>x2)?x1:x2;
       double y=(x1>x2)?x2:x1;
       m=-log(x);
       double gap=-log(y)-m;
       double A1=-(f0*y-f1*y-f1+f2)/((1.0-x)*(x-y));
       double A2=(f0*x-f1*x-f1+f2)/((x-y)*(1.0-y));
       c0=(f0*x*y-f1*x-f1*y+f2)/((1.0-y)*(1.0-x));
       A=A1*exp(double(tval)*m);
       B=(A2/A1)*exp(double(tval)*gap);
       DD=sqrt(gap);
       return true;}}
 return false;
}
*/


     //  initial guess for a two-exponential fit 
     //     A * exp( -m*t ) * (1 + B * exp(-DD^2*t) ) + c0
     //
     //  choose tfar in large time region where second
     //  exponential is negligible, then choose tnear
     //  in small time region where second exponential
     //  can be exposed;  
     //   ffar = corr(tfar), ffarp1=corr(tfar+1), ffarp2=corr(tfar+2)
     //  and fnear=corr(tnear), fnearnext=corr(tnear+1)

void TimeForwardTwoExponentialPlusConstant::eval_guess(
              int tfar, double ffar, int tfarp1, double ffarp1, int tfarp2, double ffarp2,
              int tnear, double fnear, int tnearnext, double fnearnext,
              double& A, double& m, double& B, double& DD, double& c0)
{
 if (((tfarp2-tfarp1)!=(tfarp1-tfar))||(tfar==tfarp1))
    throw(std::invalid_argument("Two exponential + const -- could not compute an initial guess"));
 double cor0=ffarp1-ffar;
 double cor1=ffarp2-ffarp1;
 double s=cor0/cor1;
 double d=double(tfarp1)-double(tfar);
 if (s<=0.0)
    throw(std::invalid_argument("Two exponential + const -- could not compute an initial guess"));
 m=log(s)/d;                         // guess for m
 double rr=exp(-m*double(tfar));
 A=cor0/(rr*(exp(-m*d)-1.0));        // guess for A
 c0=ffar-A*rr;                       // guess for c0
 double r1=exp(m*double(tnear))*(fnear-c0)/A-1.0;   
 double r2=exp(m*double(tnearnext))*(fnearnext-c0)/A-1.0; 
 s=r1/r2;
 if ((s<=1.0)||(tnearnext==tnear))
    throw(std::invalid_argument("Two exponential + const -- could not compute an initial guess"));
 double DDsq=log(s)/(double(tnearnext)-double(tnear)); // guess for DD
 DD=sqrt(DDsq);
 B=exp(DDsq*double(tnear))*r1;       // guess for B
}




 // ***********************************************************************************


      // The fitting function is a sum of two exponentials, time-symmetric + constant
      //
      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-Delta^2*t) ]
      //          + exp(-m*(T_period-t)) * [ 1 + B*exp(-Delta^2*(T_period-t)) ] }  + c0
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //         c0 = fitparams[4]
      //


void TimeSymTwoExponentialPlusConstant::setupInfos(XMLHandler& xmlm, 
                                vector<MCObsInfo>& fitparam_info, int taskcount) const
{
 TimeForwardTwoExponentialPlusConstant::setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeSymTwoExponentialPlusConstant::evaluate(
                 const vector<double>& fitparams, double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[3],fitparams[2],fitparams[4],
           tval,T_period,value);
}



void TimeSymTwoExponentialPlusConstant::evalGradient(
                             const vector<double>& fitparams, double tval, 
                             vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,
           T_period,grad[1],grad[0],grad[3],grad[2],grad[4]);
}


void TimeSymTwoExponentialPlusConstant::guessInitialParamValues(
                             const vector<double>& data, const vector<uint>& tvals,
                             vector<double>& fitparams) const
{
 TimeForwardTwoExponentialPlusConstant::eval_guess(data,tvals,fitparams);
}


void TimeSymTwoExponentialPlusConstant::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeSymTwoExponentialPlusConstant");
}



void TimeSymTwoExponentialPlusConstant::setFitInfo(
                   const std::vector<MCObsInfo>& fitparams_info,
                   const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                   uint fit_tmax, bool show_approach,
                   uint meff_step, double chisq_dof, double qual,
                   TCorrFitInfo& fitinfo) const
{
 if (show_approach)
    approachSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,meff_step,chisq_dof,qual,fitinfo,true);
 else
    simpleSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,chisq_dof,qual,fitinfo);
}


      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-DD^2*t) ]
      //          + exp(-m*(Nt-t)) * [ 1 + B*exp(-DD^2*(Nt-t)) ] }   + c0

void TimeSymTwoExponentialPlusConstant::eval_func(
                   double A, double m, double B, double DD, double c0,
                   double tf, int Nt, double& funcval) const
{
 double tb=double(Nt)-tf;
 funcval=A*(exp(-m*tf)*(1.0+B*exp(-DD*DD*tf))+exp(-m*tb)*(1.0+B*exp(-DD*DD*tb)))+c0;
}


void TimeSymTwoExponentialPlusConstant::eval_grad(
                   double A, double m, double B, double DD,
                   double tf, int Nt, double& dAval, double& dmval,
                   double& dBval, double& dDDval, double& dc0val) const
{
 double gap=DD*DD;
 double r1=exp(-m*tf); 
 double r2=exp(-gap*tf);
 dAval=r1*(1.0+B*r2);
 dmval=-tf*A*dAval;
 dBval=A*r1*r2;
 dDDval=-2.0*tf*B*DD*dBval;
 double tb=double(Nt)-tf;
 r1=exp(-m*tb); r2=exp(-gap*tb);
 double s=r1*(1.0+B*r2);
 dAval+=s;
 dmval-=tb*A*s;
 s=A*r1*r2;
 dBval+=s;
 dDDval-=2.0*tb*B*DD*s;
 dc0val=1.0;
}


 // ***********************************************************************************

      // The fitting function is a sum of a geometric series of exponentials, time-forward:
      //
      //    f(t) = A * exp(-m*t) / [ 1 - B*exp(-Delta^2*t) ]
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //


void TimeForwardGeomSeriesExponential::setupInfos(XMLHandler& xmlm, 
                                vector<MCObsInfo>& fitparam_info, int taskcount) const
{
 setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeForwardGeomSeriesExponential::setup(XMLHandler& xmlm, 
                   vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount)
{
 try{
 fitparam_info.resize(nparam);
 XMLHandler xmlen(xmlm,"FirstEnergy");
 string name; int index;
 xmlreadchild(xmlen,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for first energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlen,"IDIndex",index);
 fitparam_info[0]=MCObsInfo(name,index);

 XMLHandler xmla(xmlm,"FirstAmplitude");
 name.clear();
 xmlreadchild(xmla,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for first amplitude parameter"));
 index=taskcount;
 xmlreadifchild(xmla,"IDIndex",index);
 fitparam_info[1]=MCObsInfo(name,index);

 XMLHandler xmlg(xmlm,"SqrtGapToSecondEnergy");
 name.clear();
 xmlreadchild(xmlg,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for sqrt gap to second energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlg,"IDIndex",index);
 fitparam_info[2]=MCObsInfo(name,index);

 XMLHandler xmlb(xmlm,"SecondAmplitudeRatio");
 name.clear();
 xmlreadchild(xmlb,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for second amplitude ratio parameter"));
 index=taskcount;
 xmlreadifchild(xmlb,"IDIndex",index);
 fitparam_info[3]=MCObsInfo(name,index);

 for (uint k=0;k<nparam;k++)
 for (uint l=k+1;l<nparam;l++)
    if (fitparam_info[k]==fitparam_info[l])
        throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("GeomSeriesExponential -- ")+string(errmsg.what())));}

}


void TimeForwardGeomSeriesExponential::evaluate(
            const vector<double>& fitparams, double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,value);
}


void TimeForwardGeomSeriesExponential::evalGradient(
                const vector<double>& fitparams, double tval, 
                vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,
           grad[1],grad[0],grad[3],grad[2]);
}

/*
void TimeForwardGeomSeriesExponential::guessInitialParamValues(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams) const
{
 if (data.size()<4)
    throw(std::invalid_argument("GeomSeriesExponential -- Error: at least four data points needed! in GeomSeries exponential guess"));
 for (int nd=3;nd<int(data.size()-3);nd++)
    if (TimeForwardTwoExponential::eval_guess(tmin,data[0],data[1],data[nd-1],data[nd],
                   fitparams[1],fitparams[0],fitparams[3],fitparams[2])) return;
 fitparams[2]=0.0;
 fitparams[3]=0.0;
 TimeForwardSingleExponential::eval_guess(tmin,data[2],data[3],
                                          fitparams[1],fitparams[0]);
}
*/


void TimeForwardGeomSeriesExponential::guessInitialParamValues(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams) const
{
 TimeForwardTwoExponential::eval_guess(data,tvals,fitparams);
}



void TimeForwardGeomSeriesExponential::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardGeomSeriesExponential");
}


void TimeForwardGeomSeriesExponential::setFitInfo(
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


      //    f(t) = A * exp(-m*t) / [ 1 - B*exp(-DD^2*t) ] 


void TimeForwardGeomSeriesExponential::eval_func(
              double A, double m, double B, double DD,
              double tf, double& funcval) const
{
 funcval=A*(exp(-m*tf)/(1.0-B*exp(-DD*DD*tf)));
}


void TimeForwardGeomSeriesExponential::eval_grad(
              double A, double m, double B, double DD,
              double tf, double& dAval, double& dmval,
              double& dBval, double& dDDval) const
{
 double gap=DD*DD;
 double r1=exp(-m*tf); 
 double r2=exp(-gap*tf);
 double d1=1.0-B*r2;
 dAval=r1/d1;
 dmval=-tf*A*dAval;
 dBval=A*r1*r2/(d1*d1);
 dDDval=-2.0*tf*B*DD*dBval;
}


 // ***********************************************************************************


      // The fitting function is a sum of a geometric series of exponentials, time-symmetric:
      //
      //    f(t) = A * { exp(-m*t) / [ 1 - B*exp(-Delta^2*t) ]
      //          + exp(-m*(T_period-t)) / [ 1 - B*exp(-Delta^2*(T_period-t)) ] }
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //


void TimeSymGeomSeriesExponential::setupInfos(XMLHandler& xmlm, 
                                vector<MCObsInfo>& fitparam_info, int taskcount) const
{
 TimeForwardGeomSeriesExponential::setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeSymGeomSeriesExponential::evaluate(const vector<double>& fitparams, double tval, 
                                     double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,T_period,value);
}


void TimeSymGeomSeriesExponential::evalGradient(const vector<double>& fitparams, double tval, 
                                         vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,T_period,
           grad[1],grad[0],grad[3],grad[2]);
}

/*
void TimeSymGeomSeriesExponential::guessInitialParamValues(
                              const vector<double>& data, const vector<uint>& tvals,
                              vector<double>& fitparams) const
{
 if (data.size()<4)
    throw(std::invalid_argument("GeomSeriesExponential -- Error: at least four data points needed! in GeomSeries exponential guess"));
 for (int nd=3;nd<int(data.size()-3);nd++)
    if (TimeForwardTwoExponential::eval_guess(tmin,data[0],data[1],data[nd-1],data[nd],
                   fitparams[1],fitparams[0],fitparams[3],fitparams[2])) return;
 fitparams[2]=0.0;
 fitparams[3]=0.0;
 TimeForwardSingleExponential::eval_guess(tmin,data[2],data[3],
                                          fitparams[1],fitparams[0]);
} */

void TimeSymGeomSeriesExponential::guessInitialParamValues(
                              const vector<double>& data, const vector<uint>& tvals,
                              vector<double>& fitparams) const
{
 TimeForwardTwoExponential::eval_guess(data,tvals,fitparams);
}



void TimeSymGeomSeriesExponential::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeSymGeomSeriesExponential");
}


void TimeSymGeomSeriesExponential::setFitInfo(
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



      //    f(t) = A * { exp(-m*t) / [ 1 - B*exp(-DD^2*t) ]
      //          + exp(-m*(Nt-t)) / [ 1 - B*exp(-DD^2*(Nt-t)) ] }

void TimeSymGeomSeriesExponential::eval_func(
                double A, double m, double B, double DD,
                double tf, int Nt, double& funcval) const
{
 double tb=double(Nt)-tf;
 funcval=A*(exp(-m*tf)/(1.0-B*exp(-DD*DD*tf))+exp(-m*tb)/(1.0-B*exp(-DD*DD*tb)));
}


void TimeSymGeomSeriesExponential::eval_grad(
                double A, double m, double B, double DD,
                double tf, int Nt, double& dAval, double& dmval,
                double& dBval, double& dDDval) const
{
 double gap=DD*DD;
 double r1=exp(-m*tf); 
 double r2=exp(-gap*tf);
 double d1=1.0-B*r2;
 dAval=r1/d1;
 dmval=-tf*A*dAval;
 dBval=A*r1*r2/(d1*d1);
 dDDval=-2.0*tf*B*DD*dBval;
 double tb=double(Nt)-tf;
 r1=exp(-m*tb); r2=exp(-gap*tb);
 d1=1.0-B*r2;
 double s=r1/d1;
 dAval+=s;
 dmval-=tb*A*s;
 s=A*r1*r2/(d1*d1);
 dBval+=s;
 dDDval-=2.0*tb*B*DD*s;
}


 // ***********************************************************************************

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
 else if (modeltype=="TimeForwardGMO"){
    mptr=new TimeForwardGMO(in_Tperiod);}
 else if (modeltype=="TimeForwardDoubleExpRatio"){
    mptr=new TimeForwardDoubleExpRatio(in_Tperiod);}
 else if (modeltype=="TimeForwardTwoIndExp"){
    mptr=new TimeForwardTwoIndExp(in_Tperiod);}
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

// ******************************************************************************


      // Fitting function is single exponential time-forward only:
      //
      //       f(t) = A * exp( -m*t ) 
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1].


void TimeForwardSingleExponential::setupInfos(XMLHandler& xmlm, 
                          vector<MCObsInfo>& fitparam_info, int taskcount) 
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
 get_exp_guess(tvals,data,fitparams[0],fitparams[1]);
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



   //   Given a set of time separations in "tvals" and corresponding correlator
   //   values in "corrvals", this routine finds initial "best fit" guesses
   //   for energy0, amp0, approximating the correlator by
   //
   //       corr(t) = amp0 * exp(-energy0*t) 
   //
   //   Strategy:  fit a straight line to  ln(+/-corr(t)) vs t with least squares
   //   Discards correlator values that are opposite sign from small t corr.
   //   We assume that tvals[k] < tvals[k+1].

void TimeForwardSingleExponential::get_exp_guess(
                   const std::vector<uint>& tvals, const std::vector<double>& corrvals, 
                   double& energy0, double& amp0)
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enought time vals in get_exp_guess"));
 for (uint k=0;k<(tvals.size()-1);k++)
    if (tvals[k]>=tvals[k+1]) throw(std::invalid_argument("increasing tvalues needed in get_exp_guess"));
 double sgn=(corrvals[0]>0.0)?1.0:-1.0;
 vector<double> x,y;
 for (uint k=0;k<corrvals.size();k++){
    if (sgn*corrvals[k]>0.0){
       x.push_back(tvals[k]);
       y.push_back(log(sgn*corrvals[k]));}
    else
       break;}
 uint npts=y.size();
 if (npts<2) 
    throw(std::invalid_argument("number of valid times too small for get_exp_guess"));
 double r0=npts;
 double rx=0.0,ry=0.0,rxx=0.0,rxy=0.0;    
 for (uint k=0;k<npts;k++){
    rx+=x[k]; ry+=y[k]; rxx+=x[k]*x[k]; rxy+=x[k]*y[k];}
 double den=r0*rxx-rx*rx;
 energy0=(rx*ry-r0*rxy)/den;
 amp0=sgn*exp((ry*rxx-rx*rxy)/den);
}


/*
void TimeForwardSingleExponential::get_exp_guess(int tval, double corrt, int tnext, double corrtnext, 
                                                  double& A, double& m)
{
 double s=corrt/corrtnext;
 if ((s<=0.0)||(tval==tnext))
    throw(std::invalid_argument("SingleExponential -- could not compute a guess for exponential"));
 m=log(s)/(double(tnext)-double(tval));       // guess for m
 A=exp(m*double(tval))*corrt;                 // guess for A
}
*/

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
                        vector<MCObsInfo>& fitparam_info, int taskcount)
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
 TimeForwardSingleExponential::get_exp_guess(tvals,data,fitparams[0],fitparams[1]);
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
                                      vector<MCObsInfo>& fitparam_info, int taskcount) 
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
 get_exp_plus_const_guess(tvals,data,fitparams[0],fitparams[1],fitparams[2]);
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

   //  Evaluates C(t)-C(t+step)
   //   tvals[k+1] - tvals[k] = step is required.

void TimeForwardSingleExponentialPlusConstant::take_diff(
               const std::vector<uint>& tvals, const std::vector<double>& corrvals,
               std::vector<uint>& tdf, std::vector<double>& corrdiff, uint& step)
{
 uint npts=tvals.size();
 if (corrvals.size()!=npts) 
    throw(std::invalid_argument("vector size mismatch in take_diff"));
 if (npts<3)
    throw(std::invalid_argument("Need at least 3 points for take_diff"));
 corrdiff.clear();
 tdf.clear();
 step=tvals[1]-tvals[0];
 for (uint k=0;k<(npts-1);k++){
    if ((tvals[k+1]-tvals[k])==step){
       tdf.push_back(tvals[k]); 
       corrdiff.push_back(corrvals[k]-corrvals[k+1]);}}
}



   //   Given a set of time separations in "tvals" and corresponding correlator
   //   values in "corrvals", this routine finds initial "best fit" guesses
   //   for energy0, amp0, c0  approximating the correlator by
   //
   //       corr(t) = amp0 * exp(-energy0*t) + c0
   //
   //   Strategy: - compute D(t) = C(t) - C(t+step)
   //             - fit a straight line to  ln(+/-D(t)) vs t with least squares
   //             - extract c0 from first 1/4 of points
   //   Discards correlator values that are opposite sign from small t corr.
   //    tvals[k+1] - tvals[k] = step is required.

void TimeForwardSingleExponentialPlusConstant::get_exp_plus_const_guess(
                   const std::vector<uint>& tvals, const std::vector<double>& corrvals,
                   double& energy0, double& amp0, double& c0)
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enough time vals in get_exp_plus_const_guess"));
 try{
 vector<double> corrdiff; vector<uint> tdf;
 uint step;
 take_diff(tvals,corrvals,tdf,corrdiff,step);
 TimeForwardSingleExponential::get_exp_guess(tdf,corrdiff,energy0,amp0);
 amp0/=(1.0-exp(-double(step)*energy0));
 c0=0.0;
 uint nc0=corrvals.size()/4;
 if (nc0<2) nc0=2;
 for (uint k=0;k<nc0;k++)
   c0+=corrvals[k]-amp0*exp(-energy0*tvals[k]);
 c0/=double(nc0);}
 catch(const std::exception& xp){
    throw(std::runtime_error(string("Could not compute SingleExponentialPlusConstant initial guess: ")
          + xp.what()));}
}


/*
void TimeForwardSingleExponentialPlusConstant::get_exp_plus_const_guess(
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
*/

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
                                      vector<MCObsInfo>& fitparam_info, int taskcount) 
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
 TimeForwardSingleExponentialPlusConstant::get_exp_plus_const_guess(
              tvals,data,fitparams[0],fitparams[1],fitparams[2]);

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
                                vector<MCObsInfo>& fitparam_info, int taskcount) 
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




void TimeForwardTwoExponential::guessInitialParamValues(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams) const
{
 double tasymfrac=0.33;
 get_two_exp_guess(tvals,data,fitparams[0],fitparams[1],fitparams[2],fitparams[3],tasymfrac);
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


   //   Given a set of time separations in "tvals" and corresponding correlator
   //   values in "corrvals", this routine finds initial "best fit" guesses
   //   for energy0, amp0, gapsq, and gapamp, approximating the correlator by
   //
   //       corr(t) = amp0 * exp(-energy0*t) * (1 + gapamp*exp(-gapsqrt*gapsqrt*t))
   //
   //   Key assumption: the second exponential must be negligible after "tasymfrac"
   //   of the way from the minimum to the maximum time separation.  Also,
   //   we assume that tvals[k]<tvals[k+1].
   //
   //   Strategy: - first, fit a straight line to  ln(+/-corr(t)) vs t with least squares
   //                from "tasymfrac" the way from tmin to tmax to obtain amp0, energy0;
   //             - form f = ln(+/-corr(t)) -ln(+/-amp0) + energy0 * t
   //             - fit straight line to ln(+/-f) with least squares to get gapamp, gapsqrt
   //             - if last part fails, just try gapamp=0.5, gapsqrt=0.5


void TimeForwardTwoExponential::get_two_exp_guess(
                       const std::vector<uint>& tvals, const std::vector<double>& corrvals,
                       double& energy0, double& amp0, double& gapsqrt, double& gapamp,
                       double tasymfrac)
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enough time vals in get_two_exp_guess"));
 if (tvals.size()<4)
    throw(std::invalid_argument("Need at least 4 points for get_two_exp_guess"));
 double tmin=tvals[0];
 double tmax=tvals[0];
 for (uint k=1;k<corrvals.size();k++){
    if (tvals[k]>tmax) tmax=tvals[k];
    if (tvals[k]<tmin) tmin=tvals[k];}
 if ((tmax-tmin)<4)
    throw(std::invalid_argument("Insufficient time range for get_two_exp_guess"));
 double tasym=tmin+(tmax-tmin)*tasymfrac;
 vector<uint> tv; vector<double> cv;
 for (uint k=0;k<corrvals.size();k++){
    if (double(tvals[k])>tasym){
       tv.push_back(tvals[k]);
       cv.push_back(corrvals[k]);}}
 TimeForwardSingleExponential::get_exp_guess(tv,cv,energy0,amp0);
 tv.clear(); cv.clear();
 double sgn=(amp0>0.0)?1.0:-1.0;
 double A=log(sgn*amp0);
 for (uint k=0;k<corrvals.size();k++){
    if (double(tvals[k])<tasym){
       double f=log(sgn*corrvals[k])-A+energy0*double(tvals[k]);
       tv.push_back(tvals[k]);
       cv.push_back(f);}}
  try{
     TimeForwardSingleExponential::get_exp_guess(tv,cv,gapsqrt,gapamp);
     if (gapsqrt>0.0) gapsqrt=sqrt(gapsqrt);
     else throw(std::runtime_error("invalid exponential decay"));}
  catch(const std::exception& xp){
     gapsqrt=0.5; gapamp=0.5;}
}


/*
void TimeForwardTwoExponential::get_two_exp_guess(
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


void TimeForwardTwoExponential::get_two_exp_guess(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams)
{
 if (data.size()<4)
    throw(std::invalid_argument("TwoExponential -- Error: at least four data points needed! in two exponential guess"));
 uint kfar=2*data.size()/3;
 if (kfar<2) kfar=2;
 if (kfar==data.size()-1) kfar--;
 get_two_exp_guess(tvals[kfar],data[kfar],tvals[kfar+1],data[kfar+1], 
            tvals[0],data[0],tvals[1],data[1], 
            fitparams[1],fitparams[0],fitparams[3],fitparams[2]);
}
*/

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
                                vector<MCObsInfo>& fitparam_info, int taskcount) 
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


void TimeSymTwoExponential::guessInitialParamValues(
                              const vector<double>& data, const vector<uint>& tvals,
                              vector<double>& fitparams) const
{
 double tasymfrac=0.33;
 TimeForwardTwoExponential::get_two_exp_guess(tvals,data,fitparams[0],fitparams[1],
                                              fitparams[2],fitparams[3],tasymfrac);
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
                                vector<MCObsInfo>& fitparam_info, int taskcount) 
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


void TimeForwardTwoExponentialPlusConstant::guessInitialParamValues(
                        const vector<double>& data, const vector<uint>& tvals,
                        vector<double>& fitparams) const
{
 double tasymfrac=0.33;
 get_two_exp_plus_const_guess(tvals,data,fitparams[0],fitparams[1],fitparams[2],
                              fitparams[3],fitparams[4],tasymfrac);
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
//    approachSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,meff_step,chisq_dof,qual,fitinfo,true);
    approachSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,meff_step,chisq_dof,qual,fitinfo);
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


   //   Given a set of time separations in "tvals" and corresponding correlator
   //   values in "corrvals", this routine finds initial "best fit" guesses
   //   for energy0, amp0, gapsq, gapamp, and c0, approximating the correlator by
   //
   //       corr(t) = amp0 * exp(-energy0*t) * (1 + gapamp*exp(-gapsqrt*gapsqrt*t)) + c0
   //
   //   Key assumption: the second exponential must be negligible after "tasymfrac"
   //   of the way from the minimum to the maximum time separation.  Also,
   //    tvals[k+1] - tvals[k] = step is required.
   //
   //   Strategy: - for t=tasym ... tmax  compute D(t) = C(t)-C(t+step)
   //             - fit a straight line to  ln(+/-D(t)) vs t with least squares
   //             - extract c0 from first 1/4 of these points
   //             - form f = ln(+/-(corr(t)-c0)) -ln(+/-amp0) + energy0 * t
   //             - fit straight line to ln(+/-f) with least squares to get gapamp, gapsqrt
   //             - if last part fails, just try gapamp=0.5, gapsqrt=0.5

void TimeForwardTwoExponentialPlusConstant::get_two_exp_plus_const_guess(
                       const std::vector<uint>& tvals, 
                       const std::vector<double>& corrvals,
                       double& energy0, double& amp0, double& gapsqrt, double& gapamp,
                       double& c0, double tasymfrac)
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enough time vals in get_two_exp_plus_const_guess"));
 if (tvals.size()<5)
    throw(std::invalid_argument("Need at least 5 points for get_two_exp_plus_const_guess"));
 double tmin=tvals[0];
 double tmax=tvals[0];
 for (uint k=1;k<corrvals.size();k++){
    if (tvals[k]>tmax) tmax=tvals[k];
    if (tvals[k]<tmin) tmin=tvals[k];}
 if ((tmax-tmin)<5)
    throw(std::invalid_argument("Insufficient time range for get_two_exp_plus_const_guess"));
 double tasym=tmin+(tmax-tmin)*tasymfrac;
 vector<uint> tv; vector<double> cv;
 for (uint k=0;k<corrvals.size();k++){
    if (double(tvals[k])>tasym){
       tv.push_back(tvals[k]);
       cv.push_back(corrvals[k]);}}
 TimeForwardSingleExponentialPlusConstant::get_exp_plus_const_guess(tv,cv,energy0,amp0,c0);
 tv.clear(); cv.clear();
 double sgn=(amp0>0.0)?1.0:-1.0;
 double A=log(sgn*amp0);
 for (uint k=0;k<corrvals.size();k++){
    if (double(tvals[k])<tasym){
       double f=log(sgn*(corrvals[k]-c0))-A+energy0*double(tvals[k]);
       tv.push_back(tvals[k]);
       cv.push_back(f);}}
  try{
     TimeForwardSingleExponential::get_exp_guess(tv,cv,gapsqrt,gapamp);
     if (gapsqrt>0.0) gapsqrt=sqrt(gapsqrt);
     else throw(std::runtime_error("invalid exponential decay"));}
  catch(const std::exception& xp){
     gapsqrt=0.3; gapamp=0.5;}
}


/*
void TimeForwardTwoExponentialPlusConstant::get_two_exp_plus_const_guess(
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


void TimeForwardTwoExponentialPlusConstant::get_two_exp_plus_const_guess(
                        const vector<double>& data, const vector<uint>& tvals,
                        vector<double>& fitparams)
{
 if (data.size()<5)
    throw(std::invalid_argument("TwoExponentialPlusConst -- Error: at least five data points needed! in two exponential guess"));
 uint kfar=2*data.size()/3;
 if (kfar<2) kfar=2;
 if (kfar>=data.size()-2) kfar=data.size()-3;
 get_two_exp_plus_const_guess(tvals[kfar],data[kfar],tvals[kfar+1],data[kfar+1],tvals[kfar+2],data[kfar+2],
            tvals[0],data[0],tvals[1],data[1], 
            fitparams[1],fitparams[0],fitparams[3],fitparams[2],fitparams[4]);
}
*/


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
                                vector<MCObsInfo>& fitparam_info, int taskcount) 
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
 double tasymfrac=0.33;
 TimeForwardTwoExponentialPlusConstant::get_two_exp_plus_const_guess(
               tvals,data,fitparams[0],fitparams[1],fitparams[2],
               fitparams[3],fitparams[4],tasymfrac);
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
//    approachSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,meff_step,chisq_dof,qual,fitinfo,true);
    approachSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,meff_step,chisq_dof,qual,fitinfo);
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
                                vector<MCObsInfo>& fitparam_info, int taskcount) 
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


void TimeForwardGeomSeriesExponential::guessInitialParamValues(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams) const
{
 double tasymfrac=0.33;
 TimeForwardTwoExponential::get_two_exp_guess(tvals,data,fitparams[0],fitparams[1],
                                              fitparams[2],fitparams[3],tasymfrac);
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
                                vector<MCObsInfo>& fitparam_info, int taskcount) 
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


void TimeSymGeomSeriesExponential::guessInitialParamValues(
                              const vector<double>& data, const vector<uint>& tvals,
                              vector<double>& fitparams) const
{
 double tasymfrac=0.33;
 TimeForwardTwoExponential::get_two_exp_guess(tvals,data,fitparams[0],fitparams[1],
                                              fitparams[2],fitparams[3],tasymfrac);
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


// ******************************************************************************


      // Fitting function is GMO constructed correlator:
      //
      //       f(t) = AL * AS^(1/3) exp( -( mL+mS/3-2*mN/3-2*mX/3 )*t ) /  AN^(2/3) / AX^(2/3)
      //
      // where 
      //           mL = fitparams[0]
      //           AL = fitparams[1]
      //           mS = fitparams[2]
      //           AS = fitparams[3]
      //           mN = fitparams[4]
      //           AN = fitparams[5]
      //           mX = fitparams[6]
      //           AX = fitparams[7].
      //
      // For initial guess, D200 results have been hardcoded in model_corr.cc
      //   Need to change this.


void TimeForwardGMO::setupInfos(XMLHandler& xmlm, 
                          vector<MCObsInfo>& fitparam_info, int taskcount)
{
 setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeForwardGMO::setup(XMLHandler& xmlm, 
                 vector<MCObsInfo>& fitparam_info, uint nparam, int taskcount)
{
 try{
     fitparam_info.resize(nparam);

     std::vector<std::string> particles = {"Lambda","Sigma","Nucleon","Xi"};
     for( uint i = 0; i<particles.size(); i++){
         XMLHandler xmlen(xmlm,particles[i]+"Energy");
         string name; int index;
         xmlreadchild(xmlen,"Name",name);
         if (name.empty()) throw(std::invalid_argument("Must provide name for "+particles[i]+" energy parameter."));
         index=taskcount;
         xmlreadifchild(xmlen,"IDIndex",index);
         fitparam_info[2*i]=MCObsInfo(name,index);

         XMLHandler xmla(xmlm,particles[i]+"Amplitude");
         name.clear();
         xmlreadchild(xmla,"Name",name);
         if (name.empty()) throw(std::invalid_argument("Must provide name for "+particles[i]+" amplitude parameter."));
         index=taskcount;
         xmlreadifchild(xmla,"IDIndex",index);
         fitparam_info[2*i+1]=MCObsInfo(name,index);
     }

     XMLHandler xmlminit(xmlm,"InitialEnergies");
     xmlreadchild(xmlminit,"Lambda",m_mL_init);
     xmlreadchild(xmlminit,"Sigma",m_mS_init);
     xmlreadchild(xmlminit,"Nucleon",m_mN_init);
     xmlreadchild(xmlminit,"Xi",m_mX_init);

     XMLHandler xmlAinit(xmlm,"InitialAmplitudes");
     xmlreadchild(xmlAinit,"Lambda",m_AL_init);
     xmlreadchild(xmlAinit,"Sigma",m_AS_init);
     xmlreadchild(xmlAinit,"Nucleon",m_AN_init);
     xmlreadchild(xmlAinit,"Xi",m_AX_init);

     for( uint i = 0; i<nparam; i++)
         for ( uint ii = 0; ii<i; ii++)
             if (fitparam_info[i]==fitparam_info[ii])
                 throw(std::invalid_argument("Fit parameter infos must all differ"));


 }catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("TimeForwardGMO -- ")+string(errmsg.what())));
 }
}


void TimeForwardGMO::evaluate(const vector<double>& fitparams, 
                                            double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[3],fitparams[2],fitparams[5],fitparams[4],fitparams[7],fitparams[6],tval,value);
}


void TimeForwardGMO::evalGradient(const vector<double>& fitparams, 
                         double tval, vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],fitparams[3],fitparams[2],fitparams[5],fitparams[4],fitparams[7],fitparams[6],tval,
           grad[1],grad[0],grad[3],grad[2],grad[5],grad[4],grad[7],grad[6]);
}


void TimeForwardGMO::guessInitialParamValues(
                   const vector<double>& data, const vector<uint>& tvals,
                   vector<double>& fitparams) const
{
 if (data.size()<8)
    throw(std::invalid_argument("SingleExponential -- Error: at least two data points needed! in exponential guess"));
//  get_exp_guess(tvals,data,fitparams[0],fitparams[1]);
    // sad put hardcoding single baryon fit results here. Didn't know how else
    
 //L
 fitparams[0] = m_mL_init;
 fitparams[1] = m_AL_init;
    
 //S
 fitparams[2] = m_mS_init;
 fitparams[3] = m_AS_init;
    
 //N
 fitparams[4] = m_mN_init;
 fitparams[5] = m_AN_init;
    
 //X
 fitparams[6] = m_mX_init;
 fitparams[7] = m_AX_init;
}



void TimeForwardGMO::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardGMO");
}


void TimeForwardGMO::setFitInfo(
                   const std::vector<MCObsInfo>& fitparams_info,
                   const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                   uint fit_tmax, bool show_approach,
                   uint meff_step, double chisq_dof, double qual,
                   TCorrFitInfo& fitinfo) const
{ 
 simpleSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,chisq_dof,qual,fitinfo);
}
 
void TimeForwardGMO::eval_func(double AL, double mL, double AS, double mS, double AN, double mN, 
                               double AX, double mX, double tf, double& funcval) const
{
 double one_third = 1.0/3.0;
 double two_thirds = 2.0/3.0;
 funcval=AL*pow(AS,one_third)*exp( -(mL+(one_third*mS)-(two_thirds*mN)-(two_thirds*mX)) * tf )
     /pow(AN,two_thirds)/pow(AX,two_thirds);
}


void TimeForwardGMO::eval_grad(double AL, double mL, double AS, double mS, double AN, double mN, 
                               double AX, double mX, double tf, double& dALval, double& dmLval, 
                               double& dASval, double& dmSval, double& dANval, double& dmNval, 
                               double& dAXval, double& dmXval) const
{
 double one_third = 1.0/3.0;
 double two_thirds = 2.0/3.0;
 double five_thirds = 5.0/3.0;
 double exponent = exp( -(mL+(one_third*mS)-(two_thirds*mN)-(two_thirds*mX)) * tf );
 dALval=pow(AS,one_third)*exponent/pow(AN,two_thirds)/pow(AX,two_thirds); 
 dmLval=-tf*AL*dALval;
 dASval=one_third*AL*exponent/pow(AN,two_thirds)/pow(AX,two_thirds)/pow(AS,two_thirds); 
 dmSval=-tf*AS*dASval;
 dANval=-two_thirds*AL*pow(AS,one_third)*exponent/pow(AN,five_thirds)/pow(AX,two_thirds); 
 dmNval=-tf*AN*dANval;
 dAXval=-two_thirds*AL*pow(AS,one_third)*exponent/pow(AN,two_thirds)/pow(AX,five_thirds); 
 dmXval=-tf*AX*dAXval;
}

// ******************************************************************************


      // Fitting function is single exponential time-forward only:
      //
      //       f(t) = A * exp( -m*t ) 
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1].


void TimeForwardDoubleExpRatio::setupInfos(XMLHandler& xmlm, 
                          vector<MCObsInfo>& fitparam_info, int taskcount)
{
 setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeForwardDoubleExpRatio::setup(XMLHandler& xmlm, 
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
 fitparam_info[4]=MCObsInfo(name,index);
 fitparam_info.resize(nparam);
     
 XMLHandler xmlenN(xmlm,"NumGap");
 xmlreadchild(xmlenN,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlenN,"IDIndex",index);
 fitparam_info[1]=MCObsInfo(name,index);

 XMLHandler xmlaN(xmlm,"NumGapAmp");
 name.clear();
 xmlreadchild(xmlaN,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for amplitude parameter"));
 index=taskcount;
 xmlreadifchild(xmlaN,"IDIndex",index);
 fitparam_info[5]=MCObsInfo(name,index);
     
 XMLHandler xmlenSH1(xmlm,"SH1Gap");
 xmlreadchild(xmlenSH1,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlenSH1,"IDIndex",index);
 fitparam_info[2]=MCObsInfo(name,index);

 XMLHandler xmlaSH1(xmlm,"SH1GapAmp");
 name.clear();
 xmlreadchild(xmlaSH1,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for amplitude parameter"));
 index=taskcount;
 xmlreadifchild(xmlaSH1,"IDIndex",index);
 fitparam_info[6]=MCObsInfo(name,index);
     
 XMLHandler xmlenSH2(xmlm,"SH2Gap");
 xmlreadchild(xmlenSH2,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlenSH2,"IDIndex",index);
 fitparam_info[3]=MCObsInfo(name,index);

 XMLHandler xmlaSH2(xmlm,"SH2GapAmp");
 name.clear();
 xmlreadchild(xmlaSH2,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for amplitude parameter"));
 index=taskcount;
 xmlreadifchild(xmlaSH2,"IDIndex",index);
 fitparam_info[7]=MCObsInfo(name,index);

 for( uint i = 0; i<8; i++)
     for( uint j = 0; j<i; j++)
         if (fitparam_info[i]==fitparam_info[j])
             throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("TimeForwardDoubleExpRatio -- ")+string(errmsg.what())));}
}


void TimeForwardDoubleExpRatio::evaluate(const vector<double>& fitparams, 
                                            double tval, double& value) const
{
 eval_func(fitparams[4],fitparams[0],fitparams[5],fitparams[1],
           fitparams[6],fitparams[2],fitparams[7],fitparams[3],tval,value);
}


void TimeForwardDoubleExpRatio::evalGradient(const vector<double>& fitparams, 
                         double tval, vector<double>& grad) const
{
 eval_grad(fitparams[4],fitparams[0],fitparams[5],fitparams[1],
           fitparams[6],fitparams[2],fitparams[7],fitparams[3],tval,
           grad[4],grad[0],grad[5],grad[1],
           grad[6],grad[2],grad[7],grad[3]);
}


void TimeForwardDoubleExpRatio::guessInitialParamValues(
                   const vector<double>& data, const vector<uint>& tvals,
                   vector<double>& fitparams) const
{
 double tasymfrac=0.33;
 TimeForwardTwoExponential::get_two_exp_guess(tvals,data,fitparams[0],fitparams[4],fitparams[1],fitparams[5],tasymfrac);
 fitparams[2]=fitparams[1];
 fitparams[3]=fitparams[1];
 fitparams[6]=fitparams[5];
 fitparams[7]=fitparams[5]; //??????
}



void TimeForwardDoubleExpRatio::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardDoubleExpRatio");
}


void TimeForwardDoubleExpRatio::setFitInfo(
                   const std::vector<MCObsInfo>& fitparams_info,
                   const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                   uint fit_tmax, bool show_approach,
                   uint meff_step, double chisq_dof, double qual,
                   TCorrFitInfo& fitinfo) const
{ 
 simpleSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,chisq_dof,qual,fitinfo);
}
 
void TimeForwardDoubleExpRatio::eval_func(double A, double m, double AN, double mN, double ASH1, double mSH1, 
                                          double ASH2, double mSH2, double tf, double& funcval) const
{
 funcval=A*exp(-m*tf) * ( 1.0 + AN*exp(-mN*mN*tf) ) / ( (1.0+ASH1*exp(-mSH1*mSH1*tf)) * (1.0+ASH2*exp(-mSH2*mSH2*tf)) );
}


void TimeForwardDoubleExpRatio::eval_grad(double A, double m, double AN, double mN, double ASH1, double mSH1, 
                                          double ASH2, double mSH2, double tf, double& dAval, double& dmval,
                                         double& dANval, double& dmNval,double& dASH1val, double& dmSH1val,
                                         double& dASH2val, double& dmSH2val) const
{
 double rN = exp(-mN*mN*tf);
 double rSH1 = exp(-mSH1*mSH1*tf);
 double rSH2 = exp(-mSH2*mSH2*tf);
 dAval=exp(-m*tf) * ( 1.0 + AN*rN ) / ( (1.0+ASH1*rSH1) * (1.0+ASH2*rSH2) );
 dmval=-tf*A*dAval;
 dANval = rN/ ( (1.0+ASH1*rSH1) * (1.0+ASH2*rSH2) );
 dmNval = -2.0*mN*AN*tf*dANval;
 dASH1val = -A*rSH1*exp(-m*tf) * ( 1.0 + AN*rN ) / ( (1.0+ASH1*rSH1) * (1.0+ASH1*rSH1) * (1.0+ASH2*rSH2) );
 dmSH1val = -2.0*mSH1*ASH1*tf*dASH1val;
 dASH2val = -A*rSH2*exp(-m*tf) * ( 1.0 + AN*rN ) / ( (1.0+ASH1*rSH1) * (1.0+ASH2*rSH2) * (1.0+ASH2*rSH2) );
 dmSH2val = -2.0*mSH2*ASH2*tf*dASH2val;
}

// ******************************************************************************

      // Fitting function is independent double exponential time-forward only:
      //
      //       f(t) = A * exp( -m*t ) + A1 * exp( -m1*t )
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1]
      //           m1 = fitparams[0]
      //           A1 = fitparams[1].
      //
      // For initial guess, need corr[tmin], corr[tmin+1]


void TimeForwardTwoIndExp::setupInfos(XMLHandler& xmlm, 
                          vector<MCObsInfo>& fitparam_info, int taskcount) 
{
 setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeForwardTwoIndExp::setup(XMLHandler& xmlm, 
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

 fitparam_info.resize(nparam);
 XMLHandler xmlen1(xmlm,"Energy1");
 xmlreadchild(xmlen1,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for second energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlen1,"IDIndex",index);
 fitparam_info[2]=MCObsInfo(name,index);

 XMLHandler xmla1(xmlm,"Amplitude1");
 name.clear();
 xmlreadchild(xmla1,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for second amplitude parameter"));
 index=taskcount;
 xmlreadifchild(xmla1,"IDIndex",index);
 fitparam_info[3]=MCObsInfo(name,index);
     
 for( uint i = 0; i<4; i++)
     for( uint j = 0; j<i; j++)
         if (fitparam_info[i]==fitparam_info[j])
             throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("TimeForwardTwoIndExp -- ")+string(errmsg.what())));}
}


void TimeForwardTwoIndExp::evaluate(const vector<double>& fitparams, 
                                            double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,value);
}


void TimeForwardTwoIndExp::evalGradient(const vector<double>& fitparams, 
                         double tval, vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],fitparams[3],fitparams[2],tval,grad[1],grad[0],grad[3],grad[2]);
}


void TimeForwardTwoIndExp::guessInitialParamValues(
                   const vector<double>& data, const vector<uint>& tvals,
                   vector<double>& fitparams) const
{
 double tasymfrac=0.33;
 TimeForwardTwoExponential::get_two_exp_guess(tvals,data,fitparams[0],fitparams[1],fitparams[2],fitparams[3],tasymfrac);
 fitparams[2] = (fitparams[2]*fitparams[2])+fitparams[0]; 
 fitparams[3] = fitparams[1]*fitparams[3];
}



void TimeForwardTwoIndExp::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardTwoIndExp");
}


void TimeForwardTwoIndExp::setFitInfo(
                   const std::vector<MCObsInfo>& fitparams_info,
                   const std::vector<MCEstimate>& fitparams, uint fit_tmin,
                   uint fit_tmax, bool show_approach,
                   uint meff_step, double chisq_dof, double qual,
                   TCorrFitInfo& fitinfo) const
{ 
 simpleSetFitInfo(fitparams_info,fitparams,fit_tmin,fit_tmax,chisq_dof,qual,fitinfo);
}
 
void TimeForwardTwoIndExp::eval_func(double A, double m, double A1, double m1, double tf, double& funcval) const
{
 funcval=A*exp(-m*tf) + A1*exp(-m1*tf);
}


void TimeForwardTwoIndExp::eval_grad(double A, double m, double A1, double m1, double tf, double& dAval, 
                                             double& dmval, double& dA1val, double& dm1val) const
{
 dAval=exp(-m*tf); 
 dmval=-tf*A*dAval;
 dA1val=exp(-m1*tf); 
 dm1val=-tf*A1*dA1val;
}

// ******************************************************************************


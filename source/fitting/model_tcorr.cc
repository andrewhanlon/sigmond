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
 else if (modeltype=="TimeForwardGeomSeriesTwoExponential"){
    mptr=new TimeForwardGeomSeriesTwoExponential(in_Tperiod);}
 else if (modeltype=="TimeForwardGeomSeriesExtraFixedExp"){
    mptr=new TimeForwardGeomSeriesExtraFixedExp(in_Tperiod);}
 else if (modeltype=="TimeForwardGeomSeriesSTI"){
    mptr=new TimeForwardGeomSeriesSTI(in_Tperiod);}
 else if (modeltype=="TimeForwardTwoExponentialSTI"){
    mptr=new TimeForwardTwoExponentialSTI(in_Tperiod);}
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
 fitinfo.meff_approach.clear();
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
       fitinfo.meff_approach.push_back(XYPoint(tval+shift,yval));
       if ((fabs(yval-fitinfo.energy_mean)<=2.0*fitinfo.energy_err)
            &&(tval<tmin)) tmin=tval;}
    tval+=curvestep;}
 if (tmin<(fit_tmax-4))
    fitinfo.tmin=tmin;
 else
    fitinfo.tmin=fit_tmax-4;
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
// double tasymfrac=0.33;
// get_two_exp_guess(tvals,data,fitparams[0],fitparams[1],fitparams[2],fitparams[3],tasymfrac);
 double tasymfrac=0.75;
 get_two_exp_guess_method2(tvals,data,fitparams[0],fitparams[1],fitparams[2],fitparams[3],tasymfrac);
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


void TimeForwardTwoExponential::get_two_exp_guess_method2(
                       const std::vector<uint>& tvals, const std::vector<double>& corrvals,
                       double& energy0, double& amp0, double& gapsqrt, double& gapamp,
                       double tasymfrac)
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enough time vals in get_two_exp_guess_method2"));
 if (tvals.size()<5)
    throw(std::invalid_argument("Need at least 5 points for get_two_exp_guess_method2"));
 double tmin=tvals[0];
 double tmax=tvals[0];
 for (uint k=1;k<corrvals.size();k++){
    if (tvals[k]>tmax) tmax=tvals[k];
    if (tvals[k]<tmin) tmin=tvals[k];}
 if ((tmax-tmin)<5)
    throw(std::invalid_argument("Insufficient time range for get_two_exp_guess_method2"));
    // determine step for solving
 try{
 double tasym=tmin+(tmax-tmin)*tasymfrac;
 double step=floor((tasym-tmin)/5.0);
 if (step<1.01) step=1.0;
    // getting correlator values at t=tmin+k*step for k=0,1,2,3,4
 vector<double> cv;
 for (uint k=0;k<5;++k){
    double tget=tmin+k*step;
    int iget=-1;
    for (uint kk=0;kk<tvals.size();++kk){
       if (std::abs(tvals[kk]-tget)<1e-9){
          iget=kk; cv.push_back(corrvals[kk]); break;}}
    if (iget<0){
       throw(std::invalid_argument("Could not find appropriate correlator values for get_two_exp_guess_method2"));}}
    // now solve for the 4 parameters
 solve_two_exp(tmin,step,cv,energy0,amp0,gapsqrt,gapamp);}
 catch(const std::exception& xp){
    get_two_exp_guess(tvals,corrvals,energy0,amp0,gapsqrt,gapamp,1.0-tasymfrac);}
}


//  Given five values of the correlator  
//       C(tmin), C(tmin+step), C(tmin+2*step), C(tmin+3*step), C(tmin+4*step)
//  in "corrvals", this solves for the four parameters
//
//     C(t) = amp0 * exp( - energy * t ) * ( 1 + amp1ratio * exp( - sqrtgap^2 * t )  )
//
//  The first four correlator values are used to get one or more solutions.  Then
//  the solution that gets closest to the fifth correlator value is used. Maple
//  was used to find expressions for the solutions.

void TimeForwardTwoExponential::solve_two_exp(double tmin, double step, const std::vector<double>& corrvals,
                                              double& energy, double& amp0, double& sqrtgap, double& amp1ratio)
{
 if (corrvals.size()!=5)
    throw(std::invalid_argument("invalid size of corrvals in solve_two_exp"));
 double C0=corrvals[0];
 double C1=corrvals[1];
 double C2=corrvals[2];
 double C3=corrvals[3];
 double a2=C0*C2-C1*C1;
 double a1=C1*C2-C0*C3;
 double a0=C1*C3-C2*C2;
 double disc=a1*a1-4.0*a2*a0;
 if (disc<0)
    throw(std::invalid_argument("no solution in solve_two_exp"));
 disc=sqrt(disc);
 std::vector<double> en,DD;
 double signs[2]={1.0,-1.0};
 for (uint i=0;i<2;++i){
    double xx=(-a1+signs[i]*disc)/(2.0*a2);
    if (xx>0.0){
       double yy=a2*(C2-C1*xx)/((C1*a2+C0*a1)*xx+a0*C0);
       if (yy>0.0){
          double res=-log(yy)/step;
          if (res>0.0){
             en.push_back(-log(xx)/step);
             DD.push_back(sqrt(res));}}}}
 if (en.size()==0)
    throw(std::invalid_argument("no solution in solve_two_exp"));
 //for (uint i=0;i<en.size();++i)
 //   cout << en[i]<<"  "<<DD[i]<<endl;
 std::vector<double> aa0(en.size()),aa1(en.size());
 for (uint i=0;i<en.size();++i){
    double d0=exp(-en[i]*tmin);
    double dd0=exp(-en[i]*(tmin+step));
    double f0=exp(-DD[i]*DD[i]*tmin);
    double ff0=exp(-DD[i]*DD[i]*(tmin+step));
    aa0[i]=(C1*d0*f0-C0*dd0*ff0)/(d0*dd0*(f0-ff0)); 
    aa1[i]=(C1*d0-C0*dd0)/(C0*dd0*ff0-C1*d0*f0);}
    // now choose the solution that gets closest to C(tmin+5*step)
 double C4=corrvals[4];
 int ibest=0; 
 double closest=aa0[ibest]*exp(-en[ibest]*(tmin+5*step))
               *(1.0+aa1[ibest]*exp(-DD[ibest]*DD[ibest]*(tmin+5*step)));
 for (uint i=1;i<en.size();++i){
    double Ctry=aa0[i]*exp(-en[i]*(tmin+5*step))
               *(1.0+aa1[i]*exp(-DD[i]*DD[i]*(tmin+5*step)));
    if (std::abs(Ctry-C4)<std::abs(closest-C4)){
       ibest=i; closest=Ctry;}}
 energy=en[ibest];
 amp0=aa0[ibest]; 
 sqrtgap=DD[ibest];
 amp1ratio=aa1[ibest];
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
 TimeForwardTwoExponential::get_two_exp_guess_method2(tvals,data,fitparams[0],fitparams[1],
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
// double tasymfrac=0.33;
// TimeForwardTwoExponential::get_geomsum_exp_guess(tvals,data,fitparams[0],fitparams[1],
//                                                  fitparams[2],fitparams[3],tasymfrac);
 double tasymfrac=0.75;
 get_geomsum_guess(tvals,data,fitparams[0],fitparams[1],fitparams[2],fitparams[3],tasymfrac);
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


void TimeForwardGeomSeriesExponential::get_geomsum_guess(
                       const std::vector<uint>& tvals, const std::vector<double>& corrvals,
                       double& energy0, double& amp0, double& gapsqrt, double& gapamp,
                       double tasymfrac)
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enough time vals in get_geomsum_guess"));
 if (tvals.size()<5)
    throw(std::invalid_argument("Need at least 5 points for get_geomsum_guess"));
 double tmin=tvals[0];
 double tmax=tvals[0];
 for (uint k=1;k<corrvals.size();k++){
    if (tvals[k]>tmax) tmax=tvals[k];
    if (tvals[k]<tmin) tmin=tvals[k];}
 if ((tmax-tmin)<5)
    throw(std::invalid_argument("Insufficient time range for get_geomsum_guess"));
    // determine step for solving
 try{
 double tasym=tmin+(tmax-tmin)*tasymfrac;
 double step=floor((tasym-tmin)/5.0);
 if (step<1.01) step=1.0;
    // getting correlator values at t=tmin+k*step for k=0,1,2,3,4
 vector<double> cv;
 for (uint k=0;k<5;++k){
    double tget=tmin+k*step;
    int iget=-1;
    for (uint kk=0;kk<tvals.size();++kk){
       if (std::abs(tvals[kk]-tget)<1e-9){
          iget=kk; cv.push_back(corrvals[kk]); break;}}
    if (iget<0){
       throw(std::invalid_argument("Could not find appropriate correlator values for get_geomsum_guess"));}}
    // now solve for the 4 parameters
 solve_geomsum_exp(tmin,step,cv,energy0,amp0,gapsqrt,gapamp);}
 catch(const std::exception& xp){
    TimeForwardTwoExponential::get_two_exp_guess(tvals,corrvals,energy0,amp0,gapsqrt,gapamp,1.0-tasymfrac);}
}



//  Given five values of the correlator  
//       C(tmin), C(tmin+step), C(tmin+2*step), C(tmin+3*step), C(tmin+4*step)
//  in "corrvals", this solves for the four parameters
//
//     C(t) = amp0 * exp( - energy * t ) / ( 1 - amp1ratio * exp( - sqrtgap^2 * t )  )
//
//  The first four correlator values are used to get one or more solutions.  Then
//  the solution that gets closest to the fifth correlator value is used. Maple
//  was used to find expressions for the solutions.

void TimeForwardGeomSeriesExponential::solve_geomsum_exp(double tmin, double step, 
                       const std::vector<double>& corrvals,
                       double& energy, double& amp0, double& sqrtgap, double& amp1ratio)
{
 if (corrvals.size()!=5)
    throw(std::invalid_argument("invalid size of corrvals in solve_two_exp"));
 double C0=corrvals[0];
 double C1=corrvals[1];
 double C2=corrvals[2];
 double C3=corrvals[3];
 double a2=C0*C1*(C1*C3-C2*C2);
 double a1=-C1*C2*(C0*C3-C1*C2);
 double a0=C2*C3*(C0*C2-C1*C1);
 double C1sq=C1*C1;
 double C1cb=C1*C1sq;
 double C2sq=C2*C2;
 double C2cb=C2*C2sq;
 double C0sq=C0*C0;
 double disc=a1*a1-4.0*a2*a0;
 if (disc<0)
    throw(std::invalid_argument("no solution in solve_two_exp"));
 disc=sqrt(disc);
 std::vector<double> en,DD;
 double signs[2]={1.0,-1.0};
 for (uint i=0;i<2;++i){
    double xx=(-a1+signs[i]*disc)/(2.0*a2);
    if (xx>0.0){
       double yy=C3*(C1sq-C0*C2)*(((C1sq*C3-2.0*C1*C2sq)*C0+C1cb*C2)*xx+C0*C2cb-C1cb*C3)/
                ((C1*C3-C2sq)*((C1*C2*C3*C0sq+(-2.0*C1cb*C3+C1sq*C2sq)*C0)*xx
                -C0sq*C2sq*C3+C0*C1sq*C2*C3+C1sq*C1sq*C3-C1cb*C2sq));
       if (yy>0.0){
          double res=-log(yy)/step;
          if (res>0.0){
             en.push_back(-log(xx)/step);
             DD.push_back(sqrt(res));}}}}
 if (en.size()==0)
    throw(std::invalid_argument("no solution in solve_two_exp"));
 std::vector<double> aa0(en.size()),aa1(en.size());
 for (uint i=0;i<en.size();++i){
    double d0=exp(-en[i]*tmin);
    double dd0=exp(-en[i]*(tmin+step));
    double f0=exp(-DD[i]*DD[i]*tmin);
    double ff0=exp(-DD[i]*DD[i]*(tmin+step));
    double denom=C0*dd0*f0-C1*d0*ff0;
    aa0[i]=C0*C1*(f0-ff0)/denom;
    aa1[i]=(C0*dd0-C1*d0)/denom;}
    // now choose the solution that gets closest to C(tmin+5*step)
 double C4=corrvals[4];
 int ibest=0; 
 double closest=aa0[ibest]*exp(-en[ibest]*(tmin+5*step))
               /(1.0-aa1[ibest]*exp(-DD[ibest]*DD[ibest]*(tmin+5*step)));
 for (uint i=1;i<en.size();++i){
    double Ctry=aa0[i]*exp(-en[i]*(tmin+5*step))
               /(1.0-aa1[i]*exp(-DD[i]*DD[i]*(tmin+5*step)));
    if (std::abs(Ctry-C4)<std::abs(closest-C4)){
       ibest=i; closest=Ctry;}}
 energy=en[ibest];
 amp0=aa0[ibest]; 
 sqrtgap=DD[ibest];
 amp1ratio=aa1[ibest];
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


 // ***********************************************************************************

      // The fitting function is a sum of a geometric series of exponentials, time-forward:
      //
      //    f(t) = A * exp(-m*t) / [ 1 - B1*exp(-Delta1^2*t) * (1 - B2 * exp(-Delta2^2*t) ) ]
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //     Delta1 = fitparams[2]
      //         B1 = fitparams[3]
      //     Delta2 = fitparams[4]
      //         B2 = fitparams[5]
      //


void TimeForwardGeomSeriesTwoExponential::setupInfos(XMLHandler& xmlm, 
                                vector<MCObsInfo>& fitparam_info, int taskcount)
{
 setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeForwardGeomSeriesTwoExponential::setup(XMLHandler& xmlm, 
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

 XMLHandler xmlg2(xmlm,"SqrtGapToThirdEnergy");
 name.clear();
 xmlreadchild(xmlg2,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for sqrt gap to third energy parameter"));
 index=taskcount;
 xmlreadifchild(xmlg2,"IDIndex",index);
 fitparam_info[4]=MCObsInfo(name,index);

 XMLHandler xmlb2(xmlm,"ThirdAmplitudeRatio");
 name.clear();
 xmlreadchild(xmlb2,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for third amplitude ratio parameter"));
 index=taskcount;
 xmlreadifchild(xmlb2,"IDIndex",index);
 fitparam_info[5]=MCObsInfo(name,index);

 for (uint k=0;k<nparam;k++)
 for (uint l=k+1;l<nparam;l++)
    if (fitparam_info[k]==fitparam_info[l])
        throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("GeomSeriesTwoExponential -- ")+string(errmsg.what())));}

}


void TimeForwardGeomSeriesTwoExponential::evaluate(
            const vector<double>& fitparams, double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[3],fitparams[2],fitparams[5],fitparams[4],tval,value);
}


void TimeForwardGeomSeriesTwoExponential::evalGradient(
                const vector<double>& fitparams, double tval, 
                vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],fitparams[3],fitparams[2],fitparams[5],fitparams[4],tval,
           grad[1],grad[0],grad[3],grad[2],grad[5],grad[4]);
}


void TimeForwardGeomSeriesTwoExponential::guessInitialParamValues(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams) const
{
 double tasym1frac=0.33;
 double tasym2frac=0.75;
 get_geomsum2X_guess(tvals,data,fitparams[0],fitparams[1],fitparams[2],fitparams[3],
                     fitparams[4],fitparams[5],tasym1frac,tasym2frac);
}



void TimeForwardGeomSeriesTwoExponential::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardGeomSeriesTwoExponential");
}


void TimeForwardGeomSeriesTwoExponential::setFitInfo(
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


      //    f(t) = A * exp(-m*t) / [ 1 - B1*exp(-DD1^2*t) * (1 + B2 * exp(-DD2^2*t) )  ] 


void TimeForwardGeomSeriesTwoExponential::eval_func(
              double A, double m, double B1, double DD1, double B2, double DD2,
              double tf, double& funcval) const
{
 funcval=A*exp(-m*tf)/(1.0-B1*exp(-DD1*DD1*tf)*(1.0+B2*exp(-DD2*DD2*tf)));
}


void TimeForwardGeomSeriesTwoExponential::eval_grad(
              double A, double m, double B1, double DD1, double B2, double DD2,
              double tf, double& dAval, double& dmval,
              double& dB1val, double& dDD1val,
              double& dB2val, double& dDD2val) const
{
 double gap1=DD1*DD1;
 double gap2=DD2*DD2;
 double r1=exp(-m*tf); 
 double r2=exp(-gap1*tf);
 double r3=exp(-gap2*tf);
 double d2=1.0+B2*r3;
 double d1=1.0-B1*r2*d2;
 dAval=r1/d1;
 dmval=-tf*A*dAval;
 double dB=A*r1*r2/(d1*d1);
 dB1val=d2*dB;
 dDD1val=-2.0*tf*B1*DD1*dB1val;
 dB2val=B1*r3*dB;
 dDD2val=-2.0*tf*B2*DD2*dB2val;
}


void TimeForwardGeomSeriesTwoExponential::get_geomsum2X_guess(
                       const std::vector<uint>& tvals, const std::vector<double>& corrvals,
                       double& energy0, double& amp0, double& gapsqrt1, double& gapamp1,
                       double& gapsqrt2, double& gapamp2, double tasym1frac, double tasym2frac)
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enough time vals in get_geomsum2X_guess"));
 if (tvals.size()<7)
    throw(std::invalid_argument("Need at least 7 points for get_geomsum2X_guess"));
 double tmin=tvals[0];
 double tmax=tvals[0];
 for (uint k=1;k<corrvals.size();k++){
    if (tvals[k]<tvals[k-1]){
       throw(std::invalid_argument("time values must be in ascending order in get_geomsum2X_guess"));}
    if (tvals[k]>tmax) tmax=tvals[k];
    if (tvals[k]<tmin) tmin=tvals[k];}
 if ((tmax-tmin)<7)
    throw(std::invalid_argument("Insufficient time range for get_geomsum2X_guess"));
    // determine step for solving
 try{
 double tasym1=ceil(tmin+(tmax-tmin)*tasym1frac);
 if (tasym1<tvals[2]) tasym1=tvals[2];
 double tasym2=tasym1+(tmax-tasym1)*tasym2frac;
 double step=floor((tasym2-tasym1)/5.0);
 if (step<1.01) step=1.0;
    // getting correlator values at t=tmin+k*step for k=0,1,2,3,4
 vector<double> cv;
 for (uint k=0;k<5;++k){
    double tget=tasym1+k*step;
    int iget=-1;
    for (uint kk=0;kk<tvals.size();++kk){
       if (std::abs(tvals[kk]-tget)<1e-9){
          iget=kk; cv.push_back(corrvals[kk]); break;}}
    if (iget<0){
       throw(std::invalid_argument("Could not find appropriate correlator values for get_geomsum2X_guess"));}}
    // now solve for the 4 parameters
 TimeForwardGeomSeriesExponential::solve_geomsum_exp(tasym1,step,cv,energy0,amp0,gapsqrt1,gapamp1);
    // now solve for the last 2 parameters using first 2 correlators values
 double C0=corrvals[0]; double t0=tvals[0];
 double C1=corrvals[1]; double t1=tvals[1];
 double d0=(1.0-(amp0/C0)*exp(-energy0*t0))/(gapamp1*exp(-gapsqrt1*gapsqrt1*t0));
 double d1=(1.0-(amp0/C1)*exp(-energy0*t1))/(gapamp1*exp(-gapsqrt1*gapsqrt1*t1));
 double r=d1/d0;
 if ((r>0.0)&&(r<1.0)){
    gapsqrt2=sqrt(log(r)/(t0-t1));
    gapamp2=d0*exp(gapsqrt2*gapsqrt2*t0);}
 else{
    gapamp2=0.0; gapsqrt2=1.0;}}
 catch(const std::exception& xp){
    TimeForwardTwoExponential::get_two_exp_guess(tvals,corrvals,energy0,amp0,gapsqrt1,gapamp1,tasym1frac);
    gapamp2=0.0; gapsqrt2=1.0;}
}

 // ***********************************************************************************


      // The fitting function is a sum of a geometric series of exponentials, time-forward:
      //
      //    f(t) = A * exp(-m*t) / [ 1 - B1*exp(-Delta1^2*t) * (1 - B2 * exp(-Delta2*t) ) ]
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //     Delta1 = fitparams[2]
      //         B1 = fitparams[3]
      //         B2 = fitparams[4]
      //
      //   HERE, Delta2 is **set by the user****


void TimeForwardGeomSeriesExtraFixedExp::setupInfos(XMLHandler& xmlm, 
                                vector<MCObsInfo>& fitparam_info, int taskcount)
{
 setup(xmlm,fitparam_info,m_nparams,taskcount);
}


void TimeForwardGeomSeriesExtraFixedExp::setup(XMLHandler& xmlm, 
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

 XMLHandler xmlb2(xmlm,"ThirdAmplitudeRatio");
 name.clear();
 xmlreadchild(xmlb2,"Name",name);
 if (name.empty()) throw(std::invalid_argument("Must provide name for third amplitude ratio parameter"));
 index=taskcount;
 xmlreadifchild(xmlb2,"IDIndex",index);
 fitparam_info[4]=MCObsInfo(name,index);

    // read the value of the fixed parameter
 xmlreadchild(xmlm,"GapToThirdEnergy",m_Delta2);

 for (uint k=0;k<nparam;k++)
 for (uint l=k+1;l<nparam;l++)
    if (fitparam_info[k]==fitparam_info[l])
        throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("GeomSeriesTwoExponential -- ")+string(errmsg.what())));}

}


void TimeForwardGeomSeriesExtraFixedExp::evaluate(
            const vector<double>& fitparams, double tval, double& value) const
{
 eval_func(fitparams[1],fitparams[0],fitparams[3],fitparams[2],fitparams[4],tval,value);
}


void TimeForwardGeomSeriesExtraFixedExp::evalGradient(
                const vector<double>& fitparams, double tval, 
                vector<double>& grad) const
{
 eval_grad(fitparams[1],fitparams[0],fitparams[3],fitparams[2],fitparams[4],tval,
           grad[1],grad[0],grad[3],grad[2],grad[4]);
}


void TimeForwardGeomSeriesExtraFixedExp::guessInitialParamValues(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams) const
{
 double tasym1frac=0.33;
 double tasym2frac=0.75;
 get_geomsum2X_guess(tvals,data,fitparams[0],fitparams[1],fitparams[2],fitparams[3],
                     fitparams[4],tasym1frac,tasym2frac);
}



void TimeForwardGeomSeriesExtraFixedExp::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardGeomSeriesExtraFixedExp");
}


void TimeForwardGeomSeriesExtraFixedExp::setFitInfo(
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


      //    f(t) = A * exp(-m*t) / [ 1 - B1*exp(-DD1^2*t) * (1 + B2 * exp(-DD2^2*t) )  ] 


void TimeForwardGeomSeriesExtraFixedExp::eval_func(
              double A, double m, double B1, double DD1, double B2,
              double tf, double& funcval) const
{
 funcval=A*exp(-m*tf)/(1.0-B1*exp(-DD1*DD1*tf)*(1.0+B2*exp(-m_Delta2*tf)));
}


void TimeForwardGeomSeriesExtraFixedExp::eval_grad(
              double A, double m, double B1, double DD1, double B2,
              double tf, double& dAval, double& dmval,
              double& dB1val, double& dDD1val, double& dB2val) const
{
 double gap1=DD1*DD1;
 double r1=exp(-m*tf); 
 double r2=exp(-gap1*tf);
 double r3=exp(-m_Delta2*tf);
 double d2=1.0+B2*r3;
 double d1=1.0-B1*r2*d2;
 dAval=r1/d1;
 dmval=-tf*A*dAval;
 double dB=A*r1*r2/(d1*d1);
 dB1val=d2*dB;
 dDD1val=-2.0*tf*B1*DD1*dB1val;
 dB2val=B1*r3*dB;
}


void TimeForwardGeomSeriesExtraFixedExp::get_geomsum2X_guess(
                       const std::vector<uint>& tvals, const std::vector<double>& corrvals,
                       double& energy0, double& amp0, double& gapsqrt1, double& gapamp1,
                       double& gapamp2, double tasym1frac, double tasym2frac) const
{
 if (tvals.size()<corrvals.size()) 
    throw(std::invalid_argument("not enough time vals in get_geomsum2X_guess"));
 if (tvals.size()<6)
    throw(std::invalid_argument("Need at least 6 points for get_geomsum2X_guess"));
 double tmin=tvals[0];
 double tmax=tvals[0];
 for (uint k=1;k<corrvals.size();k++){
    if (tvals[k]<tvals[k-1]){
       throw(std::invalid_argument("time values must be in ascending order in get_geomsum2X_guess"));}
    if (tvals[k]>tmax) tmax=tvals[k];
    if (tvals[k]<tmin) tmin=tvals[k];}
 if ((tmax-tmin)<6)
    throw(std::invalid_argument("Insufficient time range for get_geomsum2X_guess"));
    // determine step for solving
 try{
 double tasym1=ceil(tmin+(tmax-tmin)*tasym1frac);
 if (tasym1<tvals[2]) tasym1=tvals[2];
 double tasym2=tasym1+(tmax-tasym1)*tasym2frac;
 double step=floor((tasym2-tasym1)/5.0);
 if (step<1.01) step=1.0;
    // getting correlator values at t=tmin+k*step for k=0,1,2,3,4
 vector<double> cv;
 for (uint k=0;k<5;++k){
    double tget=tasym1+k*step;
    int iget=-1;
    for (uint kk=0;kk<tvals.size();++kk){
       if (std::abs(tvals[kk]-tget)<1e-9){
          iget=kk; cv.push_back(corrvals[kk]); break;}}
    if (iget<0){
       throw(std::invalid_argument("Could not find appropriate correlator values for get_geomsum2X_guess"));}}
    // now solve for the 4 parameters
 TimeForwardGeomSeriesExponential::solve_geomsum_exp(tasym1,step,cv,energy0,amp0,gapsqrt1,gapamp1);
    // now solve for the last parameters using first correlator value
 double C0=corrvals[0]; double t0=tvals[0];
 double d0=(1.0-(amp0/C0)*exp(-energy0*t0))/(gapamp1*exp(-gapsqrt1*gapsqrt1*t0));
 gapamp2=d0*exp(m_Delta2*t0);}
 catch(const std::exception& xp){
    TimeForwardTwoExponential::get_two_exp_guess(tvals,corrvals,energy0,amp0,gapsqrt1,gapamp1,tasym1frac);
    gapamp2=0.0;}
}

 // ***********************************************************************************


      // Fitting function is sum of a geometric series of exponentials time-forward only
      // with short time improvement (STI)
      //
      //    f(t) = A * exp(-m*t) / [ 1 - B1*exp(-Delta1^2*t) * (1 + ff(t) ) ]
      //
      //    ff(t) = sum[ C_k exp( - (STIgap + (k-1) * STIdelta) * t ),  k = 1..nSTI]  
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //     Delta1 = fitparams[2]
      //         B1 = fitparams[3]
      //        C_1 = fitparams[4]
      //          ....
      //     C_nSTI = fitparams[3+nSTI]
      //
      //   HERE, STIgap, STIdelta, and nSTI are **set by the user****
      // For initial guess, need 4 corr values


void TimeForwardGeomSeriesSTI::setupInfos(XMLHandler& xmlm, 
                                vector<MCObsInfo>& fitparam_info, int taskcount)
{
 setup(xmlm,fitparam_info,taskcount);
}


void TimeForwardGeomSeriesSTI::setup(XMLHandler& xmlm, 
                   vector<MCObsInfo>& fitparam_info, int taskcount)
{
 try{
 int nsti;
 xmlreadchild(xmlm,"NumberOfSTIEnergies",nsti);
 if ((nsti<=0)||(nsti>24)){
    throw(std::invalid_argument("Must provide number of STI energies between 1 and 24"));}
 m_nSTI=nsti;
 m_nparams=m_nSTI+4;

 xmlreadchild(xmlm,"STIEnergyGap",m_STIgap);
 if (m_STIgap<=0.1){
    throw(std::invalid_argument("Must provide STI energy gap > 0.1"));}

 xmlreadchild(xmlm,"STIEnergyStep",m_STIdelta);
 if (m_STIdelta<=0.05){
    throw(std::invalid_argument("Must provide STI energy step > 0.05"));}

 m_geomTime=0;
 xmlreadifchild(xmlm,"GeomTime",m_geomTime);

 fitparam_info.resize(m_nparams);
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

 for (uint k=0;k<m_nSTI;++k){
    string tagname("STIAmplitude");
    tagname=tagname+int_to_string(k+1);
    XMLHandler xmlb2(xmlm,tagname);
    name.clear();
    xmlreadchild(xmlb2,"Name",name);
    if (name.empty()) throw(std::invalid_argument("Must provide name for each STI amplitude"));
    index=taskcount;
    xmlreadifchild(xmlb2,"IDIndex",index);
    fitparam_info[4+k]=MCObsInfo(name,index);}

 for (uint k=0;k<m_nparams;k++)
 for (uint l=k+1;l<m_nparams;l++)
    if (fitparam_info[k]==fitparam_info[l])
        throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("GeomSeriesSTI -- ")+string(errmsg.what())));}

}


void TimeForwardGeomSeriesSTI::evaluate(
            const vector<double>& fitparams, double tval, double& value) const
{
 double ff=0.0;
 for (uint k=0;k<m_nSTI;++k){
    ff+=fitparams[4+k]*exp(-(m_STIgap+k*m_STIdelta)*tval);}
 double A=fitparams[1];
 double m=fitparams[0];
 double B1=fitparams[3];
 double DD1=fitparams[2];
 value=A*exp(-m*tval)/(1.0-B1*exp(-DD1*DD1*tval)*(1.0+ff));
}


void TimeForwardGeomSeriesSTI::evalGradient(
                const vector<double>& fitparams, double tval, 
                vector<double>& grad) const
{
 double ff=0.0;
 vector<double> rr(m_nSTI);
 for (uint k=0;k<m_nSTI;++k){
    rr[k]=exp(-(m_STIgap+k*m_STIdelta)*tval);
    ff+=fitparams[4+k]*rr[k];}
 double A=fitparams[1];
 double m=fitparams[0];
 double B1=fitparams[3];
 double DD1=fitparams[2];
 double gap1=DD1*DD1;
 double r1=exp(-m*tval); 
 double r2=exp(-gap1*tval);
 double d1=1.0-B1*r2*(1.0+ff);
 double dAval=r1/d1;
 double dmval=-tval*A*dAval;
 double dden=A*r1*r2/(d1*d1);
 double dB1val=dden*(1.0+ff);
 double dDD1val=-2.0*tval*B1*DD1*dB1val;
 dden*=B1;
 grad.resize(4+m_nSTI);
 grad[0]=dmval;
 grad[1]=dAval;
 grad[2]=dDD1val;
 grad[3]=dB1val;
 for (uint k=0;k<m_nSTI;++k){
    grad[4+k]=dden*rr[k];}
}


void TimeForwardGeomSeriesSTI::guessInitialParamValues(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams) const
{
 double tasymfrac=0.75;
 if (m_geomTime<=0){
    TimeForwardGeomSeriesExponential::get_geomsum_guess(tvals,data,fitparams[0],
                        fitparams[1],fitparams[2],fitparams[3],tasymfrac);}
 else{
    vector<double> geomdata;
    vector<uint> geomtvals;
    for (uint k=0;k<tvals.size();++k){
       if (tvals[k]>=m_geomTime){
          geomdata.push_back(data[k]);
          geomtvals.push_back(tvals[k]);}}
       TimeForwardGeomSeriesExponential::get_geomsum_guess(geomtvals,geomdata,fitparams[0],
                           fitparams[1],fitparams[2],fitparams[3],tasymfrac);}
 for (uint k=0;k<m_nSTI;++k){
    fitparams[4+k]=0.0;}
}



void TimeForwardGeomSeriesSTI::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardGeomSeriesSTI");
}


void TimeForwardGeomSeriesSTI::setFitInfo(
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

 // ***********************************************************************************


      // Fitting function is sum of two exponentials time-forward only
      // with short time improvement (STI)
      //
      //    f(t) = A * exp(-m*t) * [ 1 + B1*exp(-Delta1^2*t) * (1 + ff(t) ) ]
      //
      //    ff(t) = sum[ C_k exp( - (STIgap + (k-1) * STIdelta) * t ),  k = 1..nSTI]  
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //     Delta1 = fitparams[2]
      //         B1 = fitparams[3]
      //        C_1 = fitparams[4]
      //          ....
      //     C_nSTI = fitparams[3+nSTI]
      //
      //   HERE, STIgap, STIdelta, and nSTI are **set by the user****
      // For initial guess, need 4 corr values


void TimeForwardTwoExponentialSTI::setupInfos(XMLHandler& xmlm, 
                                vector<MCObsInfo>& fitparam_info, int taskcount)
{
 setup(xmlm,fitparam_info,taskcount);
}


void TimeForwardTwoExponentialSTI::setup(XMLHandler& xmlm, 
                   vector<MCObsInfo>& fitparam_info, int taskcount)
{
 try{
 int nsti;
 xmlreadchild(xmlm,"NumberOfSTIEnergies",nsti);
 if ((nsti<=0)||(nsti>24)){
    throw(std::invalid_argument("Must provide number of STI energies between 1 and 24"));}
 m_nSTI=nsti;
 m_nparams=m_nSTI+4;

 xmlreadchild(xmlm,"STIEnergyGap",m_STIgap);
 if (m_STIgap<=0.1){
    throw(std::invalid_argument("Must provide STI energy gap > 0.1"));}

 xmlreadchild(xmlm,"STIEnergyStep",m_STIdelta);
 if (m_STIdelta<=0.05){
    throw(std::invalid_argument("Must provide STI energy step > 0.05"));}

 m_twoExpTime=0;
 xmlreadifchild(xmlm,"TwoExpTime",m_twoExpTime);

 fitparam_info.resize(m_nparams);
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

 for (uint k=0;k<m_nSTI;++k){
    string tagname("STIAmplitude");
    tagname=tagname+int_to_string(k+1);
    XMLHandler xmlb2(xmlm,tagname);
    name.clear();
    xmlreadchild(xmlb2,"Name",name);
    if (name.empty()) throw(std::invalid_argument("Must provide name for each STI amplitude"));
    index=taskcount;
    xmlreadifchild(xmlb2,"IDIndex",index);
    fitparam_info[4+k]=MCObsInfo(name,index);}

 for (uint k=0;k<m_nparams;k++)
 for (uint l=k+1;l<m_nparams;l++)
    if (fitparam_info[k]==fitparam_info[l])
        throw(std::invalid_argument("Fit parameter infos must all differ"));}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("TwoExponentialSTI -- ")+string(errmsg.what())));}

}


void TimeForwardTwoExponentialSTI::evaluate(
            const vector<double>& fitparams, double tval, double& value) const
{
 double ff=0.0;
 for (uint k=0;k<m_nSTI;++k){
    ff+=fitparams[4+k]*exp(-(m_STIgap+k*m_STIdelta)*tval);}
 double A=fitparams[1];
 double m=fitparams[0];
 double B1=fitparams[3];
 double DD1=fitparams[2];
 value=A*exp(-m*tval)*(1.0+B1*exp(-DD1*DD1*tval)*(1.0+ff));
}


void TimeForwardTwoExponentialSTI::evalGradient(
                const vector<double>& fitparams, double tval, 
                vector<double>& grad) const
{
 double ff=0.0;
 vector<double> rr(m_nSTI);
 for (uint k=0;k<m_nSTI;++k){
    rr[k]=exp(-(m_STIgap+k*m_STIdelta)*tval);
    ff+=fitparams[4+k]*rr[k];}
 double A=fitparams[1];
 double m=fitparams[0];
 double B1=fitparams[3];
 double DD1=fitparams[2];
 double gap1=DD1*DD1;
 double r1=exp(-m*tval); 
 double r2=exp(-gap1*tval);
 double d1=1.0+B1*r2*(1.0+ff);
 double dAval=r1*d1;
 double dmval=-tval*A*dAval;
 double dB1val=A*r1*r2*(1.0+ff);
 double dDD1val=-2.0*tval*B1*DD1*dB1val;
 double dp=A*r1*r2*B1;
 grad.resize(4+m_nSTI);
 grad[0]=dmval;
 grad[1]=dAval;
 grad[2]=dDD1val;
 grad[3]=dB1val;
 for (uint k=0;k<m_nSTI;++k){
    grad[4+k]=dp*rr[k];}
}


void TimeForwardTwoExponentialSTI::guessInitialParamValues(
                     const vector<double>& data, const vector<uint>& tvals,
                     vector<double>& fitparams) const
{
 double tasymfrac=0.75;
 if (m_twoExpTime<=0){
    TimeForwardTwoExponential::get_two_exp_guess_method2(tvals,data,
                 fitparams[0],fitparams[1],fitparams[2],fitparams[3],tasymfrac);}
 else{
    vector<double> twoexpdata;
    vector<uint> twoexptvals;
    for (uint k=0;k<tvals.size();++k){
       if (tvals[k]>=m_twoExpTime){
          twoexpdata.push_back(data[k]);
          twoexptvals.push_back(tvals[k]);}}
       TimeForwardTwoExponential::get_two_exp_guess_method2(twoexptvals,twoexpdata,fitparams[0],
                           fitparams[1],fitparams[2],fitparams[3],tasymfrac);}
 for (uint k=0;k<m_nSTI;++k){
    fitparams[4+k]=0.0;}
}



void TimeForwardTwoExponentialSTI::output_tag(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardTwoExponentialSTI");
}


void TimeForwardTwoExponentialSTI::setFitInfo(
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

 // ***********************************************************************************

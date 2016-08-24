#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include "tmatrix.h"
#include "model_funcs.h"
using namespace std;
using namespace ModelFitFunctions;



class RealTemporalCorrelatorFit 
{
    uint m_tmin, m_tmax, T_period;

 public:

    RealTemporalCorrelatorFit(int tmin, int tmax, int Nt) 
      :  m_tmin(tmin), m_tmax(tmax), T_period(Nt) {}



    void exponential_timeforward(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;
    void exponential_timeforward_grad(const std::vector<double>& fitparams,
                                      RMatrix& gradients) const;
    void exponential_guess(std::vector<double>& fitparams, 
                           const RVector& obs_mean) const;


    void exponential_timesym(const std::vector<double>& fitparams,
                             std::vector<double>& modelpoints) const;
    void exponential_timesym_grad(const std::vector<double>& fitparams,
                                  RMatrix& gradients) const;



    void exponential_timeforward_plus_const(const std::vector<double>& fitparams,
                                            std::vector<double>& modelpoints) const;
    void exponential_timeforward_plus_const_grad(const std::vector<double>& fitparams,
                                                 RMatrix& gradients) const;
    void exponential_plus_const_guess(std::vector<double>& fitparams, 
                                      const RVector& obs_mean) const;



    void exponential_timesym_plus_const(const std::vector<double>& fitparams,
                                        std::vector<double>& modelpoints) const;
    void exponential_timesym_plus_const_grad(const std::vector<double>& fitparams,
                                             RMatrix& gradients) const;



    void two_exponential_timeforward(const std::vector<double>& fitparams,
                                     std::vector<double>& modelpoints) const;
    void two_exponential_timeforward_grad(const std::vector<double>& fitparams,
                                          RMatrix& gradients) const;
    void two_exponential_guess(std::vector<double>& fitparams, 
                               const RVector& obs_mean) const;


    void two_exponential_timesym(const std::vector<double>& fitparams,
                                 std::vector<double>& modelpoints) const;
    void two_exponential_timesym_grad(const std::vector<double>& fitparams,
                                      RMatrix& gradients) const;


void func_two_exponential_timeforward_with_const(double A, double m, double B, double DD, double c0,
                                                 int t, double& funcval);

void grad_two_exponential_timeforward_with_const(double A, double m, double B, double DD,
                                                 int t, double& dAval, double& dmval,
                                                 double& dBval, double& dDDval, double& dc0val);

void guess_two_exponential_with_const(int tval, double f0, double f1, double f2, 
                                      double f3, double f4, double& A, double& m, 
                                      double& B, double& DD, double& c0);

void func_two_exponential_timesym_with_const(double A, double m, double B, double DD, double c0,
                                             int t, int Nt, double& funcval);

void grad_two_exponential_timesym_with_const(double A, double m, double B, double DD, 
                                             int t, int Nt, double& dAval, double& dmval,
                                             double& dBval, double& dDDval, double& dc0val);

};





 // ********************************************************************

      // Fitting function is single exponential time-forward only:
      //
      //       f(t) = A * exp( -m*t ) 
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1].

/*
void RealTemporalCorrelatorFit::exponential_setup(XMLHandler& xmlm)
{
 XMLHandler xmlen(xmlm,"Energy");
 string name; int index;
 xmlreadchild(xmlen,"Name",name);
 if (name.empty()) throw(string("Must provide name for energy parameter"));
 index=0;
 xmlreadifchild(xmlen,"IDIndex",index);
 m_fitparam_info[0]=MCObsInfo(string("Energy_")+name,index);

 XMLHandler xmla(xmlm,"Amplitude");
 name.clear();
 xmlreadchild(xmla,"Name",name);
 if (name.empty()) throw(string("Must provide name for amplitude parameter"));
 index=0;
 xmlreadifchild(xmla,"IDIndex",index);
 m_fitparam_info[1]=MCObsInfo(string("Amplitude_")+name,index);
}
*/

void RealTemporalCorrelatorFit::exponential_timeforward(
                          const vector<double>& fitparams,
                          vector<double>& modelpoints) const
{
 double m=fitparams[0];       // energy of lowest-lying state
 double A=fitparams[1];        // coefficient related to wavefunction
 for (int tt=m_tmin;tt<=int(m_tmax);++tt)
    func_exponential_timeforward(A,m,tt,modelpoints[tt-m_tmin]);
}

     
void RealTemporalCorrelatorFit::exponential_timeforward_grad(
                         const vector<double>& fitparams,
                         RMatrix& gradients) const
{
 double m=fitparams[0];       // energy of lowest-lying state
 double A=fitparams[1];       // coefficient related to wavefunction
 for (int k=0, tt=m_tmin; tt<=int(m_tmax);++k,++tt){
    grad_exponential_timeforward(A,m,tt,gradients(k,1),gradients(k,0));}
}

/*
void RealTemporalCorrelatorFit::exponential_timeforward_output(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardSingleExponential");
}
*/

      // Returns in "fitparam" guesses for the parameters,
      // to be used as starting values in the minimization routine.

void RealTemporalCorrelatorFit::exponential_guess(vector<double>& fitparams, 
                                                  const RVector& obs_mean) const
{
 if (obs_mean.size()<2)
    throw(string("Error: at least two data points needed! in exponential guess"));
 guess_exponential(m_tmin,obs_mean[0],obs_mean[1],fitparams[1],fitparams[0]);
}


// ******************************************************************


      // Fitting function is single exponential time-symmetric
      // (forwards and backwards):
      //
      //       f(t) = A * {exp( -m*t ) +  exp( -m*(T_period-t) )}
      //
      // where 
      //           m = fitparams[0]
      //           A = fitparams[1].

void RealTemporalCorrelatorFit::exponential_timesym(
                          const vector<double>& fitparams,
                          vector<double>& modelpoints) const
{
 double m=fitparams[0];        // energy of lowest-lying state
 double A=fitparams[1];        // amplitude
 for (int tt=m_tmin;tt<=int(m_tmax);++tt)
    func_exponential_timesym(A,m,tt,T_period,modelpoints[tt-m_tmin]);
}

     
void RealTemporalCorrelatorFit::exponential_timesym_grad(
                         const vector<double>& fitparams,
                         RMatrix& gradients) const
{
 double m=fitparams[0];       // energy of lowest-lying state
 double A=fitparams[1];       // coefficient related to wavefunction
 for (int k=0, tt=m_tmin; tt<=int(m_tmax);++k,++tt){
    grad_exponential_timesym(A,m,tt,T_period,gradients(k,1),gradients(k,0));}
}

/*
void RealTemporalCorrelatorFit::exponential_timesym_output(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeSymSingleExponential");
}
*/

// *******************************************************************


      // The fitting function is a single exponential time-forward
      // with an added constant:
      //
      //       f(t) = A * exp( -m*t ) + c0
      //
      // where 
      //           m = fitparam[0]
      //           A = fitparam[1]
      //          c0 = fitparam[2]
      //

/*
void RealTemporalCorrelatorFit::exponential_plus_const_setup(XMLHandler& xmlm)
{
 XMLHandler xmlen(xmlm,"Energy");
 string name; int index;
 xmlreadchild(xmlen,"Name",name);
 if (name.empty()) throw(string("Must provide name for energy parameter"));
 index=0;
 xmlreadifchild(xmlen,"IDIndex",index);
 m_fitparam_info[0]=MCObsInfo(string("Energy_")+name,index);

 XMLHandler xmla(xmlm,"Amplitude");
 name.clear();
 xmlreadchild(xmla,"Name",name);
 if (name.empty()) throw(string("Must provide name for amplitude parameter"));
 index=0;
 xmlreadifchild(xmla,"IDIndex",index);
 m_fitparam_info[1]=MCObsInfo(string("Amplitude_")+name,index);

 XMLHandler xmlc(xmlm,"AddedConstant");
 name.clear();
 xmlreadchild(xmlc,"Name",name);
 if (name.empty()) throw(string("Must provide name for amplitude parameter"));
 index=0;
 xmlreadifchild(xmlc,"IDIndex",index);
 m_fitparam_info[2]=MCObsInfo(string("AddConstant_")+name,index);

}
*/

void RealTemporalCorrelatorFit::exponential_timeforward_plus_const(
                          const vector<double>& fitparams,
                          vector<double>& modelpoints) const
{
 double m=fitparams[0];        // energy of lowest-lying state
 double A=fitparams[1];        // amplitude
 double c0=fitparams[2];       // added constant
 for (int tt=m_tmin;tt<=int(m_tmax);++tt)
    func_exponential_timeforward_with_const(A,m,c0,tt,modelpoints[tt-m_tmin]);
}

     
void RealTemporalCorrelatorFit::exponential_timeforward_plus_const_grad(
                         const vector<double>& fitparams,
                         RMatrix& gradients) const
{
 double m=fitparams[0];       // energy of lowest-lying state
 double A=fitparams[1];       // coefficient related to wavefunction
 for (int k=0, tt=m_tmin; tt<=int(m_tmax);++k,++tt){
    grad_exponential_timeforward_with_const(A,m,tt,gradients(k,1),
                            gradients(k,0),gradients(k,2));}
}


void RealTemporalCorrelatorFit::exponential_plus_const_guess(vector<double>& fitparams, 
                                                             const RVector& obs_mean) const
{
 if (obs_mean.size()<3)
    throw(string("Error: at least three data points needed! in exponential+const guess"));
 guess_exponential_plus_const(m_tmin,obs_mean[0],obs_mean[1],obs_mean[2],
                              fitparams[1],fitparams[0],fitparams[2]);
}

/*
void RealTemporalCorrelatorFit::exponential_timeforward_plus_const_output(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardSingleExponentialPlusConstant");
}
*/


// ******************************************************************


      // The fitting function is a single exponential time-symmetric
      // (forwards and backwards) with an added constant:
      //
      //       f(t) = A * {exp( -m*t ) +  exp( -m*(T_period-t) )} + c0
      //
      // where 
      //           m = fitparam[0]
      //           A = fitparam[1]
      //          c0 = fitparam[2]
      //


void RealTemporalCorrelatorFit::exponential_timesym_plus_const(
                          const vector<double>& fitparams,
                          vector<double>& modelpoints) const
{
 double m=fitparams[0];        // energy of lowest-lying state
 double A=fitparams[1];        // amplitude
 double c0=fitparams[2];       // added constant
 for (int tt=m_tmin;tt<=int(m_tmax);++tt)
    func_exponential_timesym_with_const(A,m,c0,tt,T_period,modelpoints[tt-m_tmin]);
}

    
void RealTemporalCorrelatorFit::exponential_timesym_plus_const_grad(
                         const vector<double>& fitparams,
                         RMatrix& gradients) const
{
 double m=fitparams[0];       // energy of lowest-lying state
 double A=fitparams[1];       // coefficient related to wavefunction
 for (int k=0, tt=m_tmin; tt<=int(m_tmax);++k,++tt){
    grad_exponential_timesym_with_const(A,m,tt,T_period,gradients(k,1),
                      gradients(k,0),gradients(k,2));}
}

/*
void RealTemporalCorrelatorFit::exponential_timesym_plus_const_output(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeSymSingleExponentialPlusConstant");
}
*/

// ******************************************************************

      // The fitting function is a sum of two exponentials, time-forward:
      //
      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-Delta^2*t) ]
      //
      //  where 
      //          m = fitparams[0]
      //          A = fitparams[1]
      //      Delta = fitparams[2]
      //          B = fitparams[3]
      //

/*
void RealTemporalCorrelatorFit::two_exponential_setup(XMLHandler& xmlm)
{
 XMLHandler xmlen(xmlm,"FirstEnergy");
 string name; int index;
 xmlreadchild(xmlen,"Name",name);
 if (name.empty()) throw(string("Must provide name for first energy parameter"));
 index=0;
 xmlreadifchild(xmlen,"IDIndex",index);
 m_fitparam_info[0]=MCObsInfo(string("FirstEnergy_")+name,index);

 XMLHandler xmla(xmlm,"FirstAmplitude");
 name.clear();
 xmlreadchild(xmla,"Name",name);
 if (name.empty()) throw(string("Must provide name for first amplitude parameter"));
 index=0;
 xmlreadifchild(xmla,"IDIndex",index);
 m_fitparam_info[1]=MCObsInfo(string("FirstAmplitude_")+name,index);

 XMLHandler xmlg(xmlm,"SqrtGapToSecondEnergy");
 name.clear();
 xmlreadchild(xmlg,"Name",name);
 if (name.empty()) throw(string("Must provide name for sqrt gap to second energy parameter"));
 index=0;
 xmlreadifchild(xmlg,"IDIndex",index);
 m_fitparam_info[2]=MCObsInfo(string("SqrtGapToSecondEnergy_")+name,index);

 XMLHandler xmlb(xmlm,"SecondAmplitudeRatio");
 name.clear();
 xmlreadchild(xmlb,"Name",name);
 if (name.empty()) throw(string("Must provide name for second amplitude ratio parameter"));
 index=0;
 xmlreadifchild(xmlb,"IDIndex",index);
 m_fitparam_info[3]=MCObsInfo(string("SecondAmplitudeRatio_")+name,index);

}
*/

void RealTemporalCorrelatorFit::two_exponential_timeforward(
                          const vector<double>& fitparams,
                          vector<double>& modelpoints) const
{
 double m=fitparams[0];        // energy of lowest-lying state
 double A=fitparams[1];        // amplitude
 double delta=fitparams[2];    // sqrt root of gap to the next exponential
 double B=fitparams[3];       // coefficient of the second exponential
 for (int tt=m_tmin;tt<=int(m_tmax);++tt)
    func_two_exponential_timeforward(A,m,B,delta,tt,modelpoints[tt-m_tmin]);
}


void RealTemporalCorrelatorFit::two_exponential_timeforward_grad(
                         const vector<double>& fitparams,
                         RMatrix& gradients) const
{
 double m=fitparams[0];        // energy of lowest-lying state
 double A=fitparams[1];        // amplitude
 double delta=fitparams[2];    // sqrt root of gap to the next exponential
 double B=fitparams[3];       // coefficient of the second exponential
 for (int k=0, tt=m_tmin; tt<=int(m_tmax);++k,++tt){
    grad_two_exponential_timeforward(A,m,B,delta,tt,gradients(k,1),
                gradients(k,0),gradients(k,3),gradients(k,2));}
}



void RealTemporalCorrelatorFit::two_exponential_guess(vector<double>& fitparams, 
                                                      const RVector& obs_mean) const
{
 if (obs_mean.size()<4)
    throw(string("Error: at least four data points needed! in two exponential guess"));
 guess_two_exponential(m_tmin,obs_mean[0],obs_mean[1],obs_mean[2],obs_mean[3],
                       fitparams[1],fitparams[0],fitparams[3],fitparams[2]);
}

/*
void RealTemporalCorrelatorFit::two_exponential_timeforward_output(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeForwardTwoExponential");
}
*/

// ******************************************************************


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


void RealTemporalCorrelatorFit::two_exponential_timesym(
                          const vector<double>& fitparams,
                          vector<double>& modelpoints) const
{
 double m=fitparams[0];        // energy of lowest-lying state
 double A=fitparams[1];        // amplitude
 double delta=fitparams[2];    // sqrt root of gap to the next exponential
 double B=fitparams[3];       // coefficient of the second exponential

 for (int tt=m_tmin;tt<=int(m_tmax);++tt)
    func_two_exponential_timesym(A,m,B,delta,tt,T_period,modelpoints[tt-m_tmin]);
}



void RealTemporalCorrelatorFit::two_exponential_timesym_grad(
                         const vector<double>& fitparams,
                         RMatrix& gradients) const
{
 double m=fitparams[0];        // energy of lowest-lying state
 double A=fitparams[1];        // amplitude
 double delta=fitparams[2];    // sqrt root of gap to the next exponential
 double B=fitparams[3];       // coefficient of the second exponential
 for (int k=0, tt=m_tmin; tt<=int(m_tmax);++k,++tt){
    grad_two_exponential_timesym(A,m,B,delta,tt,T_period,gradients(k,1),
                gradients(k,0),gradients(k,3),gradients(k,2));}
}


/*
void RealTemporalCorrelatorFit::two_exponential_timesym_output(XMLHandler& xmlout) const
{
 xmlout.set_root("Model","TimeSymTwoExponential");
}
*/

// *********************************************************************


double func(double A, double back, double m, double B, double DD, double c0,
            int t, int Nt)
{
 double tf=double(t);
 double tb=double(Nt-t);
 return      A*exp(-m*tf)*(1.0+B*exp(-DD*DD*tf))
       +back*A*exp(-m*tb)*(1.0+B*exp(-DD*DD*tb))
       +c0;
}

double dfunc_dA(double A, double back, double m, double B, double DD, double c0,
                int t, int Nt)
{
 double tf=double(t);
 double tb=double(Nt-t);
 return       exp(-m*tf)*(1.0+B*exp(-DD*DD*tf))
        +back*exp(-m*tb)*(1.0+B*exp(-DD*DD*tb));
}


double dfunc_dm(double A, double back, double m, double B, double DD, double c0,
                int t, int Nt)
{
 double tf=double(t);
 double tb=double(Nt-t);
 return      -A*tf*exp(-m*tf)*(1+B*exp(-DD*DD*tf))
        -back*A*tb*exp(-m*tb)*(1+B*exp(-DD*DD*tb));
}


double dfunc_dB(double A, double back, double m, double B, double DD, double c0,
                int t, int Nt)
{
 double tf=double(t);
 double tb=double(Nt-t);
 return       A*exp(-m*tf)*exp(-DD*DD*tf)
        +back*A*exp(-m*tb)*exp(-DD*DD*tb);
}


double dfunc_dDD(double A, double back, double m, double B, double DD, double c0,
                 int t, int Nt)
{
 double tf=double(t);
 double tb=double(Nt-t);
 return       -2*A*exp(-m*tf)*B*DD*tf*exp(-DD*DD*tf)
        +back*(-2*A*exp(-m*tb)*B*DD*tb*exp(-DD*DD*tb));
}


double dfunc_dc0(double A, double back, double m, double B, double DD, double c0,
                 int t, int Nt)
{
 return 1.0;
}


typedef double(*funcptr)(double,double,double,double,double,double,int,int);


int main(){

int tmin=3;
int tmax=62;
int Tperiod=64;
int npoints=tmax-tmin+1;
vector<double> modelpoints(npoints);

RealTemporalCorrelatorFit Tester(tmin,tmax,Tperiod);

double Avalue=567.923;
double mvalue=0.1325;

double Bvalue=6.23;
double DDvalue=0.879;

double c0value=-1232.874;

double A,m,B,DD,c0,back;

int nparam;
vector<double> fitparams;
vector<funcptr> derivs;
RMatrix grad;

cout << "Testing time forward single exp"<<endl;

nparam=2;
fitparams.resize(nparam);
fitparams[0]=m=mvalue;
fitparams[1]=A=Avalue;
back=0.0;
B=0.0; DD=0.0;
c0=0.0;
derivs.resize(nparam);
derivs[0]=dfunc_dm;
derivs[1]=dfunc_dA;
grad.resize(npoints,nparam);

Tester.exponential_timeforward(fitparams,modelpoints);
for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
Tester.exponential_timeforward_grad(fitparams,grad);
for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}



cout << "Testing time sym single exp"<<endl;
back=1.0;
Tester.exponential_timesym(fitparams,modelpoints);
for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
Tester.exponential_timesym_grad(fitparams,grad);
for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}


cout << "Testing time forward single exp with constant"<<endl;
nparam=3;
fitparams.resize(nparam);
fitparams[0]=mvalue; m=mvalue;
fitparams[1]=Avalue; A=Avalue;
fitparams[2]=c0value; c0=c0value;
back=0.0;
B=0.0; DD=0.0;
derivs.resize(nparam);
derivs[0]=dfunc_dm;
derivs[1]=dfunc_dA;
derivs[2]=dfunc_dc0;
grad.resize(npoints,nparam);


Tester.exponential_timeforward_plus_const(fitparams,modelpoints);
for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
Tester.exponential_timeforward_plus_const_grad(fitparams,grad);
for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}



cout << "Testing time symmetric single exp with constant"<<endl;
nparam=3;
fitparams.resize(nparam);
fitparams[0]=mvalue; m=mvalue;
fitparams[1]=Avalue; A=Avalue;
fitparams[2]=c0value; c0=c0value;
back=1.0;
B=0.0; DD=0.0;
derivs.resize(nparam);
derivs[0]=dfunc_dm;
derivs[1]=dfunc_dA;
derivs[2]=dfunc_dc0;
grad.resize(npoints,nparam);


Tester.exponential_timesym_plus_const(fitparams,modelpoints);
for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
Tester.exponential_timesym_plus_const_grad(fitparams,grad);
for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}


cout << "Testing time forward two exp"<<endl;

nparam=4;
fitparams.resize(nparam);
fitparams[0]=m=mvalue;
fitparams[1]=A=Avalue;
fitparams[2]=DD=DDvalue;
fitparams[3]=B=Bvalue;
back=0.0;
c0=0.0;
derivs.resize(nparam);
derivs[0]=dfunc_dm;
derivs[1]=dfunc_dA;
derivs[2]=dfunc_dDD;
derivs[3]=dfunc_dB;
grad.resize(npoints,nparam);

Tester.two_exponential_timeforward(fitparams,modelpoints);
for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
Tester.two_exponential_timeforward_grad(fitparams,grad);
for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}


cout << "Testing time sym two exp"<<endl;

nparam=4;
fitparams.resize(nparam);
fitparams[0]=m=mvalue;
fitparams[1]=A=Avalue;
fitparams[2]=DD=DDvalue;
fitparams[3]=B=Bvalue;
back=1.0;
c0=0.0;
derivs.resize(nparam);
derivs[0]=dfunc_dm;
derivs[1]=dfunc_dA;
derivs[2]=dfunc_dDD;
derivs[3]=dfunc_dB;
grad.resize(npoints,nparam);

Tester.two_exponential_timesym(fitparams,modelpoints);
for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
Tester.two_exponential_timesym_grad(fitparams,grad);
for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}


RVector obs_mean(npoints);
A=Avalue; back=0.0; m=mvalue; B=0.0; DD=0.0; c0=0.0;
for (int t=tmin, k=0; t<=tmax;++t,k++)
   obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
nparam=2;
fitparams.resize(nparam);
Tester.exponential_guess(fitparams,obs_mean);
for (int p=0;p<nparam;p++) cout << "guess["<<p<<"] = "<<fitparams[p]<<endl;

A=Avalue; back=0.0; m=mvalue; B=0.0; DD=0.0; c0=c0value;
for (int t=tmin, k=0; t<=tmax;++t,k++)
   obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
nparam=3;
fitparams.resize(nparam);
Tester.exponential_plus_const_guess(fitparams,obs_mean);
for (int p=0;p<nparam;p++) cout << "guess["<<p<<"] = "<<fitparams[p]<<endl;


A=Avalue; back=0.0; m=mvalue; B=Bvalue; DD=DDvalue; c0=0.0;
for (int t=tmin, k=0; t<=tmax;++t,k++)
   obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
nparam=4;
fitparams.resize(nparam);
Tester.two_exponential_guess(fitparams,obs_mean);
for (int p=0;p<nparam;p++) cout << "guess["<<p<<"] = "<<fitparams[p]<<endl;


A=Avalue; back=0.0; m=mvalue; B=Bvalue; DD=DDvalue; c0=22.1;
for (int t=tmin, k=0; t<=tmax;++t,k++)
   obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
nparam=5;
fitparams.resize(nparam);
Tester.two_exponential_with_const_guess(fitparams,obs_mean);
for (int p=0;p<nparam;p++) cout << "guess["<<p<<"] = "<<fitparams[p]<<endl;


  // ************************************************************************


int tvalue=14; int Nt=24;
Avalue=3.536; mvalue=0.567; Bvalue=2.313; DDvalue=0.231; c0value=23.45;
double result;
double dAval,dmval,dBval,dDDval,dc0val;

      //    f(t) = A * exp(-m*t)  

func_exponential_timeforward(Avalue,mvalue,tvalue,result);
grad_exponential_timeforward(Avalue,mvalue,tvalue,dAval,dmval);

cout << "Test 1"<<endl;
A=Avalue; back=0.0; m=mvalue; B=0.0;  DD=0.0; c0=0.0;
cout << func(A,back,m,B,DD,c0,tvalue,Nt)<<"  "<<result<<endl;
cout << dAval<<"  "<< dfunc_dA(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dmval<<"  "<< dfunc_dm(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dBval<<"  "<< dfunc_dB(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dDDval<<"  "<< dfunc_dDD(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dc0val<<"  "<< dfunc_dc0(A,back,m,B,DD,c0,tvalue,Nt)<<endl;


      //    f(t) = A * { exp(-m*t)  + exp(-m*(Nt-t))  } 

func_exponential_timesym(Avalue,mvalue,tvalue,Nt,result);
grad_exponential_timesym(Avalue,mvalue,tvalue,Nt,dAval,dmval);

cout << "Test 2"<<endl;
A=Avalue; back=1.0; m=mvalue; B=0.0;  DD=0.0; c0=0.0;
cout << func(A,back,m,B,DD,c0,tvalue,Nt)<<"  "<<result<<endl;
cout << dAval<<"  "<< dfunc_dA(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dmval<<"  "<< dfunc_dm(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dBval<<"  "<< dfunc_dB(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dDDval<<"  "<< dfunc_dDD(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dc0val<<"  "<< dfunc_dc0(A,back,m,B,DD,c0,tvalue,Nt)<<endl;


      //    f(t) = A * exp(-m*t)  + c0

func_exponential_timeforward_with_const(Avalue,mvalue,c0value,tvalue,result);
grad_exponential_timeforward_with_const(Avalue,mvalue,tvalue,
                                        dAval,dmval,dc0val);
cout << "Test 3"<<endl;
A=Avalue; back=0.0; m=mvalue; B=0.0;  DD=0.0; c0=c0value;
cout << func(A,back,m,B,DD,c0,tvalue,Nt)<<"  "<<result<<endl;
cout << dAval<<"  "<< dfunc_dA(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dmval<<"  "<< dfunc_dm(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dBval<<"  "<< dfunc_dB(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dDDval<<"  "<< dfunc_dDD(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dc0val<<"  "<< dfunc_dc0(A,back,m,B,DD,c0,tvalue,Nt)<<endl;


      //    f(t) = A * { exp(-m*t)  + exp(-m*(Nt-t))  }   + c0

func_exponential_timesym_with_const(Avalue,mvalue,c0value,tvalue,Nt,result);
grad_exponential_timesym_with_const(Avalue,mvalue,tvalue,Nt,
                                    dAval,dmval,dc0val);

cout << "Test 4"<<endl;
A=Avalue; back=1.0; m=mvalue; B=0.0;  DD=0.0; c0=c0value;
cout << func(A,back,m,B,DD,c0,tvalue,Nt)<<"  "<<result<<endl;
cout << dAval<<"  "<< dfunc_dA(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dmval<<"  "<< dfunc_dm(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dBval<<"  "<< dfunc_dB(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dDDval<<"  "<< dfunc_dDD(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dc0val<<"  "<< dfunc_dc0(A,back,m,B,DD,c0,tvalue,Nt)<<endl;


      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-DD^2*t) ]  }

func_two_exponential_timeforward(Avalue,mvalue,Bvalue,DDvalue,tvalue,result);
grad_two_exponential_timeforward(Avalue,mvalue,Bvalue,DDvalue,tvalue,
                                 dAval,dmval,dBval,dDDval);

cout << "Test 5"<<endl;
A=Avalue; back=0.0; m=mvalue; B=Bvalue;  DD=DDvalue; c0=0.0;
cout << func(A,back,m,B,DD,c0,tvalue,Nt)<<"  "<<result<<endl;
cout << dAval<<"  "<< dfunc_dA(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dmval<<"  "<< dfunc_dm(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dBval<<"  "<< dfunc_dB(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dDDval<<"  "<< dfunc_dDD(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dc0val<<"  "<< dfunc_dc0(A,back,m,B,DD,c0,tvalue,Nt)<<endl;


     //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-DD^2*t) ]
      //          + exp(-m*(Nt-t)) * [ 1 + B*exp(-DD^2*(Nt-t)) ] }

func_two_exponential_timesym(Avalue,mvalue,Bvalue,DDvalue,tvalue,Nt,result);
grad_two_exponential_timesym(Avalue,mvalue,Bvalue,DDvalue,tvalue,Nt,
                             dAval,dmval,dBval,dDDval);
cout << "Test 6"<<endl;
A=Avalue; back=1.0; m=mvalue; B=Bvalue;  DD=DDvalue; c0=0.0;
cout << func(A,back,m,B,DD,c0,tvalue,Nt)<<"  "<<result<<endl;
cout << dAval<<"  "<< dfunc_dA(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dmval<<"  "<< dfunc_dm(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dBval<<"  "<< dfunc_dB(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dDDval<<"  "<< dfunc_dDD(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
//cout << dc0val<<"  "<< dfunc_dc0(A,back,m,B,DD,c0,tvalue,Nt)<<endl;


      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-DD^2*t) ] }   + c0

func_two_exponential_timeforward_with_const(Avalue,mvalue,Bvalue,DDvalue,c0value,
                                            tvalue,result);
grad_two_exponential_timeforward_with_const(Avalue,mvalue,Bvalue,DDvalue,
                                            tvalue,dAval,dmval,dBval,dDDval,dc0val);

cout << "Test 7"<<endl;
A=Avalue; back=0.0; m=mvalue; B=Bvalue;  DD=DDvalue; c0=c0value;
cout << func(A,back,m,B,DD,c0,tvalue,Nt)<<"  "<<result<<endl;
cout << dAval<<"  "<< dfunc_dA(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dmval<<"  "<< dfunc_dm(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dBval<<"  "<< dfunc_dB(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dDDval<<"  "<< dfunc_dDD(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dc0val<<"  "<< dfunc_dc0(A,back,m,B,DD,c0,tvalue,Nt)<<endl;


      //    f(t) = A * { exp(-m*t) * [ 1 + B*exp(-DD^2*t) ]
      //          + exp(-m*(Nt-t)) * [ 1 + B*exp(-DD^2*(Nt-t)) ] }   + c0

func_two_exponential_timesym_with_const(Avalue,mvalue,Bvalue,DDvalue,c0value,
                                        tvalue,Nt,result);
grad_two_exponential_timesym_with_const(Avalue,mvalue,Bvalue,DDvalue,
                                        tvalue,Nt,dAval,dmval,dBval,dDDval,dc0val);

cout << "Test 8"<<endl;
A=Avalue; back=1.0; m=mvalue; B=Bvalue;  DD=DDvalue; c0=c0value;
cout << func(A,back,m,B,DD,c0,tvalue,Nt)<<"  "<<result<<endl;
cout << dAval<<"  "<< dfunc_dA(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dmval<<"  "<< dfunc_dm(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dBval<<"  "<< dfunc_dB(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dDDval<<"  "<< dfunc_dDD(A,back,m,B,DD,c0,tvalue,Nt)<<endl;
cout << dc0val<<"  "<< dfunc_dc0(A,back,m,B,DD,c0,tvalue,Nt)<<endl;

 // ***********************************************************************************


return 0;
}

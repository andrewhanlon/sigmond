#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include "xml_handler.h"
#include "model_tcorr.h"
#include "matrix.h"
#include "chisq_tcorr.h"
using namespace std;


namespace TcorrTester {

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

}


using namespace TcorrTester;



// ********************************************


void testChisqTcorr(XMLHandler& xml_in, int taskcount)
{
 if (xml_tag_count(xml_in,"TestChisqTcorr")==0)
 return;

 cout << endl << "Starting test_chisq_tcorr"<<endl;

 MCObsGetHandler *m_getter=0;
 MCObsHandler *m_obs=0;

 try{
    XMLHandler xmlr(xml_in);
    MCBinsInfo bininfo(xmlr);
    MCSamplingInfo sampinfo(xmlr);
    xml_tag_assert(xmlr,"MCObservables");
    m_getter=new MCObsGetHandler(xmlr,bininfo,sampinfo);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument("Bad data read"));}
 try{
    m_obs=new MCObsHandler(*m_getter);}
 catch(const std::exception& errmsg){
    delete m_getter;
    throw(std::invalid_argument("Bad MCObsHandler construction"));}

 // ***********************

 double Avalue=567.923;
 double mvalue=0.1325;
 double Bvalue=6.23;
 double DDvalue=0.879;
 double c0value=-1232.874;
 double A,m,B,DD,c0,back;
 int nobs,nparam;
 vector<funcptr> derivs;
 RMatrix grad;
 vector<double> fitparams;
 vector<double> modelpoints;
 int tmin=3;
 int tmax=125;
 int Tperiod=128;
 Vector<double> dummybins(m_obs->getNumberOfBins());
 OperatorInfo oppion("pion P=(0,0,0) A1um_1 SS_0");
 bool hermitian=true;
 CorrelatorAtTimeInfo corr(oppion,oppion,0,hermitian,false);

 // ***********************

 {cout << "Testing time forward single exp"<<endl;

 for (int tt=tmin;tt<=tmax;tt++){
    corr.resetTimeSeparation(tt);
    MCObsInfo obskey(corr);
    double res=func(Avalue,0.0,mvalue,0.0,0.0,0.0,tt,Tperiod);
    dummybins=res;
    m_obs->putBins(obskey,dummybins);}

 XMLHandler xmlt("TemporalCorrelatorFit");
 xmlt.put_child("OperatorString","pion P=(0,0,0) A1um_1 SS_0");
 xmlt.put_child("MinimumTimeSeparation",make_string(tmin));
 xmlt.put_child("MaximumTimeSeparation",make_string(tmax));
 XMLHandler xmlmod("Model");
 xmlmod.put_child("Type","TimeForwardSingleExponential");
 XMLHandler xmltmp("Energy");
 xmltmp.put_child("Name","pion");
 xmltmp.put_child("IDIndex","11");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("Amplitude");
 xmltmp.put_child("Name","A");
 xmltmp.put_child("IDIndex","5");
 xmlmod.put_child(xmltmp);
 xmlt.put_child(xmlmod);
// cout << xmlt.output()<<endl;

 RealTemporalCorrelatorFit Tester(xmlt,*m_obs,taskcount);

 nobs=Tester.getNumberOfObervables();
 nparam=Tester.getNumberOfParams();
// cout << "Number of parameters = "<<nparam<<endl;
// cout << "Number of observables = "<<nobs<<endl;

 int npoints=nobs;
 fitparams.resize(nparam);
 modelpoints.resize(nobs);

 fitparams[0]=m=mvalue;
 fitparams[1]=A=Avalue;
 back=0.0;
 B=0.0; DD=0.0;
 c0=0.0;
 derivs.resize(nparam);
 derivs[0]=dfunc_dm;
 derivs[1]=dfunc_dA;
 grad.resize(npoints,nparam);

 Tester.evalModelPoints(fitparams,modelpoints);
 for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
 Tester.evalGradients(fitparams,grad);
 for (int p=0;p<nparam;p++){
    for (int t=tmin, k=0; t<=tmax;++t,k++){
       double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
    double dg=std::fabs(grad(k,p)-g1);
    if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}

 RVector obs_mean(npoints);
 A=Avalue; back=0.0; m=mvalue; B=0.0; DD=0.0; c0=0.0;
 for (int t=tmin, k=0; t<=tmax;++t,k++)
    obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
 fitparams.resize(nparam);
 Tester.guessInitialParamValues(obs_mean,fitparams);
 cout << "guess[0] = "<<fitparams[0]<<" expect "<<mvalue<<endl;
 cout << "guess[1] = "<<fitparams[1]<<" expect "<<Avalue<<endl;
 }

 // *****************************************

 {cout << "Testing time sym single exp"<<endl;
 XMLHandler xmlt("TemporalCorrelatorFit");
 xmlt.put_child("OperatorString","pion P=(0,0,0) A1um_1 SS_0");
 xmlt.put_child("MinimumTimeSeparation",make_string(tmin));
 xmlt.put_child("MaximumTimeSeparation",make_string(tmax));
 XMLHandler xmlmod("Model");
 xmlmod.put_child("Type","TimeSymSingleExponential");
 XMLHandler xmltmp("Energy");
 xmltmp.put_child("Name","pion");
 xmltmp.put_child("IDIndex","12");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("Amplitude");
 xmltmp.put_child("Name","A");
 xmltmp.put_child("IDIndex","6");
 xmlmod.put_child(xmltmp);
 xmlt.put_child(xmlmod);
// cout << xmlt.output()<<endl;

 RealTemporalCorrelatorFit Tester(xmlt,*m_obs,taskcount);

 nobs=Tester.getNumberOfObervables();
 nparam=Tester.getNumberOfParams();
// cout << "Number of parameters = "<<nparam<<endl;
// cout << "Number of observables = "<<nobs<<endl;

 fitparams.resize(nparam);
 modelpoints.resize(nobs);

 back=1.0;
 Tester.evalModelPoints(fitparams,modelpoints);
 for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
 Tester.evalGradients(fitparams,grad);
 for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}

 RVector obs_mean(nobs);
 A=Avalue; back=0.0; m=mvalue; B=0.0; DD=0.0; c0=0.0;
 for (int t=tmin, k=0; t<=tmax;++t,k++)
    obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
 fitparams.resize(nparam);
 Tester.guessInitialParamValues(obs_mean,fitparams);
 cout << "guess[0] = "<<fitparams[0]<<" expect "<<mvalue<<endl;
 cout << "guess[1] = "<<fitparams[1]<<" expect "<<Avalue<<endl;
 }

 // *****************************************

 {cout << "Testing time forward single exp with constant"<<endl;
 XMLHandler xmlt("TemporalCorrelatorFit");
 xmlt.put_child("OperatorString","pion P=(0,0,0) A1um_1 SS_0");
 xmlt.put_child("MinimumTimeSeparation",make_string(tmin));
 xmlt.put_child("MaximumTimeSeparation",make_string(tmax));
 XMLHandler xmlmod("Model");
 xmlmod.put_child("Type","TimeForwardSingleExponentialPlusConstant");
 XMLHandler xmltmp("Energy");
 xmltmp.put_child("Name","pion");
 xmltmp.put_child("IDIndex","12");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("Amplitude");
 xmltmp.put_child("Name","A");
 xmltmp.put_child("IDIndex","6");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("AddedConstant");
 xmltmp.put_child("Name","C0");
 xmltmp.put_child("IDIndex","8");
 xmlmod.put_child(xmltmp);
 xmlt.put_child(xmlmod);
// cout << xmlt.output()<<endl;

 RealTemporalCorrelatorFit Tester(xmlt,*m_obs,taskcount);

 nobs=Tester.getNumberOfObervables();
 nparam=Tester.getNumberOfParams();
// cout << "Number of parameters = "<<nparam<<endl;
// cout << "Number of observables = "<<nobs<<endl;

 int npoints=nobs;

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

 Tester.evalModelPoints(fitparams,modelpoints);
 for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
 Tester.evalGradients(fitparams,grad);
 for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}

 RVector obs_mean(npoints);
 A=Avalue; back=0.0; m=mvalue; B=0.0; DD=0.0; c0=c0value;
 for (int t=tmin, k=0; t<=tmax;++t,k++)
    obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
 fitparams.resize(nparam);
 Tester.guessInitialParamValues(obs_mean,fitparams);
 cout << "guess[0] = "<<fitparams[0]<<" expect "<<mvalue<<endl;
 cout << "guess[1] = "<<fitparams[1]<<" expect "<<Avalue<<endl;
 cout << "guess[2] = "<<fitparams[2]<<" expect "<<c0value<<endl;
 }

 // *****************************************

 {cout << "Testing time symmetric single exp with constant"<<endl;
 XMLHandler xmlt("TemporalCorrelatorFit");
 xmlt.put_child("OperatorString","pion P=(0,0,0) A1um_1 SS_0");
 xmlt.put_child("MinimumTimeSeparation",make_string(tmin));
 xmlt.put_child("MaximumTimeSeparation",make_string(tmax));
 XMLHandler xmlmod("Model");
 xmlmod.put_child("Type","TimeSymSingleExponentialPlusConstant");
 XMLHandler xmltmp("Energy");
 xmltmp.put_child("Name","pion");
 xmltmp.put_child("IDIndex","12");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("Amplitude");
 xmltmp.put_child("Name","A");
 xmltmp.put_child("IDIndex","6");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("AddedConstant");
 xmltmp.put_child("Name","C0");
 xmltmp.put_child("IDIndex","8");
 xmlmod.put_child(xmltmp);
 xmlt.put_child(xmlmod);
// cout << xmlt.output()<<endl;

 RealTemporalCorrelatorFit Tester(xmlt,*m_obs,taskcount);

 nobs=Tester.getNumberOfObervables();
 nparam=Tester.getNumberOfParams();
// cout << "Number of parameters = "<<nparam<<endl;
// cout << "Number of observables = "<<nobs<<endl;

 int npoints=nobs;
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

 Tester.evalModelPoints(fitparams,modelpoints);
 for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
 Tester.evalGradients(fitparams,grad);
 for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}

 RVector obs_mean(npoints);
 A=Avalue; back=0.0; m=mvalue; B=0.0; DD=0.0; c0=c0value;
 for (int t=tmin, k=0; t<=tmax;++t,k++)
    obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
 fitparams.resize(nparam);
 Tester.guessInitialParamValues(obs_mean,fitparams);
 cout << "guess[0] = "<<fitparams[0]<<" expect "<<mvalue<<endl;
 cout << "guess[1] = "<<fitparams[1]<<" expect "<<Avalue<<endl;
 cout << "guess[2] = "<<fitparams[2]<<" expect "<<c0value<<endl;
 }

 // ***********************

 {cout << "Testing time forward two exp"<<endl;
 XMLHandler xmlt("TemporalCorrelatorFit");
 xmlt.put_child("OperatorString","pion P=(0,0,0) A1um_1 SS_0");
 xmlt.put_child("MinimumTimeSeparation",make_string(tmin));
 xmlt.put_child("MaximumTimeSeparation",make_string(tmax));
 XMLHandler xmlmod("Model");
 xmlmod.put_child("Type","TimeForwardTwoExponential");
 XMLHandler xmltmp("FirstEnergy");
 xmltmp.put_child("Name","pion");
 xmltmp.put_child("IDIndex","15");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("FirstAmplitude");
 xmltmp.put_child("Name","AA");
 xmltmp.put_child("IDIndex","7");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("SqrtGapToSecondEnergy");
 xmltmp.put_child("Name","DD");
 xmltmp.put_child("IDIndex","2");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("SecondAmplitudeRatio");
 xmltmp.put_child("Name","B");
 xmltmp.put_child("IDIndex","4");
 xmlmod.put_child(xmltmp);
 xmlt.put_child(xmlmod);
// cout << xmlt.output()<<endl;


 RealTemporalCorrelatorFit Tester(xmlt,*m_obs,taskcount);

 nobs=Tester.getNumberOfObervables();
 nparam=Tester.getNumberOfParams();
// cout << "Number of parameters = "<<nparam<<endl;
// cout << "Number of observables = "<<nobs<<endl;

 int npoints=nobs;

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

 Tester.evalModelPoints(fitparams,modelpoints);
 for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
 Tester.evalGradients(fitparams,grad);
 for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}

 RVector obs_mean(nobs);
 A=Avalue; back=0.0; m=mvalue; B=Bvalue; DD=DDvalue; c0=0.0;
 for (int t=tmin, k=0; t<=tmax;++t,k++)
    obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
 fitparams.resize(nparam);
 Tester.guessInitialParamValues(obs_mean,fitparams);
 cout << "guess[0] = "<<fitparams[0]<<" expect "<<mvalue<<endl;
 cout << "guess[1] = "<<fitparams[1]<<" expect "<<Avalue<<endl;
 cout << "guess[2] = "<<fitparams[2]<<" expect "<<DDvalue<<endl;
 cout << "guess[3] = "<<fitparams[3]<<" expect "<<Bvalue<<endl;
 }


 // *****************************************

 {cout << "Testing time sym two exp"<<endl;
 XMLHandler xmlt("TemporalCorrelatorFit");
 xmlt.put_child("OperatorString","pion P=(0,0,0) A1um_1 SS_0");
 xmlt.put_child("MinimumTimeSeparation",make_string(tmin));
 xmlt.put_child("MaximumTimeSeparation",make_string(tmax));
 XMLHandler xmlmod("Model");
 xmlmod.put_child("Type","TimeSymTwoExponential");
 XMLHandler xmltmp("FirstEnergy");
 xmltmp.put_child("Name","pion");
 xmltmp.put_child("IDIndex","16");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("FirstAmplitude");
 xmltmp.put_child("Name","AA");
 xmltmp.put_child("IDIndex","77");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("SqrtGapToSecondEnergy");
 xmltmp.put_child("Name","DD");
 xmltmp.put_child("IDIndex","3");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("SecondAmplitudeRatio");
 xmltmp.put_child("Name","B");
 xmltmp.put_child("IDIndex","5");
 xmlmod.put_child(xmltmp);
 xmlt.put_child(xmlmod);
// cout << xmlt.output()<<endl;


 RealTemporalCorrelatorFit Tester(xmlt,*m_obs,taskcount);

 nobs=Tester.getNumberOfObervables();
 nparam=Tester.getNumberOfParams();
// cout << "Number of parameters = "<<nparam<<endl;
// cout << "Number of observables = "<<nobs<<endl;

 int npoints=nobs;

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

 Tester.evalModelPoints(fitparams,modelpoints);
 for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
 Tester.evalGradients(fitparams,grad);
 for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}

 RVector obs_mean(npoints);
 A=Avalue; back=0.0; m=mvalue; B=Bvalue; DD=DDvalue; c0=0.0;
 for (int t=tmin, k=0; t<=tmax;++t,k++)
    obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
 fitparams.resize(nparam);
 Tester.guessInitialParamValues(obs_mean,fitparams);
 cout << "guess[0] = "<<fitparams[0]<<" expect "<<mvalue<<endl;
 cout << "guess[1] = "<<fitparams[1]<<" expect "<<Avalue<<endl;
 cout << "guess[2] = "<<fitparams[2]<<" expect "<<DDvalue<<endl;
 cout << "guess[3] = "<<fitparams[3]<<" expect "<<Bvalue<<endl;
 }

 // ***********************

 {cout << "Testing time sym two exp plus const"<<endl;
 XMLHandler xmlt("TemporalCorrelatorFit");
 xmlt.put_child("OperatorString","pion P=(0,0,0) A1um_1 SS_0");
 xmlt.put_child("MinimumTimeSeparation",make_string(tmin));
 xmlt.put_child("MaximumTimeSeparation",make_string(tmax));
 XMLHandler xmlmod("Model");
 xmlmod.put_child("Type","TimeForwardTwoExponentialPlusConstant");
 XMLHandler xmltmp("FirstEnergy");
 xmltmp.put_child("Name","pion");
 xmltmp.put_child("IDIndex","16");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("FirstAmplitude");
 xmltmp.put_child("Name","AA");
 xmltmp.put_child("IDIndex","77");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("SqrtGapToSecondEnergy");
 xmltmp.put_child("Name","DD");
 xmltmp.put_child("IDIndex","3");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("SecondAmplitudeRatio");
 xmltmp.put_child("Name","B");
 xmltmp.put_child("IDIndex","5");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("AddedConstant");
 xmltmp.put_child("Name","C0");
 xmltmp.put_child("IDIndex","9");
 xmlmod.put_child(xmltmp);
 xmlt.put_child(xmlmod);
// cout << xmlt.output()<<endl;


 RealTemporalCorrelatorFit Tester(xmlt,*m_obs,taskcount);

 nobs=Tester.getNumberOfObervables();
 nparam=Tester.getNumberOfParams();
// cout << "Number of parameters = "<<nparam<<endl;
// cout << "Number of observables = "<<nobs<<endl;

 int npoints=nobs;

 fitparams.resize(nparam);
 fitparams[0]=m=mvalue;
 fitparams[1]=A=Avalue;
 fitparams[2]=DD=DDvalue;
 fitparams[3]=B=Bvalue;
 fitparams[4]=c0value;
 back=0.0;
 c0=c0value;
 derivs.resize(nparam);
 derivs[0]=dfunc_dm;
 derivs[1]=dfunc_dA;
 derivs[2]=dfunc_dDD;
 derivs[3]=dfunc_dB;
 derivs[4]=dfunc_dc0;
 grad.resize(npoints,nparam);

 Tester.evalModelPoints(fitparams,modelpoints);
 for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
 Tester.evalGradients(fitparams,grad);
 for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}

 RVector obs_mean(npoints);
 A=Avalue; back=0.0; m=mvalue; B=Bvalue; DD=DDvalue; c0=c0value;
 for (int t=tmin, k=0; t<=tmax;++t,k++)
    obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
 fitparams.resize(nparam);
 Tester.guessInitialParamValues(obs_mean,fitparams);
 cout << "guess[0] = "<<fitparams[0]<<" expect "<<mvalue<<endl;
 cout << "guess[1] = "<<fitparams[1]<<" expect "<<Avalue<<endl;
 cout << "guess[2] = "<<fitparams[2]<<" expect "<<DDvalue<<endl;
 cout << "guess[3] = "<<fitparams[3]<<" expect "<<Bvalue<<endl;
 cout << "guess[4] = "<<fitparams[4]<<" expect "<<c0value<<endl;
 }


 {cout << "Testing time sym two exp plus const"<<endl;
 XMLHandler xmlt("TemporalCorrelatorFit");
 xmlt.put_child("OperatorString","pion P=(0,0,0) A1um_1 SS_0");
 xmlt.put_child("MinimumTimeSeparation",make_string(tmin));
 xmlt.put_child("MaximumTimeSeparation",make_string(tmax));
 XMLHandler xmlmod("Model");
 xmlmod.put_child("Type","TimeSymTwoExponentialPlusConstant");
 XMLHandler xmltmp("FirstEnergy");
 xmltmp.put_child("Name","pion");
 xmltmp.put_child("IDIndex","17");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("FirstAmplitude");
 xmltmp.put_child("Name","AA");
 xmltmp.put_child("IDIndex","9");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("SqrtGapToSecondEnergy");
 xmltmp.put_child("Name","DD");
 xmltmp.put_child("IDIndex","4");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("SecondAmplitudeRatio");
 xmltmp.put_child("Name","B");
 xmltmp.put_child("IDIndex","7");
 xmlmod.put_child(xmltmp);
 xmltmp.set_root("AddedConstant");
 xmltmp.put_child("Name","C0");
 xmltmp.put_child("IDIndex","3");
 xmlmod.put_child(xmltmp);
 xmlt.put_child(xmlmod);
// cout << xmlt.output()<<endl;


 RealTemporalCorrelatorFit Tester(xmlt,*m_obs,taskcount);

 nobs=Tester.getNumberOfObervables();
 nparam=Tester.getNumberOfParams();
// cout << "Number of parameters = "<<nparam<<endl;
// cout << "Number of observables = "<<nobs<<endl;

 int npoints=nobs;

 fitparams.resize(nparam);
 fitparams[0]=m=mvalue;
 fitparams[1]=A=Avalue;
 fitparams[2]=DD=DDvalue;
 fitparams[3]=B=Bvalue;
 fitparams[4]=c0value;
 back=1.0;
 c0=c0value;
 derivs.resize(nparam);
 derivs[0]=dfunc_dm;
 derivs[1]=dfunc_dA;
 derivs[2]=dfunc_dDD;
 derivs[3]=dfunc_dB;
 derivs[4]=dfunc_dc0;
 grad.resize(npoints,nparam);

 Tester.evalModelPoints(fitparams,modelpoints);
 for (int t=tmin, k=0; t<=tmax;++t,k++){
    double v1=func(A,back,m,B,DD,c0,t,Tperiod);
    double df=std::fabs(modelpoints[k]-v1);
    if (df>1e-10) cout << modelpoints[k]<<" "<<v1<<" " << df<<endl;}
 Tester.evalGradients(fitparams,grad);
 for (int p=0;p<nparam;p++){
   for (int t=tmin, k=0; t<=tmax;++t,k++){
      double g1=(*derivs[p])(A,back,m,B,DD,c0,t,Tperiod);
   double dg=std::fabs(grad(k,p)-g1);
   if (dg>1e-10) cout << "grad "<<p<<" diff "<< dg<<endl;}}

 RVector obs_mean(npoints);
 A=Avalue; back=0.0; m=mvalue; B=Bvalue; DD=DDvalue; c0=c0value;
 for (int t=tmin, k=0; t<=tmax;++t,k++)
    obs_mean[k]=func(A,back,m,B,DD,c0,t,Tperiod);
 fitparams.resize(nparam);
 Tester.guessInitialParamValues(obs_mean,fitparams);
 cout << "guess[0] = "<<fitparams[0]<<" expect "<<mvalue<<endl;
 cout << "guess[1] = "<<fitparams[1]<<" expect "<<Avalue<<endl;
 cout << "guess[2] = "<<fitparams[2]<<" expect "<<DDvalue<<endl;
 cout << "guess[3] = "<<fitparams[3]<<" expect "<<Bvalue<<endl;
 cout << "guess[4] = "<<fitparams[4]<<" expect "<<c0value<<endl;
 }



 delete m_obs;
 delete m_getter;

  // ************************************************************************

/*
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

int tval=4;
Avalue=3.2546; mvalue=0.5323;
double corrt=Avalue*exp(-mvalue*double(tval));
double corrtnext=Avalue*exp(-mvalue*double(tval+1));

guess_exponential(tval,corrt,corrtnext,A,m);
cout <<endl<<endl<< "Testing guess exponential"<<endl;
cout << "A = "<<A<<endl;
cout << "m = "<<m<<endl;

Avalue=7.432; mvalue=0.2912;
c0value=-21.32;
corrt=Avalue*exp(-mvalue*double(tval))+c0value;
double corrtp1=Avalue*exp(-mvalue*double(tval+1))+c0value;
double corrtp2=Avalue*exp(-mvalue*double(tval+2))+c0value;

guess_exponential_plus_const(tval,corrt,corrtp1,corrtp2,A,m,c0);
cout <<endl<<endl<< "Testing guess exponential plus cont"<<endl;
cout << "A = "<<A<<endl;
cout << "m = "<<m<<endl;
cout << "c0 = "<<c0<<endl;

Avalue=24.124; mvalue=0.2912;
Bvalue=0.9232; DDvalue=0.321;
double gap=DDvalue*DDvalue;
double f0=Avalue*exp(-mvalue*double(tval))*(1.0+Bvalue*exp(-gap*double(tval)));
double f1=Avalue*exp(-mvalue*double(tval+1))*(1.0+Bvalue*exp(-gap*double(tval+1)));
double f2=Avalue*exp(-mvalue*double(tval+2))*(1.0+Bvalue*exp(-gap*double(tval+2)));
double f3=Avalue*exp(-mvalue*double(tval+3))*(1.0+Bvalue*exp(-gap*double(tval+3)));

guess_two_exponential(tval,f0,f1,f2,f3,A,m,B,DD);
cout <<endl<<endl<< "Testing guess two exponential"<<endl;
cout << "A = "<<A<<endl;
cout << "m = "<<m<<endl;
cout << "B = "<<B<<endl;
cout << "DD = "<<DD<<endl;

Avalue=36.43; mvalue=0.111;
Bvalue=2.9232; DDvalue=0.222; c0value=55.5;
gap=DDvalue*DDvalue;
f0=Avalue*exp(-mvalue*double(tval))*(1.0+Bvalue*exp(-gap*double(tval)))+c0value;
f1=Avalue*exp(-mvalue*double(tval+1))*(1.0+Bvalue*exp(-gap*double(tval+1)))+c0value;
f2=Avalue*exp(-mvalue*double(tval+2))*(1.0+Bvalue*exp(-gap*double(tval+2)))+c0value;
f3=Avalue*exp(-mvalue*double(tval+3))*(1.0+Bvalue*exp(-gap*double(tval+3)))+c0value;
double f4=Avalue*exp(-mvalue*double(tval+4))*(1.0+Bvalue*exp(-gap*double(tval+4)))+c0value;

guess_two_exponential_with_const(tval,f0,f1,f2,f3,f4,A,m,B,DD,c0);
cout <<endl<<endl<< "Testing guess two exponential plus const"<<endl;
cout << "A = "<<A<<endl;
cout << "m = "<<m<<endl;
cout << "B = "<<B<<endl;
cout << "DD = "<<DD<<endl;
cout << "c0 = "<<c0<<endl;
*/
 // ***********************************************************************************
}

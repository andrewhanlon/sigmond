#include <iostream>
#include "xml_handler.h"
#include "obs_get_handler.h"
#include "mcobs_handler.h"
#include "chisq_tcorr.h"
#include "minimizer.h"
#include "task_utils.h"
#include "chisq_fit.h"
using namespace std;



void testChiSquare(XMLHandler& xml_in, int taskcount)
{
 if (xml_tag_count(xml_in,"TestChiSquare")==0)
 return;

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

 RealTemporalCorrelatorFit RTC(xml_in,*m_obs,taskcount);

 m_obs->begin();
 int nobs=RTC.getNumberOfObervables();
 int nparam=RTC.getNumberOfParams();
 cout << "Number of parameters = "<<nparam<<endl;
 cout << "Number of observables = "<<nobs<<endl;
 RTC.setObsMeanCov();

 vector<double> fitparams(nparam);
 vector<double> residuals(nobs);
 RTC.guessInitialFitParamValues(fitparams);
 RTC.evalResiduals(fitparams,residuals);
 cout << "Chi square = "<<RTC.evalChiSquare(residuals)<<endl;
    
 fitparams[0]=0.06918;
 fitparams[1]=149.525;
 RTC.evalResiduals(fitparams,residuals);
 cout << "Chi square = "<<RTC.evalChiSquare(residuals)<<endl;

// fitparams[0]=0.079;
// RTC.evalResiduals(fitparams,residuals);
// cout << "Chi square = "<<RTC.evalChiSquare(residuals)<<endl;
 
 RMatrix grad(nobs,nparam);
 RTC.evalResGradients(fitparams,grad);
 cout << "grad(1,0) = "<<grad(1,0)<<endl;
 cout << "grad2(1,0) = "<<grad(1,0)<<endl;

 cout << "fit param[0]="<<fitparams[0]<<endl;
 cout << "fit param[1]="<<fitparams[1]<<endl;

 const vector<MCObsInfo>& obs=RTC.getObsInfos();
 for (uint k=0;k<obs.size();++k) cout << obs[k].output()<<endl;

 const vector<MCObsInfo>& params=RTC.getFitParamInfos();
 for (uint k=0;k<params.size();++k) cout << params[k].output()<<endl;

 cout << RTC.output()<<endl;
 cout << RTC.output(3)<<endl;
 cout << RTC.str()<<endl;



   // ****************************************

 ChiSquareMinimizerInfo CSM;
 cout << "ChiSquareMinimizer Test 1:"<<endl;
 cout << "usingLMDer? "<<CSM.usingLMDer()<<endl;
 cout << "usingNL2Sol? "<<CSM.usingNL2Sol()<<endl;
 cout << "usingMinuit2? "<<CSM.usingMinuit2()<<endl;
 cout << "max iterations = "<<CSM.getMaximumIterations()<<endl;
 cout << "param rel tol = "<<CSM.getParameterRelativeTolerance()<<endl;
 cout << "chisq rel tol = "<<CSM.getChiSquareRelativeTolerance()<<endl;
 cout << "isLowVerbosity? "<<CSM.isLowVerbosity()<<endl;
 cout << "isMedVerbosity? "<<CSM.isMediumVerbosity()<<endl;
 cout << "isHiVerbosity? "<<CSM.isHighVerbosity()<<endl;
 cout << "isNotLoVerbosity? "<<CSM.isNotLowVerbosity()<<endl;
// cout << "getChiSquare output: "<<CSM.getChiSquare().output()<<endl;


 try{
    cout << "setting invalid method: should throw exception"<<endl;
    CSM.setMethod('V');
    cout << "ERROR"<<endl;}
 catch(const std::exception& xp){
    cout << "Caught exception: correct"<<endl;}

 CSM.setMethod('L');
 cout << "ChiSquareMinimizer Test 2:"<<endl;
 cout << "usingLMDer? "<<CSM.usingLMDer()<<endl;
 cout << "usingNL2Sol? "<<CSM.usingNL2Sol()<<endl;
 cout << "usingMinuit2? "<<CSM.usingMinuit2()<<endl;
 CSM.setMethod('N');
 cout << "ChiSquareMinimizer Test 3:"<<endl;
 cout << "usingLMDer? "<<CSM.usingLMDer()<<endl;
 cout << "usingNL2Sol? "<<CSM.usingNL2Sol()<<endl;
 cout << "usingMinuit2? "<<CSM.usingMinuit2()<<endl;
 try{
    CSM.setMethod('M');
    cout << "ChiSquareMinimizer Test 4:"<<endl;
    cout << "usingLMDer? "<<CSM.usingLMDer()<<endl;
    cout << "usingNL2Sol? "<<CSM.usingNL2Sol()<<endl;
    cout << "usingMinuit2? "<<CSM.usingMinuit2()<<endl;}
 catch(const std::exception& xp){
    cout << "Exception caught: no Minuit2"<<endl;}
 CSM.setLMDer();
 cout << "ChiSquareMinimizer Test 5:"<<endl;
 cout << "usingLMDer? "<<CSM.usingLMDer()<<endl;
 cout << "usingNL2Sol? "<<CSM.usingNL2Sol()<<endl;
 cout << "usingMinuit2? "<<CSM.usingMinuit2()<<endl;
 CSM.setNL2Sol();
 cout << "ChiSquareMinimizer Test 6:"<<endl;
 cout << "usingLMDer? "<<CSM.usingLMDer()<<endl;
 cout << "usingNL2Sol? "<<CSM.usingNL2Sol()<<endl;
 cout << "usingMinuit2? "<<CSM.usingMinuit2()<<endl;
 try{
    CSM.setMinuit2();
    cout << "ChiSquareMinimizer Test 7:"<<endl;
    cout << "usingLMDer? "<<CSM.usingLMDer()<<endl;
    cout << "usingNL2Sol? "<<CSM.usingNL2Sol()<<endl;
    cout << "usingMinuit2? "<<CSM.usingMinuit2()<<endl;}
 catch(const std::exception& xp){
    cout << "Exception caught: no Minuit2"<<endl;}
 CSM.setChiSquareRelativeTolerance(1.23e-4);
 CSM.setParameterRelativeTolerance(6.6e-5);
 CSM.setMaximumIterations(134);
 cout << "ChiSquareMinimizer Test 8:"<<endl;
 cout << "max iterations = "<<CSM.getMaximumIterations()<<endl;
 cout << "param rel tol = "<<CSM.getParameterRelativeTolerance()<<endl;
 cout << "chisq rel tol = "<<CSM.getChiSquareRelativeTolerance()<<endl;

 try{
    cout << "setting invalid verbosity: should throw exception"<<endl;
    CSM.setVerbosity('V');
    cout << "ERROR"<<endl;}
 catch(const std::exception& xp){
    cout << "Caught exception: correct"<<endl;}
 cout << "usingLMDer? "<<CSM.usingLMDer()<<endl;
 cout << "usingNL2Sol? "<<CSM.usingNL2Sol()<<endl;
 cout << "usingMinuit2? "<<CSM.usingMinuit2()<<endl;
 CSM.setVerbosity('L');
 cout << "ChiSquareMinimizer Test 9:"<<endl;
 cout << "isLowVerbosity? "<<CSM.isLowVerbosity()<<endl;
 cout << "isMedVerbosity? "<<CSM.isMediumVerbosity()<<endl;
 cout << "isHiVerbosity? "<<CSM.isHighVerbosity()<<endl;
 cout << "isNotLoVerbosity? "<<CSM.isNotLowVerbosity()<<endl;
 CSM.setVerbosity('M');
 cout << "ChiSquareMinimizer Test 10:"<<endl;
 cout << "isLowVerbosity? "<<CSM.isLowVerbosity()<<endl;
 cout << "isMedVerbosity? "<<CSM.isMediumVerbosity()<<endl;
 cout << "isHiVerbosity? "<<CSM.isHighVerbosity()<<endl;
 cout << "isNotLoVerbosity? "<<CSM.isNotLowVerbosity()<<endl;
 CSM.setVerbosity('H');
 cout << "ChiSquareMinimizer Test 11:"<<endl;
 cout << "isLowVerbosity? "<<CSM.isLowVerbosity()<<endl;
 cout << "isMedVerbosity? "<<CSM.isMediumVerbosity()<<endl;
 cout << "isHiVerbosity? "<<CSM.isHighVerbosity()<<endl;
 cout << "isNotLoVerbosity? "<<CSM.isNotLowVerbosity()<<endl;
 CSM.setLowVerbosity();
 cout << "ChiSquareMinimizer Test 12:"<<endl;
 cout << "isLowVerbosity? "<<CSM.isLowVerbosity()<<endl;
 cout << "isMedVerbosity? "<<CSM.isMediumVerbosity()<<endl;
 cout << "isHiVerbosity? "<<CSM.isHighVerbosity()<<endl;
 cout << "isNotLoVerbosity? "<<CSM.isNotLowVerbosity()<<endl;
 CSM.setMediumVerbosity();
 cout << "ChiSquareMinimizer Test 13:"<<endl;
 cout << "isLowVerbosity? "<<CSM.isLowVerbosity()<<endl;
 cout << "isMedVerbosity? "<<CSM.isMediumVerbosity()<<endl;
 cout << "isHiVerbosity? "<<CSM.isHighVerbosity()<<endl;
 cout << "isNotLoVerbosity? "<<CSM.isNotLowVerbosity()<<endl;
 CSM.setHighVerbosity();
 cout << "ChiSquareMinimizer Test 14:"<<endl;
 cout << "isLowVerbosity? "<<CSM.isLowVerbosity()<<endl;
 cout << "isMedVerbosity? "<<CSM.isMediumVerbosity()<<endl;
 cout << "isHiVerbosity? "<<CSM.isHighVerbosity()<<endl;
 cout << "isNotLoVerbosity? "<<CSM.isNotLowVerbosity()<<endl;

 cout << "usingLMDer? "<<CSM.usingLMDer()<<endl;
 cout << "usingNL2Sol? "<<CSM.usingNL2Sol()<<endl;
 cout << "usingMinuit2? "<<CSM.usingMinuit2()<<endl;


 cout << CSM.output()<<endl;
 cout << CSM.output(3)<<endl;
 cout << CSM.str()<<endl;

 ChiSquareMinimizerInfo CSM3;
 cout << "CSM3: "<<endl<<CSM3.output()<<endl;

 ChiSquareMinimizerInfo CSM4('N',33,22,88,'M');
 cout << "CSM4: "<<endl<<CSM4.output()<<endl;

 ChiSquareMinimizerInfo CSM5(CSM4);
 cout << "CSM5: "<<endl<<CSM5.output()<<endl;

 ChiSquareMinimizerInfo CSM6;
 CSM6=CSM4;
 cout << "CSM6: "<<endl<<CSM6.output()<<endl;

 try{
    ChiSquareMinimizerInfo CSM2(xml_in);
    cout << "CSM2: "<<CSM2.output()<<endl;}
 catch(const std::exception& xp){
    cout << "Caught exception creating CSM2"<<endl;}


 XMLHandler xmlnl;
 cout <<endl<<endl<<"Starting NL2Sol minimizations"<<endl<<endl;
 ChiSquareMinimizerInfo minfo('N',1e-6,1e-4,2048,'H');
 ChiSquareMinimizer HH(RTC,minfo);
 double chisq;
 vector<double> params_at_minimum;
 bool solflag=HH.findMinimum(chisq,params_at_minimum,xmlnl);
 cout << "solflag = "<<solflag<<endl;
 cout << "chisq = "<<chisq<<endl;
 for (uint k=0;k<params_at_minimum.size();++k)
    cout << "params_at_min["<<k<<"] = "<<params_at_minimum[k]<<endl;
 cout << xmlnl.output()<<endl;

 cout <<endl<<"Starting another NL2Sol minimization"<<endl;
 vector<double> start(params_at_minimum); 
 solflag=HH.findMinimum(start,chisq,params_at_minimum,xmlnl);
 cout << "solflag = "<<solflag<<endl;
 cout << "chisq = "<<chisq<<endl;
 for (uint k=0;k<params_at_minimum.size();++k)
    cout << "params_at_min["<<k<<"] = "<<params_at_minimum[k]<<endl;
 cout << xmlnl.output()<<endl;

 cout <<endl<<"Starting another NL2Sol minimization"<<endl;
 start[0]=1.0; start[1]=0.43;
 solflag=HH.findMinimum(start,chisq,params_at_minimum,xmlnl);
 cout << "solflag = "<<solflag<<endl;
 cout << "chisq = "<<chisq<<endl;
 for (uint k=0;k<params_at_minimum.size();++k)
    cout << "params_at_min["<<k<<"] = "<<params_at_minimum[k]<<endl;
 cout << xmlnl.output()<<endl;

 ChiSquareMinimizerInfo minfo2('L',1e-6,1e-4,2048,'H');
 HH.reset(minfo2);
 XMLHandler xmllo;
 cout <<endl<<endl<<"Starting LMDer minimizations"<<endl<<endl;

 bool ifail=HH.findMinimum(chisq,params_at_minimum,xmllo);
 cout << "ifail = "<<ifail<<endl;
 cout << "chisq = "<<chisq<<endl;
 for (uint k=0;k<params_at_minimum.size();++k)
    cout << "params_at_min["<<k<<"] = "<<params_at_minimum[k]<<endl;
 cout << xmllo.output()<<endl;

 start=params_at_minimum; 
 ifail=HH.findMinimum(start,chisq,params_at_minimum,xmllo);
 cout << endl<<"Another LMDer minimization:"<<endl;
 cout << "ifail = "<<ifail<<endl;
 cout << "chisq = "<<chisq<<endl;
 for (uint k=0;k<params_at_minimum.size();++k)
    cout << "params_at_min["<<k<<"] = "<<params_at_minimum[k]<<endl;
 cout << xmllo.output()<<endl;

 start[0]=1.0; start[1]=0.43;
 cout << endl<<"Another LMDer minimization:"<<endl;
 ifail=HH.findMinimum(start,chisq,params_at_minimum,xmllo);
 cout << "ifail = "<<ifail<<endl;
 cout << "chisq = "<<chisq<<endl;
 for (uint k=0;k<params_at_minimum.size();++k)
    cout << "params_at_min["<<k<<"] = "<<params_at_minimum[k]<<endl;
 cout << xmllo.output()<<endl;

 try{
 ChiSquareMinimizerInfo minfo3('M',1e-6,1e-4,2048,'M');
 HH.reset(minfo3);
 XMLHandler xmlout;
 bool flag=HH.findMinimum(chisq,params_at_minimum,xmlout);

 cout <<endl<<endl<<"Did a Minuit2 minimization: flag = "<<flag<<endl<<endl;
 cout << xmlout.output()<<endl;
 cout << "chisq = "<<chisq<<endl;
 for (uint k=0;k<params_at_minimum.size();++k)
    cout << "params_at_min["<<k<<"] = "<<params_at_minimum[k]<<endl;
 start[0]=1.0; start[1]=0.43;
 flag=HH.findMinimum(start,chisq,params_at_minimum,xmlout);
 
 cout <<endl<<endl<<"Did a Minuit2 minimization: flag = "<<flag<<endl<<endl;
 cout << xmlout.output()<<endl;
 cout << "chisq = "<<chisq<<endl;
 for (uint k=0;k<params_at_minimum.size();++k)
    cout << "params_at_min["<<k<<"] = "<<params_at_minimum[k]<<endl;}
 
 catch(const std::exception& xp){
   cout << "Exception caught: no Minuit2"<<endl;}

/*
 LMDerMinimizer lm(RTC);
 lm.guessInitialFitParamValues();
 int verbosity=2;
 double tol=1e-6;
 int max_its=1200;
 double chisq;
 int flag=lm.chisq_fit(tol,verbosity,max_its,chisq);
 cout << "flag = "<<flag<<endl;


 NL2SolMinimizer nl(RTC);
 nl.guessInitialFitParamValues();
 verbosity=5;
 tol=1e-6;
 max_its=1200;
 flag=nl.chisq_fit(tol,verbosity,max_its,chisq);

*/


 delete m_obs;
 delete m_getter;
}



void testGamma(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestGamma")==0)
 return;

 cout.precision(16);
 double s=100.0;
 cout << "s="<<s<<endl;
 for (double x=100.0; x<=150.0; x+=0.1){
    cout << " "<<x<<"   "<<Qgamma(s,x)<<endl;}

}

#include "chisq_fit.h"
using namespace std;


// *************************************************************************


FitResult doChiSquareFitting(ChiSquare& chisq_ref, 
                             const ChiSquareMinimizerInfo& csm_info,
                             vector<double> initial_guesses,
                             bool correlated, XMLHandler& xmlout)
{
 int nparams=chisq_ref.getNumberOfParams();
 int nobs=chisq_ref.getNumberOfObervables();
 int npriors=chisq_ref.getNumberOfPriors();
 const vector<MCObsInfo>& param_infos=chisq_ref.getFitParamInfos();
 double dof=double(nobs-nparams+npriors);
 MCObsHandler *m_obs=chisq_ref.getMCObsHandlerPtr();

 for (int p=0;p<nparams;++p)
    if (m_obs->queryFullAndSamplings(param_infos[p]))
        throw(std::invalid_argument(string("Error: samplings already available for parameter ")
             +param_infos[p].str()));

 SamplingMode mode=chisq_ref.getObsMeansSamplingMode();
 SamplingMode covmode=chisq_ref.getCovMatSamplingMode();
 if (covmode==Jackknife) xmlout.put_child("CovarianceCalculationMode","Jackknife");
 else if (covmode==Bootstrap) xmlout.put_child("CovarianceCalculationMode","Bootstrap");

 ChiSquareMinimizer CSM(chisq_ref,csm_info);
 m_obs->begin();   // start with full sample
 double chisq;
 vector<double> params_fullsample;
 RVector coveigvals;
 chisq_ref.setObsMeanCov(coveigvals,correlated);   // set means and covariance using full sample
 XMLHandler xmlcov("CovarianceMatrixEigenvalues");
 for (uint p=0;p<coveigvals.size();++p){
    xmlcov.put_child(string("Eigenvalue")+make_string(p),make_string(coveigvals[p]));}
 xmlout.put_child(xmlcov);
 xmlout.put_child("CovarianceMatrixConditionNumber",
        make_string(coveigvals[coveigvals.size()-1]/coveigvals[0]));

 XMLHandler xmlz;
    // first findMinimum guesses initial parameters
 bool flag;
 if (initial_guesses.empty()) {
   flag=CSM.findMinimum(chisq,params_fullsample,xmlz);
 }
 else if (int(initial_guesses.size()) != nparams) {
    throw(std::invalid_argument("Wrong number of initial guesses"));
 }
 else {
   flag=CSM.findMinimum(initial_guesses,chisq,params_fullsample,xmlz);
 }

 if (xmlz.good()) xmlout.put_child(xmlz);
 if (!flag){
    throw(std::invalid_argument("Fitting with full sample failed"));}
 double chisq_dof=chisq/dof;
 for (int p=0;p<nparams;++p) {
    m_obs->putCurrentSamplingValue(param_infos[p],params_fullsample[p]);
 }
 
 vector<double> start(params_fullsample);
 vector<double> params_sample;

    //   loop over the re-samplings
 for (++(*m_obs);!m_obs->end();++(*m_obs)){
    chisq_ref.setObsMean();   // reset means for this resampling, keep covariance from full
    double chisq_samp;
    bool flag=CSM.findMinimum(start,chisq_samp,params_sample);
    if (!flag){
       throw(std::invalid_argument("Fitting with one of the resamplings failed"));}
    for (int p=0;p<nparams;++p)
       m_obs->putCurrentSamplingValue(param_infos[p],params_sample[p]);}

 FitResult fit_result;
 fit_result.chisq_dof = chisq_dof;
 fit_result.quality = getChiSquareFitQuality(dof, chisq);

 fit_result.bestfit_params.resize(nparams);
 XMLHandler xmlres("BestFitResult");
 xmlres.put_child("NumberObservables",make_string(nobs));
 xmlres.put_child("NumberParameters",make_string(nparams));
 xmlres.put_child("NumberPriors",make_string(npriors));
 xmlres.put_child("DegreesOfFreedom",make_string(dof));
 xmlres.put_child("ChiSquarePerDof",make_string(chisq_dof));
 double fitqual=getChiSquareFitQuality(dof,chisq);
 xmlres.put_child("FitQuality",make_string(fitqual));
 for (int p=0;p<nparams;++p){
    XMLHandler xmlp("FitParameter"+make_string(p));
    XMLHandler xmlpi;
    param_infos[p].output(xmlpi);
    xmlp.put_child(xmlpi);
    fit_result.bestfit_params[p]=m_obs->getEstimate(param_infos[p],mode);
    XMLHandler xmlfp;
    fit_result.bestfit_params[p].output(xmlfp);
    xmlp.put_child(xmlfp);
    xmlres.put_child(xmlp);}

 xmlout.put_child(xmlres);

 return fit_result;
}




   // *************************************************
   //
   // The ratios of incomplete gamma function are defined below:
   //
   //       double Qgamma(double s, double x);    s>0, x>=0
   //       double Pgamma(double s, double x);
   //
   // If an error occurs, an exception is thrown.
   //
   // *************************************************

   // Returns the value of ln(Gamma(xx)) for xx>0

double gammln(double xx)
{
 double x,y,tmp,ser;
 static double cof[6]={76.18009172947146,-86.50532032941677,
                       24.01409824083091,-1.231739572450155,
                       0.1208650973866179e-2,-0.5395239384953e-5};
 int j;
 y=x=xx;
 tmp=x+5.5;
 tmp -= (x+0.5)*log(tmp);
 ser=1.000000000190015;
 for (j=0;j<=5;j++) ser += cof[j]/++y;
 return -tmp+log(2.5066282746310005*ser/x);
}


void gcf(double& gammcf, double a, double x, double& gln)
{
 int i;
 double an,b,c,d,del,h;
 const double eps=3.0e-12;
 const int itmax=100;
 const double fpmin=1.0e-30;

 gln=gammln(a);
 b=x+1.0-a;
 c=1.0/fpmin;
 d=1.0/b;
 h=d;
 for (i=1;i<=itmax;i++) {
   an = -i*(i-a);
   b += 2.0;
   d=an*d+b;
   if (fabs(d) < fpmin) d=fpmin;
   c=b+an/c;
   if (fabs(c) < fpmin) c=fpmin;
   d=1.0/d;
   del=d*c;
   h *= del;
   if (fabs(del-1.0) < eps) break;
   }
 if (i > itmax){
     throw(std::invalid_argument("a too large, itmax too small in gcf"));}
 gammcf=exp(-x+a*log(x)-gln)*h;
}


void gser(double& gamser, double a, double x, double& gln)
{
 int n;
 double sum,del,ap;
 const int itmax=100;
 const double eps=3.0e-12;

 gln=gammln(a);
 if (x <= 0.0) {
    if (x < 0.0) throw(std::invalid_argument("x less than 0 in routine gser"));
    gamser=0.0;
 } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=itmax;n++) {
       ++ap;
       del *= x/ap;
       sum += del;
       if (fabs(del) < fabs(sum)*eps) {
          gamser=sum*exp(-x+a*log(x)-gln);
          return;
       }
    }
    throw(std::invalid_argument("a too large, itmax too small in routine gser"));
 }
}

double Qgamma(double a, double x)
{
 double gamser,gammcf,gln;
 if (x < 0.0 || a <= 0.0){
    throw(std::invalid_argument("Invalid arguments in routine gammq"));}
 if (x < (a+1.0)) {
    gser(gamser,a,x,gln);
    return 1.0-gamser;
 } else {
    gcf(gammcf,a,x,gln);
    return gammcf;
 }
}

double Pgamma(double a, double x)
{
 double gamser,gammcf,gln;
 if (x < 0.0 || a <= 0.0){
    throw(std::invalid_argument("Invalid arguments in routine gammp"));}
 if (x < (a+1.0)) {
    gser(gamser,a,x,gln);
    return gamser;
 } else {
    gcf(gammcf,a,x,gln);
    return 1.0-gammcf;
 }
}

   //  Returns the chi-square quality of fit, given by
   //
   //       Qgamma( dof/2,  chisquare/2 )


double getChiSquareFitQuality(unsigned int dof, double chisquare)
{
 return Qgamma(0.5*double(dof), 0.5*chisquare);
}


// ****************************************************************************

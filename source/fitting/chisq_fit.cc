#include "chisq_fit.h"
#include <cmath>
using namespace std;


// *************************************************************************


void doChiSquareFitting(ChiSquare& chisq_ref, 
                        const ChiSquareMinimizerInfo& csm_info,
                        double& chisq_dof, double& fitqual, 
                        vector<MCEstimate>& bestfit_params,
                        XMLHandler& xmlout)
{
 uint nparams=chisq_ref.getNumberOfParams();
 const vector<MCObsInfo>& param_infos=chisq_ref.getFitParamInfos();
 double dof=double(chisq_ref.getNumberOfObservables()-nparams);
 MCObsHandler *m_obs=chisq_ref.getMCObsHandlerPtr();

 for (uint p=0;p<nparams;++p)
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
 chisq_ref.setObsMeanCov(coveigvals);   // set means and covariance using full sample
 XMLHandler xmlcov("CovarianceMatrixEigenvalues");
 for (uint p=0;p<coveigvals.size();++p){
    xmlcov.put_child(string("Eigenvalue")+make_string(p),make_string(coveigvals[p]));}
 xmlout.put_child(xmlcov);
 xmlout.put_child("CovarianceMatrixConditionNumber",
        make_string(coveigvals[coveigvals.size()-1]/coveigvals[0]));

 XMLHandler xmlz;
    // first findMinimum guesses initial parameters
 bool flag=CSM.findMinimum(chisq,params_fullsample,xmlz);

 if (xmlz.good()) xmlout.put_child(xmlz);
 if (!flag){
    throw(std::invalid_argument("Fitting with full sample failed"));}
 chisq_dof=chisq/dof;
 for (uint p=0;p<nparams;++p)
    m_obs->putCurrentSamplingValue(param_infos[p],params_fullsample[p]);
 
 vector<double> start(params_fullsample);
 vector<double> params_sample;

    //   loop over the re-samplings
 for (++(*m_obs);!m_obs->end();++(*m_obs)){
    chisq_ref.setObsMean();   // reset means for this resampling, keep covariance from full
    double chisq_samp;
    bool flag=CSM.findMinimum(start,chisq_samp,params_sample);
    if (!flag){
       throw(std::invalid_argument("Fitting with one of the resamplings failed"));}
    for (uint p=0;p<nparams;++p)
       m_obs->putCurrentSamplingValue(param_infos[p],params_sample[p]);}

 bestfit_params.resize(nparams);
 XMLHandler xmlres("BestFitResult");
 xmlres.put_child("NumberObservables",make_string(chisq_ref.getNumberOfObservables()));
 xmlres.put_child("NumberParameters",make_string(nparams));
 xmlres.put_child("DegreesOfFreedom",make_string(dof));
 xmlres.put_child("ChiSquarePerDof",make_string(chisq_dof));
 fitqual=getChiSquareFitQuality(dof,chisq);
 xmlres.put_child("FitQuality",make_string(fitqual));
 for (uint p=0;p<nparams;++p){
    XMLHandler xmlp("FitParameter"+make_string(p));
    XMLHandler xmlpi;
    param_infos[p].output(xmlpi);
    xmlp.put_child(xmlpi);
    bestfit_params[p]=m_obs->getEstimate(param_infos[p],mode);
    XMLHandler xmlfp;
    bestfit_params[p].output(xmlfp);
    xmlp.put_child(xmlfp);
    xmlres.put_child(xmlp);}

 xmlout.put_child(xmlres);

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
//is the first chisqr better than the second?
bool is_chi_sqr_better( const double chisqr_dof1, const double chisqr_dof2, const double tol){
    double this_chisqr_dof1 = chisqr_dof1;
    double this_chisqr_dof2 = chisqr_dof2;
    if(chisqr_dof1<1.0){this_chisqr_dof1 = 1.0/chisqr_dof1;}
    if(chisqr_dof2<1.0){this_chisqr_dof2 = 1.0/chisqr_dof2;}
//     if(this_chisqr_dof1<1.){return 0;}
    if((this_chisqr_dof1-1.0)<=tol){return 1;}
    if(this_chisqr_dof1<=this_chisqr_dof2){return 1;}
    return 0;
}

void doMultiSeriesFitting(XMLHandler& fit_xml, const int taskcount , 
                        RealMultiTemporalCorrelatorFit& chisq_ref, 
                        const ChiSquareMinimizerInfo& csm_info,
                        double& chisq_dof, double& fitqual, 
                        vector<MCEstimate>& bestfit_params,
                        XMLHandler& xmlout, uint& final_tmin)
{
        
 uint nparams=chisq_ref.getNumberOfParams();
 const vector<MCObsInfo>& param_infos=chisq_ref.getFitParamInfos();
 const vector<MCObsInfo> corr_t_infos = chisq_ref.getObsInfos();
 MCObsHandler *m_obs=chisq_ref.getMCObsHandlerPtr();
 const uint max_level = chisq_ref.getMaxLevel();

 for (uint p=0;p<nparams;++p){
    if (m_obs->queryFullAndSamplings(param_infos[p]))
        throw(std::invalid_argument(string("Error: samplings already available for parameter ")
             +param_infos[p].str()));}

 SamplingMode mode=chisq_ref.getObsMeansSamplingMode();
 SamplingMode covmode=chisq_ref.getCovMatSamplingMode();
 if (covmode==Jackknife) xmlout.put_child("CovarianceCalculationMode","Jackknife");
 else if (covmode==Bootstrap) xmlout.put_child("CovarianceCalculationMode","Bootstrap");
    
 double dof=double(chisq_ref.getNumberOfObservables()-nparams);

     
 //  chisq_ref.setObsMean();
 bool flag;
 double diff = 0;
 double std = 0;
 double this_fit = 0;
 double this_err = 0;
 double last_fit = 0;
 double last_err = 0;
 vector<uint> tvals = chisq_ref.getTvalues_tot();
 bestfit_params.resize(nparams);
 uint this_nparams, this_max_level = 1;
 
 uint Nt = tvals.size()-2;
 double chisq;
 XMLHandler::copymode cmode= XMLHandler::copy;
 uint stable_tmin[this_max_level] = {Nt,Nt,Nt,Nt,Nt};
//  uint stable_tmin[4] = {0,0,0,0};
 double max_chi_sqr = 1.45;
 double min_chi_sqr = 0.69;
 double stable_tmin_chisqr[this_max_level] = {1000.0*max_chi_sqr,1000.0*max_chi_sqr,1000.0*max_chi_sqr,
                                              1000.0*max_chi_sqr,1000.0*max_chi_sqr,1000.0*max_chi_sqr};
 uint level = 0;
    
 uint best_chi_sqr_tmin[this_max_level] = {Nt,Nt,Nt,Nt,Nt,Nt};
 double best_chisqr = 1000.0*max_chi_sqr;
 double best_chisqrs[this_max_level] = {1000.0*max_chi_sqr,1000.0*max_chi_sqr,1000.0*max_chi_sqr,
                                             1000.0*max_chi_sqr,1000.0*max_chi_sqr,1000.0*max_chi_sqr};
 double best_sym_errs[this_max_level] = {1000.0*max_chi_sqr,1000.0*max_chi_sqr,1000.0*max_chi_sqr,
                                             1000.0*max_chi_sqr,1000.0*max_chi_sqr,1000.0*max_chi_sqr};
 uint best_chisqr_level = 0;
 int success_level = -1;
 bool success = false;
 double n_start = chisq_ref.getInitialGap();
 double E1,E2,R21,n2=1.0+n_start,n3=n_start,n4=n_start,n5=n_start;
 double max_n = 6.0, n_increment = chisq_ref.getRepeatingGaps();
 double this_best_chisqr = 1000.0*max_chi_sqr;
    
 //track min chisq
 uint stable_min_level = 0;
 uint stable_level = 0;
 const double chisqr_tol = 0.175;////cm_info.getChiSquareRelativeTolerance(); //have this be an input paramter?
 double level_chisqr_tol = chisqr_tol;
 uint stable_region_depth = 3; //make into input paramter
 double all_params[this_max_level][nparams];
 double these_params[nparams];
 bool stable_region_success[this_max_level] = {0,0,0,0,0,0};
 bool any_fit_success[this_max_level] = {0,0,0,0,0,0};
 bool good_fit[this_max_level] = {0,0,0,0,0,0};
//  const vector<vector<MCObsInfo>> 
 
 //chek that any higher order eponential is constrained to the max singexp? Don't know of any other way to determine a fit fail
 double max_error = 0.0;
    
 if(Nt<=0) throw(std::invalid_argument(string("Error: not enough data for MultiSeries Fit")));
     
 //single_exponential fit
 bool full_samplings;
    
 for( level = 0; level<this_max_level; level++ ){
     if(level) if (!any_fit_success[level-1]){ 
//          std::cout<<"no fit success"<<std::endl; 
         break;}
     this_nparams = 2*(level+1);
     if(level>=2) this_nparams-=(level-1);
//      if(level==3) this_nparams-=2;
//      if(level==4) this_nparams-=3;
//      if(level==5) this_nparams-=4;
     
     Nt = tvals.size() - this_nparams;
     uint first_stable_tmin = Nt;
     last_fit = 0.0;
     last_err = 1000.0;
     this_fit = 0.0;
     this_err = 1000.0;
     full_samplings = false;
     XMLHandler xmlz0,xmlz1("IntermediateFitResults");
     this_best_chisqr = 1000.0*max_chi_sqr;
     for( uint i = 0; i<Nt; i++){
         bool this_good_fit=1; 
         if( level ){ 
             if( i>best_chi_sqr_tmin[level-1] ){ 
//                  std::cout<<"bigger than best tmin "<<best_chi_sqr_tmin[level-1]<<std::endl;
                 break; 
             }
             if( ( ( (int) level-1)==success_level) && (i>=best_chi_sqr_tmin[level-1]) ){ 
//                  std::cout<<"bigger than best tmin "<<best_chi_sqr_tmin[level-1]<<std::endl;
                 break; 
             }
         }
         
         XMLHandler this_fit_xml(fit_xml,cmode);
         this_fit_xml.seek_child("MinimumTimeSeparation");
         this_fit_xml.seek_first_child();
         this_fit_xml.set_text_content(to_string(tvals[i]));
         this_fit_xml.seek_root();
         RealMultiTemporalCorrelatorFit RTC_multi(this_fit_xml,*m_obs,taskcount);
         RTC_multi.m_model_ptr->set_fit_level(level); 
         if(level>=2) RTC_multi.m_model_ptr->set_n2(n2); 
         if(level>=3) RTC_multi.m_model_ptr->set_n3(n2+n3); 
         if(level>=4) RTC_multi.m_model_ptr->set_n4(n2+n3+n4); 
         if(level>=5) RTC_multi.m_model_ptr->set_n5(n2+n3+n4+n5); 
         chisq = 0;
         ChiSquareMinimizer CSM(RTC_multi,csm_info);
         vector<double> params_fullsample;

         if( i ){
             last_fit = this_fit;
             last_err = this_err;
         }
         
         if( full_samplings ){
             
             m_obs->begin();
             for (uint p=0;p<nparams;++p) {
                 m_obs->eraseSamplings(param_infos[p]);
             }
             for (++(*m_obs);!m_obs->end();++(*m_obs)){
                 for (uint p=0;p<nparams;++p) {
                     m_obs->eraseSamplings(param_infos[p]);
                 }
             }
             full_samplings = false;
         }
         m_obs->begin();   // start with full sample 
         RVector coveigvals;
         RTC_multi.setObsMeanCov(coveigvals);  
         vector<double> start0;
         start0.resize(nparams);
             // first findMinimum guesses initial parameters
         if( level==0 ){
             flag=CSM.findMinimum(chisq,params_fullsample,xmlz0);
         }
         else if( level==1 ){
             if(any_fit_success[level]){
                 start0[0]=all_params[level][0]; 
                 start0[1]=all_params[level][1];
                 start0[2]=all_params[level][2]; 
                 start0[3]=all_params[level][3]; 
             }
             else{
                 start0[0]=all_params[level-1][0];
                 start0[1]=sqrt(all_params[level-1][0]);
                 start0[2]=all_params[level-1][2];
                 start0[3]=1.0; 
             }
             flag=CSM.findMinimum(start0,chisq,params_fullsample,xmlz0);
         }else if(level==2){
             if(any_fit_success[level]){
                 start0[0]=all_params[level][0];
                 start0[1]=all_params[level][1];
                 start0[2]=all_params[level][2];
                 start0[3]=all_params[level][3];
                 start0[4]=all_params[level][4];
             }
             else 
             {
                 start0[0]=all_params[level-1][0];
                 start0[1]=all_params[level-1][1];
                 start0[2]=all_params[level-1][2];
                 start0[3]=all_params[level-1][3];
                 start0[4]=all_params[level-1][3];
             }
             flag=CSM.findMinimum(start0,chisq,params_fullsample,xmlz0);
         }else if(level==3){
             if(any_fit_success[level]){
                 start0[0]=all_params[level][0]; 
                 start0[1]=all_params[level][1];
                 start0[2]=all_params[level][2];
                 start0[3]=all_params[level][3]; 
                 start0[4]=all_params[level][4]; 
                 start0[5]=all_params[level][5]; 
             }
             else 
             {
                 start0[0]=all_params[level-1][0];
                 start0[1]=all_params[level-1][1];
                 start0[2]=all_params[level-1][2];
                 start0[3]=all_params[level-1][3];
                 start0[4]=all_params[level-1][4];
                 start0[5]=all_params[level-1][4];
             }
             flag=CSM.findMinimum(start0,chisq,params_fullsample,xmlz0);
         }else if(level==4){
             if(any_fit_success[level]){
                 start0[0]=all_params[level][0]; 
                 start0[1]=all_params[level][1];
                 start0[2]=all_params[level][2];
                 start0[3]=all_params[level][3]; 
                 start0[4]=all_params[level][4]; 
                 start0[5]=all_params[level][5]; 
                 start0[6]=all_params[level][6]; 
             }
             else 
             {
                 start0[0]=all_params[level-1][0];
                 start0[1]=all_params[level-1][1];
                 start0[2]=all_params[level-1][2];
                 start0[3]=all_params[level-1][3];
                 start0[4]=all_params[level-1][4];
                 start0[5]=all_params[level-1][5];
                 start0[6]=all_params[level-1][5];
             }
             flag=CSM.findMinimum(start0,chisq,params_fullsample,xmlz0);
         }else if(level==5){
             if(any_fit_success[level]){
                 start0[0]=all_params[level][0]; 
                 start0[1]=all_params[level][1];
                 start0[2]=all_params[level][2];
                 start0[3]=all_params[level][3]; 
                 start0[4]=all_params[level][4]; 
                 start0[5]=all_params[level][5]; 
                 start0[6]=all_params[level][6]; 
                 start0[7]=all_params[level][7]; 
             }
             else 
             {
                 start0[0]=all_params[level-1][0];
                 start0[1]=all_params[level-1][1];
                 start0[2]=all_params[level-1][2];
                 start0[3]=all_params[level-1][3];
                 start0[4]=all_params[level-1][4];
                 start0[5]=all_params[level-1][5];
                 start0[6]=all_params[level-1][6];
                 start0[7]=all_params[level-1][6];
             }
             flag=CSM.findMinimum(start0,chisq,params_fullsample,xmlz0);
         }

         if (!flag){
//              std::cout<<"bad fit"<<std::endl;
             continue;
            throw(std::invalid_argument("Fitting with full sample failed for single exponential fit in multifit series with tmin of "+to_string(i)));}

         vector<double> start(params_fullsample);
         vector<double> params_sample;
         bool naan_value = false;
         for (uint p=0;p<nparams;++p){
            if(std::isnan(params_fullsample[p])) naan_value=true;
         }
         if(naan_value){ 
//              std::cout<<"naan"<<std::endl; 
             continue; }
         for (uint p=0;p<nparams;++p){
            m_obs->putCurrentSamplingValue(param_infos[p],params_fullsample[p]);
         }
            //   loop over the re-samplings
         for (++(*m_obs);!m_obs->end();++(*m_obs)){
            RTC_multi.setObsMean();   // reset means for this resampling, keep covariance from full
            double chisq_samp;
            flag=CSM.findMinimum(start,chisq_samp,params_sample);
            if (!flag){
                 this_good_fit = 0; 
                 for (uint p=0;p<nparams;++p)
                   m_obs->putCurrentSamplingValue(param_infos[p],start[p]);
            }else{
                for (uint p=0;p<nparams;++p){
                   if(std::isnan(params_sample[p])) naan_value=true;
                   m_obs->putCurrentSamplingValue(param_infos[p],params_sample[p]);
                }
            }

         }   
         if(naan_value){ 
//              std::cout<<"naan"<<std::endl; 
             continue; }
         full_samplings = true;
         for (uint p=0;p<nparams;++p){
            bestfit_params[p]=m_obs->getEstimate(param_infos[p],mode);
         }
         for (uint p=0;p<nparams;++p){
            these_params[p]=params_fullsample[p];
         }
         
//          if( level==3 ){
//              xmlz0.put_child("FitLevel",make_string(level));
//              xmlz0.put_child("EnergyFitValue",make_string(bestfit_params[0].getFullEstimate()));
//              xmlz0.put_child("Tmin",make_string(tvals[i]));
//              xmlz0.put_child("ChiSqDof",make_string(chisq_dof));
//              if (xmlz0.good()) xmlout.put_child(xmlz0);
//          }
//          if(level>0) if( abs(these_params[1])>100000.0 ) this_good_fit=false;
//          if(level>1) if( abs(these_params[2])>100000.0 ) this_good_fit=false;
//          if(level>2) if( abs(these_params[3])>100000.0 ) this_good_fit=false;
//          if(level>0) if( abs(these_params[5])>100000.0 ) this_good_fit=false;
//          if(level>1) if( abs(these_params[6])>100000.0 ) this_good_fit=false;
//          if(level>2) if( abs(these_params[7])>100000.0 ) this_good_fit=false;
//          if(level>0) if( abs(these_params[5])<0.00001 ) this_good_fit=false;
//          if(level>1) if( abs(these_params[6])<0.00001 ) this_good_fit=false;
//          if(level>2) if( abs(these_params[7])<0.00001 ) this_good_fit=false;
         if(level>=0) if( these_params[2]<0.0 ) this_good_fit=false;
         if(level>=1) if( these_params[3]<0.0 ) this_good_fit=false;
         if(level>=2) if( these_params[4]<0.0 ) this_good_fit=false;
         if(level>=3) if( these_params[5]<0.0 ) this_good_fit=false;
         if(level>=4) if( these_params[6]<0.0 ) this_good_fit=false;
         if(level>=5) if( these_params[7]<0.0 ) this_good_fit=false;
//          if(level==1){
//              double extra_term_at_min_tmin = these_params[4]*these_params[5]*exp( -(these_params[0]+these_params[1])*tvals[0] );
//              if( abs(extra_term_at_min_tmin)<0.00001 ){
//                  this_good_fit=false;
//                  std::cout<<"negligible term"<<std::endl;
//              }
//          }
//          else if(level==2){
//              double extra_term_at_min_tmin = these_params[4]*these_params[5]*these_params[6];
//              extra_term_at_min_tmin *= exp( -(these_params[0]+these_params[1]+these_params[2])*tvals[0] );
//              if( abs(extra_term_at_min_tmin)<0.00001 ){ 
//                  this_good_fit=false;
//                  std::cout<<"negligible term"<<std::endl;
//              }
//          }
//          else if(level==3){
//              double extra_term_at_min_tmin = these_params[4]*these_params[5]*these_params[6]*these_params[7];
//              extra_term_at_min_tmin *= exp( -(these_params[0]+these_params[1]+these_params[2]+these_params[3])*tvals[0] );
//              std::cout<<"negligible term"<<std::endl;
//              if( abs(extra_term_at_min_tmin)<0.00001 ) this_good_fit=false;
//          }
         if(!this_good_fit){ 
//              std::cout<<"bad fit"<<std::endl; 
             continue;}
         
         this_fit=bestfit_params[0].getFullEstimate();
         this_err=bestfit_params[0].getSymmetricError();
         if(this_err>max_error){
             if (level==0) max_error=this_err;
             if (level>0){
//                  std::cout<<"energy error too large"<<std::endl;  
                 continue; }
         }
         any_fit_success[level] = 1;

         dof=double(RTC_multi.getNumberOfObservables()-this_nparams);
         chisq_dof=chisq/dof;
         
         std::cout<<level<<" "<<i<<" "<<chisq_dof<<std::endl;
         //new priority: have the first good chi_sqr determine the new fit, unless, all are bad, then use stable region
         if(is_chi_sqr_better(chisq_dof,this_best_chisqr,0.0)){
             this_best_chisqr = chisq_dof;
         }
         if(is_chi_sqr_better(chisq_dof,best_chisqrs[level],0.0)){
             best_chi_sqr_tmin[level] = i;
             best_chisqrs[level] = chisq_dof;
             for (uint p=0;p<nparams;++p){
                all_params[level][p]=params_fullsample[p];
             }
             best_sym_errs[level] = this_err;
         }
         if(is_chi_sqr_better(chisq_dof,best_chisqr,0.0)){
             best_chi_sqr_tmin[level] = i;
             best_chisqr = chisq_dof;
             best_chisqrs[level] = chisq_dof;
             best_chisqr_level = level;
             for (uint p=0;p<nparams;++p){
                all_params[level][p]=params_fullsample[p];
             }
             best_sym_errs[level] = this_err;
         }
         if( (chisq_dof<max_chi_sqr) && (chisq_dof>min_chi_sqr) ){
//              std::cout<<"best fit found"<<std::endl;
             if(success_level< (int) level) Nt=i+2;
             success_level = (int) level;
             if(i==0){ success = true;break;}
         }
     }
     if (level==0) max_error*=10.0; //make input paramter
     
     xmlz1.put_child("FitLevel",make_string(level));
     xmlz1.put_child("EnergyFitValue",make_string(all_params[level][0]));
     xmlz1.put_child("EnergyErrValue",make_string(best_sym_errs[level]));
     if (level>=2) xmlz1.put_child("N2",make_string(n2));
     if (level>=3) xmlz1.put_child("N3",make_string(n2+n3));
     if (level>=4) xmlz1.put_child("N4",make_string(n2+n3+n4));
     if (level>=5) xmlz1.put_child("N5",make_string(n2+n3+n4+n5));
     xmlz1.put_child("Tmin",make_string(tvals[best_chi_sqr_tmin[level]]));
     xmlz1.put_child("ChiSqDof",make_string(best_chisqrs[level]));
     if (xmlz1.good()) xmlout.put_child(xmlz1);
     
     if (success) break;
     
     if( (level==2) && (n2<=max_n) ){
         if( best_chisqrs[level]==this_best_chisqr ){
             if( (success_level!=(int) level) && (best_chisqr_level!=level) ){
                 n2+=n_increment;
                 any_fit_success[level] = 0;
                 level--;
             }
         }else{
             n2-=n_increment;
             any_fit_success[level] = 1;
         }
     }
     else if( (level==3) && (n3<=max_n) ){
         if( best_chisqrs[level]==this_best_chisqr ){
             if( (success_level!=(int) level) && (best_chisqr_level!=level) ){
                 n3+=n_increment;
                 any_fit_success[level] = 0;
                 level--;
             }
         }else{
             n3-=n_increment;
             any_fit_success[level] = 1;
         }
     }
     else if( (level==4) && (n4<=max_n) ){
         if( best_chisqrs[level]==this_best_chisqr ){
             if( (success_level!=(int) level) && (best_chisqr_level!=level) ){
                 n4+=n_increment;
                 any_fit_success[level] = 0;
                 level--;
             }
         }else{
             n4-=n_increment;
             any_fit_success[level] = 1;
         }
     }
     else if( (level==5) && (n5<=max_n) ){
         if( best_chisqrs[level]==this_best_chisqr ){
             if( (success_level!=(int) level) && (best_chisqr_level!=level) ){
                 n5+=n_increment;
                 any_fit_success[level] = 0;
                 level--;
             }
         }else{
             n5-=n_increment;
             any_fit_success[level] = 1;
         }
     }
 }
    
 if(success_level>=0) level = (uint) success_level;
 else level = best_chisqr_level;
    
//  if(any_fit_success[4]) level=4;
    
 final_tmin = tvals[best_chi_sqr_tmin[level]];
     
 if(full_samplings){
     m_obs->begin();
     for (uint p=0;p<nparams;++p) {
         m_obs->eraseSamplings(param_infos[p]);
     }
     for (++(*m_obs);!m_obs->end();++(*m_obs)){
         for (uint p=0;p<nparams;++p) {
             m_obs->eraseSamplings(param_infos[p]);
         }
     }
     m_obs->begin();
 }

 XMLHandler this_fit_xml(fit_xml,cmode);
 this_fit_xml.seek_child("MinimumTimeSeparation");
 this_fit_xml.seek_first_child();
 this_fit_xml.set_text_content(to_string(final_tmin));
 this_fit_xml.seek_root();
 
 RealMultiTemporalCorrelatorFit RTC_multi(this_fit_xml,*m_obs,taskcount); 
 RTC_multi.m_model_ptr->set_fit_level(level); 
 chisq_ref.m_model_ptr->set_fit_level(level); 
 RTC_multi.m_model_ptr->set_n2(n2); 
 chisq_ref.m_model_ptr->set_n2(n2); 
 RTC_multi.m_model_ptr->set_n3(n2+n3); 
 chisq_ref.m_model_ptr->set_n3(n2+n3); 
 RTC_multi.m_model_ptr->set_n4(n2+n3+n4); 
 chisq_ref.m_model_ptr->set_n4(n2+n3+n4);
 RTC_multi.m_model_ptr->set_n5(n2+n3+n4+n5); 
 chisq_ref.m_model_ptr->set_n5(n2+n3+n4+n5);
    
 RVector coveigvals;
 XMLHandler xmlz;
 chisq = 0;
 ChiSquareMinimizer CSM(RTC_multi,csm_info);
 vector<double> params_fullsample;
 vector<double> start0;

 start0.resize(nparams);
 start0[0]=all_params[level][0];
 start0[1]=all_params[level][1];
 start0[2]=all_params[level][2];
 start0[3]=all_params[level][3];
 start0[4]=all_params[level][4];
 start0[5]=all_params[level][5];
 start0[6]=all_params[level][6];
 start0[7]=all_params[level][7];

 RTC_multi.setObsMeanCov(coveigvals); 

 XMLHandler xmlcov("CovarianceMatrixEigenvalues");
 for (uint p=0;p<coveigvals.size();++p){
    xmlcov.put_child(string("Eigenvalue")+make_string(p),make_string(coveigvals[p]));}
 xmlout.put_child(xmlcov);
 xmlout.put_child("CovarianceMatrixConditionNumber",
        make_string(coveigvals[coveigvals.size()-1]/coveigvals[0]));
 flag=CSM.findMinimum(start0,chisq,params_fullsample,xmlz);
 if (xmlz.good()) xmlout.put_child(xmlz);
 if (!flag){
    throw(std::invalid_argument("Fitting with full sample failed for final exponential fit in multifit series"));}
 this_nparams = (level+1)*2;
 if(level>=2) this_nparams -= (level-1);
//  if(level==3) this_nparams -= 2;
//  if(level==4) this_nparams -= 3;
  
 dof=double(RTC_multi.getNumberOfObservables()-this_nparams);
 chisq_dof=chisq/dof;
 vector<double> start(params_fullsample);
 vector<double> params_sample;
 m_obs->begin();   // start with full sample 
 for (uint p=0;p<nparams;++p)
    m_obs->putCurrentSamplingValue(param_infos[p],params_fullsample[p]);
    //   loop over the re-samplings
 for (++(*m_obs);!m_obs->end();++(*m_obs)){
    RTC_multi.setObsMean();   // reset means for this resampling, keep covariance from full
    double chisq_samp;
    flag=CSM.findMinimum(start,chisq_samp,params_sample);
    for (uint p=0;p<nparams;++p) if(std::isnan(params_sample[p])) throw(std::invalid_argument("Fitting with one of the resamplings failed for final exponential fit in multifit series"));
    if (!flag){
        throw(std::invalid_argument("Fitting with one of the resamplings failed for final exponential fit in multifit series"));
    }
    for (uint p=0;p<nparams;++p) m_obs->putCurrentSamplingValue(param_infos[p],params_sample[p]);
 }
    
    
    
 XMLHandler xmlres("BestFitResult");
 xmlres.put_child("FitLevel",make_string(level));
 xmlres.put_child("FinalTmin",make_string(final_tmin));
 if (level>=2) xmlres.put_child("N2",make_string(n2));
 if (level>=3) xmlres.put_child("N3",make_string(n2+n3));
 if (level>=4) xmlres.put_child("N4",make_string(n2+n3+n4));
 if (level>=5) xmlres.put_child("N5",make_string(n2+n3+n4+n5));
 xmlres.put_child("NumberObservables",make_string(RTC_multi.getNumberOfObservables()));
 xmlres.put_child("NumberParameters",make_string(this_nparams));
 xmlres.put_child("DegreesOfFreedom",make_string(dof));
 xmlres.put_child("ChiSquarePerDof",make_string(chisq_dof));
 fitqual=getChiSquareFitQuality(dof,chisq);
 xmlres.put_child("FitQuality",make_string(fitqual));
 for (uint p=0;p<nparams;++p){
    XMLHandler xmlp("FitParameter"+make_string(p));
    XMLHandler xmlpi;
    param_infos[p].output(xmlpi);
    xmlp.put_child(xmlpi);
    bestfit_params[p]=m_obs->getEstimate(param_infos[p],mode);
    XMLHandler xmlfp;
    bestfit_params[p].output(xmlfp);
    xmlp.put_child(xmlfp);
    xmlres.put_child(xmlp);}

 xmlout.put_child(xmlres);

}

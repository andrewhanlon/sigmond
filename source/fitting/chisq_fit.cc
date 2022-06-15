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

void doMultiSeriesFitting(XMLHandler& fit_xml, const int taskcount , RealMultiTemporalCorrelatorFit& chisq_ref, 
                        const ChiSquareMinimizerInfo& csm_info,
                        double& chisq_dof, double& fitqual, 
                        vector<MCEstimate>& bestfit_params,
                        XMLHandler& xmlout, uint& final_tmin)
{
        
 uint nparams=chisq_ref.getNumberOfParams();
 const vector<MCObsInfo>& param_infos=chisq_ref.getFitParamInfos();
 const vector<MCObsInfo> corr_t_infos = chisq_ref.getObsInfos();
 MCObsHandler *m_obs=chisq_ref.getMCObsHandlerPtr();

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
 
 uint Nt = tvals.size()-2;
 double chisq;
 XMLHandler::copymode cmode= XMLHandler::copy;
 uint stable_tmin[4] = {Nt,Nt,Nt,Nt};
 double max_chi_sqr = 1.99;
 double stable_tmin_chisqr[4] = {100.0*max_chi_sqr,100.0*max_chi_sqr,100.0*max_chi_sqr,100.0*max_chi_sqr};
 uint level = 0;
    
 //track min chisq
//  uint min_level = 0;// set to negative and fail fit if not changed?
 uint stable_min_level = 0;
 uint stable_level = 0;
 const double chisqr_tol = 0.175;////cm_info.getChiSquareRelativeTolerance(); //have this be an input paramter?
 double level_chisqr_tol = chisqr_tol;
 uint stable_region_depth = 3; //make into input paramter
 double all_params[4][nparams];
 bool stable_region_success[4] = {0,0,0,0};
 bool any_fit_success[4] = {0,0,0,0};
 bool good_fit[4] = {0,0,0,0};
 //chek that any higher order eponential is constrained to the max singexp? Don't know of any other way to determine a fit fail
 double max_error = 0.0;
    
 if(Nt<=0) throw(std::invalid_argument(string("Error: not enough data for MultiSeries Fit")));
     
 //single_exponential fit
 bool full_samplings;
    
 for( level = 0; level<4; level++ ){
     uint this_nparams = 2*(level+1);
     Nt = tvals.size() - this_nparams;
     uint first_stable_tmin = Nt;
     last_fit = 0.0;
     last_err = 1000.0;
     this_fit = 0.0;
     this_err = 1000.0;
     full_samplings = false;
     for( uint i = 0; i<Nt; i++){
         bool this_good_fit=1;
         if( i>first_stable_tmin+stable_region_depth ) break; 
         
         XMLHandler xmlz0;
         XMLHandler this_fit_xml(fit_xml,cmode);
         this_fit_xml.seek_child("MinimumTimeSeparation");
         this_fit_xml.seek_first_child();
         this_fit_xml.set_text_content(to_string(tvals[i]));
         this_fit_xml.seek_root();
         RealMultiTemporalCorrelatorFit RTC_multi(this_fit_xml,*m_obs,taskcount);
         RTC_multi.m_model_ptr->set_fit_level(level); 

         chisq = 0;
         ChiSquareMinimizer CSM(RTC_multi,csm_info);
         vector<double> params_fullsample;

         if( i>0 ){
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
         if( level==0 ) flag=CSM.findMinimum(chisq,params_fullsample,xmlz0);
         else if( level==1 ){
             if(stable_region_success[level-1]){
                 start0[0]=all_params[level-1][0]; //use fit to contamination but how?
                 start0[1]=sqrt(all_params[level-1][0]);
                 start0[4]=all_params[level-1][4]*0.8;
                 start0[5]=0.0;
             }else{
                 start0[0]=bestfit_params[0].getFullEstimate(); //use fit to contamination but how?
                 start0[1]=sqrt(bestfit_params[0].getFullEstimate());
                 start0[4]=bestfit_params[4].getFullEstimate()*0.8;
                 start0[5]=0.0;
             }
             flag=CSM.findMinimum(start0,chisq,params_fullsample,xmlz0);
         }else if(level==2){
             if(stable_region_success[level-1]){
                 start0[0]=all_params[level-1][0];
                 start0[1]=all_params[level-1][1];
                 start0[2]=all_params[level-1][1]*2.0;
                 start0[4]=all_params[level-1][4]*0.8;
                 start0[5]=all_params[level-1][5]*0.9;
                 start0[6]=0.0;
             }else{
                 start0[0]=bestfit_params[0].getFullEstimate();
                 start0[1]=bestfit_params[1].getFullEstimate();
                 start0[2]=bestfit_params[1].getFullEstimate()*2.0;
                 start0[4]=bestfit_params[4].getFullEstimate()*0.8;
                 start0[5]=bestfit_params[5].getFullEstimate()*0.9;
                 start0[6]=0.0;
             }
             flag=CSM.findMinimum(start0,chisq,params_fullsample,xmlz0);
         }else if(level==3){
             if(stable_region_success[level-1]){
                 start0[0]=all_params[level-1][0];
                 start0[1]=all_params[level-1][1];
                 start0[2]=all_params[level-1][2];
                 start0[3]=all_params[level-1][2]*2.0;
                 start0[4]=all_params[level-1][4]*0.9;
                 start0[5]=all_params[level-1][5]*0.9;
                 start0[6]=all_params[level-1][6]*0.9;
                 start0[7]=0.0;
             }else{
                 start0[0]=bestfit_params[0].getFullEstimate();
                 start0[1]=bestfit_params[1].getFullEstimate();
                 start0[2]=bestfit_params[2].getFullEstimate();
                 start0[3]=bestfit_params[2].getFullEstimate()*2.0;
                 start0[4]=bestfit_params[4].getFullEstimate()*0.9;
                 start0[5]=bestfit_params[5].getFullEstimate()*0.9;
                 start0[6]=bestfit_params[6].getFullEstimate()*0.9;
                 start0[7]=0.0;
             }
             flag=CSM.findMinimum(start0,chisq,params_fullsample,xmlz0);
         }

         if (!flag){
             continue;
            throw(std::invalid_argument("Fitting with full sample failed for single exponential fit in multifit series with tmin of "+to_string(i)));}

         vector<double> start(params_fullsample);
         vector<double> params_sample;
         for (uint p=0;p<nparams;++p)
            m_obs->putCurrentSamplingValue(param_infos[p],params_fullsample[p]);
            //   loop over the re-samplings
         for (++(*m_obs);!m_obs->end();++(*m_obs)){
            RTC_multi.setObsMean();   // reset means for this resampling, keep covariance from full
            double chisq_samp;
            flag=CSM.findMinimum(start,chisq_samp,params_sample);
             if (!flag){
                 this_good_fit = 0;
//                for (uint p=0;p<nparams;++p) std::cout<<params_sample[p]<<std::endl;
//                throw(std::invalid_argument("Fitting with one of the resamplings failed for three exponential fit in multifit series"));
               for (uint p=0;p<nparams;++p)
                   m_obs->putCurrentSamplingValue(param_infos[p],start[p]);
            }else{
                
                for (uint p=0;p<nparams;++p)
                   m_obs->putCurrentSamplingValue(param_infos[p],params_sample[p]);
            }

         }   
         full_samplings = true;
         for (uint p=0;p<nparams;++p){
            bestfit_params[p]=m_obs->getEstimate(param_infos[p],mode);
         }
         
         //check amplitudes
         if( bestfit_params[4].getFullEstimate()<0.0 ) continue;
         if(level>0) if( bestfit_params[5].getFullEstimate()<0.0 ) continue;
         if(level>1) if( bestfit_params[6].getFullEstimate()<0.0 ) continue;
         if(level>2) if( bestfit_params[7].getFullEstimate()<0.0 ) continue;
         if(level>0) if( bestfit_params[5].getFullEstimate()>10000.0*bestfit_params[4].getFullEstimate() ) continue;
         if(level>1) if( bestfit_params[6].getFullEstimate()>10000.0*bestfit_params[4].getFullEstimate() ) continue;
         if(level>2) if( bestfit_params[7].getFullEstimate()>10000.0*bestfit_params[4].getFullEstimate() ) continue;
             
         //check energies
         if(level>1) if( bestfit_params[2].getFullEstimate()<bestfit_params[1].getFullEstimate() ) continue;
         if(level>2) if( bestfit_params[3].getFullEstimate()<bestfit_params[1].getFullEstimate() ) continue;
         if(level>2) if( bestfit_params[3].getFullEstimate()<bestfit_params[2].getFullEstimate() ) continue;
         
         this_fit=bestfit_params[0].getFullEstimate();
         this_err=bestfit_params[0].getSymmetricError();
         if(this_err>max_error){
             if (level==0) max_error=this_err;
             if (level>0) continue;
         }
         any_fit_success[level] = 1;

         dof=double(RTC_multi.getNumberOfObservables()-this_nparams);
         chisq_dof=chisq/dof;

         if( (i>0) ){
             diff = abs(last_fit-this_fit);
             std = max(last_err,this_err);

             if(std>diff){
                 if(!stable_region_success[level]){
                     first_stable_tmin = i;
                     stable_tmin[level] = i;
                     stable_tmin_chisqr[level] = chisq_dof;
                     for (uint p=0;p<nparams;++p){
                        all_params[level][p]=bestfit_params[p].getFullEstimate();
                     }
                     good_fit[level] = this_good_fit;
                     stable_region_success[level] = 1;
                     stable_level=level;
                     xmlz0.put_child("FitLevel",make_string(level));
                     xmlz0.put_child("EnergyFitValue",make_string(all_params[level][0]));
                     xmlz0.put_child("Tmin",make_string(tvals[stable_tmin[level]]));
                     xmlz0.put_child("ChiSqDof",make_string(stable_tmin_chisqr[level]));
                     if (xmlz0.good()) xmlout.put_child(xmlz0);
                 }else{
                     if( (stable_tmin_chisqr[level]>max_chi_sqr && is_chi_sqr_better(chisq_dof,stable_tmin_chisqr[level],0.0) && this_good_fit) || (this_good_fit && !good_fit[level]) ){
                         stable_tmin[level] = i;
                         stable_tmin_chisqr[level] = chisq_dof;
                         for (uint p=0;p<nparams;++p){
                            all_params[level][p]=bestfit_params[p].getFullEstimate();
                         }
                         xmlz0.put_child("FitLevel",make_string(level));
                         xmlz0.put_child("EnergyFitValue",make_string(all_params[level][0]));
                         xmlz0.put_child("Tmin",make_string(tvals[stable_tmin[level]]));
                         xmlz0.put_child("ChiSqDof",make_string(stable_tmin_chisqr[level]));
                         if (xmlz0.good()) xmlout.put_child(xmlz0);
                         good_fit[level] = this_good_fit;
                     }
                 }
             }
//              else if(stable_region_success[level]){
//                  break;
//              }
         }

     }
     if (level==0) max_error*=2.0; //make input paramter
     
     
     if (!any_fit_success[level]) break;
     if (stable_tmin[level]<2) break;
 }
    
//     std::cout<<stable_level<<std::endl;
//     std::cout<<level<<","<<good_fit[stable_level]<<","<<stable_region_success[stable_level]<<std::endl;
 for( level=stable_level; level>0; level--){
//      std::cout<<level<<","<<good_fit[level]<<","<<stable_region_success[level]<<std::endl;
     if(good_fit[level]&&stable_region_success[level]) break;
 }
 
    
 final_tmin = tvals[stable_tmin[level]];
 

 if(!stable_region_success[level]) {
     throw(std::invalid_argument("Multiseries fit failed. No single exponential fits were possible."));
 }
    
 if(!good_fit[level]) {
     throw(std::invalid_argument("Multiseries fit failed. No single exponential fits were possible."));
 }
 

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
     // first findMinimum guesses initial parameters
 flag=CSM.findMinimum(start0,chisq,params_fullsample,xmlz);
 if (xmlz.good()) xmlout.put_child(xmlz);

 if (!flag){
    throw(std::invalid_argument("Fitting with full sample failed for final exponential fit in multifit series"));}
 dof=double(RTC_multi.getNumberOfObservables()-(level+1)*2);
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
    if (!flag){
        throw(std::invalid_argument("Fitting with one of the resamplings failed for final exponential fit in multifit series"));
    }
    for (uint p=0;p<nparams;++p)
       m_obs->putCurrentSamplingValue(param_infos[p],params_sample[p]);
 }
    
    
    
 XMLHandler xmlres("BestFitResult");
 xmlres.put_child("FitLevel",make_string(level));
 xmlres.put_child("FinalTmin",make_string(final_tmin));
 xmlres.put_child("NumberObservables",make_string(RTC_multi.getNumberOfObservables()));
 xmlres.put_child("NumberParameters",make_string((level+1)*2));
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

#include "task_handler.h"
#include "chisq_anisotropy.h"
#include "chisq_disp.h"
#include "chisq_tcorr.h"
#include "chisq_fit.h"
#include "chisq_logtcorr.h"
#include "create_plots.h"
#include "task_utils.h"

using namespace std;

struct rebinResult {
    MCEstimate energy_fit;
    MCEstimate bs_error_estimate;
    double chi_square;
    double bin_stdn;
    double bin_err;
    uint nbins;
};

double error_factor(uint nbins){
   return sqrt(2.0*(double(nbins)-1.0))/double(nbins);
}

void TaskHandler::doRebinAnalysis(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 string fittype;
 xmlreadchild(xmltask,"Type",fittype,"DoRebinAnalysis");
 xmlout.set_root("DoRebinAnalysis");
 xmlout.put_child("Type",fittype); 

 ChiSquareMinimizerInfo mz_info;  // default minimizer info
 if (xmltask.count_among_children("MinimizerInfo")>0){
    ChiSquareMinimizerInfo mz_user(xmltask);
    mz_info=mz_user;}
 XMLHandler xmlmz;
 mz_info.output(xmlmz);
 xmlout.put_child(xmlmz);

 vector<int> rebin_values;
 string rebin_str;
 xmlreadchild(xmltask,"RebinValues",rebin_str);
 extract_from_string(rebin_str,rebin_values);
 xmlout.put_child("RebinValues", make_string(rebin_values));

 sort(rebin_values.begin(), rebin_values.end()); 
 if(rebin_values[0]<1) throw("Cannot rebin by less than 1. Invalid value passed to RebinValues in DoRebinAnalysis.");
 else if(rebin_values[0]!=1) rebin_values.insert( rebin_values.begin(), 1);
 
 m_obs->setToDefaultSamplingMode();
 m_obs->setCovMatToDefaultSamplingMode();
 SamplingMode mode=m_obs->getCurrentSamplingMode();
 string instr;
 if (xmlreadifchild(xmltask,"SamplingMode",instr)){
    if (instr=="Bootstrap") mode=Bootstrap;
    else if (instr=="Jackknife") mode=Jackknife;
    else throw(std::invalid_argument("Bad sampling mode"));
    if (mode==Bootstrap){
       xmlout.put_child("SamplingMode","Bootstrap");
       m_obs->setToBootstrapMode();
       m_obs->setCovMatToBootstrapMode();}
    else{
       xmlout.put_child("SamplingMode","Jackknife");
       m_obs->setToJackknifeMode();
       m_obs->setCovMatToJackknifeMode();}}
 SamplingMode covcalcmode=m_obs->getCovMatCurrentSamplingMode();
 instr.clear();
 if (xmlreadifchild(xmltask,"CovMatCalcSamplingMode",instr)){
    if (instr=="Bootstrap") covcalcmode=Bootstrap;
    else if (instr=="Jackknife") covcalcmode=Jackknife;
    else throw(std::invalid_argument("Bad cov mat calc sampling mode"));
    if (covcalcmode==Bootstrap){
       xmlout.put_child("CovMatCalcSamplingMode","Bootstrap");
       m_obs->setCovMatToBootstrapMode();}
    else{
       xmlout.put_child("CovMatCalcSamplingMode","Jackknife");
       m_obs->setCovMatToJackknifeMode();}}
 mode=m_obs->getCurrentSamplingMode();

 bool uncorrelated=(xmltask.count_among_children("Uncorrelated")>0) ? true: false;
 if (uncorrelated){
   m_obs->setToUnCorrelated();
   xmlout.put_child("Uncorrelated");}
 else
   m_obs->setToCorrelated();

 //trck rebin 1 information
 double rebin1_std2, rebin1_std2_error_factor;

 map<uint, rebinResult> rebin_results;
 for(uint i=0;i<rebin_values.size();i++){
    XMLHandler rebin_output("RebinFit");

    m_bins_info->setRebin(rebin_values[i]);
    XMLHandler xml_obs(m_getter->m_xmlin,XMLHandler::copy);
    bool precompute = false;
    if(m_obs->isBootstrapMode()) precompute = m_obs->getBootstrapper().isPrecomputeMode();
    MCObsInfo rebin_key = MCObsInfo("RebinFitErr",rebin_values[i]);

   //  for 
    delete m_obs;
    delete m_getter;
    
   //  MCSamplingInfo* new_sampling = new MCSamplingInfo(m_samp_info->getNumberOfReSamplings(m_bins_info) , m_samp_info->getRNGSeed(), m_samp_info->getSkipValue() );
    m_getter=new MCObsGetHandler(xml_obs,*m_bins_info,*m_samp_info);
    m_obs=new MCObsHandler(*m_getter,precompute);
    rebin_output.put_child("RebinValue",make_string(m_bins_info->getRebinFactor()));

    //possibly delete and update seed in sampling info redo fit a lot. 

    vector<MCEstimate> bestfit_params;

    try{
        string fittype; XMLHandler xmlof;
        //rotate + fit!!
        if (xmltask.count("RotateTask")){
            XMLHandler xmlr(xmltask,"RotateTask");
            doCorrMatrixRotation(xmlr,xmlof,taskcount);
        }
        
        if (xmltask.count("TemporalCorrelatorFit")) fittype = "TemporalCorrelatorFit";
        else if (xmltask.count("TemporalCorrelatorInteractionRatioFit")) fittype = "TemporalCorrelatorInteractionRatioFit";

        XMLHandler xmlf(xmltask,fittype);
        XMLHandler xmlfitfinal;

        if(fittype == "TemporalCorrelatorInteractionRatioFit"){
         XMLHandler xmltf(xmlf,XMLHandler::subtree_copy);
         bool erase_orig=false;
         setUpRatioFit( *m_obs, xmlf, xmltf, true, xmlout, erase_orig );
         xmlfitfinal = xmltf;
        } else xmlfitfinal = xmlf;

        RealTemporalCorrelatorFit RTC(xmlfitfinal,*m_obs,taskcount);
        RTC.output(xmlof);
        rebin_output.put_child(xmlof);


        double chisq_dof,qual;
        doChiSquareFitting(RTC,mz_info,chisq_dof,qual,bestfit_params,rebin_output);

        const vector<MCObsInfo>& fitparam_infos=RTC.getFitParamInfos(); //fit to current sample
        

        double fit_mean = bestfit_params[0].getFullEstimate();
        uint j=0;
        for (m_obs->begin();!m_obs->end();++(*m_obs),j++){
            double current_sampling = m_obs->getCurrentSamplingValue(fitparam_infos[0]);
            m_obs->putCurrentSamplingValue(fitparam_infos[0],fit_mean);
            m_obs->begin();
            m_obs->putCurrentSamplingValue(fitparam_infos[0],current_sampling);
            m_obs->setCurrentSamplingIndex(j);
            double stddev = m_obs->getEstimate(fitparam_infos[0]).getSymmetricError();
            m_obs->putCurrentSamplingValue(rebin_key,stddev*stddev);
            m_obs->putCurrentSamplingValue(fitparam_infos[0],current_sampling);
        }
        MCEstimate bs_error_estimate =  m_obs->getEstimate(rebin_key);

        rebinResult this_result;
        if(rebin_values[i]==1){
            rebin1_std2 = bestfit_params[0].getSymmetricError()*bestfit_params[0].getSymmetricError();
            rebin1_std2_error_factor = error_factor(uint(m_bins_info->getNumberOfBins()));
            this_result = rebinResult{bestfit_params[0],bs_error_estimate,chisq_dof,1.0,(2.0*rebin1_std2_error_factor+rebin1_std2_error_factor*rebin1_std2_error_factor),uint(m_bins_info->getNumberOfBins())};
        } else {
            double std2 = bestfit_params[0].getSymmetricError()*bestfit_params[0].getSymmetricError();
            double err_num = error_factor(uint(m_bins_info->getNumberOfBins()));
            double err_den = rebin1_std2_error_factor;
            this_result = rebinResult{bestfit_params[0],bs_error_estimate,chisq_dof,std2/rebin1_std2,std2/rebin1_std2*(err_num+err_den+err_num*err_den), uint(m_bins_info->getNumberOfBins())};
        }
        rebin_results[rebin_values[i]] = this_result;
      
       for (uint k=0;k<fitparam_infos.size();++k)
          m_obs->eraseSamplings(fitparam_infos[k]);

       m_obs->eraseSamplings(rebin_key);
            
    }catch(const std::exception& errmsg){
            rebin_output.put_child("Error",string("DoRebinAnalysis with type TemporalCorrelator encountered an error: ")
                    +string(errmsg.what()));
    }
    xmlout.put_child(rebin_output);
 }
 //generate csv?
 cout<<"rebin,std2n,std2n_err,chi2dof,bserr,bserr_err"<<endl;
 for(map<uint, rebinResult>::iterator it = rebin_results.begin(); it != rebin_results.end(); ++it){
   cout<<it->first<<", "<<it->second.bin_stdn<<", "<<it->second.bin_err<<", "<<it->second.chi_square<<", "<<it->second.bs_error_estimate.getFullEstimate()<<", "<<it->second.bs_error_estimate.getSymmetricError()<<endl;
 }
 cout<<endl;
 
 
 }
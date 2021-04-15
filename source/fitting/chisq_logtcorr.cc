#include "chisq_logtcorr.h"
#include "task_utils.h"
#include <string>
#include <map>
using namespace std;


// *************************************************************


LogRealTemporalCorrelatorFit::LogRealTemporalCorrelatorFit(
                  XMLHandler& xmlin, MCObsHandler& OH, int taskcount)   :  ChiSquare(OH)
{
 XMLHandler xmlf(xmlin,"LogTemporalCorrelatorFit");
 OperatorInfo COp(xmlf);
 m_op=COp;
 m_subt_vev=false;
 if (xmlf.count_among_children("SubtractVEV")>0) m_subt_vev=true;
// if ((m_subt_vev)&&(m_obs->isJackknifeMode()))
//    throw(std::invalid_argument("Must use Bootstrap mode with VEV subtractions"));

 CorrelatorInfo Cor(m_op,m_op);
 T_period=OH.getLatticeTimeExtent();
 uint tmin,tmax;
 xmlreadchild(xmlf,"MinimumTimeSeparation",tmin);
 xmlreadchild(xmlf,"MaximumTimeSeparation",tmax);
 if (tmin<0) throw(std::invalid_argument("Invalid MinimumTimeSeparation"));
 if (tmax<=tmin) throw(std::invalid_argument("Invalid MaximumTimeSeparation"));
 vector<int> texclude;
 string exstr;
 if (xmlreadifchild(xmlf,"ExcludeTimes",exstr)){
    extract_from_string(exstr,texclude);}
 m_tvalues=form_tvalues(tmin,tmax,texclude);  // ascending order

   //  Reads and computes correlator estimates, returning estimates
   //  in the "results" map.  The integer key in this map is the
   //  time separation.  Note that results are also "put" into 
   //  the memory of the MObsHandler pointed to by "moh".

 map<double,MCEstimate> corr_results;
 getCorrelatorEstimates(m_obs,Cor,true,m_subt_vev,RealPart, 
                        m_obs->getCurrentSamplingMode(),corr_results);

// for (map<double,MCEstimate>::const_iterator it=corr_results.begin();it!=corr_results.end();it++)
//    cout << "t = "<<it->first<<"  corr = "<<it->second.getFullEstimate()
//         <<" with error = "<<it->second.getSymmetricError()<<endl;

 m_noisecutoff=0.0;
 xmlreadifchild(xmlf,"LargeTimeNoiseCutoff",m_noisecutoff);

 // check data availability, determine if tmax should be lowered due to noise cutoff

 map<double,MCEstimate>::const_iterator rt;
 for (uint k=0;k<m_tvalues.size();k++){
    uint tt=m_tvalues[k];
    rt=corr_results.find(tt);
    if (rt==corr_results.end())
       throw(std::invalid_argument(string("Data not available for time = ")+make_string(tt)
                      +string(" for ")+Cor.str()));
    if (m_noisecutoff>0.0){
       double corval=std::abs(rt->second.getFullEstimate());
       double err=rt->second.getSymmetricError();
       if (corval<m_noisecutoff*err){ 
          m_tvalues.erase(m_tvalues.begin()+k,m_tvalues.end());
          break;}}}
 //if (m_tvalues.size()<4) throw(std::invalid_argument("Less than 4 points after cutoff"));

 m_nobs=m_tvalues.size();

 XMLHandler xmlm(xmlf,"LogModel");
 string modeltype;
 xmlreadchild(xmlm,"Type",modeltype);

 try{
    create_logtcorr_model(modeltype,T_period,m_model_ptr);
    m_nparams=m_model_ptr->getNumberOfParams();
    m_model_ptr->setupInfos(xmlm,m_fitparam_info,taskcount);}
 catch(const std::exception& errmsg){
    m_model_ptr=0;
    throw(std::invalid_argument(string("Invalid Model in LogRealTemporalCorrelatorFit: ")
                 +string(errmsg.what())));}

 int dof = m_nobs-m_model_ptr->getNumberOfParams();
 if (dof < 1) throw(std::invalid_argument("Degrees of Freedom must be greater than zero"));

 allocate_obs_memory();

 CorrelatorAtTimeInfo CorTime(Cor,0,true,m_subt_vev);
 for (uint k=0;k<m_tvalues.size();k++){
    uint tt=m_tvalues[k];
    CorTime.resetTimeSeparation(tt);
    MCObsInfo log_corrt("logcorr", tt);
    doLogBySamplings(*m_obs,CorTime,log_corrt);
    m_obs_info[k]=log_corrt;}

}


LogRealTemporalCorrelatorFit::~LogRealTemporalCorrelatorFit()
{
 delete m_model_ptr;
}



void LogRealTemporalCorrelatorFit::evalModelPoints(
                               const vector<double>& fitparams,
                               vector<double>& modelpoints) const
{
 for (uint k=0;k<m_tvalues.size();k++){
    uint tt=m_tvalues[k];
    m_model_ptr->evaluate(fitparams,double(tt),modelpoints[k]);}
}


void LogRealTemporalCorrelatorFit::evalGradients(
                               const vector<double>& fitparams,
                               RMatrix& gradients) const
{
 uint nparam=m_model_ptr->getNumberOfParams();
 vector<double> grad(nparam);
 for (uint k=0;k<m_tvalues.size();k++){
    uint tt=m_tvalues[k];
    m_model_ptr->evalGradient(fitparams,double(tt),grad);
    for (int p=0;p<int(nparam);p++) gradients(k,p)=grad[p];}
}


void LogRealTemporalCorrelatorFit::guessInitialParamValues(
                               const RVector& datapoints,
                               vector<double>& fitparams) const
{
 vector<double> corr(datapoints.size());
 for (uint k=0;k<m_tvalues.size();k++){
    corr[k]=datapoints[k];}
 m_model_ptr->guessInitialParamValues(corr,m_tvalues,fitparams);  
}

string LogRealTemporalCorrelatorFit::getParameterName(uint param_index) const
{
 return m_model_ptr->getParameterName(param_index);
}

uint LogRealTemporalCorrelatorFit::getParameterIndex(const string& param_name) const
{
 return m_model_ptr->getParameterIndex(param_name);
}


void LogRealTemporalCorrelatorFit::do_output(XMLHandler& xmlout) const
{
 xmlout.set_root("LogTemporalCorrelatorFit");
 XMLHandler xmlop;
 m_op.output(xmlop);
 xmlout.put_child(xmlop);
 if (m_subt_vev) xmlout.put_child("SubtractVEV");
 xmlout.put_child("TimeSeparations",make_string(m_tvalues));
 if (m_noisecutoff>0.0)
    xmlout.put_child("LargeTimeNoiseCutoff",make_string(m_noisecutoff));
 XMLHandler xmlmodel;
 m_model_ptr->output_tag(xmlmodel);
 xmlout.put_child(xmlmodel); 
}



// *********************************************************************

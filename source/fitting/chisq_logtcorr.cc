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
 if (tmax<(tmin+4)) throw(std::invalid_argument("Invalid MaximumTimeSeparation"));
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

 m_single_exp = new TimeForwardSingleExponential(T_period);

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
 if (m_tvalues.size()<4) throw(std::invalid_argument("Less than 4 points after cutoff"));

 m_nobs=m_tvalues.size();

 XMLHandler xmlmp(xmlf,"Model");
 XMLHandler xmlen(xmlmp,"Energy");
 string obsname;
 xmlreadchild(xmlen,"Name",obsname);
 uint index=taskcount;
 xmlreadifchild(xmlen,"IDIndex",index);
 MCObsInfo enkey(obsname,index);
 m_fitparam_info.push_back(enkey);
 XMLHandler xmlamp(xmlmp,"LogAmplitude");
 xmlreadchild(xmlamp,"Name",obsname);
 index=taskcount;
 xmlreadifchild(xmlamp,"IDIndex",index);
 MCObsInfo ampkey(obsname,index);
 m_fitparam_info.push_back(ampkey);
 m_nparams=m_fitparam_info.size();

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
 delete m_single_exp;
}



void LogRealTemporalCorrelatorFit::evalModelPoints(
                               const vector<double>& fitparams,
                               vector<double>& modelpoints) const
{
 for (uint k=0;k<m_tvalues.size();k++){
    uint tt=m_tvalues[k];
    modelpoints[k] = fitparams[1] - fitparams[0]*double(tt);}
}


void LogRealTemporalCorrelatorFit::evalGradients(
                               const vector<double>& fitparams,
                               RMatrix& gradients) const
{
 for (uint k=0;k<m_tvalues.size();k++){
    uint tt=m_tvalues[k];
    gradients(k,1)=1.0;
    gradients(k,0)=-double(tt);}
}


void LogRealTemporalCorrelatorFit::guessInitialParamValues(
                               const RVector& datapoints,
                               vector<double>& fitparams) const
{
 vector<double> corr(datapoints.size());
 for (uint k=0;k<m_tvalues.size();k++){
    corr[k]=exp(datapoints[k]);}
 m_single_exp->guessInitialParamValues(corr,m_tvalues,fitparams);  
 fitparams[1] = log(fitparams[1]);
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
 XMLHandler xmlmodel("Model");
 XMLHandler xmlp("LogAmplitude");
 XMLHandler xmln;
 m_fitparam_info[1].output(xmln);
 xmlp.put_child(xmln); xmlmodel.put_child(xmlp);
 xmlp.set_root("Energy");
 m_fitparam_info[0].output(xmln);
 xmlp.put_child(xmln); xmlmodel.put_child(xmlp);
 xmlout.put_child(xmlmodel);
}



// *********************************************************************

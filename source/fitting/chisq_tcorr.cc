#include "chisq_tcorr.h"
#include "task_utils.h"
#include <string>
#include <map>
using namespace std;


// *************************************************************


RealTemporalCorrelatorFit::RealTemporalCorrelatorFit(
                  XMLHandler& xmlin, MCObsHandler& OH, int taskcount)   :  ChiSquare(OH)
{
 XMLHandler xmlf(xmlin,"TemporalCorrelatorFit");
 OperatorInfo COp(xmlf);
 m_op=COp;
 m_subt_vev=false;
 if (xmlf.count_among_children("SubtractVEV")>0) m_subt_vev=true;
 if ((m_subt_vev)&&(m_obs->isJackknifeMode()))
    throw(std::invalid_argument("Must use Bootstrap mode with VEV subtractions"));

 CorrelatorInfo Cor(m_op,m_op);
 T_period=OH.getLatticeTimeExtent();
 xmlreadchild(xmlf,"MinimumTimeSeparation",m_tmin);
 xmlreadchild(xmlf,"MaximumTimeSeparation",m_tmax);
 if (m_tmin<0) throw(std::invalid_argument("Invalid MinimumTimeSeparation"));
 if (m_tmax<(m_tmin+4)) throw(std::invalid_argument("Invalid MaximumTimeSeparation"));

   //  Reads and computes correlator estimates, returning estimates
   //  in the "results" map.  The integer key in this map is the
   //  time separation.  Note that results are also "put" into 
   //  the memory of the MObsHandler pointed to by "moh".

 map<int,MCEstimate> corr_results;
 getCorrelatorEstimates(m_obs,Cor,true,m_subt_vev,RealPart, 
                        m_obs->getCurrentSamplingMode(),corr_results);

// for (map<int,MCEstimate>::const_iterator it=corr_results.begin();it!=corr_results.end();it++)
//    cout << "t = "<<it->first<<"  corr = "<<it->second.getFullEstimate()
//         <<" with error = "<<it->second.getSymmetricError()<<endl;

 m_noisecutoff=0.0;
 xmlreadifchild(xmlf,"LargeTimeNoiseCutoff",m_noisecutoff);

 // check data availability, determine if m_tmax should be lowered due to noise cutoff

 map<int,MCEstimate>::const_iterator rt;
 for (uint tt=m_tmin;tt<=m_tmax;++tt){
    rt=corr_results.find(tt);
    if (rt==corr_results.end())
       throw(std::invalid_argument((string("Data not available for time = ")+make_string(tt)
                      +string(" for ")+Cor.str()).c_str()));
    if (m_noisecutoff>0.0){
       double corval=std::abs(rt->second.getFullEstimate());
       double err=rt->second.getSymmetricError();
       if (corval<m_noisecutoff*err){ m_tmax=tt-1; break;}}}
 if (m_tmax<(m_tmin+4)) throw(std::invalid_argument("Less than 4 points after cutoff"));

 m_nobs=m_tmax-m_tmin+1;

 XMLHandler xmlm(xmlf,"Model");
 string modeltype;
 xmlreadchild(xmlm,"Type",modeltype);

 try{
    create_tcorr_model(modeltype,T_period,m_model_ptr);
    m_nparams=m_model_ptr->getNumberOfParams();
    m_model_ptr->setupInfos(xmlm,m_fitparam_info,taskcount);}
 catch(const std::exception& errmsg){
    m_model_ptr=0;
    throw(std::invalid_argument((string("Invalid Model in RealTemporalCorrelatorFit: ")
                 +string(errmsg.what())).c_str()));}

 allocate_obs_memory();

 CorrelatorAtTimeInfo CorTime(Cor,0,true,m_subt_vev);
 for (uint tt=m_tmin;tt<=m_tmax;++tt){
    CorTime.resetTimeSeparation(tt);
    m_obs_info[tt-m_tmin]=CorTime;}

}


RealTemporalCorrelatorFit::~RealTemporalCorrelatorFit()
{
 delete m_model_ptr;
}



void RealTemporalCorrelatorFit::evalModelPoints(
                               const vector<double>& fitparams,
                               vector<double>& modelpoints) const
{
 for (int tt=m_tmin;tt<=int(m_tmax);++tt)
    m_model_ptr->evaluate(fitparams,double(tt),modelpoints[tt-m_tmin]);
}


void RealTemporalCorrelatorFit::evalGradients(
                               const vector<double>& fitparams,
                               RMatrix& gradients) const
{
 uint nparam=m_model_ptr->getNumberOfParams();
 vector<double> grad(nparam);
 for (int k=0, tt=m_tmin; tt<=int(m_tmax);++k,++tt){
    m_model_ptr->evalGradient(fitparams,double(tt),grad);
    for (int p=0;p<int(nparam);p++) gradients(k,p)=grad[p];}
}


void RealTemporalCorrelatorFit::guessInitialParamValues(
                               const RVector& datapoints,
                               vector<double>& fitparams) const
{
 vector<double> corr(datapoints.size());
 for (int k=0, tt=m_tmin; tt<=int(m_tmax);++k,++tt){
    corr[k]=datapoints[k];}
 m_model_ptr->guessInitialParamValues(corr,m_tmin,fitparams);  
}


void RealTemporalCorrelatorFit::do_output(XMLHandler& xmlout) const
{
 xmlout.set_root("TemporalCorrelatorFit");
 XMLHandler xmlop;
 m_op.output(xmlop);
 xmlout.put_child(xmlop);
 if (m_subt_vev) xmlout.put_child("SubtractVEV");
 xmlout.put_child("MinimumTimeSeparation",make_string(m_tmin));
 xmlout.put_child("MaximumTimeSeparation",make_string(m_tmax));
 if (m_noisecutoff>0.0)
    xmlout.put_child("LargeTimeNoiseCutoff",make_string(m_noisecutoff));
 XMLHandler xmlmodel;
 m_model_ptr->output_tag(xmlmodel);
 xmlout.put_child(xmlmodel); 
}




// *********************************************************************



/*

TwoRealTemporalCorrelatorFit::TwoRealTemporalCorrelatorFit(
                  XMLHandler& xmlin, MCObsHandler& OH, int taskcount)   :  ChiSquare(OH)
{
 XMLHandler xmlf(xmlin,"TwoTemporalCorrelatorFit");

 XMLHandler xmlc(xmlf,"Correlator");
 OperatorInfo COp(xmlc);
 m_op=COp;
 m_subt_vev=false;
 if (xmlc.count_among_children("SubtractVEV")>0) m_subt_vev=true;
 if ((m_subt_vev)&&(m_obs->isJackknifeMode()))
    throw(std::invalid_argument("Must use Bootstrap mode with VEV subtractions"));

 CorrelatorInfo Cor(m_op,m_op);
 T_period=OH.getLatticeTimeExtent();
 xmlreadchild(xmlc,"MinimumTimeSeparation",m_tmin);
 xmlreadchild(xmlc,"MaximumTimeSeparation",m_tmax);
 if (m_tmin<0) throw(std::invalid_argument("Invalid MinimumTimeSeparation"));
 if (m_tmax<(m_tmin+4)) throw(std::invalid_argument("Invalid MaximumTimeSeparation"));

   //  Reads and computes correlator estimates, returning estimates
   //  in the "results" map.  The integer key in this map is the
   //  time separation.  Note that results are also "put" into 
   //  the memory of the MObsHandler pointed to by "moh".

 map<int,MCEstimate> corr_results;
 getCorrelatorEstimates(m_obs,Cor,true,m_subt_vev,RealPart, 
                        m_obs->getCurrentSamplingMode(),corr_results);

// for (map<int,MCEstimate>::const_iterator it=corr_results.begin();it!=corr_results.end();it++)
//    cout << "t = "<<it->first<<"  corr = "<<it->second.getFullEstimate()
//         <<" with error = "<<it->second.getSymmetricError()<<endl;

 m_noisecutoff=0.0;
 xmlreadifchild(xmlc,"LargeTimeNoiseCutoff",m_noisecutoff);

 // check data availability, determine if m_tmax should be lowered due to noise cutoff

 map<int,MCEstimate>::const_iterator rt;
 for (uint tt=m_tmin;tt<=m_tmax;++tt){
    rt=corr_results.find(tt);
    if (rt==corr_results.end())
       throw(std::invalid_argument("Data not available for time = ")+make_string(tt)+" for "+Cor.str());
    if (m_noisecutoff>0.0){
       double corval=std::abs(rt->second.getFullEstimate());
       double err=rt->second.getSymmetricError();
       if (corval<m_noisecutoff*err){ m_tmax=tt-1; break;}}}
 if (m_tmax<(m_tmin+4)) throw(std::invalid_argument("Less than 4 points after cutoff"));

 m_nobs=m_tmax-m_tmin+1;

 XMLHandler xmlm(xmlc,"Model");
 string modeltype;
 xmlreadchild(xmlm,"Type",modeltype);

 try{
    create_tcorr_model(modeltype,T_period,m_model_ptr);
    m_nparams=m_model_ptr->getNumberOfParams();
    m_model_ptr->setupInfos(xmlm,m_fitparam_info,taskcount);}
 catch(const string& errmsg){
    m_model_ptr=0;
    throw(std::invalid_argument("Invalid Model in TwoRealTemporalCorrelatorFit: ")+errmsg);}

 allocate_obs_memory();

 CorrelatorAtTimeInfo CorTime(Cor,0,true,m_subt_vev);
 for (uint tt=m_tmin;tt<=m_tmax;++tt){
    CorTime.resetTimeSeparation(tt);
    m_obs_info[tt-m_tmin]=CorTime;}

 XMLHandler xmlref(xmlf,"CorrelatorRef");
}


TwoRealTemporalCorrelatorFit::~TwoRealTemporalCorrelatorFit()
{
 delete m_model_ptr;
 delete mref_model_ptr;
}



void TwoRealTemporalCorrelatorFit::evalModelPoints(
                               const vector<double>& fitparams,
                               vector<double>& modelpoints) const
{
 for (int tt=m_tmin;tt<=int(m_tmax);++tt)
    m_model_ptr->evaluate(fitparams,double(tt),modelpoints[tt-m_tmin]);
}


void TwoRealTemporalCorrelatorFit::evalGradients(
                               const vector<double>& fitparams,
                               RMatrix& gradients) const
{
 uint nparam=m_model_ptr->getNumberOfParams();
 vector<double> grad(nparam);
 for (int k=0, tt=m_tmin; tt<=int(m_tmax);++k,++tt){
    m_model_ptr->evalGradient(fitparams,double(tt),grad);
    for (int p=0;p<int(nparam);p++) gradients(k,p)=grad[p];}
}


void TwoRealTemporalCorrelatorFit::guessInitialParamValues(
                               const RVector& datapoints,
                               vector<double>& fitparams) const
{
 vector<double> corr(datapoints.size());
 for (int k=0, tt=m_tmin; tt<=int(m_tmax);++k,++tt){
    corr[k]=datapoints[k];}
 m_model_ptr->guessInitialParamValues(corr,m_tmin,fitparams);  
}


void TwoRealTemporalCorrelatorFit::do_output(XMLHandler& xmlout) const
{
 xmlout.set_root("TwoTemporalCorrelatorFit");
 XMLHandler xmlop;
 m_op.output(xmlop);
 xmlout.put_child(xmlop);
 if (m_subt_vev) xmlout.put_child("SubtractVEV");
 xmlout.put_child("MinimumTimeSeparation",make_string(m_tmin));
 xmlout.put_child("MaximumTimeSeparation",make_string(m_tmax));
 if (m_noisecutoff>0.0)
    xmlout.put_child("LargeTimeNoiseCutoff",make_string(m_noisecutoff));
 XMLHandler xmlmodel;
 m_model_ptr->output_tag(xmlmodel);
 xmlout.put_child(xmlmodel); 
}

*/
// *********************************************************************

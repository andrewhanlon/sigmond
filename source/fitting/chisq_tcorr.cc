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

 XMLHandler xmlm(xmlf,"Model");
 string modeltype;
 xmlreadchild(xmlm,"Type",modeltype);

 try{
    create_tcorr_model(modeltype,T_period,m_model_ptr);
    m_nparams=m_model_ptr->getNumberOfParams();
    m_model_ptr->setupInfos(xmlm,m_fitparam_info,taskcount);}
 catch(const std::exception& errmsg){
    m_model_ptr=0;
    throw(std::invalid_argument(string("Invalid Model in RealTemporalCorrelatorFit: ")
                 +string(errmsg.what())));}

 int dof = m_nobs-m_model_ptr->getNumberOfParams();
 if (dof < 1) throw(std::invalid_argument("Degrees of Freedom must be greater than zero"));

 allocate_obs_memory();

 CorrelatorAtTimeInfo CorTime(Cor,0,true,m_subt_vev);
 for (uint k=0;k<m_tvalues.size();k++){
    uint tt=m_tvalues[k];
    CorTime.resetTimeSeparation(tt);
    m_obs_info[k]=CorTime;}

}


RealTemporalCorrelatorFit::~RealTemporalCorrelatorFit()
{
 delete m_model_ptr;
}



void RealTemporalCorrelatorFit::evalModelPoints(
                               const vector<double>& fitparams,
                               vector<double>& modelpoints) const
{
 for (uint k=0;k<m_tvalues.size();k++){
    uint tt=m_tvalues[k];
    m_model_ptr->evaluate(fitparams,double(tt),modelpoints[k]);}
}


void RealTemporalCorrelatorFit::evalGradients(
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


void RealTemporalCorrelatorFit::guessInitialParamValues(
                               const RVector& datapoints,
                               vector<double>& fitparams) const
{
 vector<double> corr(datapoints.size());
 for (uint k=0;k<m_tvalues.size();k++){
    corr[k]=datapoints[k];}
 m_model_ptr->guessInitialParamValues(corr,m_tvalues,fitparams);  
}


void RealTemporalCorrelatorFit::do_output(XMLHandler& xmlout) const
{
 xmlout.set_root("TemporalCorrelatorFit");
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




TwoRealTemporalCorrelatorFit::TwoRealTemporalCorrelatorFit(
                  XMLHandler& xmlin, MCObsHandler& OH, int taskcount)   :  ChiSquare(OH)
{
 XMLHandler xmlf(xmlin,"TwoTemporalCorrelatorFit");

 XMLHandler xmlen(xmlf,"EnergyRatio");
 string rationame; int ratioindex;
 xmlreadchild(xmlen,"Name",rationame);
 if (rationame.empty()) throw(std::invalid_argument("Must provide name for energy ratio parameter"));
 ratioindex=taskcount;
 xmlreadifchild(xmlen,"IDIndex",ratioindex);
 m_energyratio=MCObsInfo(rationame,ratioindex);

 XMLHandler xmlc(xmlf,"CorrelatorOne");
 OperatorInfo COp1(xmlc);
 m_op1=COp1;
 m_subt_vev1=false;
 if (xmlc.count_among_children("SubtractVEV")>0) m_subt_vev1=true;
// if ((m_subt_vev)&&(m_obs->isJackknifeMode()))
//    throw(std::invalid_argument("Must use Bootstrap mode with VEV subtractions"));

 XMLHandler xmlt(xmlf,"CorrelatorTwo");
 OperatorInfo COp2(xmlt);
 m_op2=COp2;
 m_subt_vev2=false;
 if (xmlt.count_among_children("SubtractVEV")>0) m_subt_vev2=true;
// if ((m_subt_vev)&&(m_obs->isJackknifeMode()))
//    throw(std::invalid_argument("Must use Bootstrap mode with VEV subtractions"));

 CorrelatorInfo Cor1(m_op1,m_op1);
 T_period=OH.getLatticeTimeExtent();
 uint tmin1,tmax1;
 xmlreadchild(xmlc,"MinimumTimeSeparation",tmin1);
 xmlreadchild(xmlc,"MaximumTimeSeparation",tmax1);
 if (tmin1<0) throw(std::invalid_argument("Invalid MinimumTimeSeparation"));
 if (tmax1<=tmin1) throw(std::invalid_argument("Invalid MaximumTimeSeparation"));
 vector<int> texclude;
 string exstr;
 if (xmlreadifchild(xmlc,"ExcludeTimes",exstr)){
    extract_from_string(exstr,texclude);}
 m_tvalues1=form_tvalues(tmin1,tmax1,texclude); // ascending order

 CorrelatorInfo Cor2(m_op2,m_op2);
 uint tmin2,tmax2;
 xmlreadchild(xmlt,"MinimumTimeSeparation",tmin2);
 xmlreadchild(xmlt,"MaximumTimeSeparation",tmax2);
 if (tmin2<0) throw(std::invalid_argument("Invalid MinimumTimeSeparation"));
 if (tmax2<=tmin2) throw(std::invalid_argument("Invalid MaximumTimeSeparation"));
 texclude.clear();
 exstr.clear();
 if (xmlreadifchild(xmlt,"ExcludeTimes",exstr)){
    extract_from_string(exstr,texclude);}
 m_tvalues2=form_tvalues(tmin2,tmax2,texclude); // ascending order

   //  Reads and computes correlator estimates, returning estimates
   //  in the "results" map.  The integer key in this map is the
   //  time separation.  Note that results are also "put" into 
   //  the memory of the MObsHandler pointed to by "moh".

 map<double,MCEstimate> corr1_results, corr2_results;
 getCorrelatorEstimates(m_obs,Cor1,true,m_subt_vev1,RealPart, 
                        m_obs->getCurrentSamplingMode(),corr1_results);
 getCorrelatorEstimates(m_obs,Cor2,true,m_subt_vev2,RealPart, 
                        m_obs->getCurrentSamplingMode(),corr2_results);

// for (map<double,MCEstimate>::const_iterator it=corr_results.begin();it!=corr_results.end();it++)
//    cout << "t = "<<it->first<<"  corr = "<<it->second.getFullEstimate()
//         <<" with error = "<<it->second.getSymmetricError()<<endl;

 m_noisecutoff1=0.0;
 xmlreadifchild(xmlc,"LargeTimeNoiseCutoff",m_noisecutoff1);
 m_noisecutoff2=0.0;
 xmlreadifchild(xmlt,"LargeTimeNoiseCutoff",m_noisecutoff2);

 // check data availability, determine if m_tmax should be lowered due to noise cutoff

 map<double,MCEstimate>::const_iterator rt;
 for (uint k=0;k<m_tvalues1.size();k++){
    uint tt=m_tvalues1[k];
    rt=corr1_results.find(tt);
    if (rt==corr1_results.end())
       throw(std::invalid_argument(string("Data not available for time = ")+make_string(tt)
                      +string(" for ")+Cor1.str()));
    if (m_noisecutoff1>0.0){
       double corval=std::abs(rt->second.getFullEstimate());
       double err=rt->second.getSymmetricError();
       if (corval<m_noisecutoff1*err){ 
          m_tvalues1.erase(m_tvalues1.begin()+k,m_tvalues1.end());
          break;}}}
 //if (m_tvalues1.size()<4) throw(std::invalid_argument("Less than 4 points after cutoff"));

 for (uint k=0;k<m_tvalues2.size();k++){
    uint tt=m_tvalues2[k];
    rt=corr2_results.find(tt);
    if (rt==corr2_results.end())
       throw(std::invalid_argument(string("Data not available for time = ")+make_string(tt)
                      +string(" for ")+Cor2.str()));
    if (m_noisecutoff2>0.0){
       double corval=std::abs(rt->second.getFullEstimate());
       double err=rt->second.getSymmetricError();
       if (corval<m_noisecutoff2*err){ 
          m_tvalues2.erase(m_tvalues2.begin()+k,m_tvalues2.end());
          break;}}}
 //if (m_tvalues2.size()<4) throw(std::invalid_argument("Less than 4 points after cutoff"));

 m_nobs=m_tvalues1.size()+m_tvalues2.size();

 XMLHandler xmlm(xmlc,"Model");
 string modeltype1;
 xmlreadchild(xmlm,"Type",modeltype1);
 XMLHandler xmlmm(xmlt,"Model");
 string modeltype2;
 xmlreadchild(xmlmm,"Type",modeltype2);

 try{
    create_tcorr_model(modeltype1,T_period,m_model1_ptr);
    create_tcorr_model(modeltype2,T_period,m_model2_ptr);
    m_nparams=m_model1_ptr->getNumberOfParams()+m_model2_ptr->getNumberOfParams();
    vector<MCObsInfo> fitparam_info2; 
    m_model1_ptr->setupInfos(xmlm,m_fitparam_info,taskcount);
    m_model2_ptr->setupInfos(xmlmm,fitparam_info2,taskcount);

       // check for different parameter names in the two correlators
    for (vector<MCObsInfo>::iterator it1=m_fitparam_info.begin();it1!=m_fitparam_info.end();it1++)
    for (vector<MCObsInfo>::iterator it2=fitparam_info2.begin();it2!=fitparam_info2.end();it2++)
       if (*it1==*it2) throw(std::invalid_argument("Identical infos in correlators 1 and 2 in TwoRealTemporalCorrelatorFit"));

    m_fitparam_info.insert(m_fitparam_info.end(),fitparam_info2.begin(),
                           fitparam_info2.end());}
 catch(const std::exception& errmsg){
    m_model1_ptr=m_model2_ptr=0;
    throw(std::invalid_argument(string("Invalid Model in TwoRealTemporalCorrelatorFit: ")
                 +string(errmsg.what())));}

 int dof = m_nobs-(m_model1_ptr->getNumberOfParams()+m_model2_ptr->getNumberOfParams());
 if (dof < 1) throw(std::invalid_argument("Degrees of Freedom must be greater than zero"));

 allocate_obs_memory();

 CorrelatorAtTimeInfo CorTime1(Cor1,0,true,m_subt_vev1);
 for (uint k=0;k<m_tvalues1.size();k++){
    uint tt=m_tvalues1[k];
    CorTime1.resetTimeSeparation(tt);
    m_obs_info[k]=CorTime1;}
 CorrelatorAtTimeInfo CorTime2(Cor2,0,true,m_subt_vev2);
 uint shift=m_tvalues1.size();
 for (uint k=0;k<m_tvalues2.size();k++){
    uint tt=m_tvalues2[k];
    CorTime2.resetTimeSeparation(tt);
    m_obs_info[k+shift]=CorTime2;}

}


TwoRealTemporalCorrelatorFit::~TwoRealTemporalCorrelatorFit()
{
 delete m_model1_ptr;
 delete m_model2_ptr;
}



void TwoRealTemporalCorrelatorFit::evalModelPoints(
                               const vector<double>& fitparams,
                               vector<double>& modelpoints) const
{
 uint nparam1=m_model1_ptr->getNumberOfParams();
 vector<double> fitparams1(fitparams.begin(),fitparams.begin()+nparam1);
 vector<double> fitparams2(fitparams.begin()+nparam1,fitparams.end());
 for (uint k=0;k<m_tvalues1.size();k++){
    uint tt=m_tvalues1[k];
    m_model1_ptr->evaluate(fitparams1,double(tt),modelpoints[k]);}
 uint shift=m_tvalues1.size();
 for (uint k=0;k<m_tvalues2.size();k++){
    uint tt=m_tvalues2[k];
    m_model2_ptr->evaluate(fitparams2,double(tt),modelpoints[k+shift]);}
}


void TwoRealTemporalCorrelatorFit::evalGradients(
                               const vector<double>& fitparams,
                               RMatrix& gradients) const
{
 uint nparam1=m_model1_ptr->getNumberOfParams();
 uint nparam2=m_model2_ptr->getNumberOfParams();
 vector<double> fitparams1(fitparams.begin(),fitparams.begin()+nparam1);
 vector<double> fitparams2(fitparams.begin()+nparam1,fitparams.end());
 vector<double> grad1(nparam1);
 for (uint k=0;k<m_tvalues1.size();k++){
    uint tt=m_tvalues1[k];
    m_model1_ptr->evalGradient(fitparams1,double(tt),grad1);
    for (int p=0;p<int(nparam1);p++) gradients(k,p)=grad1[p];
    for (int p=0;p<int(nparam2);p++) gradients(k,p+nparam1)=0.0;}
 vector<double> grad2(nparam2);
 uint shift=m_tvalues1.size();
 for (uint k=0;k<m_tvalues2.size();k++){
    uint tt=m_tvalues2[k];
    uint kk=k+shift;
    m_model2_ptr->evalGradient(fitparams2,double(tt),grad2);
    for (int p=0;p<int(nparam1);p++) gradients(kk,p)=0.0;
    for (int p=0;p<int(nparam2);p++) gradients(kk,p+nparam1)=grad2[p];}
}


void TwoRealTemporalCorrelatorFit::guessInitialParamValues(
                               const RVector& datapoints,
                               vector<double>& fitparams) const
{
 vector<double> corr1(m_tvalues1.size());
 for (uint k=0;k<m_tvalues1.size();k++){
    corr1[k]=datapoints[k];}
 vector<double> fitparams1(m_model1_ptr->getNumberOfParams());
 m_model1_ptr->guessInitialParamValues(corr1,m_tvalues1,fitparams1);  
 vector<double> corr2(m_tvalues2.size());
 uint shift=m_tvalues1.size();
 for (uint k=0;k<m_tvalues2.size();k++){
    corr2[k]=datapoints[k+shift];}
 vector<double> fitparams2(m_model2_ptr->getNumberOfParams());
 m_model2_ptr->guessInitialParamValues(corr2,m_tvalues2,fitparams2);
 std::copy(fitparams1.begin(),fitparams1.end(),fitparams.begin() );
 std::copy(fitparams2.begin(),fitparams2.end(),
     fitparams.begin()+m_model1_ptr->getNumberOfParams() );
}


void TwoRealTemporalCorrelatorFit::do_output(XMLHandler& xmlout) const
{
 xmlout.set_root("TwoTemporalCorrelatorFit");
 XMLHandler xmlrat("EnergyRatio");
 XMLHandler xmln; 
 m_energyratio.output(xmln);
 xmlrat.put_child(xmln);
 xmlout.put_child(xmlrat);
 XMLHandler xmlc("CorrelatorOne");
 XMLHandler xmlop;
 m_op1.output(xmlop);
 xmlc.put_child(xmlop);
 if (m_subt_vev1) xmlc.put_child("SubtractVEV");
 xmlc.put_child("TimeSeparations",make_string(m_tvalues1));
 if (m_noisecutoff1>0.0)
    xmlc.put_child("LargeTimeNoiseCutoff",make_string(m_noisecutoff1));
 XMLHandler xmlmodel;
 m_model1_ptr->output_tag(xmlmodel);
 xmlc.put_child(xmlmodel); 
 xmlout.put_child(xmlc);
 xmlc.set_root("CorrelatorTwo");
 m_op2.output(xmlop);
 xmlc.put_child(xmlop);
 if (m_subt_vev2) xmlc.put_child("SubtractVEV");
 xmlc.put_child("TimeSeparations",make_string(m_tvalues2));
 if (m_noisecutoff2>0.0)
    xmlc.put_child("LargeTimeNoiseCutoff",make_string(m_noisecutoff2));
 m_model2_ptr->output_tag(xmlmodel);
 xmlc.put_child(xmlmodel); 
 xmlout.put_child(xmlc);
}


// *********************************************************************

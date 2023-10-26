#include "chisq_disp.h"
#include "task_utils.h"
#include <string>
#include <map>
#include <stdexcept>
using namespace std;


// *************************************************************


DispersionFit::DispersionFit(
                  XMLHandler& xmlin, MCObsHandler& OH, int taskcount)   :  ChiSquare(OH)
{
 vector<MCObsInfo> temp_obs_infos;
 vector<MCObsInfo> ensq_infos;
 try{
 XMLHandler xmlf(xmlin,"DispersionFit");
 xmlreadchild(xmlf,"SpatialExtentNumSites",m_lat_spatial_extent);
 if (m_lat_spatial_extent<3) 
    throw(std::invalid_argument("Lattice spatial extent too small for dispersion fit"));
 m_momsq_quantum=6.2831853071795864770/double(m_lat_spatial_extent);
 m_momsq_quantum*=m_momsq_quantum;
 list<XMLHandler> xmlobs=xmlf.find("Energy");
 if (xmlobs.size()<3)
    throw(std::invalid_argument("At least three observations needed in DispersionFit"));
 for (list<XMLHandler>::iterator it=xmlobs.begin();it!=xmlobs.end();it++){
    string obsname;
    xmlreadchild(*it,"Name",obsname);
    uint index=taskcount;
    xmlreadifchild(*it,"IDIndex",index);
    MCObsInfo obskey(obsname,index);
    MCObsInfo sqkey(obsname+"_SQ",index);   // for the squares of the energies
    uint imomsq;
    xmlreadchild(*it,"IntMomSquared",imomsq);
    for (vector<MCObsInfo>::iterator mt=temp_obs_infos.begin();mt!=temp_obs_infos.end();mt++)
       if (*mt==obskey) 
          throw(std::invalid_argument("Keys must all be unique in DispersionFit"));
    temp_obs_infos.push_back(obskey);
    ensq_infos.push_back(sqkey);
    m_imomsq.push_back(imomsq);}
 m_nobs=temp_obs_infos.size();
 XMLHandler xmlc(xmlf,"Coefficient");
 string obsname;
 xmlreadchild(xmlc,"Name",obsname);
 uint index=taskcount;
 xmlreadifchild(xmlc,"IDIndex",index);
 MCObsInfo ckey(obsname,index);
 m_fitparam_info.push_back(ckey);
 XMLHandler xmlm0(xmlf,"RestMassSquared");
 xmlreadchild(xmlm0,"Name",obsname);
 index=taskcount;
 xmlreadifchild(xmlm0,"IDIndex",index);
 MCObsInfo m0key(obsname,index);
 m_fitparam_info.push_back(m0key);
 m_nparams=m_fitparam_info.size();

 allocate_obs_memory();
 for (uint k=0;k<m_nobs;k++)
    m_obs_info[k]=ensq_infos[k];

   // check data availability

 for (uint k=0;k<m_nobs;++k){
    if (m_obs->queryFullAndSamplings(temp_obs_infos[k])){
       doSquareBySamplings(*m_obs,temp_obs_infos[k],m_obs_info[k]);}  // reads samplings
    else{
       throw(std::invalid_argument(string("Observation: ")+temp_obs_infos[k].str()
             +string(" not available")));}}}

 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("Invalid Model in AnisotropyFromDispersionFit: ")
                 +string(errmsg.what())));}
}


DispersionFit::~DispersionFit()
{
}



void DispersionFit::evalModelPoints(
                               const vector<double>& fitparams,
                               vector<double>& modelpoints) const
{
 for (uint k=0;k<m_nobs;++k)
    modelpoints[k]=fitparams[1]+fitparams[0]*m_momsq_quantum*m_imomsq[k];
}


void DispersionFit::evalGradients(
                               const vector<double>& fitparams,
                               RMatrix& gradients) const
{
 for (uint k=0;k<m_nobs;++k){
    gradients(k,1)=1.0;
    gradients(k,0)=m_momsq_quantum*m_imomsq[k];}
}


void DispersionFit::guessInitialParamValues(
                               const RVector& datapoints,
                               vector<double>& fitparams) const
{
 uint p=0;
 for (uint k=1;k<m_nobs;++k){
    if (m_imomsq[k]!=m_imomsq[0]){
       p=k; break;}}
 if (p==0) throw(std::invalid_argument("Could not guess initial parameter values"));
 fitparams[0]=(datapoints[p]-datapoints[0])/
   (m_momsq_quantum*(double(m_imomsq[p])-double(m_imomsq[0])));
 bool flag=false;
 for (uint k=0;k<m_nobs;++k){
    if (m_imomsq[k]==0){
       flag=true; fitparams[1]=datapoints[k]; break;}}
 if (!flag){
    fitparams[1]=datapoints[0]-fitparams[0]*m_momsq_quantum*m_imomsq[0];}
}


void DispersionFit::do_output(XMLHandler& xmlout) const
{
 xmlout.set_root("DispersionFit");
 xmlout.put_child("SpatialExtentNumSites",to_string(m_lat_spatial_extent));
 XMLHandler xmlp("Coefficient");
 XMLHandler xmln;
 m_fitparam_info[0].output(xmln); 
 xmlp.put_child(xmln); xmlout.put_child(xmlp);
 xmlp.set_root("RestMassSquared");
 m_fitparam_info[1].output(xmln); 
 xmlp.put_child(xmln); xmlout.put_child(xmlp);
 for (uint k=0;k<m_nobs;k++){
    xmlp.set_root("EnergyObservable");
    m_obs_info[k].output(xmln); 
    xmlp.put_child(xmln); 
    xmlp.put_child("IntMomSquared",to_string(m_imomsq[k]));
    xmlout.put_child(xmlp);}
}




// *********************************************************************

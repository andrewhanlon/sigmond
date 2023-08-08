#include "prior.h"
using namespace std;

int Prior::seed = 0;

Prior::Prior(XMLHandler& xmlin, MCObsHandler& OH)    : m_obs(&OH)
{
 m_resampled = (xmlin.count_among_children("Name")>0);
 srand(seed);
 seed++;
 if (m_resampled) {
    string obs_name; 
    uint obs_id;
    xmlreadchild(xmlin,"Name",obs_name);
    xmlreadchild(xmlin,"IDIndex",obs_id);
    m_prior = MCObsInfo(obs_name,obs_id);
    m_error = m_obs->getEstimate(m_prior).getSymmetricError();
 }
 else {
    xmlreadchild(xmlin,"Mean",m_mean);
    xmlreadchild(xmlin,"Error",m_error);
    resample();
 }
}

double Prior::mean() const
{
 return m_obs->getCurrentSamplingValue(m_prior);
//  if (m_resampled){
//     return m_obs->getCurrentSamplingValue(m_prior);}
//  else {
//     return m_mean;}
}

double Prior::error() const
{
 return m_error;
//  if (m_resampled){
//     return m_obs->getStandardDeviation(m_prior);}
//  else {
//     return m_error;}
}


void Prior::output(XMLHandler& xmlout) const
{
 if (m_resampled){
    XMLHandler xmlop;
    m_prior.output(xmlop);
    xmlout.put_child(xmlop);
 } else {
    xmlout.put_child("Mean",to_string(m_mean));
    xmlout.put_child("Error",to_string(m_error));
 }
}

void Prior::resample(){
    m_prior = MCObsInfo("Prior"+to_string(m_mean)+"-"+to_string(m_error),0);
    m_obs->begin();
    m_obs->putCurrentSamplingValue(m_prior,m_mean);
//     srand(0); //add input for name and use better rand system
    double p0 = m_mean;
    for (++(*m_obs);!m_obs->end();++(*m_obs)){
        double delta = ((float)(rand())/(float)(RAND_MAX))*2.0*m_error-m_error;
//         double delta = ((float)(rand())/(float)(RAND_MAX))*m_error-(0.5*m_error);
        double r = (float)(rand())/(float)(RAND_MAX);
        double prob_dist = exp( -(p0+delta-m_mean)*(p0+delta-m_mean)/m_error/m_error );
        if(r<prob_dist){
            p0=p0+delta;
        }
//         std::cout<<delta<<" "<<r<<" "<<prob_dist<<" "<<p0<<std::endl;
        m_obs->putCurrentSamplingValue(m_prior,p0); //check this
    }
//     m_resampled = true;
}
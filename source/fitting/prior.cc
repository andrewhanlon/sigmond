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
    m_mean = m_obs->getEstimate(m_prior).getFullEstimate();
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
//  return m_mean;
 if (m_resampled){
    return m_obs->getEstimate(m_prior).getFullEstimate();}
   //  return m_obs->getCurrentSamplingValue(m_prior);}
 else {
    return m_mean;}
}

double Prior::error() const
{
//  return m_error;
 if (m_resampled){
    return m_obs->getEstimate(m_prior).getSymmetricError();}
 else {
    return m_error;}
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
    m_obs->eraseSamplings(m_prior);
    m_obs->begin();
    m_obs->putCurrentSamplingValue(m_prior,m_mean);
//     srand(0); //add input for name and use better rand system
    double p0 = m_mean;
    double accdelta_sum = 0.0;
    double alldelta_sum = 0.0;
    // std::cout<<m_mean<<" "<<m_error<<std::endl;
    double search_multiplier = 5.0;
    for (++(*m_obs);!m_obs->end();++(*m_obs)){
        // double delta = ((float)(rand())/(float)(RAND_MAX))*6.0*m_error-1.1*(3.0*m_error);
        // double delta = ((float)(rand())/(float)(RAND_MAX))*2.0*m_error-(1.01*m_error);
        double delta;
        //loop through this until it finds a good new value
        for( uint i=0; i<500; i++){
            delta = ((float)(rand())/(float)(RAND_MAX))*search_multiplier*m_error-1.0*(0.5*search_multiplier*m_error);
            double r = (float)(rand())/(float)(RAND_MAX);
            double prob_dist = exp( -(p0+delta-m_mean)*(p0+delta-m_mean)/m_error/m_error/2.0 );
            if(r<prob_dist){
                p0=p0+delta;
                accdelta_sum+=delta;
                break;
            } else delta = 0.0;
        }
        alldelta_sum+=delta;
        // std::cout<<delta<<","<<r<<","<<prob_dist<<","<<p0<<std::endl;
        // std::cout<<p0<<std::endl;
        m_obs->putCurrentSamplingValue(m_prior,p0); //check this
    }
    // std::cout<<std::endl;

    //check prior distribution
    // MCEstimate prior_stuff = m_obs->getEstimate(m_prior);
    // XMLHandler check_prior;
    // prior_stuff.output(check_prior);
    // std::cout<<check_prior.output()<<std::endl;
    // std::cout<<"Accepted Delta Sum: "<<accdelta_sum<<std::endl;
    // std::cout<<"All Delta Sum: "<<alldelta_sum<<std::endl;

//     m_resampled = true;
}
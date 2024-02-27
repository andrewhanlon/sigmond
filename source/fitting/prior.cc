#include <random>
#include "prior.h"
using namespace std;

int Prior::m_seed = 0;

Prior::Prior(XMLHandler& xmlin, MCObsHandler& OH)    : m_obs(&OH)
{
 m_resampled = (xmlin.count_among_children("Name")>0);
 srand(m_seed);
 m_seed++;
 //add input to keep after, otherwise erase, 
 //add lognormal, but figure out wherelse in the code needs to be changed
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
    xmlreadchild(xmlin,"Mean",m_inmean);
    xmlreadchild(xmlin,"Error",m_inerror);
    setPriorValues();
    resample();
 }
}

double Prior::mean() const
{
 return m_mean;
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
    xmlout.put_child("Mean",to_string(m_inmean)); //output input for log
    xmlout.put_child("Error",to_string(m_inerror));
 }
}

void Prior::resample(){ //filter the random numbers and then filter back, get this distribution a better way
    uint index=0;
    while(true){
//     for( uint index=0; index<100; index++){
        m_prior = MCObsInfo("Prior"+to_string(m_mean)+"-"+to_string(m_error),index);
        if(!m_obs->queryFullAndSamplings(m_prior)) break;
        index++;
    }
    m_obs->eraseSamplings(m_prior);
    m_obs->begin();

    if(m_type==1) { //log normal here ==> transform inputs here?
      // m_obs->putCurrentSamplingValue(m_prior,m_mean);
      std::default_random_engine generator;

      std::lognormal_distribution<double> distribution(/*mean=*/m_mean, /*stddev=*/m_error);
      double randomNumber;
      for (m_obs->begin();!m_obs->end();++(*m_obs)){
         randomNumber = distribution(generator); //add limitations to the number?
         m_obs->putCurrentSamplingValue(m_prior,randomNumber);
      }
    } else {
      m_obs->putCurrentSamplingValue(m_prior,m_mean);
      std::default_random_engine generator;
      std::normal_distribution<double> distribution(/*mean=*/m_mean, /*stddev=*/m_error);
      double randomNumber;
      for (m_obs->setSamplingNext();!m_obs->end();++(*m_obs)){
         randomNumber = distribution(generator);
         while( (randomNumber<(m_mean-2.0*m_error)) && (randomNumber>(m_mean+2.0*m_error)) ){
             randomNumber = distribution(generator); //add limitations to the number?
         }
         m_obs->putCurrentSamplingValue(m_prior,randomNumber);
      }
    }

    // check prior distribution
   //  std::cout<<std::endl;
   //  std::cout<<m_inmean<<" "<<m_inerror<<" "<<m_type<<std::endl;
   //  std::cout<<m_mean<<" "<<m_error<<" "<<m_type<<std::endl;
   //  MCEstimate prior_stuff = m_obs->getEstimate(m_prior);
   //  XMLHandler check_prior;
   //  prior_stuff.output(check_prior);
   //  std::cout<<check_prior.output()<<std::endl;

}

void Prior::resample_current_index() const{
    if(m_type==1) { 
      std::default_random_engine generator;
      std::lognormal_distribution<double> distribution(/*mean=*/m_mean, /*stddev=*/m_error);
      double randomNumber = distribution(generator); 
      m_obs->putCurrentSamplingValue(m_prior,randomNumber);
    } else {
      std::default_random_engine generator;
      std::normal_distribution<double> distribution(/*mean=*/m_mean, /*stddev=*/m_error);
      double randomNumber = distribution(generator);
      while( (randomNumber<(m_mean-2.0*m_error)) && (randomNumber>(m_mean+2.0*m_error)) ){
         randomNumber = distribution(generator); //add limitations to the number?
      }
      m_obs->putCurrentSamplingValue(m_prior,randomNumber);
    }

}

double Prior::evalPriorResidual(const double param_val) const{
   if(m_type==1) return lognormalResidual(param_val);
   return normalResidual(param_val);
}
double Prior::evalPriorGradient(const double param_val) const{
   if(m_type==1) return lognormalGradient(param_val);
   return normalGradient();
}

double Prior::normalResidual(const double param_val) const{
   return (param_val - m_mean)/m_error;
}
double Prior::lognormalResidual(const double param_val) const{
   return (log(abs(param_val)) - m_mean)/m_error; 
}

double Prior::normalGradient() const{
   return 1.0/m_error;
}
double Prior::lognormalGradient(const double param_val) const{
   return 1.0/(param_val*m_error);
}

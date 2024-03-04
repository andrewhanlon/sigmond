#ifndef PRIOR_H
#define PRIOR_H

#include "xml_handler.h"
#include "mcobs_info.h"
#include "mcobs_handler.h"

// *********************************************************************
// *                                                                   *
// *           <Energy>                                                *
// *               <Name>pion</Name>                                   *
// *               <IDIndex>0</IDIndex>                                *
// *           </Energy>                                               *
// *                                                                   *
// *           <Amplitude>                                             *
// *               <Mean>2.34</Mean>                                   *
// *               <Error>0.24</Error>                                 *
// *           </Amplitude>                                            *
// *                                                                   *
// *********************************************************************

class Prior
{
  MCObsHandler *m_obs;
  bool m_resampled;
  double m_mean;
  double m_error;
  double m_inmean;
  double m_inerror;
  static int m_seed;
  uint m_type = 0; //0=normal, 1=lognormal
  
  void resample();

 public:
  MCObsInfo m_prior;

  Prior(XMLHandler& xmlin, MCObsHandler& OH);

  Prior(double in_mean, double in_error)
      : m_resampled(false), m_mean(in_mean), m_error(in_error) {}

  Prior(MCObsInfo& in_prior, MCObsHandler& OH)
      : m_obs(&OH), m_resampled(true), m_prior(in_prior) {}


  ~Prior(){}
//   ~Prior(){m_obs->eraseSamplings(m_prior);}

  double mean() const;
  double error() const;
  void output(XMLHandler& xmlout) const;
  void setType( uint newtype ){
    if(newtype>1) throw(std::invalid_argument("Unknown prior type.")); 
    m_type=newtype;
    setPriorValues();
    resample();
  }

  double evalPriorResidual(const double param_val) const;
  double evalPriorGradient(const double param_val) const;

  double normalResidual(const double param_val) const;
  double lognormalResidual(const double param_val) const;

  double normalGradient() const;
  double lognormalGradient(const double param_val) const;

  void setPriorValues(){
    if(m_type==1) {
      m_mean = log(m_inmean*m_inmean/sqrt(m_inmean*m_inmean+m_inerror*m_inerror));
      m_error = sqrt(log(1.0+m_inerror*m_inerror/m_inmean/m_inmean));
    } else {
      m_mean = m_inmean;
      m_error = m_inerror;
    }
  }
  
  void resample_current_index() const;
  
};

#endif

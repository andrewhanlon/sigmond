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
  MCObsInfo m_prior;
  static int seed;
  
  void resample();

 public:

  Prior(XMLHandler& xmlin, MCObsHandler& OH);

  Prior(double in_mean, double in_error)
      : m_resampled(false), m_mean(in_mean), m_error(in_error) {}

  Prior(MCObsInfo& in_prior, MCObsHandler& OH)
      : m_obs(&OH), m_resampled(true), m_prior(in_prior) {}


  ~Prior(){}

  double mean() const;
  double error() const;
  void output(XMLHandler& xmlout) const;
  
};

#endif

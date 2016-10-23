#ifndef MC_ESTIMATE_H
#define MC_ESTIMATE_H
#include <vector>
#include <string>
#include "xml_handler.h"
#include <stdexcept>

   
enum SamplingMode { Jackknife, Bootstrap };

// **************************************************************************
// *                                                                        *
// *   The class "MCEstimate" and the struct "SimpleMCEstimate" are         *
// *   defined in this file.                                                *
// *                                                                        *
// *                                                                        *
// *   "SimpleMCEstimate" has two double-valued members:                    *
// *       SimpleMCEstimate sm;                                             *
// *       double mean=sm.mean;      // mean value                          *
// *       double stddev=sm.stddev;  // standard deviation                  *
// *                                                                        *
// *                                                                        *
// *   An object of the class "MCEstimate" holds information about the      *
// *   statistical analysis of a particular observable estimated by the     *
// *   Monte Carlo method.  The information held depends on whether the     *
// *   resampling mode is "Jackknife" or "Bootstrap".                       *
// *                                                                        *
// *       MCEstimate mc;                                                   *
// *                                                                        *
// *   (1) The mean value using the entire ensemble (no resampling).        *
// *                                                                        *
// *       mc.getFullEstimate();                                            *
// *                                                                        *
// *   (2) The average of all resampled mean values.  Will not be exactly   *
// *   the same as the full estimate above, but should be very close.       *
// *                                                                        *
// *       mc.getAverageEstimate();                                         *
// *                                                                        *
// *   (3) The symmetric error from all resampled values.                   *
// *                                                                        *
// *       mc.getSymmetricError();                                          *
// *                                                                        *
// *   Jackknife:                                                           *
// *     Let OJ[k] be the estimate from jackknife resampling "k", the       *
// *     mean value obtained after removing bin "k".  Define                *
// *     OJavg= ( sum_k OJ[k] ) / Nbins,  then the symmetric error          *
// *     is given by                                                        *
// *          symerr = sqrt( (1-1/Nbins) * [sum_k (OJ[k]-OJavg)^2]  )       *
// *                                                                        *
// *   Bootstrap:                                                           *
// *     Let OB[k] be the estimate from bootstrap resampling "k".           *
// *     Define  OBavg = ( sum_k OB[k] ) / Nboot, then the symmetric        *
// *     error is given by                                                  *
// *          symerr = sqrt( [sum_k (OB[k]-OBavg)^2] / Nboot  )             *
// *                                                                        *
// *   (4) Confidence limits in "Bootstrap" mode: 84% of all bootstrap      *
// *   resamplings lie above the lower confidence limit, 84% of all         *
// *   bootstrap resamplings lie below the upper confidence limit,          *
// *   and 50% of all bootstrap resamplings lie above the median.           *
// *                                                                        *
// *       mc.getLowerConfLimit();    // bootstrap mode only                *
// *       mc.getMedian();            // bootstrap mode only                *
// *       mc.getUpperConfLimit();    // bootstrap mode only                *
// *                                                                        *
// **************************************************************************


class MCEstimate
{

   std::vector<double> m_store;
   SamplingMode m_mode;
 
 public:

   MCEstimate();
   MCEstimate(SamplingMode inmode);
   MCEstimate(const MCEstimate& in);
   MCEstimate& operator=(const MCEstimate& in);

   double getFullEstimate() const;
   double getAverageEstimate() const;
   double getSymmetricError() const;
   double getLowerConfLimit() const;    // bootstrap mode only
   double getMedian() const;            // bootstrap mode only
   double getUpperConfLimit() const;    // bootstrap mode only

   bool isJackknifeMode() const;
   bool isBootstrapMode() const;

   std::string output(int indent=0) const;  // XML output 
   std::string str() const;                 // XML output 
   void output(XMLHandler& xmlout) const;   // XML output

   void rescale(double r);

 private:

   void jackassign(double full, double avg, double error);
   void bootassign(double full, double avg, double error,
                   double low, double med, double upp);

   friend class MCObsHandler; 

};




class SimpleMCEstimate
{
   double m_mean;
   double m_stddev;

 public:

   SimpleMCEstimate() : m_mean(0.0), m_stddev(0.1) {}

   SimpleMCEstimate(double in_mean, double in_stddev)
     : m_mean(in_mean), m_stddev(in_stddev)
    {if (in_stddev<=0.0) throw(std::invalid_argument("Cannot have negative or zero error in SimpleMCEstimate"));}

   SimpleMCEstimate(const SimpleMCEstimate& in)
     : m_mean(in.m_mean), m_stddev(in.m_stddev) {}

   SimpleMCEstimate& operator=(const SimpleMCEstimate& in)
     {m_mean=in.m_mean; m_stddev=in.m_stddev; return *this;}

   double mean() const {return m_mean;}

   double stddev() const {return m_stddev;}

        // returns string of form mean(stddev), no scientific notation
   std::string str(unsigned int nerr_digits) const; 

};


// *************************************************************
#endif

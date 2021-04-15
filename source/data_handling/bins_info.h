#ifndef BINS_INFO_H
#define BINS_INFO_H

#include "xml_handler.h"
#include "ensemble_info.h"

// ***************************************************************************
// *                                                                         *
// *  "MCBinsInfo" stores information about a known ensemble of Monte        *
// *  Carlo measurements, as well as binning and configuration omission      *
// *  information (in case of corrupted data).  The constructor takes        *
// *  XML input of the following form:                                       *
// *                                                                         *
// *      <MCBinsInfo>                                                       *
// *        <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo>      *
// *        <TweakEnsemble>  (optional)                                      *
// *           <Rebin>2</Rebin>                                              *
// *           <Omissions>2 7 11</Omissions>                                 *
// *        </TweakEnsemble>                                                 *
// *      </MCBinsInfo>                                                      *
// *                                                                         *
// *   However, an object of class "MCBinsInfo" can be constructed using     *
// *   just the <MCEnsembleInfo> tag (in which case, no omissions and no     *
// *   rebinning is assumed).                                                *
// *                                                                         *
// *   If you are using only a small fraction of the configurations, it is   *
// *   easier to use                                                         *
// *       <TweakEnsemble>                                                   *
// *         <KeepFirst>3</KeepFirst>  (index of first config to keep)       *
// *         <KeepLast>12</KeepLast>   (index of last config to keep)        *
// *       </TweakEnsemble>                                                  *
// *                                                                         *
// *   Rebinning the data can also be done as indicated below.  The          *
// *   rebinning factor is always taken with respect to the original         *
// *   measurements. In other words, using setRebin(2), then setRebin(3)     *
// *   does **not** produce a rebinning by a factor of 6.                    *
// *                                                                         *
// *       MCObsInfo MB(....);                                               *
// *       int rebin=3;  MB.setRebin(rebin);                                 *
// *                                                                         *
// *   Including one or more additional omissions (corrupted measurements,   *
// *   for example) can also be done as shown below.  Omissions are always   *
// *   referenced with respect to the original measurements.   Once the      *
// *   omissions are all taken into account, rebinning is then done using    *
// *   the remaining measurements.                                           *
// *                                                                         *
// *       int index=5;  MB.addOmission(index);                              *
// *       set<int> omissions;                                               *
// *       omissions.insert(41); omissions.insert(324);                      *
// *       MB.addOmissions(omissions);                                       *
// *                                                                         *
// *   All omissions can be removed as shown below.                          *
// *                                                                         *
// *       MB.clearOmissions();                                              *
// *                                                                         *
// *                                                                         *
// ***************************************************************************



class MCBinsInfo
{

   MCEnsembleInfo* m_ensemble;
   unsigned int m_nmeasures;                // total number of measurements in ensemble
   unsigned int m_rebin;                    // rebinning factor (number of meas in each bin)
   std::set<unsigned int> m_omit;           // set of measurements to omit
   unsigned int m_nbins;                    // number of bins

   MCBinsInfo() = delete;


 public:

   MCBinsInfo(XMLHandler& xml_in);

   MCBinsInfo(const MCEnsembleInfo& ens);

   MCBinsInfo(const MCBinsInfo& fin);

   MCBinsInfo& operator=(const MCBinsInfo& fin);

   ~MCBinsInfo(){delete m_ensemble;}


   void setRebin(int rebin);

   void addOmission(int index);

   void addOmissions(std::set<int> indices);

   void clearOmissions();



   const MCEnsembleInfo& getMCEnsembleInfo() const { return *m_ensemble; }

   unsigned int getNumberOfMeasurements() const {return m_nmeasures;}

   unsigned int getNumberOfBins() const {return m_nbins;}

   unsigned int getRebinFactor() const {return m_rebin;}

   const std::set<unsigned int>& getOmissions() const {return m_omit;}

   uint getNumberOfConfigs() const { return m_nmeasures; }

   uint getLatticeTimeExtent() const;

   uint getLatticeXExtent() const;

   uint getLatticeYExtent() const;

   uint getLatticeZExtent() const;



   void checkEqual(const MCBinsInfo& rhs) const;

   bool operator==(const MCBinsInfo& rhs) const;

   bool operator!=(const MCBinsInfo& rhs) const;

   bool isConsistentWith(const MCBinsInfo& rhs) const;  // equal, but rhs.rebin
                                                        // can be multiple of rebin
   void output(XMLHandler& xmlout) const;

   std::string output(int indent=0) const;

   std::string str() const;


 private:

   void reset_nbins();

};


// ***************************************************************
#endif  

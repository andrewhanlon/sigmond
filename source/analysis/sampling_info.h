#ifndef SAMPLING_INFO_H
#define SAMPLING_INFO_H

#include "xml_handler.h"
#include "mc_estimate.h"
#include "bootstrapper.h"
#include "bins_info.h"

// ***************************************************************************
// *                                                                         *
// *  "MCSamplingInfo" stores information about the resampling to use        *
// *  in analyzing Monte Carlo measurements.  The constructor takes          *
// *  XML input of the following form:                                       *
// *                                                                         *
// *      <MCSamplingInfo>                                                   *
// *         <Jackknife/>                                                    *
// *      </MCSamplingInfo>                                                  *
// *                                                                         *
// *      <MCSamplingInfo>                                                   *
// *         <Bootstrapper>                                                  *
// *            <NumberResamplings>2048</NumberResamplings>                  *
// *            <Seed>6754</Seed>                                            *
// *            <BootSkip>127</BootSkip>                                     *
// *         </Bootstrapper>                                                 *
// *      </MCSamplingInfo>                                                  *
// *                                                                         *
// *   Only single-point jackknifing is done since the data can be           *
// *   rebinned.  This class is mainly used for header information in        *
// *   SigMonD sampling files.  Default values for the bootstrap parameters  *
// *   are provided if any of their tags are absent.  If "Seed" is zero,     *
// *   this means set the seed using the current time.                       *
// *                                                                         *
// *                                                                         *
// ***************************************************************************



class MCSamplingInfo
{

   unsigned int icode;     // jackknife if zero, else bootstrap (12 bits skip, 20 bits nsamp)
   unsigned long rngseed;  // unused for jackknife


 public:

   MCSamplingInfo();   // default is jackknife mode

   MCSamplingInfo(XMLHandler& xml_in);

   MCSamplingInfo(const MCSamplingInfo& fin);

   MCSamplingInfo& operator=(const MCSamplingInfo& fin);

   MCSamplingInfo(uint nbootsamp, unsigned long bootseed, uint bootskip);


   ~MCSamplingInfo(){}


   bool isJackknifeMode() const;

   bool isBootstrapMode() const;

   SamplingMode getSamplingMode() const;

   void setToJackknifeMode();

   void setToBootstrapMode(const Bootstrapper& boot);


   unsigned int getNumberOfReSamplings(const MCBinsInfo& binfo) const;

   unsigned long getRNGSeed() const;   // returns zero for jackknife mode

   unsigned int getSkipValue() const;   // returns zero for jackknife mode



   void checkEqual(const MCSamplingInfo& rhs) const;

   bool operator==(const MCSamplingInfo& rhs) const;

   bool operator!=(const MCSamplingInfo& rhs) const;

   void output(XMLHandler& xmlout) const;

   std::string output(int indent=0) const;

   std::string str() const;


 private:

   void encode_bootstrap(unsigned int num_resamplings, 
                         unsigned long seed, unsigned int skip_value);

};


// ***************************************************************
#endif  

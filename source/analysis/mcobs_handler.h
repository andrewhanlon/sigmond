#ifndef MC_OBS_HANDLER_H
#define MC_OBS_HANDLER_H
#include <map>
#include <set>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include "scalar_defs.h"
#include "matrix.h"
#include "bootstrapper.h"
#include "obs_get_handler.h"
#include "mcobs_info.h"
#include "mc_estimate.h"
#include "io_handler.h"

// *********************************************************************************
// *                                                                               *
// *   The class "MCObsHandler", a container class for Monte Carlo data and        *
// *   non-simple observables, is defined in this file.  This is one of the        *
// *   most important and central classes in "sigmond".  An object of this         *
// *   class is used to get the Monte Carlo measurements, store them in            *
// *   memory, and perform statistical analysis on the data.  When calculations    *
// *   are carried out on the data or fits are performed, objects of this          *
// *   class store them so that bootstrap and jackknife errors can be              *
// *   subsequently computed.  Objects of this class also carry out reading        *
// *   and writing of bootstrap/jackknife results to and from files, as well as    *
// *   binned data of simple observables, such as rotated correlators.             *
// *                                                                               *
// *   Two important members of an "MCObsHandler" object include: a reference to   *
// *   an object of the "MCObsGetHandler" class in the member "m_in_handler"       *
// *   which is used for getting the data from files; a pointer to a               *
// *   "Bootstrapper" object in "Bptr" which handles all bootstrap resampling.     *
// *                                                                               *
// *   An "MCObsHandler" object contains the Monte Carlo data in one of three      *
// *   maps:                                                                       *
// *                                                                               *
// *     map<MCObsInfo,RVector > m_obs_simple;                                     *
// *     map<MCObsInfo,pair<RVector,uint> > m_jacksamples;                         *
// *     map<MCObsInfo,pair<RVector,uint> > m_bootsamples;                         *
// *                                                                               *
// *   All simple observables are stored as binned values in "m_obs_simple" as     *
// *   RVector data types, and are accessed with record key of class "MCObsInfo".  *
// *   The binning factor and configuration omissions from "m_in_handler" are      *
// *   used, and the RVector values in the "m_obs_simple" map contain the bins.    *
// *   Recall that a simple observable is one that corresponds to an integrand     *
// *   in a Monte Carlo estimate.                                                  *
// *                                                                               *
// *   Resampling values of both simple and non-simple observables are stored      *
// *   in either "m_jacksamples" or "m_bootsamples".  In these maps, the keys are  *
// *   of class "MCObsInfo" and the values are pairs of RVectors and unsigned      *
// *   integers. The RVectors contain the resampling values: the 0th element       *
// *   is always the value for the full sampling.  For "Jackknife" mode, the       *
// *   index "k" varies from 0 to Nbins, the number of bins, where k=0 is          *
// *   the full sample, k=1 is the first jackknife sampling (bin "k" removed),     *
// *   and so on.  The full sample, stored in the k=0 element, is used to simplify *
// *   the computation of each jackknife sample.  For "Bootstrap" mode, the        *
// *   index "k" ranges from 0 to Nboot, the number of bootstrap resamplings,      *
// *   where k=0 is the full sample, k=1 is the first bootstrap sampling, and so   *
// *   on.  The unsigned integer in each pair indicates how many resamplings       *
// *   have already been calculated and stored in memory: uncalculated values      *
// *   are stored as quiet NaNs; once this unsigned integer equals the size of     *
// *   associated RValue, then all samplings are available.                        *
// *                                                                               *
// *   The Monte Carlo measurements are stored in memory, so beware of memory      *
// *   exhaustion for a large number of observables.  No routines are provided     *
// *   for explicitly reading data, but computing means and covariances causes     *
// *   data to be read from file if not already in memory.  A "clearData" member   *
// *   is provided to release memory.                                              *
// *                                                                               *
// *   Data analysis using both Jackknife and Bootstrap resampling can be done.    *
// *   However, reading/writing sampling data from/to files can only be done       *
// *   using the "default" sampling mode, as specified in the associated           *
// *   "MCObsGetHandler" object.  This **default** sampling mode cannot be         *
// *   changed.                                                                    *
// *                                                                               *
// *   Usage:                                                                      *
// *                                                                               *
// *   (1) The constructor takes an "MCObsGetHandler" object and a boolean.  The   *
// *       ensemble info, rebinning factors, configurations to discard (if any),   *
// *       and default resampling info, including any bootstrapping details, are   *
// *       taken from this object and cannot be changed during execution.  To use  *
// *       a different resampling scheme or binning factor, rerun the program.     *
// *       Once the omissions are all taken into account, rebinning is then done   *
// *       using the remaining measurements.  The boolean value indicates if the   *
// *       Bootstrapper should pre-compute all of its index mappings and store     *
// *       them in memory, or recompute them as needed on the fly.                 *
// *                                                                               *
// *       MCObsGetHandler& in_handler;                                            *
// *       bool boot_precompute=true;                                              *
// *       MCObsHandler MH(in_handler,boot_precompute);                            *
// *                                                                               *
// *   (2) Information about current parameters can be obtained as below.          *
// *                                                                               *
// *       MH.getNumberOfMeasurements();                                           *
// *       MH.getNumberOfBins();                                                   *
// *       MH.getNumberOfBootstrapResamplings();                                   *
// *       MH.getBootstrapper();  // returns const reference                       *
// *                                                                               *
// *   (3) Erasing all data corresponding to a particular "MCObsInfo"              *
// *   or all data can be accomplished as shown below.                             *
// *                                                                               *
// *       MCObsInfo obskey;                                                       *
// *       MH.eraseData(obskey);      // erase simple or nonsimple                 *
// *       MH.eraseSamplings(obskey);  // erase only nonsimple                     *
// *       MH.clearData();         // clears all simple and nonsimple              *
// *       MH.clearSamplings();    // clears only nonsimple                        *
// *                                                                               *
// *   (4) Getting and querying simple observable data is accomplished as          *
// *   shown below.  "getBins" first checks to see if the requested                *
// *   observable is already in memory, and if so, a const reference to it         *
// *   is returned, but if not, then "m_in_handler" is called upon to read         *
// *   the data from files and store the results in memory, throwing an            *
// *   exception if it cannot be found.  "getBin" can be used for a single         *
// *   bin, but note that **all** bins are read into memory if "m_in_handler"      *
// *   is called upon.  "queryBins" returns true if the requested observable       *
// *   is already in memory or can be obtained from file; if not already           *
// *   in memory, the data is **not** read into memory.  "getHermCorrBins"         *
// *   reads the bins, but if the read fails, it then tries the time-flipped       *
// *   correlator matrix element key.                                              *
// *                                                                               *
// *       MCObsInfo obskey;                                                       *
// *       const RVector& binsref=MH.getBins(obskey);                              *
// *       int bin_index=76;                                                       *
// *       double binval=MH.getBin(obskey,bin_index);                              *
// *       bool avail=MH.queryBins(obskey);                                        *
// *                                                                               *
// *   (5) Putting simple observable data into memory is accomplished as           *
// *   shown below.  This is usually done for simple observables derived           *
// *   from the simple data in files, such as rotated correlators or               *
// *   other linear combinations of the correlators or vevs.                       *
// *                                                                               *
// *       RVector newvalues(nbins); <-- computed somehow                          *
// *       MH.putBins(obskey,newvalues);                                           *
// *                                                                               *
// *   (6) Resampling can be done using either the jackknife or the                *
// *   bootstrap method.  Use the subroutines below to set the current method to   *
// *   use for computing observable **means** or to query which method is          *
// *   currently being used.  A current method for computing **covariances** is    *
// *   also maintained, and this can be changed and queried.                       *
// *   Note that if you request a method different from the current                *
// *   method, the data currently stored from that resampling method               *
// *   is still retained.  Also, if you call "setToBootstrapMode()" and            *
// *   "setBootstrapper" has not already been called to set up the                 *
// *   bootstrapping parameters, an exception is thrown. Again, the default        *
// *   sampling mode is the mode for reading/writing samping data from/to files,   *
// *   as specified in the associated "MCObsGetHandler" object.                    *
// *                                                                               *
// *       MH.setToJackknifeMode();                                                *
// *       MH.setToBootstrapMode();                                                *
// *       SamplingMode inmode=Jackknife;                                          *
// *       MH.setSamplingMode(inmode);                                             *
// *       MH.setToDefaultSamplingMode();                                          *
// *                                                                               *
// *       MH.isJackknifeMode();                                                   *
// *       MH.isBootstrapMode();                                                   *
// *                                                                               *
// *       MH.setCovMatToJackknifeMode();                                          *
// *       MH.setCovMatToBootstrapMode();                                          *
// *       MH.setCovMatSamplingMode(inmode);                                       *
// *       MH.setCovMatToDefaultSamplingMode();                                    *
// *                                                                               *
// *       MH.isCovMatJackknifeMode();                                             *
// *       MH.isCovMatBootstrapMode();                                             *
// *                                                                               *
// *    (7) Iterating over the resamplings is often needed, either for             *
// *    getting the values or putting the values into memory.  Instead             *
// *    of looping over an integer index, the routines below should be             *
// *    used.  The "begin" here is the **full** estimate.  So in Bootstrap         *
// *    mode, this loops over 1 + MH.getNumberOfBootstrapResamplings() values,     *
// *    and in Jackknife mode, loops over 1 + MH.getNumberOfBins() values.         *
// *                                                                               *
// *       for (MH.setSamplingBegin();!MH.isSamplingEnd();MH.setSamplingNext()){   *
// *           .... }                                                              *
// *       for (MH.begin();!MH.end();++MH){                                        *
// *           .... }                                                              *
// *                                                                               *
// *    (8) There are member routines for getting expected values from the entire  *
// *    ensemble for an observable "obskey" or from a particular resampling.  For  *
// *    a nonsimple observable (except a vev-subtracted correlator), it is         *
// *    assumed that these expected values were previously computed and "put"      *
// *    in memory, and these routines simple retrieve the values, or an exception  *
// *    is thrown.  For a VEV-subtracted correlator matrix element at one time,    *
// *    these routines call "getBins" for the correlator as well as the VEVs       *
// *    and compute all samplings.  Similarly, for a simple observable, "getBins"  *
// *    is called and all samplings are computed and put into memory.              *
// *                                                                               *
// *         // get mean or expected value using entire ensemble                   *
// *       double fullmean=MH.getFullSampleValue(obskey);                          *
// *                                                                               *
// *         // get observable value using current resampling                      *
// *       double thismean=MH.getCurrentSamplingValue(obskey);                     *
// *                                                                               *
// *         // get full and all samplings                                         *
// *       SamplingMode mode=Jackknife;                                            *
// *       const RVector& sampvec=getFullAndSamplingValues(obskey,mode);           *
// *                                                                               *
// *         // query if all samplings are available (include **full** estimate)   *
// *       bool flag=MH.queryAllSamplings(obskey);                                 *
// *                                                                               *
// *    (9) The expected value of a nonsimple observable from a particular         *
// *    resampling or the entire ensemble must be computed outside of this         *
// *    class, but then the result must be "put" into this class so that           *
// *    statistical analysis using all resampling values can be done.              *
// *    For example, if a fit is done using a chi-square based on the current      *
// *    resampling, the fit values should be "put" so that after all fits are      *
// *    done for all resamplings, errors on the fit values can be determined.      *
// *                                                                               *
// *       double fitvalue=....                                                    *
// *       bool overwrite=false; // default value                                  *
// *       MH.putCurrentSamplingValue(obskey,fitvalue,overwrite);                  *
// *                                                                               *
// *    (10) Covariances for two observables are calculated from the entire        *
// *    ensemble using either a jackknife estimate or a bootstrap estimate.        *
// *    The square root of the jackknife covariance for obskey1==obskey2 is the    *
// *    standard deviation.                                                        *
// *                                                                               *
// *       double thiscov12=MH.getCovariance(obskey1,obskey2,Bootstrap);           *
// *       double thiscov12=MH.getCovariance(obskey1,obskey2); //current samp mode *
// *       double stddev=MH.getStandardDeviation(obskey);                          *
// *                                                                               *
// *    (11) Autocorrelation of a simple observable for a particular               *
// *    Markov "time" separation.                                                  *
// *                                                                               *
// *       uint markovtime=3;                                                      *
// *       double autocorr=MH.getAutoCorrelation(obskey,markovtime);               *
// *                                                                               *
// *    (12) Estimates for all resamplings (excluding **full** estimates). Can be  *
// *    used regardless of the current resampling mode.   Results are also stored  *
// *    in memory.  Notice that to exclude the full estimate, the return type      *
// *    is an "ORVectorRef" and not an "RVector".                                  *
// *                                                                               *
// *       ORVectorRef Vsamps=MH.getJackknifeSamplingValues(obskey);               *
// *       ORVectorRef& vsamps;                                                    *
// *       MH.getJackknifeSamplingValues(obskey,vsamps);                           *
// *                                                                               *
// *       ORVectorRef Vsamps=MH.getBootstrapSamplingValues(obskey);               *
// *       ORVectorRef& vsamps;                                                    *
// *       MH.getBootstrapSamplingValues(obskey,vsamps);                           *
// *                                                                               *
// *       SamplingMode mode=Jackknife;                                            *
// *       ORVectorRef Vsamps=MH.getSamplingValues(obskey,mode);                   *
// *       ORVectorRef& vsamps;                                                    *
// *       MH.getSamplingValues(obskey,samp,mode);                                 *
// *                                                                               *
// *     To include the full estimates in samp[0], use                             *
// *                                                                               *
// *       MH.getFullAndSamplingValues(obskey,samp,mode);                          *
// *                                                                               *
// *                                                                               *
// *    (13) J-point jackknife errors (estimates removing J neighboring bins).     *
// *    Can be used regardless of the current resampling mode.  Results are        *
// *    **not** stored in memory.                                                  *
// *                                                                               *
// *       uint jacksize=4;                                                        *
// *       double err=MH.getJackKnifeError(obskey,jacksize);                       *
// *                                                                               *
// *    (14) Statistical analysis:  These routines first call the routines         *
// *    in (12) above, then return the information as an object of the             *
// *    class "MCEstimate".  See the file "mc_estimate.h" for further              *
// *    information of this class.                                                 *
// *                                                                               *
// *       MCEstimate est=MH.getEstimate(obskey);                                  *
// *       SamplingMode mode=Bootstrap;                                            *
// *       MCEstimate est=MH.getEstimate(obskey,mode);                             *
// *       MCEstimate est=MH.getJackknifeEstimate(obskey);                         *
// *       MCEstimate est=MH.getBootstrapEstimate(obskey);                         *
// *                                                                               *
// *    (15) Input/output of all samplings (including full estimates):             *
// *    To write samplings (the **default** mode only) that are already in memory  *
// *    (all samplings include **full* esimates must be available) to file, use    *
// *                                                                               *
// *       set<MCObsInfo> obskeys; ....                                            *
// *       string filename;                                                        *
// *       XMLHandler xmlout;   (for output)                                       *
// *       WriteMode wmode=Protect; // default value or Update or Overwrite        *
// *       MH.writeSamplingValuesToFile(obskeys,filename,xmlout,putmode);          *
// *                                                                               *
// *    If "filename" does not exist, it will be created.  If "filename" exists    *
// *    and mode is "Overwrite", the old file will be destroyed and completely     *
// *    overwritten.  If "filename" exists but mode is "Protect", then the         *
// *    header is checked for consistency and these bins are added to the          *
// *    file, as long as the observables do not already exist in the file.         *
// *    To enable updating observables that are already in the file, use mode      *
// *    "Update".                                                                  *
// *                                                                               *
// *    To read from file and put into memory,                                     *
// *                                                                               *
// *       MH.readSamplingValuesFromFile(filename,mode,xmlout);                    *
// *       MH.readSamplingValuesFromFile(obskeys,filename,mode,xmlout);            *
// *                                                                               *
// *    The first version puts the samplings of all observables encountered        *
// *    in the file into memory; the second version only puts samplings for        *
// *    those observables in the set "obskey" and outputs an error message         *
// *    if any of these observables cannot be found.  The header information       *
// *    in the file must be consistent with the current MCObsHandler object.       *
// *                                                                               *
// *    The format of the files is that for an IOMap.  The header for the          *
// *    file contains the information                                              *
// *                                                                               *
// *        <SigmondSamplingsFile>                                                 *
// *           <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo>         *
// *           <NumberOfMeasurements>551</NumberOfMeasurements>                    *
// *           <NumberOfBins>274</NumberOfBins>                                    *
// *           <TweakEnsemble>    (optional)                                       *
// *              <Rebin>2</Rebin>                                                 *
// *              <OmitConfig>0</OmitConfig>                                       *
// *              <OmitConfig>1</OmitConfig>                                       *
// *              <OmitConfig>78</OmitConfig>                                      *
// *           </TweakEnsemble>                                                    *
// *           <JackKnifeMode/>  or  <Bootstrapper>                                *
// *                             <NumberResamplings>2048</NumberResamplings>       *
// *                                   <Seed>6754</Seed>                           *
// *                                   <BootSkip>127</BootSkip>                    *
// *                                 </Bootstrapper>                               *
// *        </SigmondSamplingsFile>                                                *
// *                                                                               *
// *    The record keys are unsigned integers, and the values contain              *
// *    the MCObservable XML string, followed by the vector of samplings.          *
// *    An IOMap requires the keys to all have the same size in bytes, so          *
// *    an MCObsInfo cannot be used as a key since its encoding can be             *
// *    different sizes.                                                           *
// *                                                                               *
// *    (16) Input/output of bins of computed values of simple observables,        *
// *    such as rotated correlators.  To write bins that are already in memory     *
// *    to file, use                                                               *
// *                                                                               *
// *       set<MCObsInfo> obskeys; ....                                            *
// *       string filename;                                                        *
// *       XMLHandler xmlout;   (for output)                                       *
// *       PutMode putmode = Protect; // default value or Update or Overwrite      *                          *
// *       MH.writeBinsToFile(obskeys,filename,xmlout,putmode);                    *
// *                                                                               *
// *    If "filename" does not exist, it will be created.  If "filename" exists    *
// *    and mode is "Overwrite", the old file will be destroyed and completely     *
// *    overwritten.  If "filename" exists but mode is "Protect", then the         *
// *    header is checked for consistency and these bins are added to the          *
// *    file, as long as the observables do not already exist in the file.         *
// *    To enable updating observables that are already in the file, use mode      *
// *    "Update".                                                                  *
// *                                                                               *
// *    To read from file and put into memory,                                     *
// *                                                                               *
// *       MH.readBinsFromFile(filename,xmlout);                                   *
// *       MH.readBinsFromFile(obskeys,filename,xmlout);                           *
// *                                                                               *
// *    The first version puts the bins of all observables encountered             *
// *    in the file into memory; the second version only puts samplings for        *
// *    those observables in the set "obskey" and outputs an error message         *
// *    if any of these observables cannot be found.  The header information       *
// *    in the file must be consistent with the current MCObsHandler object.       *
// *                                                                               *
// *    The format of the files is described in the documentation for the          *
// *    classes "BinsGetHandler" and "SamplingsGetHandler".  The record keys are   *
// *    of class "MCObsInfo".  The data types are Vector<double>.  For bin files,  *
// *    the Vector size should be the number of bins; for sampling files, the      *
// *    Vector size should be the number of resamplings PLUS one (for the full     *
// *    sample).  The headers for these files contain have the forms               *
// *                                                                               *
// *       <SigmondBinsFile>                                                       *
// *         <MCBinsInfo> ... </MCBinsInfo>                                        *
// *       </SigmondBinsFile>                                                      *
// *                                                                               *
// *       <SigmondSamplingsFile>                                                  *
// *         <MCBinsInfo> ... </MCBinsInfo>                                        *
// *         <MCSamplingInfo> ... </MCSamplingInfo>                                *
// *       </SigmondSamplingsFile>                                                 *
// *                                                                               *
// *                                                                               *
// *********************************************************************************


class MCObsHandler
{

   MCObsGetHandler& m_in_handler;           // the input handler class
   Bootstrapper* Bptr;                      // pointer to the boot strapper

   std::map<MCObsInfo,RVector > m_obs_simple;        // contains the bins of simple observables
   std::map<MCObsInfo,std::pair<RVector,uint> > m_jacksamples;   // for storing resamplings
   std::map<MCObsInfo,std::pair<RVector,uint> > m_bootsamples;   // for storing resamplings

   SamplingMode m_curr_sampling_mode;   //  0 = Jackknife, 1 = Bootstrap
   uint m_curr_sampling_index;          //  0 = full sample, 1..m_max are boot/jack samplings
   uint m_curr_sampling_max;
   std::map<MCObsInfo,std::pair<RVector,uint> > *m_curr_samples;

   SamplingMode m_curr_covmat_sampling_mode;   // current mode to use when computing covariances
   bool m_is_weighted;

            // prevent copying
#ifndef NO_CXX11
   MCObsHandler() = delete;
   MCObsHandler(const MCObsHandler& indata) = delete;
   MCObsHandler& operator=(const MCObsHandler& indata) = delete;
#else
   MCObsHandler();
   MCObsHandler(const MCObsHandler& indata);
   MCObsHandler& operator=(const MCObsHandler& indata);
#endif


 public:

   MCObsHandler(MCObsGetHandler& in_handler, bool boot_precompute=true);

   ~MCObsHandler();


   unsigned int getNumberOfMeasurements() const;

   unsigned int getNumberOfBins() const;

   SamplingMode getDefaultSamplingMode() const;

   const MCBinsInfo& getBinsInfo() const;

   const MCSamplingInfo& getSamplingInfo() const;

   unsigned int getNumberOfBootstrapResamplings() const;

   unsigned int getRebinFactor() const;

   const std::set<unsigned int>& getOmissions() const;

   uint getLatticeTimeExtent() const;

   uint getLatticeXExtent() const;

   uint getLatticeYExtent() const;

   uint getLatticeZExtent() const;

   void getFileMap(XMLHandler& xmlout) const;

   uint getBinsInMemorySize() const;

   uint getJackknifeSamplingsInMemorySize() const;

   uint getBootstrapSamplingsInMemorySize() const;

   uint getCurrentModeSamplingsInMemorySize() const;

   uint getSamplingsInMemorySize(SamplingMode mode) const;



   const Bootstrapper& getBootstrapper() const;

   const Vector<uint>& getBootstrapperResampling(uint bootindex) const;

   void clearData();

   void eraseData(const MCObsInfo& obskey);

   void clearSamplings();

   void eraseSamplings(const MCObsInfo& obskey);


   const RVector& getBins(const MCObsInfo& obskey);

   void getBins(const MCObsInfo& obskey, RVector& values);

   bool queryBins(const MCObsInfo& obskey);

   bool queryBinsInMemory(const MCObsInfo& obskey);

   double getBin(const MCObsInfo& obskey, int bin_index);

   const RVector& putBins(const MCObsInfo& obskey, const RVector& values);



   void setToJackknifeMode();

   void setToBootstrapMode();

   void setSamplingMode(SamplingMode inmode);

   void setToDefaultSamplingMode();

   bool isJackknifeMode() const;

   bool isBootstrapMode() const;

   bool isDefaultSamplingMode() const;

   SamplingMode getCurrentSamplingMode() const {return m_curr_sampling_mode;}
 

   MCObsHandler& setSamplingBegin();

   MCObsHandler& begin();

   MCObsHandler& setSamplingNext();

   MCObsHandler& operator++();

   bool isSamplingEnd() const;

   bool end() const;


   void setCovMatToJackknifeMode();

   void setCovMatToBootstrapMode();

   void setCovMatSamplingMode(SamplingMode inmode);

   void setCovMatToDefaultSamplingMode();

   bool isCovMatJackknifeMode() const;

   bool isCovMatBootstrapMode() const;

   SamplingMode getCovMatCurrentSamplingMode() const {return m_curr_covmat_sampling_mode;}
 


   bool queryFullAndSamplings(const MCObsInfo& obskey);

   bool queryFullAndSamplings(const MCObsInfo& obskey, SamplingMode mode);

   const RVector& getFullAndSamplingValues(const MCObsInfo& obskey, SamplingMode mode);

   void getFullAndSamplingValues(const MCObsInfo& obskey, 
                                 RVector& samples, SamplingMode mode);

   double getFullSampleValue(const MCObsInfo& obskey);

   double getFullSampleValue(const MCObsInfo& obskey, SamplingMode mode);

   double getCurrentSamplingValue(const MCObsInfo& obskey);

   bool getCurrentSamplingValueMaybe(const MCObsInfo& obskey, double& result);

   void putCurrentSamplingValue(const MCObsInfo& obskey, 
                                double value, bool overwrite=true);


   RVector getJackknifeSamplingValues(const MCObsInfo& obskey);

   void getJackknifeSamplingValues(const MCObsInfo& obskey, RVector& samp);

   RVector getBootstrapSamplingValues(const MCObsInfo& obskey);

   void getBootstrapSamplingValues(const MCObsInfo& obskey, RVector& samp);

   RVector getSamplingValues(const MCObsInfo& obskey, SamplingMode mode);

   void getSamplingValues(const MCObsInfo& obskey, RVector& samp,
                          SamplingMode mode);



   MCEstimate getEstimate(const MCObsInfo& obskey);

   MCEstimate getEstimate(const MCObsInfo& obskey, SamplingMode inmode);

   MCEstimate getJackknifeEstimate(const MCObsInfo& obskey);

   MCEstimate getBootstrapEstimate(const MCObsInfo& obskey);



   double getCovariance(const MCObsInfo& obskey1,
                        const MCObsInfo& obskey2, SamplingMode inmode);

   double getCovariance(const MCObsInfo& obskey1,     // current sampling mode
                        const MCObsInfo& obskey2);

   double getStandardDeviation(const MCObsInfo& obskey);

   double getJackKnifeError(const MCObsInfo& obskey, uint jacksize);


   double getAutoCorrelation(const MCObsInfo& obskey, uint markovtime);



             // read all samplings from file and put into memory (second version
             // only reads those records matching the MCObsInfo objects in "obskeys")
             // NOTE: only the default sampling method can be used.

   void readSamplingValuesFromFile(const std::string& filename, 
                                   XMLHandler& xmlout);

   void readSamplingValuesFromFile(const std::set<MCObsInfo>& obskeys, 
                                   const std::string& filename,
                                   XMLHandler& xmlout);

             // write from memory into file (only if all samplings available)

   void writeSamplingValuesToFile(const std::set<MCObsInfo>& obskeys, 
                                  const std::string& filename,
                                  XMLHandler& xmlout,
                                  WriteMode wmode = Protect);

             // read all bins from file and put into memory (second version
             // only reads those records matching the MCObsInfo objects in "obskeys")

   void readBinsFromFile(const std::string& filename, 
                         XMLHandler& xmlout);

   void readBinsFromFile(const std::set<MCObsInfo>& obskeys, 
                         const std::string& filename,
                         XMLHandler& xmlout);

             // write all bins from memory into file 

   void writeBinsToFile(const std::set<MCObsInfo>& obskeys, 
                        const std::string& filename,
                        XMLHandler& xmlout, WriteMode = Protect);


 private:

   void assert_simple(const MCObsInfo& obskey, const std::string& name);

   const RVector& get_bins(const MCObsInfo& obskey);

   const RVector& get_full_and_sampling_values(const MCObsInfo& obskey, 
                        std::map<MCObsInfo,std::pair<RVector,uint> > *samp_ptr,
                        SamplingMode mode, bool allow_not_all_available=false);

   double get_a_sampling_value(const MCObsInfo& obskey, 
                        std::map<MCObsInfo,std::pair<RVector,uint> > *samp_ptr,
                        SamplingMode mode, uint sampindex);

   const RVector* get_full_and_sampling_values_maybe(const MCObsInfo& obskey, 
                        std::map<MCObsInfo,std::pair<RVector,uint> > *samp_ptr,
                        SamplingMode mode);

   bool get_a_sampling_value_maybe(const MCObsInfo& obskey, 
                        std::map<MCObsInfo,std::pair<RVector,uint> > *samp_ptr,
                        SamplingMode mode, uint sampindex, double& result);

   const RVector& put_samplings_in_memory(const MCObsInfo& obskey,
                        const RVector& samplings,
                        std::map<MCObsInfo,std::pair<RVector,uint> > *samp_ptr);

   void put_a_sampling_in_memory(const MCObsInfo& obskey, uint sampling_index,
                        double value, bool overwrite, uint sampling_max,
                        std::map<MCObsInfo,std::pair<RVector,uint> > *samp_ptr);

   void calc_simple_jack_samples(const RVector& bins, RVector& samplings);

   void calc_simple_boot_samples(const RVector& bins, RVector& samplings);

#ifdef COMPLEXNUMBERS
   void calc_corr_subvev(const RVector& corr_re_samples, const RVector& corr_im_samples,
                         const RVector& snkvev_re_samples, const RVector& snkvev_im_samples,
                         const RVector& srcvev_re_samples, const RVector& srcvev_im_samples,
                         RVector& corrsubvev_re_samples, RVector& corrsubvev_im_samples);
#else
   void calc_corr_subvev(const RVector& corr_samples, const RVector& snkvev_samples, 
                         const RVector& srcvev_samples, RVector& corrsubvev_samples);
#endif

   const RVector& calc_samplings_from_bins(const MCObsInfo& obskey,
                      std::map<MCObsInfo,std::pair<RVector,uint> > *samp_ptr,
                      void (MCObsHandler::*simpcalc_ptr)(const RVector&,RVector&));

   const RVector* calc_corrsubvev_from_samplings(const MCObsInfo& obskey,
                      std::map<MCObsInfo,std::pair<RVector,uint> > *samp_ptr);

   bool query_samplings_from_bins(const MCObsInfo& obskey);

   bool query_from_samplings_file(const MCObsInfo& obskey);

   void jack_analyze(const RVector& sampvals, MCEstimate& result);

   void boot_analyze(RVector& sampvals, MCEstimate& result);

   double jack_covariance(const RVector& sampvals1, 
                          const RVector& sampvals2);

   double boot_covariance(const RVector& sampvals1, 
                          const RVector& sampvals2);


};

// *************************************************************
#endif

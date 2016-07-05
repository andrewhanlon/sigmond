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
// *   An "MCObsHandler" object contains the Monte Carlo data in one of three      *
// *   maps:                                                                       *
// *                                                                               *
// *     map<MCObsInfo,Vector<double> > m_obs_simple;                              *
// *     map<pair<MCObsInfo,uint>,double> m_jacksamples;                           *
// *     map<pair<MCObsInfo,uint>,double> m_bootsamples;                           *
// *                                                                               *
// *   All simple observables are stored in "m_obs_simple" as Vector<double>       *
// *   data types, and are accessed with record key of class "MCObsInfo".          *
// *   Rebinning can be done, and individual Monte Carlo measurements can be       *
// *   omitted.  Each Vector<double> data contains the bins of the simple          *
// *   observable.  Recall that a simple observable is one that corresponds        *
// *   to an integrand in a Monte Carlo estimate.                                  *
// *                                                                               *
// *   Non-simple observables are stored in either "m_jacksamples" or              *
// *   "m_bootsamples" depending on whether the current resampling mode is         *
// *   set to "Jackknife" or "Bootstrap".  In these maps, the "MCObsInfo"          *
// *   is combined with an integer index to form the record key, and the           *
// *   data values are individual doubles.  The meaning of the integer             *
// *   index depends on the resampling mode.   For "Jackknife" mode, the           *
// *   index "k" varies from 0 to Nbins, the number of bins, where k=0 is          *
// *   the full sample, k=1 is the first jackknife sampling (bin "k" removed),     *
// *   and so on.  The full sample, stored in the k=0 element, is used to simplify *
// *   the computation of each jackknife sample.  For "Bootstrap" mode, the        *
// *   index "k" ranges from 0 to Nboot, the number of bootstrap resamplings,      *
// *   where k=0 is the full sample, k=1 is the first bootstrap sampling, and so   *
// *   on.  The full sample stored in k=0 is not really needed, but is used        *
// *   to be consistent with jackknife mode.                                       *
// *                                                                               *
// *   The Monte Carlo measurements are stored in memory, so beware of memory      *
// *   exhaustion for a large number of observables.  No routines are provided     *
// *   for explicitly reading data, but computing means and covariances causes     *
// *   data to be read from file if not already in memory.  A "clearData" member   *
// *   is provided to release memory.                                              *
// *                                                                               *
// *   Other members of an "MCObsHandler" object include: a reference to an        *
// *   object of the "MCObsGetHandler" class in the member "m_in_handler"          *
// *   which is used for getting the data from files; a pointer to a               *
// *   "Bootstrapper" object in "Bptr" which handles all bootstrap resampling.     *
// *                                                                               *
// *   Usage:                                                                      *
// *                                                                               *
// *   (1) There are three constructors:                                           *
// *                                                                               *
// *       LaphEnv::MCObsGetHandler& in_handler;                                   *
// *                 // use all data, no rebinning                                 *
// *       MCObsHandler MH(in_handler);                                            *
// *                                                                               *
// *                 // use all data, but rebin by 2                               *
// *       int rebin=2;                                                            *
// *       MCObsHandler MH(in_handler,rebin);                                      *
// *                                                                               *
// *             // discard original measurements 4 and 32 since corrupted         *
// *       set<int> omissions;                                                     *
// *       omissions.insert(4); omissions.insert(32);                              *
// *             // then rebin by 2                                                *
// *       MCObsHandler MH(in_handler, rebin, omissions);                          *
// *                                                                               *
// *                                                                               *
// *   (2) Rebinning the data can also be done as indicated below.  Note that      *
// *   this clears **all** data currently stored in memory.  The rebinning         *
// *   factor is always taken with respect to the original measurements.           *
// *   In other words, using setRebin(2), then setRebin(3) does **not**            *
// *   produce a rebinning by a factor of 6.                                       *
// *                                                                               *
// *       int rebin=3;  MH.setRebin(rebin);                                       *
// *                                                                               *
// *   (3) Including one or more additional omissions (corrupted measurements,     *
// *   for example) can also be done as shown below.  Omissions are always         *
// *   referenced with respect to the original measurements.  Note that            *
// *   this clears **all** data currently stored in memory.  Once the              *
// *   omissions are all taken into account, rebinning is then done using          *
// *   the remaining measurements.                                                 *
// *                                                                               *
// *       int index=5;  MH.addOmission(index);                                    *
// *       set<int> omissions;                                                     *
// *       omissions.insert(41); omissions.insert(324);                            *
// *       MH.addOmissions(omissions);                                             *
// *                                                                               *
// *   (4) All omissions can be removed as shown below.  This also clears          *
// *   **all** data currently stored in memory.                                    *
// *                                                                               *
// *       MH.clearOmissions();                                                    *
// *                                                                               *
// *   (5) Bootstrapping details are set as shown below.  Any current              *
// *   bootstrap data is erased first.  This must be called before any             *
// *   bootstrapping can be done.                                                  *
// *                                                                               *
// *       int num_resamplings=1024;                                               *
// *       unsigned long bootseed=3215467;  // 0 means set by time                 *
// *       unsigned int bootskip=64;                                               *
// *       bool precompute=false;                                                  *
// *       MH.setBootstrapper(num_resamplings,bootseed,                            *
// *                          bootskip,precompute);                                *
// *                                                                               *
// *   (6) Information about current parameters can be obtained as below.          *
// *                                                                               *
// *       MH.getNumberOfMeasurements();                                           *
// *       MH.getNumberOfBins();                                                   *
// *       MH.getNumberOfBootstrapResamplings();                                   *
// *       MH.getBootstrapper();  // returns const reference                       *
// *                                                                               *
// *   (7) Erasing all data corresponding to a particular "MCObsInfo"              *
// *   or all data can be accomplished as shown below.                             *
// *                                                                               *
// *       MCObsInfo obskey;                                                       *
// *       MH.eraseData(obskey);      // erase simple or nonsimple                 *
// *       MH.eraseSamplings(obskey);  // erase only nonsimple                     *
// *       MH.clearData();         // clears all simple and nonsimple              *
// *       MH.clearSamplings();    // clears only nonsimple                        *
// *                                                                               *
// *   (8) Getting and querying simple observable data is accomplished as          *
// *   shown below.  "getBins" first checks to see if the requested                *
// *   observable is already in memory, and if so, a const reference to it         *
// *   is returned, but if not, then "m_in_handler" is called upon to read         *
// *   the data from files and store the results in memory, throwing an            *
// *   exception if it cannot be found.  "getBin" can be used for a single         *
// *   bin, but note that **all** bins are read into memory if "m_in_handler"      *
// *   is called upon.  "queryBins" returns true if the requested observable       *
// *   is already in memory or can be obtained from file; if not already           *
// *   in memory, the data is **not** read into memory.                            *
// *                                                                               *
// *       MCObsInfo obskey;                                                       *
// *       const Vector<double>& binsref=MH.getBins(obskey);                       *
// *       int bin_index=76;                                                       *
// *       double binval=MH.getBin(obskey,bin_index);                              *
// *       bool avail=MH.queryBins(obskey);                                        *
// *                                                                               *
// *   (9) Putting simple observable data into memory is accomplished as           *
// *   shown below.  This is usually done for simple observables derived           *
// *   from the simple data in files, such as rotated correlators or               *
// *   other linear combinations of the correlators or vevs.                       *
// *                                                                               *
// *       Vector<double> newvalues(nbins); <-- computed somehow                   *
// *       MH.putBins(obskey,newvalues);                                           *
// *                                                                               *
// *                                                                               *
// *                                                                               *
// *   (10) Resampling can be done using either the jackknife or the               *
// *   bootstrap method.  Use the subroutines below to set the current             *
// *   method to use or to query which method is currently being used.             *
// *   Note that if you request a method different from the current                *
// *   method, the data currently stored from that resampling method               *
// *   is still retained.  Also, if you call "setToBootstrapMode()" and            *
// *   "setBootstrapper" has not already been called to set up the                 *
// *   bootstrapping parameters, an exception is thrown.                           *
// *                                                                               *
// *       MH.setToJackknifeMode();                                                *
// *       MH.setToBootstrapMode();                                                *
// *       SamplingMode inmode=Jackknife;                                          *
// *       MH.setSamplingMode(inmode);                                             *
// *                                                                               *
// *       MH.isJackknifeMode();                                                   *
// *       MH.isBootstrapMode();                                                   *
// *                                                                               *
// *                                                                               *
// *    (11) Iterating over the resamplings is often needed, either for            *
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
// *    (12) Getting expected values from the entire ensemble for an observable    *
// *    "obskey" or from a particular resampling.  For a simple observable, these  *
// *    routines call "getBins".  For a nonsimple observable, it is expected that  *
// *    these expected values were previously computed and "put" in memory, and    *
// *    these routines simple retrieve the values, or an exception is thrown       *
// *    (exception: for a VEV-subtracted correlator matrix element at one time,    *
// *    these routines call "getBins" for the correlator and well as the VEVs).    *
// *                                                                               *
// *         // get mean or expected value using entire ensemble                   *
// *       double fullmean=MH.getFullSampleValue(obskey);                          *
// *                                                                               *
// *         // get observable value using current resampling                      *
// *       double thismean=MH.getCurrentSamplingValue(obskey);                     *
// *                                                                               *
// *         // query if all samplings are available (include **full** estimate)   *
// *       bool flag=MH.queryAllSamplings(obskey);                                 *
// *                                                                               *
// *    (13) The expected value of a nonsimple observable from a particular        *
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
// *    (14) Covariances for two observables from entire ensemble or the           *
// *    current resampling.  These are often needed for correlated chi-square      *
// *    fitting.  If both observables are simple, the standard covariance          *
// *    calculation is used; if both observables are vev-subtracted correlators    *
// *    at one time slice, then a single-point jackknife is used to compute        *
// *    the covariance.  For all other situations, the standard jackknife or       *
// *    bootstrap covariance is computed for "getFullSampleCovariance", but        *
// *    an exception is thrown for "getCurrentSamplingCovariance".                 *
// *    The square root of the covariance for obskey1==obskey2 is the standard     *
// *    deviation.                                                                 *
// *                                                                               *
// *       double fullcov12=MH.getFullSampleCovariance(obskey1,obskey2);           *
// *       double thiscov12=MH.getCurrentSamplingCovariance(obskey1,obskey2);      *
// *       double stddev=MH.getStandardDeviation(obskey);                          *
// *                                                                               *
// *    (15) Autocorrelation of a simple observable for a particular               *
// *    Markov "time" separation.                                                  *
// *                                                                               *
// *       uint markovtime=3;                                                      *
// *       double autocorr=MH.getAutoCorrelation(obskey,markovtime);               *
// *                                                                               *
// *    (16) Estimates for all resamplings (excluding **full** estimates). Can be  *
// *    used regardless of the current resampling mode.   Results are also stored  *
// *    in memory.                                                                 *
// *                                                                               *
// *       Vector<double> Vsamps=MH.getJackknifeSamplingValues(obskey);            *
// *       std::vector<double>& vsamps;                                            *
// *       MH.getJackknifeSamplingValues(obskey,vsamps);                           *
// *                                                                               *
// *       Vector<double> Vsamps=MH.getBootstrapSamplingValues(obskey);            *
// *       std::vector<double>& vsamps;                                            *
// *       MH.getBootstrapSamplingValues(obskey,vsamps);                           *
// *                                                                               *
// *       SamplingMode mode=Jackknife;                                            *
// *       Vector<double> Vsamps=MH.getSamplingValues(obskey,mode);                *
// *       std::vector<double>& vsamps;                                            *
// *       MH.getSamplingValues(obskey,samp,mode);                                 *
// *                                                                               *
// *     To include the full estimates in samp[0], use                             *
// *                                                                               *
// *       MH.getFullAndSamplingValues(obskey,samp,mode);                          *
// *                                                                               *
// *                                                                               *
// *    (17) J-point jackknife errors (estimates removing J neighboring bins).     *
// *    Can be used regardless of the current resampling mode.  Results are        *
// *    **not** stored in memory.                                                  *
// *                                                                               *
// *       uint jacksize=4;                                                        *
// *       double err=MH.getJackKnifeError(obskey,jacksize);                       *
// *                                                                               *
// *    (18) Statistical analysis:  These routines first call the routines         *
// *    in (16) above, then return the information as an object of the             *
// *    class "MCEstimate".  See the file "mc_estimate.h" for further              *
// *    information of this class.                                                 *
// *                                                                               *
// *       MCEstimate est=MH.getEstimate(obskey);                                  *
// *       SamplingMode mode=Bootstrap;                                            *
// *       MCEstimate est=MH.getEstimate(obskey,mode);                             *
// *       MCEstimate est=MH.getJackknifeEstimate(obskey);                         *
// *       MCEstimate est=MH.getBootstrapEstimate(obskey);                         *
// *                                                                               *
// *    (19) Input/output of all samplings (including full estimates):             *
// *    To write samplings (bootstrap or jackknife) that are already in memory     *
// *    (all samplings include **full* esimates must be available) to file, use    *
// *                                                                               *
// *       set<MCObsInfo> obskeys; ....                                            *
// *       string filename;                                                        *
// *       SamplingMode mode=Bootstrap;                                            *
// *       XMLHandler xmlout;   (for output)                                       *
// *       bool overwrite=false;  (file overwriting)                               *
// *       MH.writeSamplingValuesToFile(obskeys,filename,mode,xmlout,overwrite);   *
// *                                                                               *
// *    If "filename" does not exist, it will be created.  If "filename" exists    *
// *    and "overwrite" is true, the old file will be destroyed and completely     *
// *    overwritten.  If "filename" exists but "overwrite" is false, then the      *
// *    header is checked for consistency and these samplings are added to the     *
// *    file, as long as the observables do not already exist in the file.         *
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
// *    (20) Input/output of bins of computed values of simple observables,        *
// *    such as rotated correlators.  To write bins that are already in memory     *
// *    to file, use                                                               *
// *                                                                               *
// *       set<MCObsInfo> obskeys; ....                                            *
// *       string filename;                                                        *
// *       XMLHandler xmlout;   (for output)                                       *
// *       bool overwrite=false;  (file overwriting)                               *
// *       MH.writeBinsToFile(obskeys,filename,xmlout,overwrite);                  *
// *                                                                               *
// *    If "filename" does not exist, it will be created.  If "filename" exists    *
// *    and "overwrite" is true, the old file will be destroyed and completely     *
// *    overwritten.  If "filename" exists but "overwrite" is false, then the      *
// *    header is checked for consistency and these bins are added to the          *
// *    file, as long as the observables do not already exist in the file.         *
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
// *    The format of the files is that for an IOMap.  The header for the          *
// *    file contains the information                                              *
// *                                                                               *
// *        <SigmondBinsFile>                                                      *
// *           <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo>         *
// *           <NumberOfMeasurements>551</NumberOfMeasurements>                    *
// *           <NumberOfBins>274</NumberOfBins>                                    *
// *           <TweakEnsemble>    (optional)                                       *
// *              <Rebin>2</Rebin>                                                 *
// *              <OmitConfig>0</OmitConfig>                                       *
// *              <OmitConfig>1</OmitConfig>                                       *
// *              <OmitConfig>78</OmitConfig>                                      *
// *           </TweakEnsemble>                                                    *
// *        </SigmondBinsFile>                                                     *
// *                                                                               *
// *    The record keys are unsigned integers, and the values contain              *
// *    the MCObservable XML string, followed by the vector of samplings.          *
// *    An IOMap requires the keys to all have the same size in bytes, so          *
// *    an MCObsInfo cannot be used as a key since its encoding can be             *
// *    different sizes.                                                           *
// *                                                                               *
// *                                                                               *
// *********************************************************************************


class MCObsHandler
{

   LaphEnv::MCObsGetHandler& m_in_handler;  // the input handler class
   unsigned int m_nmeasures;                // total number of measurements in ensemble
   unsigned int m_rebin;                    // rebinning factor (number of meas in each bin)
   std::set<unsigned int> m_omit;           // set of measurements to omit
   unsigned int m_nbins;                    // number of bins
   Bootstrapper* Bptr;                      // pointer to the boot strapper

   std::map<MCObsInfo,Vector<double> > m_obs_simple;        // contains the bins of simple observables
   std::map<std::pair<MCObsInfo,uint>,double> m_jacksamples;   // for storing resamplings
   std::map<std::pair<MCObsInfo,uint>,double> m_bootsamples;   // for storing resamplings

   SamplingMode m_curr_sampling_mode;   //  0 = Jackknife, 1 = Bootstrap
   uint m_curr_sampling_index;          //  0 = full sample, 1..m_max are boot/jack samplings
   uint m_curr_sampling_max;
   std::map<std::pair<MCObsInfo,uint>,double> *m_curr_samples;

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

   MCObsHandler(LaphEnv::MCObsGetHandler& in_handler);

   MCObsHandler(LaphEnv::MCObsGetHandler& in_handler, int rebin);

   MCObsHandler(LaphEnv::MCObsGetHandler& in_handler, int rebin, 
                const std::set<int> omissions);

   void setRebin(int rebin);

   void addOmission(int index);

   void addOmissions(std::set<int> indices);

   void clearOmissions();

   void setBootstrapper(int num_resamplings=1024, unsigned long bootseed=0, // 0 means set by time
                        unsigned int bootskip=64, bool precompute=false);

   ~MCObsHandler();

   unsigned int getNumberOfMeasurements() const {return m_nmeasures;}

   unsigned int getNumberOfBins() const {return m_nbins;}

   unsigned int getNumberOfBootstrapResamplings() const 
    {return (Bptr)?Bptr->getNumberOfResamplings():0;}

   unsigned int getRebinFactor() const {return m_rebin;}

   const std::set<unsigned int>& getOmissions() const {return m_omit;}

   uint getLatticeTimeExtent() const;

   uint getLatticeXExtent() const;

   uint getLatticeYExtent() const;

   uint getLatticeZExtent() const;


   Bootstrapper& getBootstrapper() const;

   void clearData();

   void eraseData(const MCObsInfo& obskey);

   void clearSamplings();

   void eraseSamplings(const MCObsInfo& obskey);


   const Vector<double>& getBins(const MCObsInfo& obskey);

   void getBins(const MCObsInfo& obskey, std::vector<double>& values);

   bool queryBins(const MCObsInfo& obskey);

   bool queryBinsInMemory(const MCObsInfo& obskey);

   double getBin(const MCObsInfo& obskey, int bin_index);

   void putBins(const MCObsInfo& obskey, const Vector<double>& values);


   void setToJackknifeMode();

   void setToBootstrapMode();

   void setSamplingMode(SamplingMode inmode);

   bool isJackknifeMode() const;

   bool isBootstrapMode() const;

   SamplingMode getCurrentSamplingMode() const {return m_curr_sampling_mode;}
 
   MCObsHandler& setSamplingBegin();

   MCObsHandler& begin();

   MCObsHandler& setSamplingNext();

   MCObsHandler& operator++();

   bool isSamplingEnd() const;

   bool end() const;


   bool queryAllSamplings(const MCObsInfo& obskey);

   double getFullSampleValue(const MCObsInfo& obskey);

   double getFullSampleValue(const MCObsInfo& obskey, SamplingMode mode);

   double getCurrentSamplingValue(const MCObsInfo& obskey);

   void putCurrentSamplingValue(const MCObsInfo& obskey, 
                                double value, bool overwrite=true);

   double getFullSampleCovariance(const MCObsInfo& obskey1,   // for simple observables
                                  const MCObsInfo& obskey2);
   double getCurrentSamplingCovariance(const MCObsInfo& obskey1,   // for simple observables
                                       const MCObsInfo& obskey2);


   Vector<double> getJackknifeSamplingValues(const MCObsInfo& obskey);

   void getJackknifeSamplingValues(const MCObsInfo& obskey, std::vector<double>& samp);

   Vector<double> getBootstrapSamplingValues(const MCObsInfo& obskey);

   void getBootstrapSamplingValues(const MCObsInfo& obskey, std::vector<double>& samp);

   Vector<double> getSamplingValues(const MCObsInfo& obskey, SamplingMode mode);

   void getSamplingValues(const MCObsInfo& obskey, std::vector<double>& samp,
                          SamplingMode mode);

   void getFullAndSamplingValues(const MCObsInfo& obskey, 
                                 std::vector<double>& samples, SamplingMode mode);

   double getStandardDeviation(const MCObsInfo& obskey);

   double getJackKnifeError(const MCObsInfo& obskey, uint jacksize);

   MCEstimate getEstimate(const MCObsInfo& obskey);

   MCEstimate getEstimate(const MCObsInfo& obskey, SamplingMode inmode);

   MCEstimate getJackknifeEstimate(const MCObsInfo& obskey);

   MCEstimate getBootstrapEstimate(const MCObsInfo& obskey);


   double getAutoCorrelation(const MCObsInfo& obskey, uint markovtime);

             // read all samplings from file and put into memory (second version
             // only reads those records matching the MCObsInfo objects in "obskeys")

   void readSamplingValuesFromFile(const std::string& filename, 
                                   SamplingMode mode, XMLHandler& xmlout);

   void readSamplingValuesFromFile(const std::set<MCObsInfo>& obskeys, 
                                   const std::string& filename,
                                   SamplingMode mode, XMLHandler& xmlout);

             // write from memory into file (only if all samplings available)

   void writeSamplingValuesToFile(const std::set<MCObsInfo>& obskeys, 
                                  const std::string& filename,
                                  SamplingMode mode, XMLHandler& xmlout,
                                  bool overwrite=false);

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
                        XMLHandler& xmlout, bool overwrite=false);


 private:

   void assert_simple(const MCObsInfo& obskey, const std::string& name);
   void assert_simple(const MCObsInfo& obskey1, const MCObsInfo& obskey2, const std::string& name);
   const Vector<double>& get_data(const MCObsInfo& obskey);
   bool query_data(const MCObsInfo& obskey);
   void reset_nbins();
   void read_data(const MCObsInfo& obskey, Vector<Scalar>& data);  // "Scalar" could be complex

   double calc_mean(const Vector<double>& dvals);
   double calc_mean(const Vector<double>& dvals, double full, uint binremove);
   double calc_mean(const Vector<double>& dvals, const Vector<uint>& indmap);
   double calc_jack_mean(const MCObsInfo& obskey, uint jindex);
   double calc_boot_mean(const MCObsInfo& obskey, uint bindex);

   double get_full_sample_value(const MCObsInfo& obskey,
                 std::map<std::pair<MCObsInfo,uint>,double> *samp_ptr);

   double calc_simple_covariance(const Vector<double>& dvals1, const Vector<double>& dvals2);
   double calc_simple_covariance(const Vector<double>& dvals1, const Vector<double>& dvals2,
                                 uint binremove);
   double calc_simple_covariance(const Vector<double>& dvals1, const Vector<double>& dvals2,
                                 const Vector<uint>& indmap);

   void jack_analyze(double full, std::vector<double>& sampvals, MCEstimate& result);
   void boot_analyze(double full, std::vector<double>& sampvals, MCEstimate& result);

#ifdef COMPLEXNUMBERS
   void get_jack_samples(const Vector<double>& corrbins, bool realpart,
                         const Vector<double>& snkvevrebins, const Vector<double>& snkvevimbins,
                         const Vector<double>& srcvevrebins, const Vector<double>& srcvevimbins,
                         std::vector<double>& sampvals);
   void get_jack_samples(const Vector<double>& corrbins, bool realpart,
                         const Vector<double>& snkvevrebins, const Vector<double>& snkvevimbins,
                         const Vector<double>& srcvevrebins, const Vector<double>& srcvevimbins,
                         const Vector<uint>& indmap, std::vector<double>& sampvals);
   void set_up_corvev_bins(const MCObsInfo& obskey, const Vector<double>* &corrbins, bool& realpart,
                        const Vector<double>* &snkvevrebins, const Vector<double>* &snkvevimbins,
                        const Vector<double>* &srcvevrebins, const Vector<double>* &srcvevimbins);
#else
   void get_jack_samples(const Vector<double>& corrbins,
                         const Vector<double>& snkvevrebins, const Vector<double>& srcvevrebins, 
                         std::vector<double>& sampvals);
   void get_jack_samples(const Vector<double>& corrbins,
                         const Vector<double>& snkvevrebins, const Vector<double>& srcvevrebins, 
                         const Vector<uint>& indmap, std::vector<double>& sampvals);
   void set_up_corvev_bins(const MCObsInfo& obskey, const Vector<double>* &corrbins, 
                        const Vector<double>* &snkvevrebins, const Vector<double>* &srcvevrebins);
#endif
   double jack_covariance(const std::vector<double>& sampvals1, 
                          const std::vector<double>& sampvals2);

   void output_samplings_header(XMLHandler& xmlout, SamplingMode mode);
   bool check_samplings_header(XMLHandler& xmlin, SamplingMode mode);
   void output_bins_header(XMLHandler& xmlout);
   bool check_bins_header(XMLHandler& xmlin);

   void getJackknifeSamplingValues(const MCObsInfo& obskey, double* samp, bool exclude_full);
   void getBootstrapSamplingValues(const MCObsInfo& obskey, double* samp, bool exclude_full);

 public:

   class DataType
   {
      std::string header;
      std::vector<double> data;
   
    public:

      DataType() {}

      DataType(const DataType& in) 
          : header(in.header), data(in.data) {}

      DataType& operator=(const DataType& in) 
       {header=in.header; data=in.data; return *this;}

      void read(LaphEnv::IOHandler& ioh)
       { ioh.read(header); ioh.read(data); }

      void write(LaphEnv::IOHandler& ioh) const
       { ioh.write(header); ioh.write(data); }

      size_t numbytes(LaphEnv::IOHandler& ioh) const
       { return ioh.numbytes(header)+ioh.numbytes(data); }

      friend class MCObsHandler;

   };


};


inline void read(LaphEnv::IOHandler& ioh, MCObsHandler::DataType& samps)
 { samps.read(ioh); }

inline void write(LaphEnv::IOHandler& ioh, const MCObsHandler::DataType& samps)
 { samps.write(ioh); }

inline size_t numbytes(LaphEnv::IOHandler& ioh, const MCObsHandler::DataType& samps)
 { return samps.numbytes(ioh); }


// *************************************************************
#endif

#include "mcobs_handler.h"
#include <algorithm>
#include <limits>

using namespace LaphEnv;
using namespace std;

// *********************************************************************************
// *                                                                               *
// *   The class "MCObsHandler", a container class for Monte Carlo data and        *
// *   non-simple observables, is defined in this file.  This is one of the        *
// *   most important and central classes in "sigmond".  An object of this         *
// *   class is used to get the Monte Carlo measurements, store them in            *
// *   memory, and perform statistical analysis on the data.  When calculations    *
// *   are carried out on the data or fits are performed, objects of this          *
// *   class store them so that bootstrap and jackknife errors can be              *
// *   subsequently computed.                                                      *
// *                                                                               *
// *********************************************************************************



MCObsHandler::MCObsHandler(MCObsGetHandler& in_handler, bool bootprecompute)  
   : m_in_handler(in_handler), Bptr(0), 
     m_curr_sampling_mode(in_handler.getDefaultSamplingMode()), m_curr_sampling_index(0),
     m_curr_sampling_max(in_handler.getNumberOfDefaultResamplings()), 
     m_curr_samples(in_handler.getSamplingInfo().isJackknifeMode() ? &m_jacksamples : &m_bootsamples),
     m_curr_covmat_sampling_mode(in_handler.getDefaultSamplingMode())
{
 /*
 if (getNumberOfMeasurements()<24){
    //cout << "Number of measurements is too small"<<endl;
    throw(std::invalid_argument("Invalid MCObsHandler: Number of measurements is too small"));}*/
 if (in_handler.getSamplingInfo().isBootstrapMode()){
    const MCSamplingInfo& samp=getSamplingInfo();
    Bptr=new Bootstrapper(getNumberOfBins(),samp.getNumberOfReSamplings(getBinsInfo()),
                          samp.getRNGSeed(),samp.getSkipValue(),bootprecompute);}
 m_is_weighted=m_in_handler.getEnsembleInfo().isWeighted();
}


MCObsHandler::~MCObsHandler()
{
 if (Bptr) delete Bptr;
}


unsigned int MCObsHandler::getNumberOfMeasurements() const
{
 return m_in_handler.getNumberOfMeasurements();
}


unsigned int MCObsHandler::getNumberOfBins() const
{
 return m_in_handler.getBinsInfo().getNumberOfBins();
}


SamplingMode MCObsHandler::getDefaultSamplingMode() const
{
 return m_in_handler.getSamplingInfo().getSamplingMode();
}


const MCBinsInfo& MCObsHandler::getBinsInfo() const
{
 return m_in_handler.getBinsInfo();
}


const MCSamplingInfo& MCObsHandler::getSamplingInfo() const
{
 return m_in_handler.getSamplingInfo();
}


unsigned int MCObsHandler::getNumberOfBootstrapResamplings() const 
{
 return (Bptr)?Bptr->getNumberOfResamplings():0;
}


unsigned int MCObsHandler::getRebinFactor() const
{
 return m_in_handler.getBinsInfo().getRebinFactor();
}


const std::set<unsigned int>& MCObsHandler::getOmissions() const
{
 return m_in_handler.getBinsInfo().getOmissions();
}


uint MCObsHandler::getLatticeTimeExtent() const
{
 return m_in_handler.getBinsInfo().getLatticeTimeExtent();
}


uint MCObsHandler::getLatticeXExtent() const
{
 return m_in_handler.getBinsInfo().getLatticeXExtent();
}


uint MCObsHandler::getLatticeYExtent() const
{
 return m_in_handler.getBinsInfo().getLatticeYExtent();
}


uint MCObsHandler::getLatticeZExtent() const
{
 return m_in_handler.getBinsInfo().getLatticeZExtent();
}



const Bootstrapper& MCObsHandler::getBootstrapper() const
{
 if (Bptr) return *Bptr;
 cout << "Fatal error: requested reference to nonexistence Bootstrapper"<<endl;
 throw(std::invalid_argument("Nonexistent Bootstrapper"));
}


const Vector<uint>& MCObsHandler::getBootstrapperResampling(uint bootindex) const
{
 if (Bptr) return Bptr->getResampling(bootindex);
 cout << "Fatal error: requested reference to nonexistence Bootstrapper"<<endl;
 throw(std::invalid_argument("Nonexistent Bootstrapper"));
}


void MCObsHandler::clearData()
{
 m_obs_simple.clear();
 clearSamplings();
}


void MCObsHandler::eraseData(const MCObsInfo& obskey)
{
 m_obs_simple.erase(obskey);
 eraseSamplings(obskey);
}


void MCObsHandler::clearSamplings()
{
 m_jacksamples.clear();
 m_bootsamples.clear();
}

void MCObsHandler::eraseSamplings(const MCObsInfo& obskey)
{
 m_jacksamples.erase(obskey);
 m_bootsamples.erase(obskey);
}


void MCObsHandler::getFileMap(XMLHandler& xmlout) const
{
 m_in_handler.getFileMap(xmlout);
}

uint MCObsHandler::getBinsInMemorySize() const
{
 return m_obs_simple.size();
}


uint MCObsHandler::getJackknifeSamplingsInMemorySize() const
{
 return m_jacksamples.size();
}


uint MCObsHandler::getBootstrapSamplingsInMemorySize() const
{
 return m_bootsamples.size();
}


uint MCObsHandler::getCurrentModeSamplingsInMemorySize() const
{
 return m_curr_samples->size();
}


uint MCObsHandler::getSamplingsInMemorySize(SamplingMode mode) const
{
 return (mode==Jackknife) ? m_jacksamples.size() : m_bootsamples.size();
}


// ************************************************************


void MCObsHandler::assert_simple(const MCObsInfo& obskey, const string& name)
{
 if (obskey.isNonSimple())
    throw(std::invalid_argument(string("Error in ")+name
          +string(":  NonSimple observable where simple required")));
}


const RVector& MCObsHandler::get_bins(const MCObsInfo& obskey)
{
 assert_simple(obskey,"getBins");
 map<MCObsInfo,RVector>::const_iterator dt=m_obs_simple.find(obskey);
 if (dt!=m_obs_simple.end()){ return (dt->second);}

 if (!obskey.hasNoRelatedFlip()){
    MCObsInfo tkey(obskey.getTimeFlipped());
    map<MCObsInfo,RVector>::const_iterator dt=m_obs_simple.find(tkey);
    if (dt!=m_obs_simple.end()){
       RVector buf(dt->second);
#ifdef COMPLEXNUMBERS
       if (tkey.isImaginaryPart())
          for (uint k=0;k<buf.size();++k) 
             buf[k]=-buf[k];
#endif
       return putBins(obskey,buf);}}

 if (m_obs_simple.size()>262144){
    cout << "Data exhaustion!"<<endl;
    m_obs_simple.clear();
    exit(1);}

#ifdef COMPLEXNUMBERS
 try{
    RVector bins_re, bins_im;
    m_in_handler.getBinsComplex(obskey,bins_re,bins_im);
    pair< map<MCObsInfo,RVector>::iterator,bool> ret,ret2;
    MCObsInfo key2(obskey);
    if (obskey.isRealPart()){
       ret=m_obs_simple.insert(make_pair(obskey,bins_re));
       key2.setToImaginaryPart();
       ret2=m_obs_simple.insert(make_pair(key2,bins_im));}
    else{
       ret=m_obs_simple.insert(make_pair(obskey,bins_im));
       key2.setToRealPart();
       ret2=m_obs_simple.insert(make_pair(key2,bins_re));}
    if ((!ret.second)||(!ret2.second)){
       //cout << "Error inserting in MCObsHandler map"<<endl;
       throw(std::invalid_argument("Insertion error: Error inserting in MCObsHandler map"));}
    return (ret.first)->second;}

#else
 try{
    RVector bins;
    m_in_handler.getBins(obskey,bins);
    pair<map<MCObsInfo,RVector>::iterator,bool> ret;
    ret=m_obs_simple.insert(make_pair(obskey,bins));
    if (!ret.second){
       //cout << "Error inserting in MCObsHandler map"<<endl;
       throw(std::invalid_argument("Insertion error: Error inserting in MCObsHandler map"));}
    return (ret.first)->second;}
#endif

 catch(const std::exception& errmsg){
    throw(std::runtime_error(string("Error in MCObsHandler::getBins: ")+errmsg.what()));}
}

  //  Returns bins for "obskey".  If not found, then it deduces the bins
  //  from those of the time-flipped correlator.

const RVector& MCObsHandler::getBins(const MCObsInfo& obskey)
{
 return get_bins(obskey);
}


void MCObsHandler::getBins(const MCObsInfo& obskey, RVector& values)
{
 try{
    values=getBins(obskey).c_vector();}
 catch(const std::exception& errmsg){
    //cout << "Error in MCObsHandler::getBins: "<<errmsg<<endl;
    throw(std::runtime_error(string("Error in MCObsHandler::getBins: ")+errmsg.what()));}
}


bool MCObsHandler::queryBins(const MCObsInfo& obskey)
{
 if (obskey.isNonSimple()) return false;
 map<MCObsInfo,RVector>::const_iterator dt=m_obs_simple.find(obskey);
 if (dt!=m_obs_simple.end()) return true;
 return m_in_handler.queryBins(obskey);
}



bool MCObsHandler::queryBinsInMemory(const MCObsInfo& obskey)
{
 if (obskey.isNonSimple()) return false;
 map<MCObsInfo,RVector>::const_iterator dt=m_obs_simple.find(obskey);
 return (dt!=m_obs_simple.end());
}



double MCObsHandler::getBin(const MCObsInfo& obskey, int bin_index)
{
 try{
    const RVector& d=getBins(obskey);
    return d[bin_index];}
 catch(const std::exception& errmsg){
    //cout << "Error in MCObsHandler::getBin: "<<errmsg<<endl;
    throw(std::runtime_error(string("Error in MCObsHandler::getBin: ")+errmsg.what()));}
}


const RVector& MCObsHandler::putBins(const MCObsInfo& obskey, const RVector& values)
{
 assert_simple(obskey,"putBins");
 if (values.size()!=getNumberOfBins())
    throw(std::invalid_argument("Invalid Vector size in putBins"));
// if (m_in_handler.queryBins(obskey))
//    throw(std::invalid_argument("Cannot put Bins for data contained in the files"));
 m_obs_simple.erase(obskey);
 try{
    pair<map<MCObsInfo,RVector>::iterator,bool> flag=m_obs_simple.insert(make_pair(obskey,values));
    if (!(flag.second)) throw(std::runtime_error("Could not putBins"));
    return (flag.first)->second;}
 catch(std::exception& xp){
    cout << "FATAL: could not do putBins"<<endl; exit(1);}
}


// *****************************************************************


void MCObsHandler::setToJackknifeMode()
{
 if (m_curr_sampling_mode!=Jackknife){
    m_curr_sampling_mode=Jackknife;
    m_curr_samples=&m_jacksamples;
    m_curr_sampling_index=0;
    m_curr_sampling_max=getNumberOfBins();}
}


void MCObsHandler::setToBootstrapMode()
{
 if (m_curr_sampling_mode!=Bootstrap){
    if (Bptr==0) throw(std::invalid_argument("Must initialize bootstrapper to go to bootstrap mode"));
    m_curr_sampling_mode=Bootstrap;
    m_curr_samples=&m_bootsamples;
    m_curr_sampling_index=0;
    m_curr_sampling_max=getNumberOfBootstrapResamplings();}
}


void MCObsHandler::setSamplingMode(SamplingMode inmode)
{
 if (inmode==Jackknife)
    setToJackknifeMode();
 else
    setToBootstrapMode();
}


void MCObsHandler::setToDefaultSamplingMode()
{
 setSamplingMode(m_in_handler.getDefaultSamplingMode());
}



bool MCObsHandler::isJackknifeMode() const
{
 return (m_curr_sampling_mode==Jackknife);
}


bool MCObsHandler::isBootstrapMode() const
{
 return (m_curr_sampling_mode==Bootstrap);
}


bool MCObsHandler::isDefaultSamplingMode() const
{
 return (m_curr_sampling_mode==m_in_handler.getDefaultSamplingMode());
}





void MCObsHandler::setCovMatToJackknifeMode()
{
 m_curr_covmat_sampling_mode=Jackknife;
}


void MCObsHandler::setCovMatToBootstrapMode()
{
 if (m_curr_covmat_sampling_mode!=Bootstrap){
    if (Bptr==0) throw(std::invalid_argument("Must initialize bootstrapper to go to bootstrap mode"));
    m_curr_covmat_sampling_mode=Bootstrap;}
}


void MCObsHandler::setCovMatSamplingMode(SamplingMode inmode)
{
 if (inmode==Jackknife)
    setCovMatToJackknifeMode();
 else
    setCovMatToBootstrapMode();
}


void MCObsHandler::setCovMatToDefaultSamplingMode()
{
 setCovMatSamplingMode(m_in_handler.getDefaultSamplingMode());
}



bool MCObsHandler::isCovMatJackknifeMode() const
{
 return (m_curr_covmat_sampling_mode==Jackknife);
}


bool MCObsHandler::isCovMatBootstrapMode() const
{
 return (m_curr_covmat_sampling_mode==Bootstrap);
}



MCObsHandler& MCObsHandler::setSamplingBegin()
{
 m_curr_sampling_index=0;
 return *this;
}


MCObsHandler& MCObsHandler::begin()
{
 m_curr_sampling_index=0;
 return *this;
}


MCObsHandler& MCObsHandler::setSamplingNext()
{
 m_curr_sampling_index++;
 return *this;
}


MCObsHandler& MCObsHandler::operator++()
{
 m_curr_sampling_index++;
 return *this;
}


bool MCObsHandler::isSamplingEnd() const
{
 return (m_curr_sampling_index>m_curr_sampling_max);
}


bool MCObsHandler::end() const
{
 return (m_curr_sampling_index>m_curr_sampling_max);
}



// *****************************************************************


bool MCObsHandler::queryFullAndSamplings(const MCObsInfo& obskey)
{
 map<MCObsInfo,pair<RVector,uint> >::const_iterator dt=m_curr_samples->find(obskey);
 if (dt!=m_curr_samples->end()) 
    if ((dt->second).second==(dt->second).first.size()) return true;
 if (!obskey.hasNoRelatedFlip()){
    MCObsInfo tkey(obskey.getTimeFlipped());
    map<MCObsInfo,pair<RVector,uint> >::const_iterator dt=m_curr_samples->find(tkey);
    if (dt!=m_curr_samples->end()){
       if ((dt->second).second==(dt->second).first.size()) return true;}}
 if (query_from_samplings_file(obskey)) return true;
 return query_samplings_from_bins(obskey);
}


bool MCObsHandler::queryFullAndSamplings(const MCObsInfo& obskey, SamplingMode mode)
{
 SamplingMode keep=getCurrentSamplingMode();
 setSamplingMode(mode);
 bool result=queryFullAndSamplings(obskey);
 setSamplingMode(keep);
 return result;
}


const RVector& MCObsHandler::get_full_and_sampling_values(const MCObsInfo& obskey, 
                      map<MCObsInfo,pair<RVector,uint> > *samp_ptr,
                      SamplingMode mode, bool allow_not_all_available)
{
 map<MCObsInfo,pair<RVector,uint> >::const_iterator dt=samp_ptr->find(obskey);
 if (dt!=samp_ptr->end()){
    if (((dt->second).second==(dt->second).first.size())||(allow_not_all_available))
       return (dt->second).first;}
 if (!obskey.hasNoRelatedFlip()){
    MCObsInfo tkey(obskey.getTimeFlipped());
    map<MCObsInfo,pair<RVector,uint> >::const_iterator dt=samp_ptr->find(tkey);
    if (dt!=samp_ptr->end()){
       if (((dt->second).second==(dt->second).first.size())||(allow_not_all_available)){
          RVector samples(dt->second.first);
          if (obskey.isImaginaryPart()){
             samples*=-1.0;}
          return put_samplings_in_memory(obskey,samples,samp_ptr);}}}
 if (mode==m_in_handler.getDefaultSamplingMode()){
    RVector samples;
    if (m_in_handler.getSamplingsMaybe(obskey,samples)){
       return put_samplings_in_memory(obskey,samples,samp_ptr);}
    const RVector* res=calc_corrsubvev_from_samplings(obskey,samp_ptr);
    if (res!=0){return *res;}}
 try{
    void (MCObsHandler::*simpcalc_ptr)(const RVector&,RVector&)
      =(mode==Jackknife) ? &MCObsHandler::calc_simple_jack_samples
                         : &MCObsHandler::calc_simple_boot_samples;
    return calc_samplings_from_bins(obskey,samp_ptr,simpcalc_ptr);}
 catch(std::exception& xp){
    throw(std::runtime_error(string("getSamplings failed: ")+xp.what()));}
}


double MCObsHandler::get_a_sampling_value(const MCObsInfo& obskey, 
                        std::map<MCObsInfo,std::pair<RVector,uint> > *samp_ptr,
                        SamplingMode mode, uint sampindex)
{
 const RVector& samps=get_full_and_sampling_values(obskey,m_curr_samples,m_curr_sampling_mode,true);
 if (std::isnan(samps[sampindex]))
    throw(std::runtime_error(string("getSampling failed for ")+obskey.str()+string(" for index = ")
          +make_string(sampindex)));
 return samps[sampindex];
}


const RVector* MCObsHandler::get_full_and_sampling_values_maybe(const MCObsInfo& obskey, 
                      map<MCObsInfo,pair<RVector,uint> > *samp_ptr,
                      SamplingMode mode)
{
 map<MCObsInfo,pair<RVector,uint> >::const_iterator dt=samp_ptr->find(obskey);
 if (dt!=samp_ptr->end()){
    return &((dt->second).first);}
 if (!obskey.hasNoRelatedFlip()){
    MCObsInfo tkey(obskey.getTimeFlipped());
    map<MCObsInfo,pair<RVector,uint> >::const_iterator dt=samp_ptr->find(tkey);
    if (dt!=samp_ptr->end()){
       RVector samples(dt->second.first);
       if (obskey.isImaginaryPart()){
          samples*=-1.0;}
       return &put_samplings_in_memory(obskey,samples,samp_ptr);}}
 if (mode==m_in_handler.getDefaultSamplingMode()){
    RVector samples;
    if (m_in_handler.getSamplingsMaybe(obskey,samples)){
       return &(put_samplings_in_memory(obskey,samples,samp_ptr));}
    const RVector* res=calc_corrsubvev_from_samplings(obskey,samp_ptr);
    if (res!=0) return res;}
 try{
    void (MCObsHandler::*simpcalc_ptr)(const RVector&,RVector&)
      =(mode==Jackknife) ? &MCObsHandler::calc_simple_jack_samples
                         : &MCObsHandler::calc_simple_boot_samples;
    return &(calc_samplings_from_bins(obskey,samp_ptr,simpcalc_ptr));}
 catch(std::exception& xp){
    return 0;}
}


bool MCObsHandler::get_a_sampling_value_maybe(const MCObsInfo& obskey, 
                        std::map<MCObsInfo,std::pair<RVector,uint> > *samp_ptr,
                        SamplingMode mode, uint sampindex, double& result)
{
 result=0.0;
 const RVector* samps=get_full_and_sampling_values_maybe(obskey,m_curr_samples,m_curr_sampling_mode);
 if (samps==0) return false;
 if (std::isnan((*samps)[sampindex])) return false;
 result=(*samps)[sampindex];
 return true;
}


const RVector& MCObsHandler::put_samplings_in_memory(const MCObsInfo& obskey,
                      const RVector& samplings,
                      map<MCObsInfo,pair<RVector,uint> > *samp_ptr)
{
 samp_ptr->erase(obskey);
 pair<map<MCObsInfo,pair<RVector,uint> >::iterator,bool> ret;
 ret=samp_ptr->insert(make_pair(obskey,make_pair(samplings,samplings.size())));
 if (ret.second==false){
    throw(std::runtime_error("put samplings into memory failed"));}
 return ((ret.first)->second).first;
}


const RVector& MCObsHandler::getFullAndSamplingValues(const MCObsInfo& obskey, 
                                                      SamplingMode mode)
{
 map<MCObsInfo,pair<RVector,uint> > *samp_ptr
     =(mode==Jackknife) ? &m_jacksamples : &m_bootsamples;
 return get_full_and_sampling_values(obskey,samp_ptr,mode);
}


void MCObsHandler::getFullAndSamplingValues(const MCObsInfo& obskey, 
                          RVector& samples, SamplingMode mode)
{
 samples=getFullAndSamplingValues(obskey,mode);
}


double MCObsHandler::getFullSampleValue(const MCObsInfo& obskey)
{
 return get_a_sampling_value(obskey,m_curr_samples,m_curr_sampling_mode,0);
}


double MCObsHandler::getFullSampleValue(const MCObsInfo& obskey, SamplingMode mode)
{
 map<MCObsInfo,pair<RVector,uint> > *samp_ptr
     =(mode==Jackknife) ? &m_jacksamples : &m_bootsamples;
 return get_a_sampling_value(obskey,samp_ptr,mode,0);
}


double MCObsHandler::getCurrentSamplingValue(const MCObsInfo& obskey)
{
 return get_a_sampling_value(obskey,m_curr_samples,m_curr_sampling_mode,m_curr_sampling_index);
}


bool MCObsHandler::getCurrentSamplingValueMaybe(const MCObsInfo& obskey, double& result)
{
 return get_a_sampling_value_maybe(obskey,m_curr_samples,m_curr_sampling_mode,m_curr_sampling_index,result);
}


void MCObsHandler::putCurrentSamplingValue(const MCObsInfo& obskey, 
                                           double value, bool overwrite)
{
 put_a_sampling_in_memory(obskey,m_curr_sampling_index,value,overwrite, 
                          m_curr_sampling_max,m_curr_samples);
}


void MCObsHandler::put_a_sampling_in_memory(
                        const MCObsInfo& obskey, uint sampling_index, 
                        double value, bool overwrite, uint sampling_max,
                        std::map<MCObsInfo,std::pair<RVector,uint> > *samp_ptr)
{
 if (sampling_index>sampling_max)
    throw(std::invalid_argument("invalid index in put_a_sampling_in_memory"));
 map<MCObsInfo,pair<RVector,uint> >::iterator dt=samp_ptr->find(obskey);
 if (dt!=samp_ptr->end()){
    double& entry=(dt->second.first)[sampling_index];
    if (std::isnan(entry)){
       entry=value; (dt->second.second)++;}
    else if (overwrite){
       entry=value;}
    else{
       throw(std::invalid_argument("cannot putCurrentSamplingValue since no overwrite"));}}
 RVector buffer(sampling_max+1,std::numeric_limits<double>::quiet_NaN());
 buffer[sampling_index]=value;
 samp_ptr->insert(make_pair(obskey,make_pair(buffer,1)));
}


void MCObsHandler::calc_simple_jack_samples(const RVector& bins, RVector& samplings)
{
 uint nbins=bins.size(); 
 samplings.resize(nbins+1);  // 0 = full, 1..nbins are the jackknife samplings
 if (!m_is_weighted){
   double dm=0.0;
   for (uint k=0;k<nbins;k++) 
      dm+=bins[k];
   samplings[0]=dm/double(nbins);
   double rj=1.0/double(nbins-1);
   for (uint k=1;k<=nbins;k++)
      samplings[k]=(dm-bins[k-1])*rj;}
 else{
   const vector<double>& wts=m_in_handler.getWeights();
   double dm=0.0;
   double den=0.0;
   for (uint k=0;k<nbins;k++){
      dm+=wts[k]*bins[k];
      den+=wts[k];}
   samplings[0]=dm/den;
   for (uint k=0;k<nbins;k++)
      samplings[k+1]=(dm-wts[k]*bins[k])/(den-wts[k]);}
}


void MCObsHandler::calc_simple_boot_samples(const RVector& bins, RVector& samplings)
{
 if (Bptr==0)
    throw(std::runtime_error("No Bootstrapper object!"));
 uint nsamps=Bptr->getNumberOfResamplings();
 uint nbins=Bptr->getNumberOfObjects();
 if (nbins!=bins.size())
    throw(std::runtime_error("Mismatch in number of bins in bootstrapper"));
 samplings.resize(nsamps+1);   // 0 = full, 1..nsamps are the bootstrap samplings
 if (!m_is_weighted){
   double dm=0.0;
   for (uint k=0;k<nbins;k++) 
      dm+=bins[k];
   samplings[0]=dm/double(nbins);
   for (uint bootindex=0;bootindex<nsamps;bootindex++){
      const Vector<uint>& indmap=Bptr->getResampling(bootindex);
      dm=0.0;
      for (uint k=0;k<nbins;k++) 
         dm+=bins[indmap[k]];
      samplings[bootindex+1]=dm/double(nbins);}}
 else{
   const vector<double>& wts=m_in_handler.getWeights();
   double dm=0.0;
   double den=0.0;
   for (uint k=0;k<nbins;k++){
      dm+=wts[k]*bins[k];
      den+=wts[k];}
   samplings[0]=dm/den;
   for (uint bootindex=0;bootindex<nsamps;bootindex++){
      const Vector<uint>& indmap=Bptr->getResampling(bootindex);
      dm=den=0.0;
      for (uint k=0;k<nbins;k++){
         uint kk=indmap[k];
         dm+=wts[kk]*bins[kk];
         den+=wts[kk];}
      samplings[bootindex+1]=dm/den;}}
}


#ifdef COMPLEXNUMBERS

void MCObsHandler::calc_corr_subvev(const RVector& corr_re_samples, const RVector& corr_im_samples,
                                    const RVector& snkvev_re_samples, const RVector& snkvev_im_samples,
                                    const RVector& srcvev_re_samples, const RVector& srcvev_im_samples,
                                    RVector& corrsubvev_re_samples, RVector& corrsubvev_im_samples)
{
 uint n=corr_re_samples.size();
 if  ((corr_im_samples.size()!=n)||(snkvev_re_samples.size()!=n)
    ||(snkvev_im_samples.size()!=n)||(srcvev_re_samples.size()!=n)
    ||(srcvev_im_samples.size()!=n))
    throw(std::runtime_error("Number of resamplings do not match in calc_corr_subvev"));
 corrsubvev_re_samples.resize(n);
 corrsubvev_im_samples.resize(n);
 for (uint k=0;k<n;k++){
    corrsubvev_re_samples[k]=corr_re_samples[k]-(snkvev_re_samples[k]*srcvev_re_samples[k]
                            +snkvev_im_samples[k]*srcvev_im_samples[k]);
    corrsubvev_im_samples[k]=corr_im_samples[k]-(snkvev_im_samples[k]*srcvev_re_samples[k]
                            -snkvev_re_samples[k]*srcvev_im_samples[k]);}
}

#else

void MCObsHandler::calc_corr_subvev(const RVector& corr_samples, const RVector& snkvev_samples, 
                                    const RVector& srcvev_samples, RVector& corrsubvev_samples)
{
 uint n=corr_samples.size();
 if  ((snkvev_samples.size()!=n)||(srcvev_samples.size()!=n))
    throw(std::runtime_error("Number of resamplings do not match in calc_corr_subvev"));
 corrsubvev_samples.resize(n);
 for (uint k=0;k<n;k++){
    corrsubvev_samples[k]=corr_samples[k]-snkvev_samples[k]*srcvev_samples[k];}
}

#endif


const RVector& MCObsHandler::calc_samplings_from_bins(const MCObsInfo& obskey,
                      map<MCObsInfo,pair<RVector,uint> > *samp_ptr,
                      void (MCObsHandler::*simpcalc_ptr)(const RVector&,RVector&))
{
 if (obskey.isSimple()){
    const RVector& bins=getBins(obskey);
    RVector samplings;
    (this->*simpcalc_ptr)(bins,samplings);
    return put_samplings_in_memory(obskey,samplings,samp_ptr);}
 else if (obskey.isCorrelatorAtTime()){   // since nonsimple, must have vev subtraction
    bool realpart=obskey.isRealPart();
#ifdef REALNUMBERS
    if (!realpart) throw(std::runtime_error("Invalid"));
#endif
    bool herm=obskey.isHermitianCorrelatorAtTime();
    unsigned int tval=obskey.getCorrelatorTimeIndex();
    OperatorInfo src(obskey.getCorrelatorSourceInfo());
    OperatorInfo snk(obskey.getCorrelatorSinkInfo());
    MCObsInfo corr_re_info(snk,src,tval,herm,RealPart,false);   // no vev subtraction
    const RVector& corr_re_bins=getBins(corr_re_info);
    RVector corr_re_samplings;
    (this->*simpcalc_ptr)(corr_re_bins,corr_re_samplings);
    put_samplings_in_memory(corr_re_info,corr_re_samplings,samp_ptr);
    MCObsInfo src_re_info(src,RealPart);
    MCObsInfo snk_re_info(snk,RealPart);
    const RVector& src_re_bins=getBins(src_re_info);
    const RVector& snk_re_bins=getBins(snk_re_info);
    RVector src_re_samplings, snk_re_samplings;
    (this->*simpcalc_ptr)(src_re_bins,src_re_samplings);
    (this->*simpcalc_ptr)(snk_re_bins,snk_re_samplings);
    put_samplings_in_memory(src_re_info,src_re_samplings,samp_ptr);
    put_samplings_in_memory(snk_re_info,snk_re_samplings,samp_ptr);
    RVector corrsubvev_re_samplings;
    const RVector *ci=0;
    const RVector *cr=0;
#ifdef COMPLEXNUMBERS
    MCObsInfo corr_im_info(snk,src,tval,herm,ImaginaryPart,false);   // no vev subtraction
    const RVector& corr_im_bins=getBins(corr_im_info);
    RVector corr_im_samplings;
    (this->*simpcalc_ptr)(corr_im_bins,corr_im_samplings);
    put_samplings_in_memory(corr_im_info,corr_im_samplings,samp_ptr);
    MCObsInfo src_im_info(src,ImaginaryPart);
    MCObsInfo snk_im_info(snk,ImaginaryPart);
    const RVector& src_im_bins=getBins(src_im_info);
    const RVector& snk_im_bins=getBins(snk_im_info);
    RVector src_im_samplings, snk_im_samplings;
    (this->*simpcalc_ptr)(src_im_bins,src_im_samplings);
    (this->*simpcalc_ptr)(snk_im_bins,snk_im_samplings);
    put_samplings_in_memory(src_im_info,src_im_samplings,samp_ptr);
    put_samplings_in_memory(snk_im_info,snk_im_samplings,samp_ptr);
    RVector corrsubvev_im_samplings;
    calc_corr_subvev(corr_re_samplings,corr_im_samplings,snk_re_samplings,
                     snk_im_samplings,src_re_samplings,src_im_samplings,
                     corrsubvev_re_samplings,corrsubvev_im_samplings);
    MCObsInfo corrsubvev_im_info(snk,src,tval,herm,ImaginaryPart,true);   // no vev subtraction
    ci=&put_samplings_in_memory(corrsubvev_im_info,corrsubvev_im_samplings,samp_ptr);
#else
    calc_corr_subvev(corr_re_samplings,snk_re_samplings,src_re_samplings, 
                     corrsubvev_re_samplings);
#endif
    MCObsInfo corrsubvev_re_info(snk,src,tval,herm,RealPart,true);   // no vev subtraction
    cr=&put_samplings_in_memory(corrsubvev_re_info,corrsubvev_re_samplings,samp_ptr);
    return (realpart) ? *cr : *ci;}
 else
    throw(std::invalid_argument("Unable to get all samplings in MCObsHandler"));
}


const RVector* MCObsHandler::calc_corrsubvev_from_samplings(const MCObsInfo& obskey,
                                         map<MCObsInfo,pair<RVector,uint> > *samp_ptr)
{
 if (!obskey.isCorrelatorAtTime()) return 0;
 if (!(obskey.getCorrelatorAtTimeInfo().subtractVEV())) return 0;
 bool realpart=obskey.isRealPart();
#ifdef REALNUMBERS
 if (!realpart) return 0;
#endif
 try{
 bool herm=obskey.isHermitianCorrelatorAtTime();
 unsigned int tval=obskey.getCorrelatorTimeIndex();
 OperatorInfo src(obskey.getCorrelatorSourceInfo());
 OperatorInfo snk(obskey.getCorrelatorSinkInfo());
 MCObsInfo corr_re_info(snk,src,tval,herm,RealPart,false);   // no vev subtraction
 RVector corr_re_samplings;
 m_in_handler.getSamplings(corr_re_info,corr_re_samplings);
 put_samplings_in_memory(corr_re_info,corr_re_samplings,samp_ptr);
 MCObsInfo src_re_info(src,RealPart);
 MCObsInfo snk_re_info(snk,RealPart);
 RVector src_re_samplings, snk_re_samplings;
 m_in_handler.getSamplings(src_re_info,src_re_samplings);
 m_in_handler.getSamplings(snk_re_info,snk_re_samplings);
 put_samplings_in_memory(src_re_info,src_re_samplings,samp_ptr);
 put_samplings_in_memory(snk_re_info,snk_re_samplings,samp_ptr);
 RVector corrsubvev_re_samplings;
 const RVector *ci=0;
 const RVector *cr=0;
#ifdef COMPLEXNUMBERS
 MCObsInfo corr_im_info(snk,src,tval,herm,ImaginaryPart,false);   // no vev subtraction
 RVector corr_im_samplings;
 m_in_handler.getSamplings(corr_im_info,corr_im_samplings);
 put_samplings_in_memory(corr_im_info,corr_im_samplings,samp_ptr);
 MCObsInfo src_im_info(src,ImaginaryPart);
 MCObsInfo snk_im_info(snk,ImaginaryPart);
 RVector src_im_samplings, snk_im_samplings;
 m_in_handler.getSamplings(src_im_info,src_im_samplings);
 m_in_handler.getSamplings(snk_im_info,snk_im_samplings);
 put_samplings_in_memory(src_im_info,src_im_samplings,samp_ptr);
 put_samplings_in_memory(snk_im_info,snk_im_samplings,samp_ptr);
 RVector corrsubvev_im_samplings;
 calc_corr_subvev(corr_re_samplings,corr_im_samplings,snk_re_samplings,
                  snk_im_samplings,src_re_samplings,src_im_samplings,
                  corrsubvev_re_samplings,corrsubvev_im_samplings);
 MCObsInfo corrsubvev_im_info(snk,src,tval,herm,ImaginaryPart,true);   // no vev subtraction
 ci=&put_samplings_in_memory(corrsubvev_im_info,corrsubvev_im_samplings,samp_ptr);
#else
 calc_corr_subvev(corr_re_samplings,snk_re_samplings,src_re_samplings, 
                  corrsubvev_re_samplings);
#endif
 MCObsInfo corrsubvev_re_info(snk,src,tval,herm,RealPart,true);   // no vev subtraction
 cr=&put_samplings_in_memory(corrsubvev_re_info,corrsubvev_re_samplings,samp_ptr);
 return (realpart) ? cr : ci;}
 catch(const std::exception& xp){}
 return 0;
}



bool MCObsHandler::query_samplings_from_bins(const MCObsInfo& obskey)
{
 if (obskey.isSimple()){
    return queryBins(obskey);}
 else if (obskey.isCorrelatorAtTime()){   // since nonsimple, must have vev subtraction
    bool herm=obskey.isHermitianCorrelatorAtTime();
    unsigned int tval=obskey.getCorrelatorTimeIndex();
    OperatorInfo src(obskey.getCorrelatorSourceInfo());
    OperatorInfo snk(obskey.getCorrelatorSinkInfo());
    ComplexArg arg=(obskey.isRealPart())?RealPart:ImaginaryPart;
    MCObsInfo corr_info(snk,src,tval,herm,arg,false);   // no vev subtraction
    MCObsInfo src_re_info(src,RealPart);
    MCObsInfo snk_re_info(snk,RealPart);
    bool flag=queryBins(corr_info) && queryBins(src_re_info) && queryBins(snk_re_info);
#ifdef COMPLEXNUMBERS
    if (!flag) return false;
    MCObsInfo src_im_info(src,ImaginaryPart);
    MCObsInfo snk_im_info(snk,ImaginaryPart);
    flag=queryBins(src_im_info) && queryBins(snk_im_info);
#endif
    return flag;}
 else
    return false;
}


bool MCObsHandler::query_from_samplings_file(const MCObsInfo& obskey)
{
 if (!isDefaultSamplingMode()) return false;
 if (obskey.isSimple()){
    return m_in_handler.querySamplings(obskey);}
      // even if nonsimple, try the query first
 if (m_in_handler.querySamplings(obskey)) return true;
 if (!(obskey.isCorrelatorAtTime())) return false;
       // last chance: if vev-subtracted correlator, then
       // can still evaluate it from samplings of vevs and non-vev-subtracted
       // correlator
 bool herm=obskey.isHermitianCorrelatorAtTime();
 unsigned int tval=obskey.getCorrelatorTimeIndex();
 OperatorInfo src(obskey.getCorrelatorSourceInfo());
 OperatorInfo snk(obskey.getCorrelatorSinkInfo());
 ComplexArg arg=(obskey.isRealPart())?RealPart:ImaginaryPart;
 MCObsInfo corr_info(snk,src,tval,herm,arg,false);   // no vev subtraction
 MCObsInfo src_re_info(src,RealPart);
 MCObsInfo snk_re_info(snk,RealPart);
 bool flag=   m_in_handler.querySamplings(corr_info) 
           && m_in_handler.querySamplings(src_re_info)
           && m_in_handler.querySamplings(snk_re_info);
#ifdef COMPLEXNUMBERS
 if (!flag) return false;
 MCObsInfo src_im_info(src,ImaginaryPart);
 MCObsInfo snk_im_info(snk,ImaginaryPart);
 flag=  m_in_handler.querySamplings(src_im_info) 
     && m_in_handler.querySamplings(snk_im_info);
#endif
 return flag;
}


RVector MCObsHandler::getJackknifeSamplingValues(const MCObsInfo& obskey)
{
 return copyRVectorWithoutZerothElement(getFullAndSamplingValues(obskey,Jackknife));
}


void MCObsHandler::getJackknifeSamplingValues(const MCObsInfo& obskey, RVector& jacksamples)
{
 jacksamples=copyRVectorWithoutZerothElement(getFullAndSamplingValues(obskey,Jackknife));
}


RVector MCObsHandler::getBootstrapSamplingValues(const MCObsInfo& obskey)
{
 return copyRVectorWithoutZerothElement(getFullAndSamplingValues(obskey,Bootstrap));
}


void MCObsHandler::getBootstrapSamplingValues(const MCObsInfo& obskey, RVector& bootsamples)
{
 bootsamples=copyRVectorWithoutZerothElement(getFullAndSamplingValues(obskey,Bootstrap));
}



RVector MCObsHandler::getSamplingValues(const MCObsInfo& obskey,
                                        SamplingMode mode)
{
 return copyRVectorWithoutZerothElement(getFullAndSamplingValues(obskey,mode));
}


void MCObsHandler::getSamplingValues(const MCObsInfo& obskey, 
                                     RVector& samples, SamplingMode mode)
{
 samples=copyRVectorWithoutZerothElement(getFullAndSamplingValues(obskey,mode));
}



// ***************************************************************************


MCEstimate MCObsHandler::getEstimate(const MCObsInfo& obskey)
{
 if (isJackknifeMode()){
    MCEstimate result(Jackknife);
    const RVector& sampvals=getFullAndSamplingValues(obskey,Jackknife);
    jack_analyze(sampvals,result);
    return result;}
 else{
    MCEstimate result(Bootstrap);
    RVector sampvals;
    getFullAndSamplingValues(obskey,sampvals,Bootstrap); // sampvals will be sorted
    boot_analyze(sampvals,result);
    return result;}
}


MCEstimate MCObsHandler::getEstimate(const MCObsInfo& obskey, SamplingMode inmode)
{
 MCEstimate result(inmode);
 if (inmode==Jackknife){
    const RVector& sampvals=getFullAndSamplingValues(obskey,Jackknife);
    jack_analyze(sampvals,result);
    return result;}
 else{
    RVector sampvals;
    getFullAndSamplingValues(obskey,sampvals,Bootstrap); // sampvals will be sorted
    boot_analyze(sampvals,result);
    return result;}
}


MCEstimate MCObsHandler::getJackknifeEstimate(const MCObsInfo& obskey)
{
 MCEstimate result(Jackknife);
 const RVector& sampvals=getFullAndSamplingValues(obskey,Jackknife);
 jack_analyze(sampvals,result);
 return result;
}


MCEstimate MCObsHandler::getBootstrapEstimate(const MCObsInfo& obskey)
{
 MCEstimate result(Bootstrap);
 RVector sampvals;
 getFullAndSamplingValues(obskey,sampvals,Bootstrap); // sampvals will be sorted
 boot_analyze(sampvals,result);
 return result;
}


// **********************************************************************


double MCObsHandler::getCovariance(const MCObsInfo& obskey1,
                                   const MCObsInfo& obskey2, SamplingMode inmode)
{
 const RVector& sampvals1=getFullAndSamplingValues(obskey1,inmode);
 const RVector& sampvals2=getFullAndSamplingValues(obskey2,inmode);
 if (inmode==Jackknife) 
    return jack_covariance(sampvals1,sampvals2);
 else
    return boot_covariance(sampvals1,sampvals2);
}


double MCObsHandler::getCovariance(const MCObsInfo& obskey1,     // current sampling mode
                                   const MCObsInfo& obskey2)
{
 return getCovariance(obskey1,obskey2,m_curr_covmat_sampling_mode);
}




double MCObsHandler::getStandardDeviation(const MCObsInfo& obskey)
{
 return sqrt(getCovariance(obskey,obskey,Jackknife));
}



double MCObsHandler::getJackKnifeError(const MCObsInfo& obskey, uint jacksize)
{
 uint nbins=getNumberOfBins();
 if ((jacksize<1)||(jacksize>(nbins/12)))
    throw(std::invalid_argument("Invalid jack size in getJackKnifeError"));
 try{
    const RVector& buffer=getBins(obskey);
    uint N=nbins - (nbins % jacksize);   // N is divisible by jacksize
    uint NJ=N/jacksize;
    double zJ=((double)(N-jacksize))/((double) N);
    double sum=buffer[0];
    for (unsigned int k=1;k<N;k++)
       sum+=buffer[k];
    double avg=sum/double(N);
    double cov=0.0;
    uint count=0;
    for (uint jack=0;jack<NJ;jack++){
       double avgJ=sum;
       for (uint i=0;i<jacksize;i++){
          avgJ-=buffer[count++];}
       avgJ/=double(N-jacksize);
       cov+=(avg-avgJ)*(avg-avgJ);}
    cov*=zJ;
    return sqrt(cov);}
 catch(const std::exception& errmsg){
    cout << "Error in MCObsHandler::getJackKnifeError: "<<errmsg.what()<<endl;
    throw(std::invalid_argument("getJackKnifeError failed"));}
}

//          This function returns the autocorrelation
//          function in rho(markovtime).  By definition,
//          rho(0)=1.  Note that the JLQCD definition of the 
//          autocorrelation function is used.
//
//                                  < (O(i)-OA)*(O(i+t)-OB) >
//                 rho(t) =   --------------------------------------
//                            [<(O(i)-OA)**2><(O(i+t)-OB)**2>]^(1/2)
//
//          where OA = sum( O(i), i=0..N-1 ) / N
//                OB = sum( O(i+markovtime), i=0..N-1 ) / N   where N=Nmeas-markovtime
//          <O> means sum( O(i), i=0..N-1 ) / N
//


double MCObsHandler::getAutoCorrelation(const MCObsInfo& obskey, uint markovtime)
{
 assert_simple(obskey,"getAutoCorrelation");
 uint nbins=getNumberOfBins();
 if (markovtime>(nbins/12))
    throw(std::invalid_argument("Markov time too large in getAutoCorrelation"));
 const RVector& buffer=getBins(obskey);
 uint N=nbins-markovtime;
 double ninv=1.0/double(N);

 double OA=0.0; for (uint i=0;i<N;i++) OA+=buffer[i]; OA*=ninv;
 double OB=0.0; for (uint i=0;i<N;i++) OB+=buffer[i+markovtime]; OB*=ninv;
 double W1(0.0),W2(0.0),W3(0.0);
 for (uint i=0;i<N;i++){
    double Wt2=buffer[i]-OA;
    double Wt3=buffer[i+markovtime]-OB;
    W1+=Wt2*Wt3;
    W2+=Wt2*Wt2;
    W3+=Wt3*Wt3;}
 W1*=ninv; W2*=ninv; W3*=ninv;
 return W1/sqrt(W2*W3);
}


    //  compute jackknife error using samplings (sampvals[0] contains full sample)

void MCObsHandler::jack_analyze(const RVector& sampvals, MCEstimate& result)
{
 uint n=sampvals.size()-1;
 double avg=0.0;
 for (uint k=1;k<=n;k++)
    avg+=sampvals[k];
 avg/=double(n);   // should be the same as sampvals[0]
 double var=0.0;
 for (uint k=1;k<=n;k++){
    double x=sampvals[k]-avg;
    var+=x*x;}
 var*=(1.0-1.0/double(n));  // jackknife
 result.jackassign(sampvals[0],avg,sqrt(var));
}

    //  compute jackknife covariance using samplings 
    //  (sampvals1[0], sampvals2[0] contain full samples)

double MCObsHandler::jack_covariance(const RVector& sampvals1, 
                                     const RVector& sampvals2)
{
 uint n=sampvals1.size()-1;
 double avg1=0.0;
 double avg2=0.0;
 for (uint k=1;k<=n;k++){
    avg1+=sampvals1[k];
    avg2+=sampvals2[k];}
 avg1/=double(n);
 avg2/=double(n);
 double jackcov=0.0;
 for (uint k=1;k<=n;k++){
    jackcov+=(sampvals1[k]-avg1)*(sampvals2[k]-avg2);}
 jackcov*=(1.0-1.0/double(n));  // jackknife
 return jackcov;
}


    //  compute bootstrap covariance using samplings 
    //  (sampvals1[0], sampvals2[0] contain full samples)

double MCObsHandler::boot_covariance(const RVector& sampvals1, 
                                     const RVector& sampvals2)
{
 uint n=sampvals1.size()-1;
 double avg1=0.0;
 double avg2=0.0;
 for (uint k=1;k<=n;k++){
    avg1+=sampvals1[k];
    avg2+=sampvals2[k];}
 avg1/=double(n);
 avg2/=double(n);
 double bootcov=0.0;
 for (uint k=1;k<=n;k++){
    bootcov+=(sampvals1[k]-avg1)*(sampvals2[k]-avg2);}
 bootcov/=double(n-1);  // bootstrap
 return bootcov;
}


   // This takes an array of "Nboot" bootstrap samples
   // and returns their average (ans.boot_avg), the
   // value which 84% of the samples lie above (ans.boot_low),
   // the value which 84% of the samples lie below (ans.boot_upp),
   // and the value which 50% of the samples lie above and
   // 50% lie below (ans.boot_med).   CAUTION: this changes "sampvals"

void MCObsHandler::boot_analyze(RVector& sampvals, MCEstimate& result)
{
 uint nb=sampvals.size()-1;
 double avg=0.0;
 for (uint k=1;k<=nb;k++)
    avg+=sampvals[k];
 avg/=double(nb);
 double var=0.0;
 for (uint k=1;k<=nb;k++){
    double x=sampvals[k]-avg;
    var+=x*x;}
 var/=double(nb-1);  // bootstrap  

 const double conf_level=0.68;  // one standard deviation
 const double lowconf=0.5*(1.0-conf_level);
 const double uppconf=0.5*(1.0+conf_level);
 const double eps=1e-10;
 std::sort(sampvals.c_vector_ref().begin()+1,sampvals.c_vector_ref().end());

    // get the low and high confidence bounds and the median

 uint i=(unsigned int) floor( lowconf*nb+eps); 
 double low=0.5*(sampvals[i]+sampvals[i+1]);
 i=(unsigned int) floor( uppconf*nb+eps); 
 double upp=0.5*(sampvals[i]+sampvals[i+1]);
 i=(unsigned int) floor( 0.5*nb+eps); 
 double med=0.5*(sampvals[i]+sampvals[i+1]);

 result.bootassign(sampvals[0],avg,sqrt(var),low,med,upp);
}

// *************************************************************


             // read all samplings from file and put into memory

void MCObsHandler::readSamplingValuesFromFile(const string& filename, 
                                              XMLHandler& xmlout)
{
 xmlout.set_root("ReadSamplingsFromFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 try{
    m_in_handler.connectSamplingsFile(filename);
    xmlout.put_child("Status","Success");}
 catch(std::exception& xp){
    xmlout.put_child("Error",xp.what());}
}


             // read all samplings from file and put into memory (this version
             // only reads those records matching the MCObsInfo objects in "obskeys")

void MCObsHandler::readSamplingValuesFromFile(const set<MCObsInfo>& obskeys, 
                                              const string& filename,
                                              XMLHandler& xmlout)
{
 xmlout.set_root("ReadSamplingsFromFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 try{
    m_in_handler.connectSamplingsFile(filename,obskeys);
    xmlout.put_child("Status","Success");
    for (set<MCObsInfo>::const_iterator it=obskeys.begin();it!=obskeys.end();it++){
       XMLHandler xmlobs; it->output(xmlobs); xmlout.put_child(xmlobs);}}
 catch(std::exception& xp){
    xmlout.put_child("Error",xp.what());}
}

      // Write from memory into file (only if all samplings available).
      // If "filename" does not exist, it will be created.  If "filename"
      // exists and "overwrite" is true, the old file will be destroyed
      // and completely overwritten.  If "filename" exists but "overwrite"
      // is false, then the header is checked for consistency and these
      // samplings are added to the file, as long as the observables do
      // not already exist in the file.

void MCObsHandler::writeSamplingValuesToFile(const set<MCObsInfo>& obskeys, 
                                             const string& filename,
                                             XMLHandler& xmlout,
                                             WriteMode wmode)
{
 xmlout.set_root("WriteSamplingsToFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 xmlout.put_child("FileName",fname);
 try{
    SamplingsPutHandler SP(m_in_handler.getBinsInfo(),m_in_handler.getSamplingInfo(),
                           filename,wmode, m_in_handler.useCheckSums());
    for (set<MCObsInfo>::const_iterator it=obskeys.begin();it!=obskeys.end();it++){
       XMLHandler xmlo; it->output(xmlo);
       XMLHandler xmle("Write");
       xmle.put_child(xmlo);
       if (!queryFullAndSamplings(*it)){
          xmle.put_child("Error","Not all samplings are in memory");}
       else{
       try{
          const RVector& buffer=getFullAndSamplingValues(*it,
                          m_in_handler.getSamplingInfo().getSamplingMode());
          SP.putData(*it,buffer);
          xmle.put_child("Success");}
       catch(std::exception& xp){
          xmle.put_child("Error",string(xp.what()));}}
       xmlout.put_child(xmle);}}
 catch(const std::exception& errmsg){
    xmlout.put_child("Error",string(errmsg.what()));}
}

// ************************************************************************

             // read all bins from file (such as from rotated correlators)
             // and put into memory

void MCObsHandler::readBinsFromFile(const string& filename, XMLHandler& xmlout)
{
 xmlout.set_root("ReadBinsFromFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 try{
    m_in_handler.connectBinsFile(filename);
    xmlout.put_child("Status","Success");}
 catch(std::exception& xp){
    xmlout.put_child("Error",xp.what());}
}

             // read all bins from file and put into memory (this version
             // only reads those records matching the MCObsInfo objects in "obskeys")

void MCObsHandler::readBinsFromFile(const set<MCObsInfo>& obskeys, 
                                    const string& filename, XMLHandler& xmlout)
{
 xmlout.set_root("ReadBinsFromFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 try{
    m_in_handler.connectBinsFile(filename,obskeys);
    xmlout.put_child("Status","Success");
    for (set<MCObsInfo>::const_iterator it=obskeys.begin();it!=obskeys.end();it++){
       XMLHandler xmlobs; it->output(xmlobs); xmlout.put_child(xmlobs);}}
 catch(std::exception& xp){
    xmlout.put_child("Error",xp.what());}
}

      // Write bins from memory into file (simple observables only).
      // If "filename" does not exist, it will be created.  If "filename"
      // exists and "overwrite" is true, the old file will be destroyed
      // and completely overwritten.  If "filename" exists but "overwrite"
      // is false, then the header is checked for consistency and these
      // bins are added to the file, as long as the observables do
      // not already exist in the file.

void MCObsHandler::writeBinsToFile(const set<MCObsInfo>& obskeys, 
                                   const string& filename,
                                   XMLHandler& xmlout, WriteMode wmode)
{
 xmlout.set_root("WriteBinsToFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 xmlout.put_child("FileName",fname);
 xmlout.put_child("NumberObservablesToWrite",make_string(obskeys.size()));
 try{
    BinsPutHandler BP(m_in_handler.getBinsInfo(),filename, 
                      wmode, m_in_handler.useCheckSums());
    uint success=0;
    for (set<MCObsInfo>::const_iterator it=obskeys.begin();it!=obskeys.end();it++){
       XMLHandler xmlo; it->output(xmlo);
       XMLHandler xmle("Write");
       xmle.put_child(xmlo);
       if (!queryBins(*it)){
          xmle.put_child("Error","Not all bins are not in memory");}
       else{
       try{
          const RVector& buffer=getBins(*it);
          BP.putData(*it,buffer);
          xmle.put_child("Success"); success++;}
       catch(std::exception& xp){
          xmle.put_child("Error",string(xp.what()));}}
       xmlout.put_child(xmle);}
       xmlout.put_child("NumberObservablesSuccessfullyWrittenToFile",make_string(success));}
 catch(const std::exception& errmsg){
    xmlout.put_child("Error",string(errmsg.what()));}
}


// ************************************************************************

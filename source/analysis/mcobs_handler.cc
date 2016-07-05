#include "mcobs_handler.h"
#include <algorithm>

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



MCObsHandler::MCObsHandler(MCObsGetHandler& in_handler)  
   : m_in_handler(in_handler), m_nmeasures(in_handler.getNumberOfMeasurements()), 
     m_rebin(1), m_nbins(m_nmeasures), Bptr(0), 
     m_curr_sampling_mode(Jackknife), m_curr_sampling_index(0),
     m_curr_sampling_max(m_nbins), m_curr_samples(&m_jacksamples)
{
 if (m_nmeasures<24){
    cout << "Number of measurements is too small"<<endl;
    throw(std::invalid_argument("Invalid MCObsHandler: Number of measurements is too small"));}
}


MCObsHandler::MCObsHandler(MCObsGetHandler& in_handler, int rebin)  
   : m_in_handler(in_handler), m_nmeasures(in_handler.getNumberOfMeasurements()), 
     m_rebin(1), m_nbins(m_nmeasures), Bptr(0), 
     m_curr_sampling_mode(Jackknife), m_curr_sampling_index(0),
     m_curr_sampling_max(m_nbins), m_curr_samples(&m_jacksamples)
{
 if (m_nmeasures<24){
    cout << "Number of measurements is too small"<<endl;
    throw(std::invalid_argument("Invalid MCObsHandler: Number of measurements is too small"));}
 setRebin(rebin);
}


MCObsHandler::MCObsHandler(MCObsGetHandler& in_handler, int rebin, 
                           const set<int> omissions)  
   : m_in_handler(in_handler), m_nmeasures(in_handler.getNumberOfMeasurements()), 
     m_rebin(1), m_nbins(m_nmeasures), Bptr(0),
     m_curr_sampling_mode(Jackknife), m_curr_sampling_index(0),
     m_curr_sampling_max(m_nbins), m_curr_samples(&m_jacksamples)
{
 if (m_nmeasures<24){
    cout << "Number of measurements is too small"<<endl;
    throw(std::invalid_argument("Invalid MCObsHandler: Number of measurements is too small"));}
 addOmissions(omissions);
 setRebin(rebin);
}


void MCObsHandler::setRebin(int rebin)
{
 if ((rebin>0)&&(rebin<int(m_nmeasures/24))){
    m_rebin=rebin;
    reset_nbins();}
 else{
    cout << "Invalid value for rebin"<<endl;
    throw(std::invalid_argument("Invalid rebin value"));}
}


void MCObsHandler::addOmission(int index)
{
 if ((index>=0)&&(index<int(m_nmeasures))){
    m_omit.insert(index);
    reset_nbins();}
}


void MCObsHandler::addOmissions(set<int> indices)
{
 for (set<int>::const_iterator it=indices.begin();it!=indices.end();it++){
    if ((*it>=0)&&(*it<int(m_nmeasures))) m_omit.insert(*it);}
 reset_nbins();
}


void MCObsHandler::clearOmissions()
{
 m_omit.clear();
 reset_nbins();
}


void MCObsHandler::setBootstrapper(int num_resamplings, unsigned long bootseed,
                                   unsigned int bootskip, bool precompute)
{
 if (Bptr){
    delete Bptr;
    m_bootsamples.clear();}
 Bptr=new Bootstrapper(m_nbins,num_resamplings,bootseed,bootskip,precompute);
 if (m_curr_sampling_mode==Bootstrap){
    m_curr_samples=&m_bootsamples;
    m_curr_sampling_index=0;
    m_curr_sampling_max=Bptr->getNumberOfResamplings();}
}


void MCObsHandler::reset_nbins()
{
 m_obs_simple.clear();
 m_jacksamples.clear();
 m_bootsamples.clear();
 unsigned int new_nbins=(m_nmeasures-m_omit.size())/m_rebin;
 if (new_nbins!=m_nbins){
    m_nbins=new_nbins;
    if (Bptr){
       Bootstrapper *newboot=new Bootstrapper(m_nbins,Bptr->getNumberOfResamplings(), 
                 Bptr->getRNGSeed(),Bptr->getSkipValue(),Bptr->isPrecomputeMode());
       delete Bptr;
       Bptr=newboot;}
    if (m_curr_sampling_mode==Bootstrap){
       m_curr_samples=&m_bootsamples;
       m_curr_sampling_index=0;
       m_curr_sampling_max=Bptr->getNumberOfResamplings();}
    else{
       m_curr_samples=&m_jacksamples;
       m_curr_sampling_index=0;
       m_curr_sampling_max=m_nbins;}}
}
 

MCObsHandler::~MCObsHandler()
{
 if (Bptr) delete Bptr;
}


Bootstrapper& MCObsHandler::getBootstrapper() const
{
 if (Bptr) return *Bptr;
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
 for (uint j=0;j<=m_nbins;j++)
    m_jacksamples.erase(make_pair(obskey,j));
 for (uint j=0;j<=getNumberOfBootstrapResamplings();j++)
    m_bootsamples.erase(make_pair(obskey,j));
}




uint MCObsHandler::getLatticeTimeExtent() const
{
 return m_in_handler.getEnsembleInfo().getLatticeTimeExtent();
}

uint MCObsHandler::getLatticeXExtent() const
{
 return m_in_handler.getEnsembleInfo().getLatticeXExtent();
}

uint MCObsHandler::getLatticeYExtent() const
{
 return m_in_handler.getEnsembleInfo().getLatticeYExtent();
}

uint MCObsHandler::getLatticeZExtent() const
{
 return m_in_handler.getEnsembleInfo().getLatticeZExtent();
}




void MCObsHandler::read_data(const MCObsInfo& obskey, Vector<Scalar>& result)
{
 try{
   result.resize(m_nbins);
   if ((m_rebin==1)&&(m_omit.empty())){
      for (unsigned int k=0;k<m_nbins;k++){
         m_in_handler.getData(obskey,k,result[k]);}}
   else if ((m_rebin>1)&&(m_omit.empty())){
      unsigned int count=0;
      double r=1.0/double(m_rebin);
      Scalar buffer,buffer2;
      for (unsigned int k=0;k<m_nbins;k++){
         m_in_handler.getData(obskey,count++,buffer);
         for (unsigned int j=1;j<m_rebin;j++){
            m_in_handler.getData(obskey,count++,buffer2);
            buffer+=buffer2;}
         result[k]=buffer*r;}}
   else if ((m_rebin==1)&&(!m_omit.empty())){
      set<unsigned int>::const_iterator om=m_omit.begin();
      unsigned int count=0;
      while ((om!=m_omit.end())&&(count==*om)){om++; count++;}
      for (unsigned int k=0;k<m_nbins;k++){
         m_in_handler.getData(obskey,count,result[k]);
         count++; while ((om!=m_omit.end())&&(count==*om)){om++; count++;}}}
   else{
      double r=1.0/double(m_rebin);
      Scalar buffer,buffer2;
      set<unsigned int>::const_iterator om=m_omit.begin();
      unsigned int count=0;
      while ((om!=m_omit.end())&&(count==*om)){om++; count++;}
      for (unsigned int k=0;k<m_nbins;k++){
         m_in_handler.getData(obskey,count,buffer);
         count++; while ((om!=m_omit.end())&&(count==*om)){om++; count++;}
         for (unsigned int j=1;j<m_rebin;j++){
            m_in_handler.getData(obskey,count,buffer2);
            count++; while ((om!=m_omit.end())&&(count==*om)){om++; count++;}
            buffer+=buffer2;}
         result[k]=buffer*r;}}}
 catch(const std::exception& errmsg){
    //cout << errmsg <<endl;
    throw(std::invalid_argument((string("read_data failed  ")+string(errmsg.what())).c_str()));}
}


bool MCObsHandler::query_data(const MCObsInfo& obskey)
{
 map<MCObsInfo,Vector<double> >::const_iterator dt=m_obs_simple.find(obskey);
 if (dt!=m_obs_simple.end()) return true;

 if ((m_rebin==1)&&(m_omit.empty())){
    for (unsigned int k=0;k<m_nbins;k++){
       if (!m_in_handler.queryData(obskey,k)) return false;}}
 else if ((m_rebin>1)&&(m_omit.empty())){
    unsigned int count=0;
    for (unsigned int k=0;k<m_nbins;k++){
       if (!m_in_handler.queryData(obskey,count++)) return false;
       for (unsigned int j=1;j<m_rebin;j++)
          if (!m_in_handler.queryData(obskey,count++)) return false;}}
 else if ((m_rebin==1)&&(!m_omit.empty())){
    set<unsigned int>::const_iterator om=m_omit.begin();
    unsigned int count=0;
    while ((om!=m_omit.end())&&(count==*om)){om++; count++;}
    for (unsigned int k=0;k<m_nbins;k++){
       if (!m_in_handler.queryData(obskey,count)) return false;
       count++; while ((om!=m_omit.end())&&(count==*om)){om++; count++;}}}
 else{
    set<unsigned int>::const_iterator om=m_omit.begin();
    unsigned int count=0;
    while ((om!=m_omit.end())&&(count==*om)){om++; count++;}
    for (unsigned int k=0;k<m_nbins;k++){
       if (!m_in_handler.queryData(obskey,count)) return false;
       count++; while ((om!=m_omit.end())&&(count==*om)){om++; count++;}
       for (unsigned int j=1;j<m_rebin;j++){
          if (!m_in_handler.queryData(obskey,count)) return false;
          count++; while ((om!=m_omit.end())&&(count==*om)){om++; count++;}}}}
 return true;
}


#ifdef COMPLEXNUMBERS

    //  If numbers are stored as complex in the data files, then might as
    //  well read and store both the real and imaginary parts, since most
    //  likely, both will eventually be needed.


const Vector<double>& MCObsHandler::get_data(const MCObsInfo& obskey)
{
 map<MCObsInfo,Vector<double> >::const_iterator dt=m_obs_simple.find(obskey);
 if (dt!=m_obs_simple.end()) return (dt->second);

 if (m_obs_simple.size()>65536){
    cout << "Data exhaustion!"<<endl;
    m_obs_simple.clear();
    throw(std::invalid_argument("Data exhaustion"));}

 try{
   Vector<Scalar> result(m_nbins);
   read_data(obskey,result);
   uint n=result.size();
   Vector<double> result_re(n);
   Vector<double> result_im(n);
   for (uint k=0;k<n;k++){
      result_re[k]=realpart(result[k]);
      result_im[k]=imaginarypart(result[k]);}
   pair< map<MCObsInfo,Vector<double> >::iterator,bool> ret,ret2;
   MCObsInfo key2(obskey);
   if (obskey.isRealPart()){
      ret=m_obs_simple.insert(make_pair(obskey,result_re));
      key2.setToImaginaryPart();
      ret2=m_obs_simple.insert(make_pair(key2,result_im));}
   else{
      ret=m_obs_simple.insert(make_pair(obskey,result_im));
      key2.setToRealPart();
      ret2=m_obs_simple.insert(make_pair(key2,result_re));}
   if ((!ret.second)||(!ret2.second)){
       cout << "Error inserting in MCObsHandler map"<<endl;
      throw(std::invalid_argument("Insertion error: Error inserting in MCObsHandler map"));}
   return (ret.first)->second;}
 catch(const std::exception& errmsg){
    //cout << errmsg <<endl;
    throw(std::invalid_argument((string("get_data failed")+string(errmsg.what())).c_str()));}
}


#else


const Vector<double>& MCObsHandler::get_data(const MCObsInfo& obskey)
{
 map<MCObsInfo,Vector<double> >::const_iterator dt=m_obs_simple.find(obskey);
 if (dt!=m_obs_simple.end()) return (dt->second);

 if (m_obs_simple.size()>65536){
    cout << "Data exhaustion!"<<endl;
    m_obs_simple.clear();
    throw(std::invalid_argument("Data exhaustion"));}

 try{
   Vector<double> result(m_nbins);
   read_data(obskey,result);
   pair<map<MCObsInfo,Vector<double> >::iterator,bool> ret;
   ret=m_obs_simple.insert(pair<MCObsInfo,Vector<double> >(obskey,result));
   if (!ret.second){
       cout << "Error inserting in MCObsHandler map"<<endl;
      throw(std::invalid_argument("Insertion error: Error inserting in MCObsHandler map"));}
   return (ret.first)->second;}
 catch(const std::exception& errmsg){
    //cout << errmsg <<endl;
    throw(std::invalid_argument((string("get_data failed: ")+string(errmsg.what())).c_str()));}
}

#endif


void MCObsHandler::assert_simple(const MCObsInfo& obskey, const string& name)
{
 if (obskey.isNonSimple())
    throw(std::invalid_argument((string("Error in ")+name
          +string(":  NonSimple observable where simple required")).c_str()));
}


void MCObsHandler::assert_simple(const MCObsInfo& obskey1, const MCObsInfo& obskey2, 
                                 const string& name)
{
 if ((obskey1.isNonSimple())||(obskey2.isNonSimple()))
    throw(std::invalid_argument((string("Error in ")+name
           +string(":  NonSimple observable(s) where simple required")).c_str()));
}


const Vector<double>& MCObsHandler::getBins(const MCObsInfo& obskey)
{
 assert_simple(obskey,"getBins");
 try{
    return get_data(obskey);}
 catch(const std::exception& errmsg){
    cout << "Error in MCObsHandler::getBins: "<<errmsg.what()<<endl;
    throw;}
}


void MCObsHandler::getBins(const MCObsInfo& obskey, vector<double>& values)
{
 assert_simple(obskey,"getBins");
 try{
    values=get_data(obskey).c_vector();}
 catch(const std::exception& errmsg){
    //cout << "Error in MCObsHandler::getBins: "<<errmsg<<endl;
    throw;}
}



bool MCObsHandler::queryBins(const MCObsInfo& obskey)
{
 assert_simple(obskey,"queryBins");
 return query_data(obskey);
}



bool MCObsHandler::queryBinsInMemory(const MCObsInfo& obskey)
{
 assert_simple(obskey,"queryBins");
 map<MCObsInfo,Vector<double> >::const_iterator dt=m_obs_simple.find(obskey);
 return (dt!=m_obs_simple.end());
}



double MCObsHandler::getBin(const MCObsInfo& obskey, int bin_index)
{
 assert_simple(obskey,"getBin");
 try{
    const Vector<double>& d=get_data(obskey);
    return d[bin_index];}
 catch(const std::exception& errmsg){
    //cout << "Error in MCObsHandler::getBin: "<<errmsg<<endl;
    throw;}
}


void MCObsHandler::putBins(const MCObsInfo& obskey, const Vector<double>& values)
{
 assert_simple(obskey,"putBins");
 if (values.size()!=m_nbins)
    throw(std::invalid_argument("Invalid Vector size in putBins"));
// set<unsigned int>::const_iterator om=m_omit.begin();
// unsigned int count=0;
// while ((om!=m_omit.end())&&(count==*om)){om++; count++;}
// for (unsigned int k=0;k<(m_nbins*m_rebin);k++){
//    if (m_in_handler.queryData(obskey,count))
//       throw(std::invalid_argument("Cannot put Bins for data contained in the files"));
//    count++; while ((om!=m_omit.end())&&(count==*om)){om++; count++;}}
 m_obs_simple.erase(obskey);
 m_obs_simple.insert(make_pair(obskey,values));
}



void MCObsHandler::setToJackknifeMode()
{
 if (m_curr_sampling_mode!=Jackknife){
    m_curr_sampling_mode=Jackknife;
    m_curr_samples=&m_jacksamples;
    m_curr_sampling_index=0;
    m_curr_sampling_max=m_nbins;}
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


bool MCObsHandler::isJackknifeMode() const
{
 return (m_curr_sampling_mode==Jackknife);
}


bool MCObsHandler::isBootstrapMode() const
{
 return (m_curr_sampling_mode==Bootstrap);
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



bool MCObsHandler::queryAllSamplings(const MCObsInfo& obskey)
{
 bool flag=true;
 for (uint index=0;index<=m_curr_sampling_max;index++){
    if (m_curr_samples->find(make_pair(obskey,index))==m_curr_samples->end())
       flag=false; break;}
 if (flag) return true;
 if (obskey.isNonSimple()) return false;
 return query_data(obskey);
}


double MCObsHandler::getFullSampleValue(const MCObsInfo& obskey)
{
 return get_full_sample_value(obskey,m_curr_samples);
}


double MCObsHandler::getFullSampleValue(const MCObsInfo& obskey, SamplingMode mode)
{
 if (mode==Jackknife)
    return get_full_sample_value(obskey,&m_jacksamples);
 else
    return get_full_sample_value(obskey,&m_bootsamples);
}



double MCObsHandler::get_full_sample_value(const MCObsInfo& obskey,
                 std::map<std::pair<MCObsInfo,uint>,double> *samp_ptr)
{
 map<pair<MCObsInfo,uint>,double>::iterator dt
            =samp_ptr->find(make_pair(obskey,0));
 if (dt!=samp_ptr->end()) 
    return (dt->second);
 if (obskey.isSimple()){
    const Vector<double>& buffer=getBins(obskey);
    double mean=calc_mean(buffer);
    samp_ptr->insert(make_pair(make_pair(obskey,0),mean));
    return mean;}
 else if (obskey.isCorrelatorAtTime()){   // since nonsimple, must have vev subtraction
    bool herm=obskey.isHermitianCorrelatorAtTime();
    unsigned int tval=obskey.getCorrelatorTimeIndex();
    OperatorInfo src(obskey.getCorrelatorSourceInfo());
    OperatorInfo snk(obskey.getCorrelatorSinkInfo());
    double vev=0.0;
    bool realpart=obskey.isRealPart();
    ComplexArg arg=(realpart)?RealPart:ImaginaryPart;
    MCObsInfo corr(snk,src,tval,herm,arg,false);   // no vev subtraction
    MCObsInfo src_re_info(src,RealPart);
    MCObsInfo snk_re_info(snk,RealPart);
    double src_re=get_full_sample_value(src_re_info,samp_ptr);
    double snk_re=get_full_sample_value(snk_re_info,samp_ptr);
#ifdef COMPLEXNUMBERS
    MCObsInfo src_im_info(src,ImaginaryPart);
    MCObsInfo snk_im_info(snk,ImaginaryPart);
    double src_im=get_full_sample_value(src_im_info,samp_ptr);
    double snk_im=get_full_sample_value(snk_im_info,samp_ptr);
    if (realpart) vev=snk_re*src_re+snk_im*src_im;
    else vev=snk_im*src_re-snk_re*src_im;
#else
    vev=snk_re*src_re;
#endif
    double corrval=get_full_sample_value(corr,samp_ptr)-vev;
    m_curr_samples->insert(make_pair(make_pair(obskey,0),corrval));
    return corrval;}
 else
    throw(std::invalid_argument("Unable to getFullSampleValue in MCObsHandler"));
 return 0.0;
}


double MCObsHandler::getCurrentSamplingValue(const MCObsInfo& obskey)
{
 if (m_curr_sampling_index==0)
    return getFullSampleValue(obskey);
 if ((obskey.isSimple())||(obskey.isNonStandard())){
    if (m_curr_sampling_mode==Jackknife)
       return calc_jack_mean(obskey,m_curr_sampling_index);
    else
       return calc_boot_mean(obskey,m_curr_sampling_index);}
 else if (obskey.isCorrelatorAtTime()){   // since nonsimple, must have vev subtraction
    bool herm=obskey.isHermitianCorrelatorAtTime();
    unsigned int tval=obskey.getCorrelatorTimeIndex();
    OperatorInfo src(obskey.getCorrelatorSourceInfo());
    OperatorInfo snk(obskey.getCorrelatorSinkInfo());
    double vev=0.0;
    bool realpart=obskey.isRealPart();
    ComplexArg arg=(realpart)?RealPart:ImaginaryPart;
    MCObsInfo corr(snk,src,tval,herm,arg,false);   // no vev subtraction
    MCObsInfo src_re_info(src,RealPart);
    MCObsInfo snk_re_info(snk,RealPart);
    double src_re=getCurrentSamplingValue(src_re_info);
    double snk_re=getCurrentSamplingValue(snk_re_info);
#ifdef COMPLEXNUMBERS
    MCObsInfo src_im_info(src,ImaginaryPart);
    MCObsInfo snk_im_info(snk,ImaginaryPart);
    double src_im=getCurrentSamplingValue(src_im_info);
    double snk_im=getCurrentSamplingValue(snk_im_info);
    if (realpart) vev=snk_re*src_re+snk_im*src_im;
    else vev=snk_im*src_re-snk_re*src_im;
#else
    vev=snk_re*src_re;
#endif
    double corrval=getCurrentSamplingValue(corr)-vev;
    m_curr_samples->insert(make_pair(make_pair(obskey,m_curr_sampling_index),corrval));
    return corrval;}
 else
    throw(std::invalid_argument("Unable to getFullSampleValue in MCObsHandler"));
 return 0.0;
}


void MCObsHandler::putCurrentSamplingValue(const MCObsInfo& obskey, 
                                           double value, bool overwrite)
{
 if (!overwrite){
    map<pair<MCObsInfo,uint>,double>::iterator dt
               =m_curr_samples->find(make_pair(obskey,m_curr_sampling_index));
    if (dt!=m_curr_samples->end()) 
       throw(std::invalid_argument("cannot putCurrentSamplingValue since no overwrite"));}
 m_curr_samples->insert(make_pair(make_pair(obskey,m_curr_sampling_index),value));
}




double MCObsHandler::getFullSampleCovariance(const MCObsInfo& obskey1,
                                             const MCObsInfo& obskey2)
{
 if ((obskey1.isSimple())&&(obskey2.isSimple())){
    const Vector<double>& dvals1=getBins(obskey1);
    const Vector<double>& dvals2=getBins(obskey2);
    return calc_simple_covariance(dvals1,dvals2);}
 else if ((obskey1.isCorrelatorAtTime())&&(obskey2.isCorrelatorAtTime())){
    vector<double> sampvals1, sampvals2;
#ifdef COMPLEXNUMBERS
    {bool realpart;
     const Vector<double>* corrbins=0;
     const Vector<double>* snkvevrebins=0; 
     const Vector<double>* snkvevimbins=0;
     const Vector<double>* srcvevrebins=0; 
     const Vector<double>* srcvevimbins=0;
     set_up_corvev_bins(obskey1,corrbins,realpart,snkvevrebins,snkvevimbins,
                        srcvevrebins,srcvevimbins);
     get_jack_samples(*corrbins,realpart,*snkvevrebins,*snkvevimbins,
                      *srcvevrebins,*srcvevimbins,sampvals1);}
    {bool realpart;
     const Vector<double>* corrbins=0;
     const Vector<double>* snkvevrebins=0; 
     const Vector<double>* snkvevimbins=0;
     const Vector<double>* srcvevrebins=0; 
     const Vector<double>* srcvevimbins=0;
     set_up_corvev_bins(obskey2,corrbins,realpart,snkvevrebins,snkvevimbins,
                        srcvevrebins,srcvevimbins);
     get_jack_samples(*corrbins,realpart,*snkvevrebins,*snkvevimbins,
                      *srcvevrebins,*srcvevimbins,sampvals2);}
#else
    {const Vector<double>* corrbins=0;
     const Vector<double>* snkvevrebins=0; 
     const Vector<double>* srcvevrebins=0; 
     set_up_corvev_bins(obskey1,corrbins,snkvevrebins,srcvevrebins);
     get_jack_samples(*corrbins,*snkvevrebins,*srcvevrebins,sampvals1);}
    {const Vector<double>* corrbins=0;
     const Vector<double>* snkvevrebins=0; 
     const Vector<double>* srcvevrebins=0; 
     set_up_corvev_bins(obskey2,corrbins,snkvevrebins,srcvevrebins);
     get_jack_samples(*corrbins,*snkvevrebins,*srcvevrebins,sampvals2);}
#endif
    return jack_covariance(sampvals1,sampvals2);}
 else{
    throw(std::invalid_argument("Could not compute FullSampleCovariance"));}
 return 0.0;}




double MCObsHandler::getCurrentSamplingCovariance(const MCObsInfo& obskey1,
                                                  const MCObsInfo& obskey2)
{
 if (m_curr_sampling_index==0)
    return getFullSampleCovariance(obskey1,obskey2);
 if ((obskey1.isSimple())&&(obskey2.isSimple())){
    const Vector<double>& dvals1=getBins(obskey1);
    const Vector<double>& dvals2=getBins(obskey2);
    double val;
    if (m_curr_sampling_mode==Jackknife){
       val=calc_simple_covariance(dvals1,dvals2,m_curr_sampling_index-1);}
    else{
       const Vector<uint>& indexmapper=Bptr->getResampling(m_curr_sampling_index-1);
       val=calc_simple_covariance(dvals1,dvals2,indexmapper);}
    return val;}
 else if ((obskey1.isCorrelatorAtTime())&&(obskey2.isCorrelatorAtTime())){
    if (!isBootstrapMode())
       throw(std::invalid_argument("getCurrentSamplingCovariance with VEVs requires Bootstrap mode"));
    vector<double> sampvals1, sampvals2;
    const Vector<uint>& indexmapper=Bptr->getResampling(m_curr_sampling_index-1);
#ifdef COMPLEXNUMBERS
    {bool realpart;
     const Vector<double>* corrbins=0;
     const Vector<double>* snkvevrebins=0; 
     const Vector<double>* snkvevimbins=0;
     const Vector<double>* srcvevrebins=0; 
     const Vector<double>* srcvevimbins=0;
     set_up_corvev_bins(obskey1,corrbins,realpart,snkvevrebins,snkvevimbins,
                        srcvevrebins,srcvevimbins);
     get_jack_samples(*corrbins,realpart,*snkvevrebins,*snkvevimbins,
                      *srcvevrebins,*srcvevimbins,indexmapper,sampvals1);}
    {bool realpart;
     const Vector<double>* corrbins=0;
     const Vector<double>* snkvevrebins=0; 
     const Vector<double>* snkvevimbins=0;
     const Vector<double>* srcvevrebins=0; 
     const Vector<double>* srcvevimbins=0;
     set_up_corvev_bins(obskey2,corrbins,realpart,snkvevrebins,snkvevimbins,
                        srcvevrebins,srcvevimbins);
     get_jack_samples(*corrbins,realpart,*snkvevrebins,*snkvevimbins,
                      *srcvevrebins,*srcvevimbins,indexmapper,sampvals2);}
#else
    {const Vector<double>* corrbins=0;
     const Vector<double>* snkvevrebins=0; 
     const Vector<double>* srcvevrebins=0; 
     set_up_corvev_bins(obskey1,corrbins,snkvevrebins,srcvevrebins);
     get_jack_samples(*corrbins,*snkvevrebins,*srcvevrebins,indexmapper,sampvals1);}
    {const Vector<double>* corrbins=0;
     const Vector<double>* snkvevrebins=0; 
     const Vector<double>* srcvevrebins=0; 
     set_up_corvev_bins(obskey2,corrbins,snkvevrebins,srcvevrebins);
     get_jack_samples(*corrbins,*snkvevrebins,*srcvevrebins,indexmapper,sampvals2);}
#endif
    return jack_covariance(sampvals1,sampvals2);}
 else{
    throw(std::invalid_argument("Could not compute CurrentSamplingCovariance"));}
 return 0.0;
}



 


Vector<double> MCObsHandler::getJackknifeSamplingValues(const MCObsInfo& obskey)
{
 Vector<double> jacksamples(m_nbins);
 getJackknifeSamplingValues(obskey,&jacksamples[0],true);
 return jacksamples;
}


void MCObsHandler::getJackknifeSamplingValues(const MCObsInfo& obskey, vector<double>& jacksamples)
{
 jacksamples.resize(m_nbins);
 getJackknifeSamplingValues(obskey,&jacksamples[0],true);
}


void MCObsHandler::getJackknifeSamplingValues(const MCObsInfo& obskey, 
                                              double* jacksamples, bool exclude_full)
{
 SamplingMode currmode=getCurrentSamplingMode();
 setToJackknifeMode();
 setSamplingBegin();
 if (exclude_full) setSamplingNext();
 uint k=0;
 while (!isSamplingEnd()){
    jacksamples[k++]=getCurrentSamplingValue(obskey);
    setSamplingNext();}
 setSamplingMode(currmode);
}



Vector<double> MCObsHandler::getBootstrapSamplingValues(const MCObsInfo& obskey)
{
 uint nboot=getNumberOfBootstrapResamplings();
 Vector<double> bootsamples(nboot);
 getBootstrapSamplingValues(obskey,&bootsamples[0],true);
 return bootsamples;
}


void MCObsHandler::getBootstrapSamplingValues(const MCObsInfo& obskey, vector<double>& bootsamples)
{
 uint nboot=getNumberOfBootstrapResamplings();
 bootsamples.resize(nboot);
 getBootstrapSamplingValues(obskey,&bootsamples[0],true);
}


void MCObsHandler::getBootstrapSamplingValues(const MCObsInfo& obskey, 
                                              double* bootsamples, bool exclude_full)
{
 if (!Bptr){
    cout << "Bootstrapper NOT set"<<endl;
    throw(std::invalid_argument("Bootstrapper not set"));}
 SamplingMode currmode=getCurrentSamplingMode();
 setToBootstrapMode();
 setSamplingBegin();
 if (exclude_full) setSamplingNext();
 uint k=0;
 while (!isSamplingEnd()){
    bootsamples[k++]=getCurrentSamplingValue(obskey);
    setSamplingNext();}
 setSamplingMode(currmode);
}



Vector<double> MCObsHandler::getSamplingValues(const MCObsInfo& obskey,
                                               SamplingMode mode)
{
 if (mode==Jackknife) return getJackknifeSamplingValues(obskey);
 else return getBootstrapSamplingValues(obskey);
}


void MCObsHandler::getSamplingValues(const MCObsInfo& obskey, 
                          vector<double>& samples, SamplingMode mode)
{
 if (mode==Jackknife) getJackknifeSamplingValues(obskey,samples);
 else getBootstrapSamplingValues(obskey,samples);
}




void MCObsHandler::getFullAndSamplingValues(const MCObsInfo& obskey, 
                          vector<double>& samples, SamplingMode mode)
{
 if (mode==Jackknife){
    samples.resize(m_nbins+1);
    getJackknifeSamplingValues(obskey,&samples[0],false);}
 else{
    uint nboot=getNumberOfBootstrapResamplings();
    samples.resize(nboot+1);
    getBootstrapSamplingValues(obskey,&samples[0],false);}
}



double MCObsHandler::getStandardDeviation(const MCObsInfo& obskey)
{
 return sqrt(getFullSampleCovariance(obskey,obskey));
}



double MCObsHandler::getJackKnifeError(const MCObsInfo& obskey, uint jacksize)
{
 if ((jacksize<1)||(jacksize>(m_nbins/12)))
    throw(std::invalid_argument("Invalid jack size in getJackKnifeError"));
 try{
    const Vector<double>& buffer=get_data(obskey);
    uint N=m_nbins - (m_nbins % jacksize);   // N is divisible by jacksize
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
    throw;}
}


MCEstimate MCObsHandler::getEstimate(const MCObsInfo& obskey)
{
 double full=getFullSampleValue(obskey);
 vector<double> sampvals;
 if (isJackknifeMode()){
    MCEstimate result(Jackknife);
    getJackknifeSamplingValues(obskey,sampvals);
    jack_analyze(full,sampvals,result);
    return result;}
 else{
    MCEstimate result(Bootstrap);
    getBootstrapSamplingValues(obskey,sampvals);
    boot_analyze(full,sampvals,result);
    return result;}
}


MCEstimate MCObsHandler::getEstimate(const MCObsInfo& obskey, SamplingMode inmode)
{
 double full=getFullSampleValue(obskey,inmode);
 vector<double> sampvals;
 MCEstimate result(inmode);
 if (inmode==Jackknife){
    getJackknifeSamplingValues(obskey,sampvals);
    jack_analyze(full,sampvals,result);
    return result;}
 else{
    getBootstrapSamplingValues(obskey,sampvals);
    boot_analyze(full,sampvals,result);
    return result;}
}



MCEstimate MCObsHandler::getJackknifeEstimate(const MCObsInfo& obskey)
{
 double full=getFullSampleValue(obskey,Jackknife);
 vector<double> sampvals;
 MCEstimate result(Jackknife);
 getJackknifeSamplingValues(obskey,sampvals);
 jack_analyze(full,sampvals,result);
 return result;
}


MCEstimate MCObsHandler::getBootstrapEstimate(const MCObsInfo& obskey)
{
 double full=getFullSampleValue(obskey,Bootstrap);
 vector<double> sampvals;
 MCEstimate result(Bootstrap);
 getBootstrapSamplingValues(obskey,sampvals);
 boot_analyze(full,sampvals,result);
 return result;
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
 if (markovtime>(m_nbins/12))
    throw(std::invalid_argument("Markov time too large in getAutoCorrelation"));
 const Vector<double>& buffer=get_data(obskey);
 uint N=m_nbins-markovtime;
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


   //     private routines


double MCObsHandler::calc_mean(const Vector<double>& dvals)
{
 unsigned int n=dvals.size();
 double dm=0.0;
 for (unsigned int k=0;k<n;k++) 
    dm+=dvals[k];
 return dm/double(n);
}


double MCObsHandler::calc_mean(const Vector<double>& dvals,
                               double full, uint binremove)
{
 unsigned int n=dvals.size();
 return (full*double(n)-dvals[binremove])/double(n-1);
}


double MCObsHandler::calc_mean(const Vector<double>& dvals,
                               const Vector<uint>& indmap)
{
 unsigned int n=indmap.size();
 double dm=0.0;
 for (unsigned int k=0;k<n;k++) 
    dm+=dvals[indmap[k]];
 return dm/double(n);
}


   // jindex=0 is full sample, jindex=1 is first jackknife sampling

double MCObsHandler::calc_jack_mean(const MCObsInfo& obskey, uint jindex)
{
 map<pair<MCObsInfo,uint>,double>::iterator dt
            =m_jacksamples.find(make_pair(obskey,jindex));
 if (dt!=m_jacksamples.end()) 
    return (dt->second);
 assert_simple(obskey,"getCurrentSamplingValue");
 const Vector<double>& buffer=getBins(obskey);
 double full=getFullSampleValue(obskey,Jackknife);
 double val=(jindex>0)?calc_mean(buffer,full,jindex-1):full;
 m_jacksamples.insert(make_pair(make_pair(obskey,jindex),val));
 return val;
}


   // bindex=0 is full sample, bindex=1 is first bootstrap sampling

double MCObsHandler::calc_boot_mean(const MCObsInfo& obskey, uint bindex)
{
 map<pair<MCObsInfo,uint>,double>::iterator dt
            =m_bootsamples.find(make_pair(obskey,bindex));
 if (dt!=m_bootsamples.end()) 
    return (dt->second);
 assert_simple(obskey,"getCurrentSamplingValue");
 const Vector<double>& buffer=getBins(obskey);
 double val;
 if (bindex>0){
    const Vector<uint>& indexmapper=Bptr->getResampling(bindex-1);
    val=calc_mean(buffer,indexmapper);}
 else
    val=getFullSampleValue(obskey,Bootstrap);
 m_bootsamples.insert(make_pair(make_pair(obskey,bindex),val));
 return val;
}





double MCObsHandler::calc_simple_covariance(const Vector<double>& dvals1,
                                            const Vector<double>& dvals2)
{
 unsigned int n=dvals1.size();
 double r1=1.0/double(n);
 double r2=1.0/(double(n-1));
 double dm1=0.0, dm2=0.0, dm12=0.0;
 for (unsigned int k=0;k<n;k++){
    dm1+=dvals1[k];
    dm2+=dvals2[k];
    dm12+=dvals1[k]*dvals2[k];}
 dm1*=r1; dm2*=r1; dm12*=r1;
 return (dm12-dm1*dm2)*r2;
}


double MCObsHandler::calc_simple_covariance(const Vector<double>& dvals1,
                                            const Vector<double>& dvals2,
                                            uint binremove)
{
 unsigned int n=dvals1.size()-1;
 double r1=1.0/double(n);
 double r2=1.0/(double(n-1));
 double dm1=0.0, dm2=0.0, dm12=0.0;
 for (unsigned int k=0;k<binremove;k++){
    dm1+=dvals1[k];
    dm2+=dvals2[k];
    dm12+=dvals1[k]*dvals2[k];}
 for (unsigned int k=binremove+1;k<=n;k++){
    dm1+=dvals1[k];
    dm2+=dvals2[k];
    dm12+=dvals1[k]*dvals2[k];}
 dm1*=r1; dm2*=r1; dm12*=r1;
 return (dm12-dm1*dm2)*r2;
}


double MCObsHandler::calc_simple_covariance(const Vector<double>& dvals1,
                                            const Vector<double>& dvals2,
                                            const Vector<uint>& indmap)
{
 unsigned int n=indmap.size();
 double r1=1.0/double(n);
 double r2=1.0/(double(n-1));
 double dm1=0.0, dm2=0.0, dm12=0.0;
 uint kk;
 for (unsigned int k=0;k<n;k++){
    kk=indmap[k];
    dm1+=dvals1[kk];
    dm2+=dvals2[kk];
    dm12+=dvals1[kk]*dvals2[kk];}
 dm1*=r1; dm2*=r1; dm12*=r1;
 return (dm12-dm1*dm2)*r2;
}




void MCObsHandler::jack_analyze(double full, vector<double>& sampvals, MCEstimate& result)
{
 uint n=sampvals.size();
 double avg=0.0;
 for (uint k=0;k<n;k++)
    avg+=sampvals[k];
 avg/=double(n);
 double var=0.0;
 for (uint k=0;k<n;k++){
    double x=sampvals[k]-avg;
    var+=x*x;}
 var*=(1.0-1.0/double(n));  // jackknife
 result.jackassign(full,avg,sqrt(var));
}

/*
double MCObsHandler::jack_covariance(const vector<double>& sampvals1, 
                                     const vector<double>& sampvals2)
{
 uint n=sampvals1.size();
 double avg1=0.0;
 double avg2=0.0;
 for (uint k=0;k<n;k++){
    avg1+=sampvals1[k];
    avg2+=sampvals2[k];}
 avg1/=double(n);
 avg2/=double(n);
 double jackcov=0.0;
 for (uint k=0;k<n;k++){
    jackcov+=(sampvals1[k]-avg1)*(sampvals2[k]-avg2);}
 jackcov*=(1.0-1.0/double(n));  // jackknife
 return jackcov;
}



double MCObsHandler::jack_covariance(const vector<double>& sampvals1, 
                                     const vector<double>& sampvals2,
                                     const Vector<uint>& indmap)
{
 uint n=sampvals.size();
 double avg1=0.0;
 double avg2=0.0;
 uint kk;
 for (uint k=0;k<n;k++){
    kk=indmap[k];
    avg1+=sampvals1[kk];
    avg2+=sampvals2[kk];}
 avg1/=double(n);
 avg2/=double(n);
 double jackcov=0.0;
 for (uint k=0;k<n;k++){
    kk=indmap[k];
    jackcov+=(sampvals1[kk]-avg1)*(sampvals2[kk]-avg2);}
 jackcov*=(1.0-1.0/double(n));  // jackknife
 return jackcov;
}
*/


   // This takes an array of "Nboot" bootstrap samples
   // and returns their average (ans.boot_avg), the
   // value which 84% of the samples lie above (ans.boot_low),
   // the value which 84% of the samples lie below (ans.boot_upp),
   // and the value which 50% of the samples lie above and
   // 50% lie below (ans.boot_med).

void MCObsHandler::boot_analyze(double full, vector<double>& sampvals, MCEstimate& result)
{
 uint nb=sampvals.size();
 double avg=0.0;
 for (uint k=0;k<nb;k++)
    avg+=sampvals[k];
 avg/=double(nb);
 double var=0.0;
 for (uint k=0;k<nb;k++){
    double x=sampvals[k]-avg;
    var+=x*x;}
 var/=double(nb-1);  // bootstrap

 const double conf_level=0.68;  // one standard deviation
 const double lowconf=0.5*(1.0-conf_level);
 const double uppconf=0.5*(1.0+conf_level);
 const double eps=1e-10;

 std::sort(sampvals.begin(),sampvals.end());

    // get the low and high confidence bounds and the median

 uint i=(unsigned int) floor( lowconf*nb+eps); 
 double low=0.5*(sampvals[i-1]+sampvals[i]);
 i=(unsigned int) floor( uppconf*nb+eps); 
 double upp=0.5*(sampvals[i-1]+sampvals[i]);
 i=(unsigned int) floor( 0.5*nb+eps); 
 double med=0.5*(sampvals[i-1]+sampvals[i]);

 result.bootassign(full,avg,sqrt(var),low,med,upp);
}

// *************************************************************

#ifdef COMPLEXNUMBERS

void MCObsHandler::set_up_corvev_bins(
                        const MCObsInfo& obskey, const Vector<double>* &corrbins, bool& realpart,
                        const Vector<double>* &snkvevrebins, const Vector<double>* &snkvevimbins,
                        const Vector<double>* &srcvevrebins, const Vector<double>* &srcvevimbins)
{
 if (!(obskey.getCorrelatorAtTimeInfo().isVEVsubtracted()))
    throw(std::invalid_argument("Could not get VEV info"));
 bool herm=obskey.isHermitianCorrelatorAtTime();
 unsigned int tval=obskey.getCorrelatorTimeIndex();
 OperatorInfo src(obskey.getCorrelatorSourceInfo());
 OperatorInfo snk(obskey.getCorrelatorSinkInfo());
 realpart=obskey.isRealPart();
 ComplexArg arg=(realpart)?RealPart:ImaginaryPart;
 MCObsInfo corr(snk,src,tval,herm,arg,false);   // no vev subtraction
 MCObsInfo src_re_info(src,RealPart);
 MCObsInfo snk_re_info(snk,RealPart);
 MCObsInfo src_im_info(src,ImaginaryPart);
 MCObsInfo snk_im_info(snk,ImaginaryPart);
 corrbins=&getBins(corr);
 srcvevrebins=&getBins(src_re_info);
 snkvevrebins=&getBins(snk_re_info);
 srcvevimbins=&getBins(src_im_info);
 snkvevimbins=&getBins(snk_im_info);
}


void MCObsHandler::get_jack_samples(const Vector<double>& corrbins, bool realpart,
                         const Vector<double>& snkvevrebins, const Vector<double>& snkvevimbins,
                         const Vector<double>& srcvevrebins, const Vector<double>& srcvevimbins,
                         std::vector<double>& sampvals)
{
 uint nbins=corrbins.size();
 double r1=1.0/double(nbins-1);
 double corrsum,snkvevresum,snkvevimsum,srcvevresum,srcvevimsum;
 corrsum=snkvevresum=snkvevimsum=srcvevresum=srcvevimsum=0.0;
 for (uint k=0;k<nbins;++k){
    corrsum+=corrbins[k];
    snkvevresum+=snkvevrebins[k];
    snkvevimsum+=snkvevimbins[k];
    srcvevresum+=srcvevrebins[k];
    srcvevimsum+=srcvevimbins[k];}
 sampvals.resize(nbins);
 for (uint k=0;k<nbins;++k){
    double corval=r1*(corrsum-corrbins[k]);
    double snk_re=r1*(snkvevresum-snkvevrebins[k]);
    double snk_im=r1*(snkvevimsum-snkvevimbins[k]);
    double src_re=r1*(srcvevresum-srcvevrebins[k]);
    double src_im=r1*(srcvevimsum-srcvevimbins[k]);
    if (realpart) sampvals[k]=corval-(snk_re*src_re+snk_im*src_im);
    else sampvals[k]=corval-(snk_im*src_re-snk_re*src_im);}
}

void MCObsHandler::get_jack_samples(const Vector<double>& corrbins, bool realpart,
                         const Vector<double>& snkvevrebins, const Vector<double>& snkvevimbins,
                         const Vector<double>& srcvevrebins, const Vector<double>& srcvevimbins,
                         const Vector<uint>& indmap, std::vector<double>& sampvals)
{
 uint nbins=corrbins.size();
 double r1=1.0/double(nbins-1);
 double corrsum,snkvevresum,snkvevimsum,srcvevresum,srcvevimsum;
 corrsum=snkvevresum=snkvevimsum=srcvevresum=srcvevimsum=0.0;
 for (uint k=0;k<nbins;++k){
    uint kk=indmap[k];
    corrsum+=corrbins[kk];
    snkvevresum+=snkvevrebins[kk];
    snkvevimsum+=snkvevimbins[kk];
    srcvevresum+=srcvevrebins[kk];
    srcvevimsum+=srcvevimbins[kk];}
 sampvals.resize(nbins);
 for (uint k=0;k<nbins;++k){
    uint kk=indmap[k];
    double corval=r1*(corrsum-corrbins[kk]);
    double snk_re=r1*(snkvevresum-snkvevrebins[kk]);
    double snk_im=r1*(snkvevimsum-snkvevimbins[kk]);
    double src_re=r1*(srcvevresum-srcvevrebins[kk]);
    double src_im=r1*(srcvevimsum-srcvevimbins[kk]);
    if (realpart) sampvals[k]=corval-(snk_re*src_re+snk_im*src_im);
    else sampvals[k]=corval-(snk_im*src_re-snk_re*src_im);}
}

#else

void MCObsHandler::set_up_corvev_bins(
                        const MCObsInfo& obskey, const Vector<double>* &corrbins, 
                        const Vector<double>* &snkvevrebins, const Vector<double>* &srcvevrebins)
{
 if (!(obskey.getCorrelatorAtTimeInfo().isVEVsubtracted()))
    throw(std::invalid_argument("Could not get VEV info"));
 bool herm=obskey.isHermitianCorrelatorAtTime();
 unsigned int tval=obskey.getCorrelatorTimeIndex();
 OperatorInfo src(obskey.getCorrelatorSourceInfo());
 OperatorInfo snk(obskey.getCorrelatorSinkInfo());
 MCObsInfo corr(snk,src,tval,herm,RealPart,false);   // no vev subtraction
 MCObsInfo src_re_info(src,RealPart);
 MCObsInfo snk_re_info(snk,RealPart);
 corrbins=&getBins(corr);
 srcvevrebins=&getBins(src_re_info);
 snkvevrebins=&getBins(snk_re_info);
}


void MCObsHandler::get_jack_samples(const Vector<double>& corrbins,
                         const Vector<double>& snkvevrebins, const Vector<double>& srcvevrebins, 
                         std::vector<double>& sampvals)
{
 uint nbins=corrbins.size();
 double r1=1.0/double(nbins-1);
 double corrsum,snkvevresum,srcvevresum;
 corrsum=snkvevresum=srcvevresum=0.0;
 for (uint k=0;k<nbins;++k){
    corrsum+=corrbins[k];
    snkvevresum+=snkvevrebins[k];
    srcvevresum+=srcvevrebins[k];}
 sampvals.resize(nbins);
 for (uint k=0;k<nbins;++k){
    double corval=r1*(corrsum-corrbins[k]);
    double snk_re=r1*(snkvevresum-snkvevrebins[k]);
    double src_re=r1*(srcvevresum-srcvevrebins[k]);
    sampvals[k]=corval-snk_re*src_re;}
}


void MCObsHandler::get_jack_samples(const Vector<double>& corrbins,
                         const Vector<double>& snkvevrebins, const Vector<double>& srcvevrebins, 
                         const Vector<uint>& indmap, std::vector<double>& sampvals)
{
 uint nbins=corrbins.size();
 double r1=1.0/double(nbins-1);
 double corrsum,snkvevresum,srcvevresum;
 corrsum=snkvevresum=srcvevresum=0.0;
 for (uint k=0;k<nbins;++k){
    uint kk=indmap[k];
    corrsum+=corrbins[kk];
    snkvevresum+=snkvevrebins[kk];
    srcvevresum+=srcvevrebins[kk];}
 sampvals.resize(nbins);
 for (uint k=0;k<nbins;++k){
    uint kk=indmap[k];
    double corval=r1*(corrsum-corrbins[kk]);
    double snk_re=r1*(snkvevresum-snkvevrebins[kk]);
    double src_re=r1*(srcvevresum-srcvevrebins[kk]);
    sampvals[k]=corval-snk_re*src_re;}
}

#endif



double MCObsHandler::jack_covariance(const std::vector<double>& sampvals1, 
                                     const std::vector<double>& sampvals2)
{
 double a1,a2;
 a1=a2=0.0;
 uint nsamp=sampvals1.size();
 double r=1.0/double(nsamp);
 for (uint k=0;k<nsamp;++k){
    a1+=sampvals1[k];
    a2+=sampvals2[k];}
 a1*=r; a2*=r;
 double cov=0.0;
 for (uint k=0;k<nsamp;++k){
    cov+=(sampvals1[k]-a1)*(sampvals2[k]-a2);}
 return (1.0-r)*cov;
}


// ************************************************************************

             // read all samplings from file and put into memory

void MCObsHandler::readSamplingValuesFromFile(const string& filename, 
                                              SamplingMode mode, XMLHandler& xmlout)
{
 xmlout.set_root("ReadFromFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 xmlout.put_child("FileName",fname);
 SamplingMode currentmode=getCurrentSamplingMode();
 string filetypeid("Sigmond--SamplingsFile");   
 string xmlheader;
 IOMap<UIntKey,DataType> iom;
 iom.openReadOnly(fname,filetypeid,xmlheader);
 XMLHandler xmlh;
 xmlh.set_from_string(xmlheader);
 if (!check_samplings_header(xmlh,mode)){
    xmlout.put_child("Error",string("Mismatch in header for file ")+fname);
    XMLHandler xmli("InFile"); xmli.put_child(xmlh);
    xmlout.put_child(xmli);
    return;}

 setSamplingMode(mode);
 set<UIntKey> ikeys;
 iom.getKeys(ikeys);
 DataType buffer;
 for (set<UIntKey>::const_iterator it=ikeys.begin();it!=ikeys.end();it++){
    XMLHandler xmle("Read");
    try{
       iom.get(*it,buffer);
       XMLHandler xmlk; xmlk.set_from_string(buffer.header);
       MCObsInfo obskey(xmlk);
       xmle.put_child(xmlk);
       uint count=0;
       for (setSamplingBegin();!isSamplingEnd();setSamplingNext()){
          putCurrentSamplingValue(obskey,buffer.data[count++],false);}}
    catch(const std::exception& errmsg){
       xmle.put_child("Error",string(errmsg.what()));}
    xmlout.put_child(xmle);}

 setSamplingMode(currentmode);
}


             // read all samplings from file and put into memory (this version
             // only reads those records matching the MCObsInfo objects in "obskeys")

void MCObsHandler::readSamplingValuesFromFile(const set<MCObsInfo>& obskeys, 
                                              const string& filename,
                                              SamplingMode mode, XMLHandler& xmlout)
{
 xmlout.set_root("ReadFromFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 xmlout.put_child("FileName",fname);
 SamplingMode currentmode=getCurrentSamplingMode();
 string filetypeid("Sigmond--SamplingsFile");   
 string xmlheader;
 IOMap<UIntKey,DataType> iom;
 iom.openReadOnly(fname,filetypeid,xmlheader);
 XMLHandler xmlh;
 xmlh.set_from_string(xmlheader);
 if (!check_samplings_header(xmlh,mode)){
    xmlout.put_child("Error",string("Mismatch in header for file ")+fname);
    XMLHandler xmli("InFile"); xmli.put_child(xmlh);
    xmlout.put_child(xmli);
    return;}

 set<MCObsInfo> seekkeys(obskeys); 
 setSamplingMode(mode);
 set<UIntKey> ikeys;
 iom.getKeys(ikeys);
 DataType buffer;
 for (set<UIntKey>::const_iterator it=ikeys.begin();it!=ikeys.end();it++){
    XMLHandler xmle("Read");
    try{
       iom.get(*it,buffer);
       XMLHandler xmlk; xmlk.set_from_string(buffer.header);
       MCObsInfo obskey(xmlk);
       if (seekkeys.find(obskey)!=seekkeys.end()){
          xmle.put_child(xmlk);
          uint count=0;
          try{
             for (setSamplingBegin();!isSamplingEnd();setSamplingNext()){
                putCurrentSamplingValue(obskey,buffer.data[count++],false);}
             seekkeys.erase(obskey);}
          catch(const std::exception& errmsgb){
             xmle.put_child("Error",string(errmsgb.what()));}
          xmlout.put_child(xmle);}}
    catch(const std::exception& errmsg){
       xmle.put_child("Error",string(errmsg.what())+" for ikey "+make_string(it->getValue()));}}

 if (!seekkeys.empty()){
    XMLHandler xmle("ReadFailures");
    xmle.put_child("Error","MCObservable keys missing or already in memory");
    for (set<MCObsInfo>::const_iterator it=seekkeys.begin();it!=seekkeys.end();it++){
       XMLHandler xmlo; it->output(xmlo);
       xmle.put_child(xmlo);}
    xmlout.put_child(xmle);}

 setSamplingMode(currentmode);
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
                                             SamplingMode mode, XMLHandler& xmlout,
                                             bool overwrite)
{
 xmlout.set_root("WriteToFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 xmlout.put_child("FileName",fname);
 SamplingMode currentmode=getCurrentSamplingMode();

 IOMap<UIntKey,DataType> iom;
 string filetypeid("Sigmond--SamplingsFile");   
 string xmlheader;
 bool fileopen=false;
 uint next_ikey=0;
 set<MCObsInfo> obskeysinfile;
 DataType buffer;

 if ((fileExists(fname))&&(!overwrite)){
    iom.openUpdate(fname,filetypeid,xmlheader);
    XMLHandler xmlr; xmlr.set_from_string(xmlheader);
    if (!check_samplings_header(xmlr,mode)){
       xmlout.put_child("Error",string("Mismatch in header for file ")+fname);
       XMLHandler xmli("InFile"); xmli.put_child(xmlr);
       xmlout.put_child(xmli);
       iom.close(); return;}
    fileopen=iom.isOpen();
    if (!fileopen) throw(std::invalid_argument("Could not open file in writeSamplingValuesToFile"));
    set<UIntKey> ikeys;
    iom.getKeys(ikeys);
    for (set<UIntKey>::const_iterator it=ikeys.begin();it!=ikeys.end();it++){
       if (it->getValue()>=next_ikey) next_ikey=it->getValue()+1;
       iom.get(*it,buffer);
       XMLHandler xmlk; xmlk.set_from_string(buffer.header);
       obskeysinfile.insert(MCObsInfo(xmlk));}}

 setSamplingMode(mode);
 for (set<MCObsInfo>::const_iterator it=obskeys.begin();it!=obskeys.end();it++){
    XMLHandler xmlo; it->output(xmlo);
    XMLHandler xmle("Write");
    if (obskeysinfile.find(*it)!=obskeysinfile.end()){
       xmle.put_child("Error","MCObservable already in file");
       xmle.put_child(xmlo);}
    else if (!queryAllSamplings(*it)){
       xmle.put_child("Error","Not all samplings are in memory");
       xmle.put_child(xmlo);}
    else{
       try{
          buffer.header=xmlo.str();
          getFullAndSamplingValues(*it,buffer.data,mode);
          if (!fileopen){   // must be new file if not already open
             XMLHandler xmlh;
             output_samplings_header(xmlh,mode);
             xmlheader=xmlh.str();
             iom.openNew(fname,filetypeid,xmlheader,false);
             if (!iom.isOpen())
                throw(std::invalid_argument("File could not be opened for output"));
             fileopen=true;}
          iom.put(UIntKey(next_ikey++),buffer);}
       catch(const std::exception& errmsg){
          xmle.put_child("Error",string(errmsg.what()));}
       xmle.put_child(xmlo);}
    xmlout.put_child(xmle);}

 setSamplingMode(currentmode);
}


void MCObsHandler::output_samplings_header(XMLHandler& xmlout, SamplingMode mode)
{
 xmlout.set_root("SigmondSamplingsFile");
 MCEnsembleInfo mcens=m_in_handler.getEnsembleInfo();
 XMLHandler xmltmp;
 mcens.output(xmltmp);
 xmlout.put_child(xmltmp);
 xmlout.put_child("NumberOfMeasurements",make_string(m_nmeasures));
 xmlout.put_child("NumberOfBins",make_string(m_nbins));
 if ((m_rebin>1)||(!m_omit.empty())){
    XMLHandler xmltw("TweakEnsemble");
    if (m_rebin>1) xmltw.put_child("Rebin",make_string(m_rebin));
    for (set<unsigned int>::const_iterator it=m_omit.begin();it!=m_omit.end();it++)
       xmltw.put_child("OmitConfig",make_string(*it));
    xmlout.put_child(xmltw);}
 if (mode==Bootstrap){
    const Bootstrapper& boot=getBootstrapper();
    XMLHandler xmlb("Bootstrapper");
    xmlb.put_child("NumberResamplings",make_string(boot.getNumberOfResamplings()));
    xmlb.put_child("Seed",make_string(boot.getRNGSeed()));
    xmlb.put_child("BootSkip",make_string(boot.getSkipValue()));
    xmlout.put_child(xmlb);}
 else{
    xmlout.put_child("JackKnifeMode");}
}


bool MCObsHandler::check_samplings_header(XMLHandler& xmlin, SamplingMode mode)
{
 try{
 XMLHandler xmlr(xmlin,"SigmondSamplingsFile");
 MCEnsembleInfo mcens(xmlr);
 if (mcens!=m_in_handler.getEnsembleInfo()) return false;
 uint nmeasread;
 xmlreadchild(xmlr,"NumberOfMeasurements",nmeasread);
 if (nmeasread!=m_nmeasures) return false;
 uint nbinread;
 xmlreadchild(xmlr,"NumberOfBins",nbinread);
 if (nbinread!=m_nbins) return false;
 uint rebin=1;
 set<uint> omit;
 if (xmlr.count_among_children("TweakEnsemble")==1){
    XMLHandler xmlk(xmlr,"TweakEnsemble");
    xmlreadifchild(xmlk,"Rebin",rebin);
    list<XMLHandler> xmlo=xmlk.find("OmitConfig");
    for (list<XMLHandler>::iterator tt=xmlo.begin();tt!=xmlo.end();tt++){
       uint anomission;
       xmlread(*tt,"OmitConfig",anomission,"dummy");
       omit.insert(anomission);}}
 if ((rebin!=m_rebin)||(omit!=m_omit)) return false;
 int bootcount=xmlr.count_among_children("Bootstrapper");
 int jackcount=xmlr.count_among_children("JackKnifeMode");
 if ((bootcount+jackcount)!=1) return false;
 if (bootcount==1){
    if (mode!=Bootstrap) return false;
    XMLHandler xmlb(xmlr,"Bootstrapper");
    uint num_resamplings=1024;
    unsigned long bootseed=0, bootskip=64;
    xmlreadifchild(xmlb,"NumberResamplings",num_resamplings);
    xmlreadifchild(xmlb,"Seed",bootseed);
    xmlreadifchild(xmlb,"BootSkip",bootskip);
    const Bootstrapper& boot=getBootstrapper();
    if ((num_resamplings!=boot.getNumberOfResamplings())
      ||(bootseed!=boot.getRNGSeed())
      ||(bootskip!=boot.getSkipValue())) return false;}
 else if ((mode!=Jackknife)||(jackcount!=1))
    return false;}
 catch(const std::exception& xp){
    return false;}
 return true;
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
 xmlout.put_child("FileName",fname);
 string filetypeid("Sigmond--BinsFile");   
 string xmlheader;
 IOMap<UIntKey,DataType> iom;
 iom.openReadOnly(fname,filetypeid,xmlheader);
 XMLHandler xmlh;
 xmlh.set_from_string(xmlheader);
 if (!check_bins_header(xmlh)){
    xmlout.put_child("Error",string("Mismatch in header for file ")+fname);
    XMLHandler xmli("InFile"); xmli.put_child(xmlh);
    xmlout.put_child(xmli);
    return;}

 set<UIntKey> ikeys;
 iom.getKeys(ikeys);
 DataType buffer;
 for (set<UIntKey>::const_iterator it=ikeys.begin();it!=ikeys.end();it++){
    XMLHandler xmle("Read");
    try{
       iom.get(*it,buffer);
       XMLHandler xmlk; xmlk.set_from_string(buffer.header);
       MCObsInfo obskey(xmlk);
       xmle.put_child(xmlk);
       putBins(obskey,Vector<double>(buffer.data));}
    catch(const std::exception& errmsg){
       xmle.put_child("Error",string(errmsg.what()));}
    xmlout.put_child(xmle);}
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
 xmlout.put_child("FileName",fname);
 string filetypeid("Sigmond--BinsFile");   
 string xmlheader;
 IOMap<UIntKey,DataType> iom;
 iom.openReadOnly(fname,filetypeid,xmlheader);
 XMLHandler xmlh;
 xmlh.set_from_string(xmlheader);
 if (!check_bins_header(xmlh)){
    xmlout.put_child("Error",string("Mismatch in header for file ")+fname);
    XMLHandler xmli("InFile"); xmli.put_child(xmlh);
    xmlout.put_child(xmli);
    return;}

 set<MCObsInfo> seekkeys(obskeys); 
 set<UIntKey> ikeys;
 iom.getKeys(ikeys);
 DataType buffer;
 for (set<UIntKey>::const_iterator it=ikeys.begin();it!=ikeys.end();it++){
    XMLHandler xmle("Read");
    try{
       iom.get(*it,buffer);
       XMLHandler xmlk; xmlk.set_from_string(buffer.header);
       MCObsInfo obskey(xmlk);
       if (seekkeys.find(obskey)!=seekkeys.end()){
          xmle.put_child(xmlk);
          try{
             putBins(obskey,Vector<double>(buffer.data));
             seekkeys.erase(obskey);}
          catch(const std::exception& errmsgb){
             xmle.put_child("Error",string(errmsgb.what()));}
          xmlout.put_child(xmle);}}
    catch(const std::exception& errmsg){
       xmle.put_child("Error",string(errmsg.what())+" for ikey "+make_string(it->getValue()));}}

 if (!seekkeys.empty()){
    XMLHandler xmle("ReadFailures");
    xmle.put_child("Error","MCObservable keys missing or already in memory");
    for (set<MCObsInfo>::const_iterator it=seekkeys.begin();it!=seekkeys.end();it++){
       XMLHandler xmlo; it->output(xmlo);
       xmle.put_child(xmlo);}
    xmlout.put_child(xmle);}
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
                                   XMLHandler& xmlout, bool overwrite)
{
 xmlout.set_root("WriteBinsToFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 xmlout.put_child("FileName",fname);

 IOMap<UIntKey,DataType> iom;
 string filetypeid("Sigmond--BinsFile");   
 string xmlheader;
 bool fileopen=false;
 uint next_ikey=0;
 set<MCObsInfo> obskeysinfile;
 DataType buffer;

 if ((fileExists(fname))&&(!overwrite)){
    iom.openUpdate(fname,filetypeid,xmlheader);
    XMLHandler xmlr; xmlr.set_from_string(xmlheader);
    if (!check_bins_header(xmlr)){
       xmlout.put_child("Error",string("Mismatch in header for file ")+fname);
       XMLHandler xmli("InFile"); xmli.put_child(xmlr);
       xmlout.put_child(xmli);
       iom.close(); return;}
    fileopen=iom.isOpen();
    if (!fileopen) throw(std::invalid_argument("Could not open file in writeBinsToFile"));
    set<UIntKey> ikeys;
    iom.getKeys(ikeys);
    for (set<UIntKey>::const_iterator it=ikeys.begin();it!=ikeys.end();it++){
       if (it->getValue()>=next_ikey) next_ikey=it->getValue()+1;
       iom.get(*it,buffer);
       XMLHandler xmlk; xmlk.set_from_string(buffer.header);
       obskeysinfile.insert(MCObsInfo(xmlk));}}

 for (set<MCObsInfo>::const_iterator it=obskeys.begin();it!=obskeys.end();it++){
    XMLHandler xmlo; it->output(xmlo);
    XMLHandler xmle("Write");
    if (obskeysinfile.find(*it)!=obskeysinfile.end()){
       xmle.put_child("Error","MCObservable already in file");
       xmle.put_child(xmlo);}
    else if (!queryBins(*it)){
       xmle.put_child("Error","Not all bins are not in memory");
       xmle.put_child(xmlo);}
    else{
       try{
       buffer.header=xmlo.str();
       getBins(*it,buffer.data);
       if (!fileopen){   // must be new file if not already open
          XMLHandler xmlh;
          output_bins_header(xmlh);
          xmlheader=xmlh.str();
          iom.openNew(fname,filetypeid,xmlheader,false);
          if (!iom.isOpen())
             throw(std::invalid_argument("File could not be opened for output"));
          fileopen=true;}
       iom.put(UIntKey(next_ikey++),buffer);}
       catch(const std::exception& errmsg){
          xmle.put_child("Error",string(errmsg.what()));}
       xmle.put_child(xmlo);}
    xmlout.put_child(xmle);}
}


void MCObsHandler::output_bins_header(XMLHandler& xmlout)
{
 xmlout.set_root("SigmondBinsFile");
 MCEnsembleInfo mcens=m_in_handler.getEnsembleInfo();
 XMLHandler xmltmp;
 mcens.output(xmltmp);
 xmlout.put_child(xmltmp);
 xmlout.put_child("NumberOfMeasurements",make_string(m_nmeasures));
 xmlout.put_child("NumberOfBins",make_string(m_nbins));
 if ((m_rebin>1)||(!m_omit.empty())){
    XMLHandler xmltw("TweakEnsemble");
    if (m_rebin>1) xmltw.put_child("Rebin",make_string(m_rebin));
    for (set<unsigned int>::const_iterator it=m_omit.begin();it!=m_omit.end();it++)
       xmltw.put_child("OmitConfig",make_string(*it));
    xmlout.put_child(xmltw);}
}


bool MCObsHandler::check_bins_header(XMLHandler& xmlin)
{
 try{
 XMLHandler xmlr(xmlin,"SigmondBinsFile");
 MCEnsembleInfo mcens(xmlr);
 if (mcens!=m_in_handler.getEnsembleInfo()) return false;
 uint nmeasread;
 xmlreadchild(xmlr,"NumberOfMeasurements",nmeasread);
 if (nmeasread!=m_nmeasures) return false;
 uint nbinread;
 xmlreadchild(xmlr,"NumberOfBins",nbinread);
 if (nbinread!=m_nbins) return false;
 uint rebin=1;
 set<uint> omit;
 if (xmlr.count_among_children("TweakEnsemble")==1){
    XMLHandler xmlk(xmlr,"TweakEnsemble");
    xmlreadifchild(xmlk,"Rebin",rebin);
    list<XMLHandler> xmlo=xmlk.find("OmitConfig");
    for (list<XMLHandler>::iterator tt=xmlo.begin();tt!=xmlo.end();tt++){
       uint anomission;
       xmlread(*tt,"OmitConfig",anomission,"dummy");
       omit.insert(anomission);}}
 if ((rebin!=m_rebin)||(omit!=m_omit)) return false;}
 catch(const std::exception& xp){
    return false;}
 return true;
}


// ************************************************************************

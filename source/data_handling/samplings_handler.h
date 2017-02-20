#ifndef SAMP_HANDLER_H
#define SAMP_HANDLER_H

#include "data_io_handler.h"
#include "bins_info.h"
#include "mcobs_info.h"
#include "sampling_info.h"
#include "matrix.h"
#include <set>


 // *********************************************************************************
 // *                                                                               *
 // *   The classes defined in this file that are important are                     *
 // *                                                                               *
 // *         SamplingsGetHandler                                                   *
 // *         SamplingsPutHandler                                                   *
 // *                                                                               *
 // *   Theses classes handle inserting and retrieving data to/from SigMonD         *
 // *   sampling files.  The "put" objects just deal with one file, whereas         *
 // *   the "get" objects deal with a set of files.  These classes do NOT store     *
 // *   the data in memory; they just handle reading and writing to files.          *
 // *   The record keys for all of these classes are of class "MCObsInfo".          *
 // *   The data types are Vector<double>. For sampling files, the Vector size      *
 // *   should be the number of resamplings PLUS one (for the full sample).         *
 // *                                                                               *
 // *   The header in the samplings files must have the form                        *
 // *                                                                               *
 // *     <SigmondSamplingsFile>                                                    *
 // *       <MCBinsInfo> ... </MCBinsInfo>                                          *
 // *       <MCSamplingInfo> ... </MCSamplingInfo>                                  *
 // *     </SigmondSamplingsFile>                                                   *
 // *                                                                               *
 // *   Recall that the above tags have the forms                                   *
 // *                                                                               *
 // *      <MCBinsInfo>                                                             *
 // *        <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo>            *
 // *        <TweakEnsemble>  (optional)                                            *
 // *           <Rebin>2</Rebin>                                                    *
 // *           <Omissions>2 7 11</Omissions>                                       *
 // *        </TweakEnsemble>                                                       *
 // *      </MCBinsInfo>                                                            *
 // *                                                                               *
 // *      <MCSamplingInfo>                                                         *
 // *         <Jackknife/>                                                          *
 // *      </MCSamplingInfo>                                                        *
 // *                                                                               *
 // *      <MCSamplingInfo>                                                         *
 // *         <Bootstrapper>                                                        *
 // *            <NumberResamplings>2048</NumberResamplings>                        *
 // *            <Seed>6754</Seed>                                                  *
 // *            <BootSkip>127</BootSkip>                                           *
 // *         </Bootstrapper>                                                       *
 // *      </MCSamplingInfo>                                                        *
 // *                                                                               *
 // *                                                                               *
 // *   The "get" handlers take a set of file names, whereas the "put" handlers     *
 // *   require only one file name.  Errors encountered in the "get" handlers are   *
 // *   fatal and cause program execution abort (if you are getting data, you       *
 // *   presumably really need it, so any error should be fatal).  The "put"        *
 // *   handlers throw exceptions if errors are encountered.                        *
 // *                                                                               *
 // *   Objects of these "put" handlers always assume an "updating" mode. Existing  *
 // *   files are never erased, and new files are created as needed.  New records   *
 // *   are added to the files.  If the key of a record to be put already exists    *
 // *   in a file, the put will only occur if "overwrite" is specified AND the      *
 // *   size of the data to be put does not exceed the size of the data already in  *
 // *   the file for that key.                                                      *
 // *                                                                               *
 // *   The sampling mode in each file must MATCH the sampling mode of the          *
 // *   handler as stored in "m_samp_info".                                         *
 // *                                                                               *
 // *                                                                               *
 // *********************************************************************************



   // **************************************************************
   // *                                                            *
   // *                    SamplingsGetHandler                     *
   // *                                                            *
   // **************************************************************


class SamplingsGetHandler
{

    MCBinsInfo m_bin_info;
    MCSamplingInfo m_samp_info;
    DataGetHandlerMF<SamplingsGetHandler,MCObsInfo,std::vector<double> > *m_get;

 public:

    SamplingsGetHandler(const MCBinsInfo& binfo, const MCSamplingInfo& sampinfo, 
                        const std::set<std::string>& file_names, bool use_checksums=false)
           : m_bin_info(binfo),  m_samp_info(sampinfo), m_get(0)
     {m_get=new DataGetHandlerMF<SamplingsGetHandler,MCObsInfo,std::vector<double> >(*this,
            file_names,std::string("Sigmond--SamplingsFile"),use_checksums);}

    void addFile(const std::string& file_name)
     {m_get->addFile(file_name);}

    bool addFile(const std::string& file_name, const std::set<MCObsInfo>& keys_to_keep)
     {return m_get->addFile(file_name,keys_to_keep);}

    void removeFile(const std::string& file_name)
     {m_get->removeFile(file_name);}

    void clear()
     {m_get->clear();}
 
    void close()
     {m_get->close();}
 
    ~SamplingsGetHandler() {delete m_get;}

    bool keepKeys(const std::set<MCObsInfo>& keys_to_keep) 
     {return m_get->keepKeys(keys_to_keep);}


    bool queryData(const MCObsInfo& rkey)
     {return m_get->queryData(rkey);}

    void getData(const MCObsInfo& rkey, Vector<double>& result)
     {std::vector<double> buffer; m_get->getData(rkey,buffer);
      result=Vector<double>(buffer);}

    bool getDataMaybe(const MCObsInfo& rkey, Vector<double>& result)
     {result.clear(); std::vector<double> buffer; 
      bool info=m_get->getDataMaybe(rkey,buffer);
      if (info) result=Vector<double>(buffer);
      return info;}


    std::set<MCObsInfo> getKeys()
     {return m_get->getKeys();}

    void outputKeys(XMLHandler& xmlout)
     {m_get->outputKeys(xmlout);}
    
    unsigned int size() const 
     {return m_get->size();}



    bool checkHeader(XMLHandler& xmlin)
     {try{XMLHandler xmlb(xmlin,"SigmondSamplingsFile");
      MCBinsInfo bchk(xmlb); MCSamplingInfo schk(xmlb);
      return ((bchk==m_bin_info)&&(schk==m_samp_info));}
      catch(std::exception& xp){ return false;}}

    
          // disallow copies and default
    SamplingsGetHandler();
    SamplingsGetHandler(const SamplingsGetHandler& in);
    SamplingsGetHandler& operator=(const SamplingsGetHandler& in);
};


   // **************************************************************
   // *                                                            *
   // *                    SamplingsPutHandler                     *
   // *                                                            *
   // **************************************************************


class SamplingsPutHandler
{

    MCBinsInfo m_bin_info;
    MCSamplingInfo m_samp_info;
    DataPutHandlerSF<SamplingsPutHandler,MCObsInfo,std::vector<double> > *m_put;

 public:

    SamplingsPutHandler(const MCBinsInfo& binfo, const MCSamplingInfo& sampinfo, 
                        const std::string& file_name, bool overwrite=false, 
                        bool use_checksums=false)
           : m_bin_info(binfo),  m_samp_info(sampinfo), m_put(0)
     {m_put=new DataPutHandlerSF<SamplingsPutHandler,MCObsInfo,std::vector<double> >(*this,
            file_name,std::string("Sigmond--SamplingsFile"),overwrite,use_checksums);}

    ~SamplingsPutHandler() {delete m_put;}

    void putData(const MCObsInfo& rkey, const Vector<double>& data)
     {if (data.size()!=m_samp_info.getNumberOfReSamplings(m_bin_info)+1)
         throw(std::runtime_error("Wrong number of Samplings in SamplingsPutHandler::putData"));
      try{
         m_put->putData(rkey,data.c_vector());}
      catch(std::exception& xp){
         throw(std::runtime_error(std::string("putData failed: permission problem or record already in file and no overwrite -- ")
           +xp.what()));}}


    void flush()
     {m_put->flush();}
  
    bool queryData(const MCObsInfo& rkey)  // already exists?
     {return m_put->queryData(rkey);}


    bool checkHeader(XMLHandler& xmlin)
     {try{XMLHandler xmlb(xmlin,"SigmondSamplingsFile");
      MCBinsInfo bchk(xmlb); MCSamplingInfo schk(xmlb);
      return ((bchk==m_bin_info)&&(schk==m_samp_info));}
      catch(std::exception& xp){ return false;}}

    void writeHeader(XMLHandler& xmlout)
     {xmlout.set_root("SigmondSamplingsFile");
      XMLHandler xmlb; m_bin_info.output(xmlb);
      xmlout.put_child(xmlb);
      m_samp_info.output(xmlb);
      xmlout.put_child(xmlb);}
    
          // disallow copies and default
    SamplingsPutHandler();
    SamplingsPutHandler(const SamplingsPutHandler& in);
    SamplingsPutHandler& operator=(const SamplingsPutHandler& in);
};


// ************************************************************************************
#endif

#ifndef BINS_HANDLER_H
#define BINS_HANDLER_H

#include "data_io_handler.h"
#include "bins_info.h"
#include "mcobs_info.h"
#include "matrix.h"
#include <set>


 // *********************************************************************************
 // *                                                                               *
 // *   The classes defined in this file that are important are                     *
 // *                                                                               *
 // *         BinsGetHandler                                                        *
 // *         BinsPutHandler                                                        *
 // *                                                                               *
 // *   Theses classes handle inserting and retrieving data to/from SigMonD bin     *
 // *   files.  The "put" objects just deal with one file, whereas the "get"        *
 // *   objects deal with a set of files.  These classes do NOT store               *
 // *   the data in memory; they just handle reading and writing to files.          *
 // *   The record keys for all of these classes are of class "MCObsInfo".          *
 // *   The data types are Vector<double>.  For bin files, the Vector size should   *
 // *   be the number of bins.                                                      *
 // *                                                                               *
 // *   The header in the bin files must have the form                              *
 // *                                                                               *
 // *     <SigmondBinsFile>                                                         *
 // *       <MCBinsInfo> ... </MCBinsInfo>                                          *
 // *     </SigmondBinsFile>                                                        *
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
 // *   The "get" handlers take a set of file names, whereas the "put" handlers     *
 // *   require only one file name.  Errors encountered in the "get" handlers are   *
 // *   fatal and cause program execution abort (if you are getting data, you       *
 // *   presumably really need it, so any error should be fatal).  The "put"        *
 // *   handlers throw exceptions if errors are encountered.                        *
 // *                                                                               *
 // *                                                                               *
 // *********************************************************************************



   // **************************************************************
   // *                                                            *
   // *                       BinsGetHandler                       *
   // *                                                            *
   // **************************************************************


class BinsGetHandler
{

    MCBinsInfo m_bins_info;
    DataGetHandlerMF<BinsGetHandler,MCObsInfo,std::vector<double> > *m_get;

 public:

    BinsGetHandler(const MCBinsInfo& binfo, const std::set<std::string>& file_names, 
                   bool use_checksums=false)
           : m_bins_info(binfo), m_get(0)
     {m_get=new DataGetHandlerMF<BinsGetHandler,MCObsInfo,std::vector<double> >(*this,
            file_names,std::string("Sigmond--BinsFile"),use_checksums);}

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
 
    ~BinsGetHandler() {delete m_get;}

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
     {try{XMLHandler xmlb(xmlin,"SigmondBinsFile");
      MCBinsInfo chk(xmlb); return (chk==m_bins_info);}
      catch(std::exception& xp){ return false;}}

    
          // disallow copies and default
    BinsGetHandler();
    BinsGetHandler(const BinsGetHandler& in);
    BinsGetHandler& operator=(const BinsGetHandler& in);
};


   // **************************************************************
   // *                                                            *
   // *                       BinsPutHandler                       *
   // *                                                            *
   // **************************************************************


class BinsPutHandler
{

    MCBinsInfo m_bins_info;
    DataPutHandlerSF<BinsPutHandler,MCObsInfo,std::vector<double> > *m_put;

 public:

    BinsPutHandler(const MCBinsInfo& binfo, const std::string& file_name, 
                   bool overwrite=false, bool use_checksums=false)
           : m_bins_info(binfo), m_put(0)
     {m_put=new DataPutHandlerSF<BinsPutHandler,MCObsInfo,std::vector<double> >(*this,
            file_name,std::string("Sigmond--BinsFile"),overwrite,use_checksums);}

    ~BinsPutHandler() {delete m_put;}

    void putData(const MCObsInfo& rkey, const Vector<double>& data)
     {if (data.size()!=m_bins_info.getNumberOfBins())
         throw(std::runtime_error("Wrong number of bins in BinsPutHandler::putData"));
      if (!rkey.isSimple())
         throw(std::runtime_error("Only simple observable allowed for BinsPutHandler::putData"));
      try{
         m_put->putData(rkey,data.c_vector());}
      catch(std::exception& xp){
         throw(std::runtime_error((std::string("putData failed: permission problem or record already in file and no overwrite -- ")
           +xp.what()).c_str()));}}

    void flush()
     {m_put->flush();}
  
    void close()
     {m_put->close();}
 
    bool queryData(const MCObsInfo& rkey)  // already exists?
     {return m_put->queryData(rkey);}


    bool checkHeader(XMLHandler& xmlin)
     {try{XMLHandler xmlb(xmlin,"SigmondBinsFile");
      MCBinsInfo chk(xmlb); return (chk==m_bins_info);}
      catch(std::exception& xp){ return false;}}

    void writeHeader(XMLHandler& xmlout)
     {xmlout.set_root("SigmondBinsFile");
      XMLHandler xmlb; m_bins_info.output(xmlb);
      xmlout.put_child(xmlb);}
    
          // disallow copies and default
    BinsPutHandler();
    BinsPutHandler(const BinsPutHandler& in);
    BinsPutHandler& operator=(const BinsPutHandler& in);
};


// ************************************************************************************
#endif

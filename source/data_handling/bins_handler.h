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
 // *   be the number of retained bins.  Note: these classes are merely for         *
 // *   reading and writing bins; there is no storage, and no rebinning is done.    *
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
 // *   handlers throw exceptions if errors are encountered.  The bins info stored  *
 // *   in the file must be consistent with that requested: same omissions, but     *
 // *   rebin factor must be integer multiple of rebin factor in file.              *
 // *                                                                               *
 // *   Objects of these "put" handlers open file in one of three modes:            *
 // *         Protect, Update, Overwrite                                            *
 // *   "Protect" will either open an existing file or create a new file if it      *
 // *   does not exist, and it will add new records to the file, but it will NOT    *
 // *   overwrite existing records.  "Update" does the same as "Protect" except     *
 // *   that it will overwrite existing records.  "Overwrite" creates a new file,   *
 // *   deleting the old one if its exists. If the key of a record to be put        *
 // *   already exists in a file, the put will only occur in "Update" or            *
 // *   "Overwrite" mode AND the size of the data to be put does not exceed the     *
 // *   size of the data already in the file for that key.                          *
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
           : m_bins_info(binfo), m_get(0)     // sets m_bins_info temporarily
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

    MCBinsInfo getBinsInfo() const
     {return m_bins_info;}
 
    bool keepKeys(const std::set<MCObsInfo>& keys_to_keep) 
     {return m_get->keepKeys(keys_to_keep);}


    bool queryData(const MCObsInfo& rkey) const
     {return m_get->queryData(rkey);}

    void getData(const MCObsInfo& rkey, Vector<double>& result) const
     {std::vector<double> buffer; m_get->getData(rkey,buffer);
      result=Vector<double>(buffer);}

    bool getDataMaybe(const MCObsInfo& rkey, Vector<double>& result) const
     {result.clear(); std::vector<double> buffer; 
      bool info=m_get->getDataMaybe(rkey,buffer);
      if (info) result=Vector<double>(buffer);
      return info;}

          // sym averages with Herm transpose if isHermitianCorrelatorAtTime

    bool querySymData(const MCObsInfo& rkey) const
     {if ((rkey.isImagDiagOfHermCorr())||(m_get->queryData(rkey))) return true;
      if (!(rkey.isHermitianCorrelatorAtTime())) return false;
      return m_get->queryData(rkey.getTimeFlipped());}

    void getSymData(const MCObsInfo& rkey, Vector<double>& result) const
     {bool info=getSymDataMaybe(rkey,result);
      if (!info) throw(std::invalid_argument("getSymData failed in BinsGetHandler"));}

    bool getSymDataMaybe(const MCObsInfo& rkey, Vector<double>& result) const
     {if (rkey.isImagDiagOfHermCorr()){
         result=Vector<double>(m_bins_info.getNumberOfBins(),0.0); return true;}
      result.clear();
      if (rkey.hasNoRelatedFlip()){
         std::vector<double> buffer; 
         bool info=m_get->getDataMaybe(rkey,buffer);
         if (info) result=Vector<double>(buffer);
         return info;}
      std::vector<double> buffer;
      bool info1=m_get->getDataMaybe(rkey,buffer);
      if (info1) result=Vector<double>(buffer);
      Vector<double> res2;
      bool info2=m_get->getDataMaybe(rkey.getTimeFlipped(),buffer);
      if (info2){
         res2=Vector<double>(buffer);
         if (rkey.isImaginaryPart()) res2*=-1.0;}
      if ((!info1)&&(!info2)) return false;
      else if ((!info1)&&(info2)){ result=res2;}
      else if (info2){ result+=res2; result*=0.5;}
      return true;}


    std::set<MCObsInfo> getKeys() const
     {return m_get->getKeys();}

    void outputKeys(XMLHandler& xmlout) const
     {m_get->outputKeys(xmlout);}
    
    void getFileMap(XMLHandler& xmlout) const
     {m_get->getFileMap(xmlout);}

    std::set<std::string> getFileNames() const
     {return m_get->getFileNames();}

    unsigned int size() const 
     {return m_get->size();}


        // check that BinsInfo in header is consistent with
        // requested BinsInfo (omissions must match, and
        // rebin must be multiple of rebin in file)

    bool checkHeader(XMLHandler& xmlin)
     {try{XMLHandler xmlb(xmlin,"SigmondBinsFile");
      MCBinsInfo chk(xmlb);
      if (chk.isConsistentWith(m_bins_info)){
         return true;}
      else
         return false;}
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
                   WriteMode wmode=Protect, bool use_checksums=false)
           : m_bins_info(binfo), m_put(0)
     {m_put=new DataPutHandlerSF<BinsPutHandler,MCObsInfo,std::vector<double> >(*this,
            file_name,std::string("Sigmond--BinsFile"),wmode,use_checksums);}

    ~BinsPutHandler() {delete m_put;}

    void putData(const MCObsInfo& rkey, const Vector<double>& data)
     {if (data.size()!=m_bins_info.getNumberOfBins())
         throw(std::runtime_error("Wrong number of bins in BinsPutHandler::putData"));
      if (!rkey.isSimple())
         throw(std::runtime_error("Only simple observable allowed for BinsPutHandler::putData"));
      try{
         m_put->putData(rkey,data.c_vector());}
      catch(std::exception& xp){
         throw(std::runtime_error(std::string("putData failed: permission problem or record already in file and no overwrite -- ")
           +xp.what()));}}

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

#ifndef VEV_DATA_HANDLER_H
#define VEV_DATA_HANDLER_H

#include "scalar_defs.h"
#include "data_io_handler.h"
#include "multi_compare.h"
#include "operator_info.h"
#include "ensemble_info.h"
#include "filelist_info.h"

namespace LaphEnv {

// *****************************************************************
// *                                                               *
// *  "VEVDataHandler" handles access to the LapH hadron-operator  *
// *  vacuum expectation values.  This is a LOW-LEVEL handler      *
// *  used by a mid-level object of class "MCObsGetHandler".       *
// *  The main routine to be used is the "getData" member.         *
// *                                                               *
// *  The constructor has the form                                 *
// *                                                               *
// *   VEVDataHandler(                                             *
// *          const std::list<FileListInfo>& inputfiles,           *
// *          const std::set<OperatorInfo>& opSet,                 *
// *          const MCEnsembleInfo *ensemble,                      *
// *          bool use_checksums=false);                           *
// *                                                               *
// *  Files containing data for the correlators to be analyzed     *
// *  are specified in "inputfiles" above.  If the pointer         *
// *  "ensemble" is not null, the constructor ensures that all     *
// *  data read corresponds to the requested ensemble; if this     *
// *  pointer is null, the ensemble information is obtained from   *
// *  the data input files, but all input files must correspond    *
// *  to the same ensemble. "use_checksums" is used to instruct    *
// *  the object to do check sum testing or not.                   *
// *                                                               *
// *  If the set "opSet" is empty, then all VEVs found while       *
// *  reading the files will be used.  If any VEVs ARE specified   *
// *  in "opSet", then only those VEVs will be considered for      *
// *  input: files corresponding to other VEVs will be ignored,    *
// *  and an exception is thrown if the file containing the data   *
// *  for any VEV cannot be found.                                 *
// *                                                               *
// *                                                               *
// *****************************************************************


// *****************************************************************
// *                                                               *
// *  Implementation notes:                                        *
// *                                                               *
// *   -  m_getter is a pointer to a "DataGetHandlerMF" object     *
// *      which handles all data retrieval.  The variable members  *
// *      "maxgetopen" and "cleanfrac" control memory limitations  *
// *      and garbage collection.                                  *
// *                                                               *
// *****************************************************************



class VEVDataHandler
{

 public:

   class RecordKey
   {
      int config_serial_index;
   
    public:
   
      RecordKey(int index) 
           : config_serial_index(index) {}

      RecordKey(const RecordKey& in) 
           : config_serial_index(in.config_serial_index) {}

      RecordKey& operator=(const RecordKey& in) 
       {config_serial_index=in.config_serial_index;
        return *this;}

      ~RecordKey() {}
   
      bool operator<(const RecordKey& rhs) const 
       {return (config_serial_index<rhs.config_serial_index);}

      bool operator==(const RecordKey& rhs) const 
       {return (config_serial_index==rhs.config_serial_index);}

      bool operator!=(const RecordKey& rhs) const 
       {return (config_serial_index!=rhs.config_serial_index);}
   
      int getConfigSerialIndex() const {return config_serial_index;}
  
      void output(XMLHandler& xmlw) const 
      {xmlw.set_root("Key");
       xmlw.put_child("ConfigSerialIndex",make_string(config_serial_index));}

      std::string output(int indent=0) const
       {XMLHandler xmlout; output(xmlout);
        return xmlout.output(indent);}

      std::string str() const
       {XMLHandler xmlout; output(xmlout);
        return xmlout.str();}

      explicit RecordKey(const unsigned int* buf) 
       { config_serial_index=buf[0];}

      int numints() const {return 1;} 

      size_t numbytes() const {return sizeof(unsigned int);}

      void copyTo(unsigned int* buf) const 
       { buf[0]=config_serial_index; }


   };
   

   typedef OperatorInfo    FileKey;
   typedef InScalar        DataType;

   DataGetHandlerMF<FileKey,RecordKey,DataType> *m_getter;
   const MCEnsembleInfo *m_ensembleptr;

   static const unsigned int maxgetopen=512;
#ifdef NO_CXX11
   static const double cleanfrac=0.33;
#else
   static constexpr double cleanfrac=0.33;
#endif

       // Prevent copying ... handler might contain large
       // amounts of data

#ifndef NO_CXX11
   VEVDataHandler() = delete;
   VEVDataHandler(const VEVDataHandler&) = delete;
   VEVDataHandler& operator=(const VEVDataHandler&) = delete;
#else
   VEVDataHandler();
   VEVDataHandler(const VEVDataHandler&);
   VEVDataHandler& operator=(const VEVDataHandler&);
#endif



 public:


   VEVDataHandler(const std::list<FileListInfo>& inputfiles,
                  const std::set<OperatorInfo>& opSet,
                  const MCEnsembleInfo *ensemble,
                  bool use_checksums=false);

   ~VEVDataHandler();

   std::set<OperatorInfo> getOperatorSet() const;

   MCEnsembleInfo getEnsemble() const;

   std::string getEnsembleId() const;

   unsigned int getNumberOfMeasurements();

   void getData(const OperatorInfo& opkey, int serial_index, 
                Scalar& data);

   bool queryData(const OperatorInfo& opkey, int serial_index);

   void close();


   void getFileMap(XMLHandler& xmlout) const;

   std::set<OperatorInfo> getFileKeys() const;

   std::set<RecordKey> getKeys(const OperatorInfo& fkey);

   void outputKeys(XMLHandler& xmlout);

};



// ***************************************************************

inline size_t numbytes(IOHandler& ioh, const VEVDataHandler::RecordKey& rkey)
{
 return rkey.numbytes();
}

// ***************************************************************
}
#endif  

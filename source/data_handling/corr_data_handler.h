#ifndef CORR_DATA_HANDLER_H
#define CORR_DATA_HANDLER_H

#include "scalar_defs.h"
#include "laph_data_io_handler.h"
#include "multi_compare.h"
#include "operator_info.h"
#include "correlator_info.h"
#include "ensemble_info.h"
#include "filelist_info.h"

namespace LaphEnv {

// *****************************************************************
// *                                                               *
// *  "BLCorrelatorDataHandler" handles access to the basic LapH   *
// *  hadronic-operator correlators.  This is a LOW-LEVEL handler  *
// *  used by a mid-level object of class "MCObsGetHandler".       *
// *  The main routines to be used are the "getData" and           *
// *  "getSymData" members.                                        *
// *                                                               *
// *  The constructor has the form                                 *
// *                                                               *
// *   BLCorrelatorDataHandler(                                    *
// *          const std::list<FileListInfo>& inputfiles,           *
// *          const std::set<CorrelatorInfo>& corrSetNoSym,        *
// *          const std::set<CorrelatorInfo>& corrSetSym,          *
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
// *  If the sets "corrSetNoSym" and "corrSetSym" are empty, then  *
// *  all correlators found while reading the files will be used.  *
// *  If any correlators ARE specified in these sets, then only    *
// *  those correlators will be considered for input: files        *
// *  corresponding to other correlators will be ignored, and an   *
// *  exception is thrown if the file containing the data for any  *
// *  correlator cannot be found.   The correlators in             *
// *  "corrSetNoSym" must all be found; for each correlator in     *
// *  "corrSetSym", either the correlator OR its time flipped      *
// *  correlator must be found.                                    *
// *                                                               *
// *                                                               *
// *****************************************************************


// *****************************************************************
// *                                                               *
// *  Implementation notes:                                        *
// *                                                               *
// *   -  m_getter is a pointer to a "LapHDataGetHandlerMF" object *
// *      which handles all data retrieval.  The variable members  *
// *      "maxgetopen" and "cleanfrac" control memory limitations  *
// *      and garbage collection.                                  *
// *                                                               *
// *****************************************************************




class BLCorrelatorDataHandler
{

 public:

   class RecordKey
   {
      int time_index;
      int config_serial_index;
   
    public:
   
      RecordKey(int intime, int index) 
           : time_index(intime), config_serial_index(index) {}

      RecordKey(const RecordKey& in) 
           : time_index(in.time_index), 
             config_serial_index(in.config_serial_index) {}

      RecordKey& operator=(const RecordKey& in) 
       {time_index=in.time_index; 
        config_serial_index=in.config_serial_index;
        return *this;}

      ~RecordKey() {}
   
      bool operator<(const RecordKey& rhs) const 
       {return multiLessThan(time_index,rhs.time_index,
               config_serial_index,rhs.config_serial_index);}

      bool operator==(const RecordKey& rhs) const 
       {return multiEqual(time_index,rhs.time_index,
               config_serial_index,rhs.config_serial_index);}

      bool operator!=(const RecordKey& rhs) const 
       {return multiNotEqual(time_index,rhs.time_index,
               config_serial_index,rhs.config_serial_index);}
   
      int getTimeIndex() const {return time_index;}

      int getConfigSerialIndex() const {return config_serial_index;}
  
      void output(XMLHandler& xmlw) const 
      {xmlw.set_root("Key");
       xmlw.put_child("TimeIndex",make_string(time_index));
       xmlw.put_child("ConfigSerialIndex",make_string(config_serial_index));}

      std::string output(int indent=0) const
       {XMLHandler xmlout; output(xmlout);
        return xmlout.output(indent);}

      std::string str() const
       {XMLHandler xmlout; output(xmlout);
        return xmlout.str();}

      explicit RecordKey(const unsigned int* buf) 
       { time_index=buf[0]; config_serial_index=buf[1];}

      static int numints() {return 2;} 

      size_t numbytes() const {return 2*sizeof(unsigned int);}

      void copyTo(unsigned int* buf) const 
       { buf[0]=time_index; buf[1]=config_serial_index; }


   };


   typedef CorrelatorInfo        FileKey;
   typedef CorrelatorAtTimeInfo  MCDataKey;
   typedef Array<InScalar>       DataType;

   LapHDataGetHandlerMF<FileKey,RecordKey,DataType> *m_getter;
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
   BLCorrelatorDataHandler() = delete;
   BLCorrelatorDataHandler(const BLCorrelatorDataHandler&) = delete;
   BLCorrelatorDataHandler& operator=(const BLCorrelatorDataHandler&) = delete;
#else
   BLCorrelatorDataHandler();
   BLCorrelatorDataHandler(const BLCorrelatorDataHandler&);
   BLCorrelatorDataHandler& operator=(const BLCorrelatorDataHandler&);
#endif


 public:


   BLCorrelatorDataHandler(const std::list<FileListInfo>& inputfiles,
                           const std::set<CorrelatorInfo>& corrSetNoSym,
                           const std::set<CorrelatorInfo>& corrSetSym,
                           const MCEnsembleInfo *ensemble,
                           bool use_checksums=false);

   ~BLCorrelatorDataHandler();

   std::set<CorrelatorInfo> getCorrelatorSet() const;

   MCEnsembleInfo getEnsemble() const;

   std::string getEnsembleId() const;

   unsigned int getNumberOfMeasurements();

   void getData(const MCDataKey& mckey, int serial_index, 
                Scalar& data);

   void getSymData(const MCDataKey& mckey, int serial_index, 
                   Scalar& data);


   bool queryData(const MCDataKey& mckey, int serial_index);

   bool querySymData(const MCDataKey& mckey, int serial_index);

   void close();


   void getFileMap(XMLHandler& xmlout) const;

   std::set<CorrelatorInfo> getFileKeys() const;

   std::set<RecordKey> getKeys(const CorrelatorInfo& fkey);

   std::string getFileName(const CorrelatorInfo& fkey);

   void outputKeys(XMLHandler& xmlout);

};


// ***************************************************************

inline size_t numbytes(IOHandler& ioh, const BLCorrelatorDataHandler::RecordKey& rkey)
{
 return rkey.numbytes();
}

// ***************************************************************
}
#endif  

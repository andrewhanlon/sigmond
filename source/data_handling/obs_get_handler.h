#ifndef OBS_GET_HANDLER_H
#define OBS_GET_HANDLER_H

#include "mcobs_info.h"
#include "ensemble_info.h"
#include "corr_data_handler.h"
#include "vev_data_handler.h"
#include "data_handler.h"
#include "bins_info.h"
#include "sampling_info.h"

// ******************************************************************
// *                                                                *
// *  "MCObsGetHandler" handles access to all Monte Carlo           *
// *  observables stored in files, including the LapH correlators   *
// *  and vacuum expectation values of hadronic operators, as well  *
// *  as fit parameters and bin files.  It is a mid-level routine:  *
// *  it is meant to be used by the higher level "MCObsHandler"     *
// *  and it uses the lower-level "BLCorrelatorDataHandler",        *
// *  "BLVEVDataHandler", "BinsGetHandler", and                     *
// *  "SamplingsGetHandler" handlers.  This handler keeps track     *
// *  of which Monte Carlo observables are available and reads      *
// *  them from file when requested.  Results are not stored in     *
// *  memory: the higher level "MCObsHandler" is responsible for    *
// *  storing bins and samplings in memory.                         *
// *                                                                *
// *  Four important member routines are                            *
// *                                                                *
// *     void getBins(const MCObsInfo& obsinfo, RVector& bins);     *
// *     void getSamplings(const MCObsInfo& obsinfo, RVector& samp);*
// *     bool queryBins(const MCObsInfo& obsinfo);                  *
// *     bool querySamplings(const MCObsInfo& obsinfo);             *
// *                                                                *
// *  The "get" routines above throw an exception if the observable *
// *  is not available.  Basic LapH data is stored as complex       *
// *  numbers, so another useful routine is                         *
// *                                                                *
// *     void getBinsComplex(const MCObsInfo& obsinfo,              *
// *               RVector& bins_re, RVector& bins_im);             *
// *                                                                *
// *  The <Arg> tag is ignored in "obsinfo", and the real and       *
// *  imaginary parts are returned in "bins_re" and "bins_im",      *
// *  respectively.  If not available, the results are returned     *
// *  in empty structures, and no exception is thrown.              *
// *                                                                *
// *  For correlator matrix elements of a Hermitian matrix, a       *
// *  symmetrized "get" is done.                                    *
// *                                                                *
// *                                                                *
// *  The constructor takes an XMLHandler as input parameter, as    *
// *  well as an "MCBinsInfo" object and an "MCSamplingInfo"        *
// *  object.  The rebinning and omission of configurations in      *
// *  "MCBinsInfo" must match the header information in the files.  *
// *  Similarly for sampling data, the precise resampling used      *
// *  as specified in the "MCSamplingInfo" object must match the    *
// *  header information in the sampling files.  The XML needed     *
// *  for the constructor must have the form                        *
// *                                                                *
// *     <MCObservables>                                            *
// *                                                                *
// *       <BLCorrelatorData>                                       *
// *         <FileListInfo>...</FileListInfo>                       *
// *         <FileListInfo>...</FileListInfo>     <-- data files    *
// *             ....                                               *
// *       </BLCorrelatorData>                                      *
// *                                                                *
// *       <BLVEVData>                                              *
// *         <FileListInfo>...</FileListInfo>                       *
// *         <FileListInfo>...</FileListInfo>     <-- data files    *
// *             ....                                               *
// *       </BLVEVData>                                             *
// *                                                                *
// *       <BinData>                                                *
// *          <FileName>...</FileName>                              *
// *             ....                                               *
// *       </BinData>                                               *
// *                                                                *
// *       <SamplingData>                                           *
// *          <FileName>...</FileName>                              *
// *             ....                                               *
// *       </SamplingData>                                          *
// *                                                                *
// *       <UseCheckSums/>                         (optional)       *
// *                                                                *
// *       <Specifications>                                         *
// *             ...                                                *
// *           specifications of observables (optional)             *
// *             ...                                                *
// *       </Specifications>                                        *
// *     </MCObservables>                                           *
// *                                                                *
// *  Files containing data for the correlators to be analyzed      *
// *  are specified in terms of <FileListInfo> tags inside a        *
// *  <CorrelatorData>.  Similarly, files containing data for       *
// *  any vacuum expectation values must be specified in terms      *
// *  <FileListInfo> tags inside a <VEVData> tag.                   *
// *                                                                *
// *  If no observables are specified, then all observables found   *
// *  while reading the files will be used.  If any observables     *
// *  ARE specified, then only those observables will be            *
// *  considered for input: files corresponding to other            *
// *  observables will be ignored, and an exception is thrown if    *
// *  the file containing the data for any requested observable     *
// *  cannot be found.                                              *
// *                                                                *
// *  Observables can be specified inside a <Specifications>        *
// *  tag as follows:                                               *
// *                                                                *
// *   (a) a single correlator (no VEV subtraction)                 *
// *                                                                *
// *       <Correlator>...</Correlator>                             *
// *                                                                *
// *   (b) a single correlator with VEV subtraction                 *
// *                                                                *
// *       <CorrelatorWithVEV>...</CorrelatorWithVEV>               *
// *                                                                *
// *   (c) a single VEV                                             *
// *                                                                *
// *       <VEV>...</VEV>                                           *
// *                                                                *
// *   (d) a Hermitian matrix of correlators (no VEV subtraction)   *
// *                                                                *
// *       <HermitianCorrelationMatrix>                             *
// *           <Operator>...</Operator>                             *
// *           <Operator>...</Operator>                             *
// *                ...                                             *
// *         <AssignName>cormat1</AssignName> (optional)            *
// *       </HermitianCorrelationMatrix>                            *
// *                                                                *
// *   (e) a Hermitian matrix of correlators with VEV subtraction   *
// *                                                                *
// *       <HermitianCorrelationMatrixWithVEVs>                     *
// *           <Operator>...</Operator>                             *
// *           <OperatorString>...</OperatorString>                 *
// *                ...                                             *
// *         <AssignName>cormat2</AssignName> (optional)            *
// *       </HermitianCorrelationMatrixWithVEVs>                    *
// *                                                                *
// *   (f) matrix of correlators (no Hermiticity, no VEV subtract)  *
// *                                                                *
// *       <CorrelationMatrix>                                      *
// *           <Operator>...</Operator>                             *
// *           <Operator>...</Operator>                             *
// *                ...                                             *
// *         <AssignName>cormat3</AssignName> (optional)            *
// *       </CorrelationMatrix>                                     *
// *                                                                *
// *   (g) matrix of correlators with VEV subtract (no Hermiticity) *
// *                                                                *
// *       <CorrelationMatrixWithVEVs>                              *
// *           <Operator>...</Operator>                             *
// *           <Operator>...</Operator>                             *
// *                ...                                             *
// *         <AssignName>cormat4</AssignName> (optional)            *
// *       </CorrelationMatrixWithVEVs>                             *
// *                                                                *
// *   (h) set of MCObsInfo in bin files (non basic LapH)           *
// *                                                                *
// *       <ObsBins>                                                *
// *           <MCObservable>...</MCObservable>                     *
// *           <MCObservable>...</MCObservable>                     *
// *                ...                                             *
// *       </ObsBins>                                               *
// *                                                                *
// *   (i) set of MCObsInfo in sampling files                       *
// *                                                                *
// *       <ObsSamplings>                                           *
// *           <MCObservable>...</MCObservable>                     *
// *           <MCObservable>...</MCObservable>                     *
// *                ...                                             *
// *       </ObsSamplings>                                          *
// *                                                                *
// *                                                                *
// *  If a tag <AllHermitian/> appears inside the <Specifications>  *
// *  tag, then all correlation matrices are treated as hermitian.  *
// *                                                                *
// *  Both bin and sampling files can be added and removed later    *
// *  using the "connect" and "disconnect" members below.  However, *
// *  basic LapH data files cannot; they must be fully specified    *
// *  in the constructor.                                           *
// *                                                                *
// *                                                                *
// ******************************************************************



// *******************************************************************
// *                                                                 *
// *  Implementation notes:                                          *
// *                                                                 *
// *  An object of this "MCObsGetHandler" class maintains a pointer  *
// *  to a "BLCorrelatorDataHandler" object, a pointer to a          *
// *  "BLVEVDataHandler" object, a pointer to a "BinsGetHandler"     *
// *  object, and a pointer to a "SamplingsGetHandler".  These       *
// *  objects perform the actual data reading from files but they do *
// *  NOT store data in memory.  "MCObsGetHandler" objects also do   *
// *  NOT store any data in memory.  Instead, the higher level       *
// *  "MCObsHandler" object maintains the input data in memory, in a *
// *  format that facilitates evaluating mean values, covariances,   *
// *  bootstrap errors, and so on.  The class "MCObsGetHandler" is   *
// *  meant to be used as a data handler for an object of the        *
// *  "MCObsHandler" class, with "MCObsInfo" as the record key type. *
// *                                                                 *
// *  Note that basic LapH data are stored as complex numbers, but   *
// *  objects of class "MCObsInfo" store the complex argument,       *
// *  that is, it includes the real or imaginary part separately.    *
// *  Hence, a "getBinsComplex" is provided for efficiency since     *
// *  the real and imaginary parts must be read anyways.             *
// *                                                                 *
// *                                                                 *
// *******************************************************************


class MCObsGetHandler
{

   LaphEnv::BLCorrelatorDataHandler *m_corrdh;
   LaphEnv::BLVEVDataHandler *m_vevdh;
   BinsGetHandler *m_binsdh;
   SamplingsGetHandler *m_sampsdh;

   MCBinsInfo m_bins_info;
   MCSamplingInfo m_sampling_info;
   bool m_use_checksums;

       // Prevent copying ... handler might contain large
       // amounts of data

#ifndef NO_CXX11
   MCObsGetHandler() = delete;
   MCObsGetHandler(const MCObsGetHandler&) = delete;
   MCObsGetHandler& operator=(const MCObsGetHandler&) = delete;
#else
   MCObsGetHandler();
   MCObsGetHandler(const MCObsGetHandler&);
   MCObsGetHandler& operator=(const MCObsGetHandler&);
#endif


 public:


   MCObsGetHandler(XMLHandler& xml_in, const MCBinsInfo& bins_info, 
                   const MCSamplingInfo& samp_info);

   ~MCObsGetHandler();


   void connectBinsFile(const std::string& file_name);

   void connectBinsFile(const std::string& file_name, 
                       const std::set<MCObsInfo>& keys_to_keep);

   void disconnectBinsFile(const std::string& file_name);


   void connectSamplingsFile(const std::string& file_name);

   void connectSamplingsFile(const std::string& file_name, 
                       const std::set<MCObsInfo>& keys_to_keep);

   void disconnectSamplingsFile(const std::string& file_name);




   std::string getEnsembleId() const;

   MCEnsembleInfo getEnsembleInfo() const;

   unsigned int getNumberOfMeasurements() const;

   unsigned int getNumberOfDefaultResamplings() const;

   SamplingMode getDefaultSamplingMode() const;

   const MCBinsInfo& getBinsInfo() const;

   const MCSamplingInfo& getSamplingInfo() const;

   bool useCheckSums() const;


           // ignores real vs imag part of "obsinfo",
           // but if Hermitian, does a symmetrized get

   void getData(const MCObsInfo& obsinfo, int serial_index, 
                Scalar& data);

           // ignores real vs imag part of "obsinfo"
           // but if Hermitian, does a symmetrized query

   bool queryData(const MCObsInfo& obsinfo, int serial_index);

   void close();


   void getFileMap(XMLHandler& xmlout) const;

   std::set<OperatorInfo> getVEVInfos() const;

   std::set<CorrelatorInfo> getCorrelatorInfos() const;



   void getBins(const MCObsInfo& obsinfo, RVector& bins);

   void getSamplings(const MCObsInfo& obsinfo, RVector& samp);

   bool getSamplingsMaybe(const MCObsInfo& obsinfo, RVector& samp);

   bool queryBins(const MCObsInfo& obsinfo);

   bool querySamplings(const MCObsInfo& obsinfo);


#ifdef COMPLEXNUMBERS
   void getBinsComplex(const MCObsInfo& obsinfo, RVector& bins_re, 
                       RVector& bins_im);
#endif


 private:

 
   void clear();

   void setup_correlators(XMLHandler& xmlin, 
                          const std::string& tagname, bool vevs, 
                          std::set<CorrelatorInfo>& corrSet,
                          std::set<OperatorInfo>& vevSet);
 
   void setup_vevs(XMLHandler& xmls, const std::string& tagname,
                   std::set<OperatorInfo>& vevSet);

   void setup_correlator_matrices(XMLHandler& xmlin,
                          const std::string& tagname,
                          bool hermitian, bool vevs, 
                          std::set<CorrelatorInfo>& corrSet,
                          std::set<OperatorInfo>& vevSet);

   void setup_obsset(XMLHandler& xmlin, const std::string& tagname, 
                     std::set<MCObsInfo>& obsSet);



   class BasicLapHGetter
   {
    BasicLapHGetter(){}
    BasicLapHGetter(const BasicLapHGetter&);
    BasicLapHGetter& operator=(const BasicLapHGetter&);
    virtual ~BasicLapHGetter(){}
    virtual void getData(int serial_index, Scalar& result) = 0;
    virtual bool queryData(int serial_index) = 0;
    friend class MCObsGetHandler;
   };

   class BasicLapHCorrGetter: private BasicLapHGetter
   {
    CorrelatorAtTimeInfo m_corr_info;
    LaphEnv::BLCorrelatorDataHandler *m_get;
    BasicLapHCorrGetter(const MCObsInfo& obskey, 
            LaphEnv::BLCorrelatorDataHandler *getptr) 
          : m_corr_info(obskey.getCorrelatorAtTimeInfo()), m_get(getptr) {}
    BasicLapHCorrGetter(const BasicLapHCorrGetter&);
    BasicLapHCorrGetter& operator=(const BasicLapHCorrGetter&);
    virtual ~BasicLapHCorrGetter(){}
    virtual void getData(int serial_index, Scalar& result)
     {m_get->getData(m_corr_info,serial_index,result);}
    virtual bool queryData(int serial_index)
     {return m_get->queryData(m_corr_info,serial_index);}
    friend class MCObsGetHandler;
   };

   class BasicLapHCorrSymGetter: private BasicLapHGetter
   {
    CorrelatorAtTimeInfo m_corr_info;
    LaphEnv::BLCorrelatorDataHandler *m_get;
    BasicLapHCorrSymGetter(const MCObsInfo& obskey, 
            LaphEnv::BLCorrelatorDataHandler *getptr) 
          : m_corr_info(obskey.getCorrelatorAtTimeInfo()), m_get(getptr) {}
    BasicLapHCorrSymGetter(const BasicLapHCorrSymGetter&);
    BasicLapHCorrSymGetter& operator=(const BasicLapHCorrSymGetter&);
    virtual ~BasicLapHCorrSymGetter(){}
    virtual void getData(int serial_index, Scalar& result)
     {m_get->getSymData(m_corr_info,serial_index,result);}
    virtual bool queryData(int serial_index)
     {return m_get->querySymData(m_corr_info,serial_index);}
    friend class MCObsGetHandler;
   };

   class BasicLapHVEVGetter: private BasicLapHGetter
   {
    OperatorInfo m_op_info;
    LaphEnv::BLVEVDataHandler *m_get;
    BasicLapHVEVGetter(const MCObsInfo& obskey, 
            LaphEnv::BLVEVDataHandler *getptr) 
          : m_op_info(obskey.getVEVInfo()), m_get(getptr) {}
    BasicLapHVEVGetter(const BasicLapHVEVGetter&);
    BasicLapHVEVGetter& operator=(const BasicLapHVEVGetter&);
    virtual ~BasicLapHVEVGetter(){}
    virtual void getData(int serial_index, Scalar& result)
     {m_get->getData(m_op_info,serial_index,result);}
    virtual bool queryData(int serial_index)
     {return m_get->queryData(m_op_info,serial_index);}
    friend class MCObsGetHandler;
   };



   void read_data(BasicLapHGetter& getter, Vector<Scalar>& result);

   bool query_data(BasicLapHGetter& getter);


#ifdef COMPLEXNUMBERS
   void get_data(MCObsGetHandler::BasicLapHGetter& getter,
                 RVector& results_re, RVector& results_im);
#else
   void get_data(MCObsGetHandler::BasicLapHGetter& getter,
                 RVector& results);
#endif
};


// ***************************************************************
#endif  

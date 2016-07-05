#ifndef OBS_GET_HANDLER_H
#define OBS_GET_HANDLER_H

#include "mcobs_info.h"
#include "ensemble_info.h"
#include "corr_data_handler.h"
#include "vev_data_handler.h"

namespace LaphEnv {


// ******************************************************************
// *                                                                *
// *  "MCObsGetHandler" handles access to the LapH correlators and  *
// *  vacuum expectation values of hadronic operators. It is meant  *
// *  to handle all simple observables, defined as any MC           *
// *  integrand function, a quantity that can be defined on a       *
// *  single gauge configuration.  It is a mid-level routine:       *
// *  it is meant to be used by the higher level "MCObsHandler"     *
// *  and it uses the lower-level "CorrelatorDataHandler" and       *
// *  "VEVDataHandler".  The "getData" member function is the       *
// *  main useful function. This handles only "standard"            *
// *  observables (see "mcobs_info.h").                             *
// *                                                                *
// *  The constructor takes an XMLHandler as input parameter.       *
// *  The XML needed for the constructor must have the form         *
// *                                                                *
// *     <MCObservables>                                            *
// *       <CorrelatorData>                                         *
// *         <FileListInfo>...</FileListInfo>                       *
// *         <FileListInfo>...</FileListInfo>     <-- data files    *
// *             ....                                               *
// *       </CorrelatorData>                                        *
// *       <VEVData>                                                *
// *         <FileListInfo>...</FileListInfo>                       *
// *         <FileListInfo>...</FileListInfo>     <-- data files    *
// *             ....                                               *
// *       </VEVData>                                               *
// *       <MCEnsembleInfo>....</MCEnsembleInfo>(optional with data)*
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
// *  If the tag <MCEnsembleInfo> is present, the constructor       *
// *  ensures that all data read corresponds to the requested       *
// *  ensemble.  If this tag is absent, the ensemble information    *
// *  is obtained from the data input files, but all input files    *
// *  must correspond to the same ensemble.                         *
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
// *  If a tag <AllHermitian/> appears inside the <Specifications>  *
// *  tag, then all correlation matrices are treated as hermitian.  *
// *                                                                *
// *                                                                *
// ******************************************************************



// *******************************************************************
// *                                                                 *
// *  Implementation notes:                                          *
// *                                                                 *
// *  An object of this "MCObsGetHandler" class maintains a pointer  *
// *  to a "CorrelatorDataHandler" object and a pointer to a         *
// *  "VEVDataHandler" object.  These objects perform the actual     *
// *  data reading from files but they do NOT store data in memory.  *
// *  "MCObsGetHandler" objects also do NOT store any data in        *
// *  memory.  Instead, the higher level "MCObsHandler" object       *
// *  maintains the input data in memory, in a format that           *
// *  facilitates evaluating mean values, covariances, bootstrap     *
// *  errors, and so on.  The class "MCObsGetHandler" is meant to    *
// *  be used as a data handler for an object of the "MCObsHandler"  *
// *  class, with "MCObsInfo" as the record key type.  Recall that   *
// *  objects of class "MCObsInfo" store the complex argument,       *
// *  that is, it includes the real or imaginary part separately.    *
// *  Here, the data returned is complex, so in "getData", the       *
// *  argument part of the "MCObsInfo" key is ignored.  Instead,     *
// *  this is handled by the "MCObsHandler" class.                   *
// *                                                                 *
// *                                                                 *
// *******************************************************************




class MCObsGetHandler
{

   CorrelatorDataHandler *m_corrdh;
   VEVDataHandler *m_vevdh;
   const MCEnsembleInfo *m_ensembleptr;

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


   MCObsGetHandler(XMLHandler& xml_in);

   ~MCObsGetHandler();

   std::string getEnsembleId() const;

   MCEnsembleInfo getEnsembleInfo() const;

   unsigned int getNumberOfMeasurements();

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


 private:


   void read_correlators(XMLHandler& xmlin, 
                         const std::string& tagname, bool vevs, 
                         std::set<CorrelatorInfo>& corrSet,
                         std::set<OperatorInfo>& vevSet);

   void read_vevs(XMLHandler& xmls, const std::string& tagname,
                  std::set<OperatorInfo>& vevSet);

   void read_correlator_matrices(XMLHandler& xmlin,
                          const std::string& tagname,
                          bool hermitian, bool vevs, 
                          std::set<CorrelatorInfo>& corrSet,
                          std::set<OperatorInfo>& vevSet);

};


// ***************************************************************
}
#endif  

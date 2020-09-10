#ifndef MCENSEMBLE_INFO_H
#define MCENSEMBLE_INFO_H

#include "xml_handler.h"

typedef unsigned int  uint;

// *******************************************************************
// *                                                                 *
// *  "MCEnsembleInfo" stores information about an ensemble of       *
// *  Monte Carlo measurements.   If the ensemble is known to        *
// *  SigMonD, the required XML input for constructing the info      *
// *  is very simple:                                                *
// *                                                                 *
// *   <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo>   *
// *                                                                 *
// *  Known ensembles are stored in a file, usually named            *
// *  "ensembles.xml".  A default name for this file is compiled     *
// *  into SigMonD using the variable                                *
// *                                                                 *
// *      std::string MCEnsembleInfo::m_known_ensembles_filename     *
// *                      =DEFAULTENSFILE;                           *
// *                                                                 *
// *  This file must have information specified in the               *
// *  following XML format:                                          *
// *                                                                 *
// *   <KnownEnsembles>                                              *
// *     <Infos>                                                     *
// *       <EnsembleInfo>...</EnsembleInfo>                          *
// *       <EnsembleInfo>...</EnsembleInfo>                          *
// *        ....                                                     *
// *     </Infos>                                                    *
// *     <CLSEnsembleWeights>                                        *
// *       <Ensemble>...</Ensemble>                                  *
// *        ....                                                     *
// *     </CLSEnsembleWeights>                                       *
// *   </KnownEnsembles>                                             *
// *                                                                 *
// *  with each ensemble in the <Infos> tags specified by            *
// *                                                                 *
// *   <EnsembleInfo>                                                *
// *      <Id>clover_s24_t128_ud840_s743</Id>                        *
// *      <NStreams>4</NStreams>                                     *
// *      <NMeas>551</NMeas>                                         *
// *      <NSpace>24</NSpace>                                        *
// *      <NTime>128</NTime>                                         *
// *      <Weighted/>  (if has CLS weights; omit otherwise)          *
// *   </EnsembleInfo>                                               *
// *                                                                 *
// *  The entries in the <CLSEnsembleWeights> tag must have the      *
// *  form:                                                          *
// *                                                                 *
// *   <Ensemble>                                                    *
// *      <Id>cls21_D200_r000</Id>                                   *
// *      <Weights> 0.999 0.998 ... </Weights>                       *
// *   </Ensemble>                                                   *
// *                                                                 *
// *  From "ensembles.xml", this info class knows how many           *
// *  Markov-chain streams are available, how many RHMC trajectory   *
// *  numbers are valid, and so on.                                  *
// *                                                                 *
// *  To specify an ensemble not known to SigMonD, you must provide  *
// *  the following information:                                     *
// *                                                                 *
// *       n_meas   n_streams  n_x  n_y  n_z  n_t                    *
// *                                                                 *
// *  This is done by appending to the id string fields that give    *
// *  this information in the order above, separated by "|".         *
// *  Example:                                                       *
// *                                                                 *
// *  <MCEnsembleInfo>idname|800|1|24|24|24|36</MCEnsembleInfo>      *
// *                                                                 *
// *  With this method, weights cannot be used.                      *
// *                                                                 *
// *******************************************************************


class MCEnsembleInfo
{

   std::string m_id;
   uint n_meas, n_streams, n_x, n_y, n_z, n_t;
   bool m_is_weighted;

   static std::string m_known_ensembles_filename;

#ifndef NO_CXX11
   MCEnsembleInfo() = delete;
#else
   MCEnsembleInfo();
#endif

 public:


   MCEnsembleInfo(XMLHandler& xml_in);

   MCEnsembleInfo(const std::string& id);

   MCEnsembleInfo(const std::string& id, const std::string& in_ensembles_filename);

   MCEnsembleInfo(const std::string& id, uint num_meas, uint num_streams, 
                  uint nx, uint ny, uint nz, uint nt);

   MCEnsembleInfo(const MCEnsembleInfo& fin);

   MCEnsembleInfo& operator=(const MCEnsembleInfo& fin);

   ~MCEnsembleInfo(){}



   std::string getId() const { return m_id; }

   uint getNumberOfMeasurements() const { return n_meas; }

   uint getNumberOfStreams() const { return n_streams; }

   bool isWeighted() const { return m_is_weighted; }   

   void getWeights(std::vector<double>& weights) const;
  
   void checkEqual(const MCEnsembleInfo& rhs) const;

   bool operator==(const MCEnsembleInfo& rhs) const {return (m_id==rhs.m_id); }

   bool operator!=(const MCEnsembleInfo& rhs) const {return (m_id!=rhs.m_id); }

   void output(XMLHandler& xmlout) const;

   std::string output(int indent=0) const;

   std::string str() const;

   uint getLatticeTimeExtent() const {return n_t;}

   uint getLatticeXExtent() const {return n_x;}

   uint getLatticeYExtent() const {return n_y;}

   uint getLatticeZExtent() const {return n_z;}

 private:

   void initialize();

   bool parse(const std::string& idstr);

   friend class TaskHandler;

};


// ***************************************************************
#endif  

#ifndef MCENSEMBLE_INFO_H
#define MCENSEMBLE_INFO_H

#include "xml_handler.h"

// *******************************************************************
// *                                                                 *
// *  "MCEnsembleInfo" stores information about a known ensemble     *
// *  of Monte Carlo measurements.   The required XML input for      *
// *  constructing the info is very simple:                          *
// *                                                                 *
// *   <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo>   *
// *                                                                 *
// *  The following identifying strings are currently supported:     *
// *     clover_s24_t128_ud840_s743                                  *
// *     clover_s24_t128_ud860_s743                                  *
// *     clover_s32_t256_ud860_s743                                  *
// *     clover_s16_t128_ud840_s743                                  *
// *                                                                 *
// *  Given the above identifying string, this info class knows      *
// *  how many Markov-chain streams are available, how many          *
// *  RHMC trajectory numbers are valid, and so on.                  *
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
// *  <MCEnsembleInfo>idname|800|1|24|24|24|36</MCEnsembleInfo>      *                                                         *
// *                                                                 *
// *******************************************************************


class MCEnsembleInfo
{

   std::string m_id;
   uint n_meas, n_streams, n_x, n_y, n_z, n_t;

#ifndef NO_CXX11
   MCEnsembleInfo() = delete;
#else
   MCEnsembleInfo();
#endif

 public:


   MCEnsembleInfo(XMLHandler& xml_in);

   MCEnsembleInfo(const std::string& id);

   MCEnsembleInfo(const std::string& id, uint num_meas, uint num_streams, 
                  uint nx, uint ny, uint nz, uint nt);

   MCEnsembleInfo(const MCEnsembleInfo& fin);

   MCEnsembleInfo& operator=(const MCEnsembleInfo& fin);

   ~MCEnsembleInfo(){}



   std::string getId() const { return m_id; }

   uint getNumberOfMeasurements() const { return n_meas; }

   uint getNumberOfStreams() const { return n_streams; }
  
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

};


// ***************************************************************
#endif  

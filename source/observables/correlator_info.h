#ifndef CORRELATOR_INFO_H
#define CORRELATOR_INFO_H

#include <map>
#include <vector>
#include "operator_info.h"


class MCObsInfo;

// *******************************************************************
// *                                                                 *
// *   The classes "CorrelatorInfo" and "CorrelatorAtTimeInfo" are   *
// *   defined in this file.  Objects of class "CorrelatorInfo"      *
// *   store identifying information about one particular temporal   *
// *   correlation between a source and a sink operator, and         *
// *   objects of class "CorrelatorAtTimeInfo" store information     *
// *   about a temporal correlator for one time separation.          *
// *   "CorrelatorInfo" construction by XML content requires XML in  *
// *   the following format:                                         *
// *                                                                 *
// *       <Correlator>                                              *
// *         <Source>                                                *
// *            <Operator>..</Operator>                              *
// *                or <OperatorString>..</OperatorString>           *
// *         </Source>                                               *
// *         <Sink>                                                  *
// *            <Operator>..</Operator>                              *
// *                or <OperatorString>..</OperatorString>           *
// *         </Sink>                                                 *
// *       </Correlator>                                             *
// *                                                                 *
// *   "CorrelatorAtTimeInfo" construction by XML content requires   *
// *   XML in the following format:                                  *
// *                                                                 *
// *       <Correlator>                                              *
// *         <Source>                                                *
// *            <Operator>..</Operator>                              *
// *                or <OperatorString>..</OperatorString>           *
// *         </Source>                                               *
// *         <Sink>                                                  *
// *            <Operator>..</Operator>                              *
// *                or <OperatorString>..</OperatorString>           *
// *         </Sink>                                                 *
// *         <TimeIndex>..</TimeIndex>                               *
// *         <HermitianMatrix/>    (optional)                        *
// *         <SubtractVEV/>    (optional)                            *
// *         <Reweight/>    (optional)                               *
// *       </Correlator>                                             *
// *                                                                 *
// *   If the <HermitianMatrix\> tag is given, this facilitates      *
// *   input which automatically averages using the complex          *
// *   conjugate elements, if available.  If the <SubtractVEV/>      *
// *   tag is present, then the observable is no longer simple.      *
// *   Of course, a correlator is a VEV of a product of two          *
// *   operators at possibly different times; however, subtracting   *
// *   the VEV here refers to the VEVs of the individual operators.  *
// *                                                                 *
// *   See the file "operator_info.h" for more information about     *
// *   the XML format needed inside the <Operator> or                *
// *   <OperatorString> tags.                                        *
// *                                                                 * 
// *******************************************************************


class CorrelatorInfo
{

   std::vector<unsigned int> icode;

#ifndef NO_CXX11
    CorrelatorInfo() = delete;
#else
    CorrelatorInfo();
#endif


 public:

   CorrelatorInfo(XMLHandler& xml_in);

   CorrelatorInfo(const OperatorInfo& sink, const OperatorInfo& source);

   CorrelatorInfo(const CorrelatorInfo& Cor) : icode(Cor.icode) {}

   CorrelatorInfo& operator=(const CorrelatorInfo& Cor)
    {icode=Cor.icode; return *this;}

   ~CorrelatorInfo(){}

    // output functions
    
   OperatorInfo getSource() const;

   OperatorInfo getSink() const;

   CorrelatorInfo getTimeFlipped() const;

   bool isSinkSourceSame() const;

   
   std::string output(bool longform=false, int indent=0) const;  // XML output 

   std::string str() const;  // XML output 

   void output(XMLHandler& xmlout, bool longform=false) const;  // XML output


   bool operator==(const CorrelatorInfo& rhs) const;

   bool operator!=(const CorrelatorInfo& rhs) const;

   bool operator<(const CorrelatorInfo& rhs) const;


 private:

   void assign(const OperatorInfo& sink, const OperatorInfo& source);

   CorrelatorInfo(const std::vector<unsigned int>& incode) 
          : icode(incode) {}

       // for construction from a CorrelatorAtTimeInfo
   CorrelatorInfo(const std::vector<unsigned int>& incode, uint srcsize) 
          : icode(incode) 
    {icode.back()=srcsize;}

   void interchange_ends(std::vector<unsigned int>& outcode,
                         const std::vector<unsigned int>& incode) const;

   friend class MCObsInfo;
   friend class OperatorInfo;
   friend class CorrelatorAtTimeInfo;
   friend class CorrelatorDataHandler;

};



// **************************************************************


class CorrelatorAtTimeInfo
{

   std::vector<unsigned int> icode;

 public:

   CorrelatorAtTimeInfo(XMLHandler& xml_in);

   CorrelatorAtTimeInfo(const OperatorInfo& sink, const OperatorInfo& source,
                        int timeval, bool hermitianmatrix,
                        bool subtractvev, bool reweight);

   CorrelatorAtTimeInfo(const CorrelatorInfo& corr, int timeval, 
                        bool hermitianmatrix, bool subtractvev,
                        bool reweight);

   CorrelatorAtTimeInfo(const MCObsInfo& obsinfo);

   CorrelatorAtTimeInfo(const CorrelatorAtTimeInfo& Cor) : icode(Cor.icode) {}

   CorrelatorAtTimeInfo& operator=(const CorrelatorAtTimeInfo& Cor)
    {icode=Cor.icode; return *this;}

   CorrelatorAtTimeInfo& resetTimeSeparation(int timeval);

   CorrelatorAtTimeInfo& resetSubtractVEV(bool subvev);

   CorrelatorAtTimeInfo& resetReweight(bool reweight);

   ~CorrelatorAtTimeInfo(){}


    // output functions
    
   CorrelatorInfo getCorrelator() const;

   OperatorInfo getSource() const;

   OperatorInfo getSink() const;

   unsigned int getTimeSeparation() const
    {return icode.back()>>9;}
    
   bool isHermitianMatrix() const
    {return isHermitian(icode);}

   bool subtractVEV() const
    {return (icode.back()&4u);}

   bool reweight() const
    {return (icode.back()&1u);}



   std::string output(bool longform=false, int indent=0) const;  // XML output 

   std::string str() const;  // XML output 

   void output(XMLHandler& xmlout, bool longform=false) const;  // XML output



   bool operator==(const CorrelatorAtTimeInfo& rhs) const;

   bool operator!=(const CorrelatorAtTimeInfo& rhs) const;

   bool operator<(const CorrelatorAtTimeInfo& rhs) const;


 private:

   void assign(const OperatorInfo& sink, const OperatorInfo& source,
               int timeval, bool hermitianmatrix, bool subvev,
               bool reweight);

   void assign(const CorrelatorInfo& corr, int timeval,
               bool hermitianmatrix, bool subvev, bool reweight);

   CorrelatorAtTimeInfo(std::vector<unsigned int>::const_iterator inbegin,
                        std::vector<unsigned int>::const_iterator inend) 
          : icode(inbegin,inend) {}

   void set_time_herm_vev_rew(uint srcsize, int timeval,
                              bool hermitianmatrix, bool subvev,
                              bool reweight);

   static bool isHermitian(const std::vector<unsigned int>& incode)
    {return (incode.back()&2u);}

   uint get_source_size() const
    {return (icode.back()>>3)&63u;}

   friend class MCObsInfo;
   friend class OperatorInfo;
   friend class CorrelatorInfo;
   friend class CorrelatorDataHandler;

};


// **************************************************
#endif  

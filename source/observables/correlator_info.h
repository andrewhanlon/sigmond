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
// *         <HermitianMatrix\>    (optional)                        *
// *         <SubtractVEV\>    (optional)                            *
// *       </Correlator>                                             *
// *                                                                 *
// *   If the <HermitianMatrix\> tag is given, this facilitates      *
// *   input which automatically averages using the complex          *
// *   conjugate elements, if available.  If the <SubtractVEV/>      *
// *   tag is present, then the observable is no longer simple.      *
// *                                                                 *
// *   See the file "operator_info.h" for more information about     *
// *   the XML format needed inside the <Operator> or                *
// *   <OperatorString> tags.                                        *
// *                                                                 * 
// *******************************************************************


class CorrelatorInfo
{

   typedef std::vector<unsigned int> ivec;
   ivec icode;

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

   
   std::string output(bool longform=true, int indent=0) const;  // XML output 

   std::string str() const;  // XML output 

   void output(XMLHandler& xmlout, bool longform=true) const;  // XML output


   bool operator==(const CorrelatorInfo& rhs) const;

   bool operator!=(const CorrelatorInfo& rhs) const;

   bool operator<(const CorrelatorInfo& rhs) const;


 private:

   static void assign(std::vector<unsigned int>& outcode,
                      const OperatorInfo& sink, const OperatorInfo& source);

   CorrelatorInfo(std::vector<unsigned int>::const_iterator inbegin,
                  std::vector<unsigned int>::const_iterator inend) 
          : icode(inbegin,inend) {}

   // CorrelatorInfo() {}   


   static ivec::iterator sourcebegin(ivec& incode)
    {return incode.begin();}

   static ivec::iterator sourceend(ivec& incode)
    {return incode.begin()+OperatorInfo::codesize(incode[0]);}

   static ivec::iterator sinkbegin(ivec& incode)
    {return incode.begin()+OperatorInfo::codesize(incode[0]);}

   static ivec::iterator sinkend(ivec& incode)
    {return incode.end();}

   static ivec::const_iterator sourcebegin(const ivec& incode)
    {return incode.begin();}

   static ivec::const_iterator sourceend(const ivec& incode)
    {return incode.begin()+OperatorInfo::codesize(incode[0]);}

   static ivec::const_iterator sinkbegin(const ivec& incode)
    {return incode.begin()+OperatorInfo::codesize(incode[0]);}

   static ivec::const_iterator sinkend(const ivec& incode)
    {return incode.end();}

   static void interchange_ends(std::vector<unsigned int>& outcode,
                                const std::vector<unsigned int>& incode);

   friend class MCObsInfo;
   friend class OperatorInfo;
   friend class CorrelatorAtTimeInfo;
   friend class CorrelatorDataHandler;

};



// **************************************************************


class CorrelatorAtTimeInfo
{

   typedef std::vector<unsigned int>  ivec;
   ivec icode;

 public:

   CorrelatorAtTimeInfo(XMLHandler& xml_in);

   CorrelatorAtTimeInfo(const OperatorInfo& sink, const OperatorInfo& source,
                        int timeval, bool hermitianmatrix=true,
                        bool subtractvev=false);

   CorrelatorAtTimeInfo(const CorrelatorInfo& corr, int timeval, 
                        bool hermitianmatrix=true, bool subtractvev=false);

   CorrelatorAtTimeInfo(const MCObsInfo& obsinfo);

   CorrelatorAtTimeInfo(const CorrelatorAtTimeInfo& Cor) : icode(Cor.icode) {}

   CorrelatorAtTimeInfo& operator=(const CorrelatorAtTimeInfo& Cor)
    {icode=Cor.icode; return *this;}

   CorrelatorAtTimeInfo& resetTimeSeparation(int timeval);

   CorrelatorAtTimeInfo& resetVEVSubtracted(bool subvev);

   ~CorrelatorAtTimeInfo(){}


    // output functions
    
   CorrelatorInfo getCorrelator() const;

   OperatorInfo getSource() const;

   OperatorInfo getSink() const;

   unsigned int getTimeSeparation() const;
 
   bool isHermitianMatrix() const;

   bool isVEVsubtracted() const;



   std::string output(bool longform=true, int indent=0) const;  // XML output 

   std::string str() const;  // XML output 

   void output(XMLHandler& xmlout, bool longform=true) const;  // XML output



   bool operator==(const CorrelatorAtTimeInfo& rhs) const;

   bool operator!=(const CorrelatorAtTimeInfo& rhs) const;

   bool operator<(const CorrelatorAtTimeInfo& rhs) const;


 private:

   static void assign(std::vector<unsigned int>& outcode,
                      const OperatorInfo& sink, const OperatorInfo& source,
                      int timeval, bool hermitianmatrix, bool subvev);

   static void assign(std::vector<unsigned int>& outcode,
                      const CorrelatorInfo& corr, int timeval,
                      bool hermitianmatrix, bool subvev);

   CorrelatorAtTimeInfo(std::vector<unsigned int>::const_iterator inbegin,
                        std::vector<unsigned int>::const_iterator inend) 
          : icode(inbegin,inend) {}

   static ivec::iterator sourcebegin(ivec& incode)
    {return incode.begin();}

   static ivec::iterator sourceend(ivec& incode)
    {return incode.begin()+OperatorInfo::codesize(incode[0]);}

   static ivec::iterator sinkbegin(ivec& incode)
    {return incode.begin()+OperatorInfo::codesize(incode[0]);}

   static ivec::iterator sinkend(ivec& incode)
    {return incode.end()-1;}

   static ivec::const_iterator sourcebegin(const ivec& incode)
    {return incode.begin();}

   static ivec::const_iterator sourceend(const ivec& incode)
    {return incode.begin()+OperatorInfo::codesize(incode[0]);}

   static ivec::const_iterator sinkbegin(const ivec& incode)
    {return incode.begin()+OperatorInfo::codesize(incode[0]);}

   static ivec::const_iterator sinkend(const ivec& incode)
    {return incode.end()-1;}

   static unsigned int extract_time(const std::vector<unsigned int>& incode)
    {return incode.back()>>2;}

   static bool isHermitian(const std::vector<unsigned int>& incode)
    {return (incode.back()&1u);}

   static bool isVEVsubtracted(const std::vector<unsigned int>& incode)
    {return (incode.back()&2u);}

   static void set_time_herm_vev(std::vector<unsigned int>& outcode, int timeval,
                                 bool hermitianmatrix, bool subvev)
    {if (timeval<0){
       throw(std::invalid_argument("Nonnegative time separation required in CorrelatorAtTimeInfo"));}
     uint tcode=timeval; tcode<<=1;
     if (subvev) tcode|=1u; tcode<<=1;
     if (hermitianmatrix) tcode|=1u;
     outcode.back()=tcode;}

   static void reset_time(std::vector<unsigned int>& outcode, int timeval)
    {if (timeval<0){
       throw(std::invalid_argument("Nonnegative time separation required in CorrelatorAtTimeInfo"));}
     uint tcode=timeval; tcode<<=2;
     tcode|= (outcode.back()&3u);
     outcode.back()=tcode;}

   static void reset_vev(std::vector<unsigned int>& outcode, bool subvev)
    {uint tcode=outcode.back()>>2; tcode<<=1;
     uint hcode=(outcode.back()&1u);
     if (subvev) tcode|=1u; 
     tcode<<=1; tcode|=hcode;
     outcode.back()=tcode;}

   friend class MCObsInfo;
   friend class OperatorInfo;
   friend class CorrelatorInfo;
   friend class CorrelatorDataHandler;

};


// **************************************************
#endif  

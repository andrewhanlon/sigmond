#ifndef MCOBS_INFO_H
#define MCOBS_INFO_H

#include "xml_handler.h"
#include "operator_info.h"
#include "correlator_info.h"
#include "scalar_defs.h"


// ********************************************************************
// *                                                                  *
// *   Objects of class "MCObsInfo" store identifying information     *
// *   about one particular Monte Carlo observable.  Each observable  *
// *   must be associated with a REAL-VALUED quantity that can be     *
// *   estimated by our Monte Carlo path integrals.  These can be     *
// *   simple quantities, such as the real or imaginary part of       *
// *   a temporal correlator for one time separation, that can be     *
// *   defined on a single gauge configuration, or much more          *
// *   complicated quantities, such as a fit parameter yielding a     *
// *   stationary-state energy, determined by fitting a decaying      *
// *   exponential to a temporal correlation function.                *
// *                                                                  *
// *   An observable is termed "simple" if it can be associated with  *
// *   the integrand of a single path integral.  Simple observables   *
// *   include                                                        *
// *    (1) the real or imaginary part of a temporal correlator       *
// *                for one time separation                           *
// *    (2) the real or imaginary part of a vacuum expectation value  *
// *   Other observables are referred to as "nonsimple".              *
// *                                                                  *
// *   The class "MCObsInfo" is meant to encompass all observables    *
// *   of interest.  Observables can be classified as "standard"      *
// *   or "nonstandard":  "standard" refers to VEVs and correlators   *
// *   of LapH operators whose descriptions must match that stored    *
// *   in the files containing their sources and sinks; "nonstandard" *
// *   refers to fit parameters, rotated correlators, and other       *
// *   user-defined observables.                                      *
// *                                                                  *
// *   Standard observables:                                          *
// *                                                                  *
// *     For "standard" observables, the constructor can take         *
// *     an XMLHandler as its single argument, or there are           *
// *     other constructor versions that use other classes, such as   *
// *     "OperatorInfo" and "CorrelatorAtTimeInfo", as arguments.     *
// *     Construction by XML content requires XML in the following    *
// *     format:                                                      *
// *                                                                  *
// *     <MCObservable>                                               *
// *       <VEV>                                                      *
// *          <Operator>..</Operator>                                 *
// *              or <OperatorString>..</OperatorString>              *
// *       </VEV>                                                     *
// *       <Arg>RealPart</Arg> or <Arg>Re</Arg>                       *
// *           or <Arg>ImaginaryPart</Arg> or <Arg>Im</Arg>           *
// *     </MCObservable>                                              *
// *                                                                  *
// *     <MCObservable>                                               *
// *       <Correlator>                                               *
// *         <Source>                                                 *
// *            <Operator>..</Operator>                               *
// *                or <OperatorString>..</OperatorString>            *
// *         </Source>                                                *
// *         <Sink>                                                   *
// *            <Operator>..</Operator>                               *
// *                or <OperatorString>..</OperatorString>            *
// *         </Sink>                                                  *
// *         <TimeIndex>..</TimeIndex>                                *
// *         <HermitianMatrix\>    (optional)                         *
// *       </Correlator>                                              *
// *       <Arg>RealPart</Arg> or <Arg>Re</Arg>                       *
// *           or <Arg>ImaginaryPart</Arg> or <Arg>Im</Arg>           *
// *     </MCObservable>                                              *
// *                                                                  *
// *     If the <HermitianMatrix\> tag is given, this facilitates     *
// *     input which automatically averages using the complex         *
// *     conjugate elements, if available.  If <Arg> is omitted, then *
// *     the real part is assumed.                                    *
// *                                                                  *
// *     See the file "operator_info.h" for more information about    *
// *     the XML format needed inside the <Operator> or               *
// *     <OperatorString> tags.                                       *
// *                                                                  * 
// *                                                                  * 
// *   Nonstandard observables:                                       *
// *                                                                  *
// *     Nonstandard observables are specified using a <ObsName> tag, *
// *     an unsigned integer <Index>, and a possible <Description>    *
// *     tag.  The <Description> tag is used only for outputting      *
// *     XML about the observable, and is NOT used for internally     *
// *     representing the observable.  Only the <ObsName> and <Index> *
// *     tags are important, as well as an optional <Arg> tag         *
// *     (assumed RealPart if absent) and an optional <Simple/>       *
// *     tag (if absent, the observable is assumed to be nonsimple.   *
// *     Once an MCObsInfo is created for a particular <Name>, the    *
// *     <Description> cannot be changed and need not be present in   *
// *     subsequent constructors.                                     *
// *                                                                  *
// *     For "nonstandard" observables, the constructor can take      *
// *     an XMLHandler as its single argument, or a version using     *
// *     the <Name> string, <Index> integer, and so on, is available. *
// *     Construction by XML content requires XML in the              *
// *     following format:                                            *
// *                                                                  *
// *     <MCObservable>                                               *
// *       <ObsName>T1upEnergy</ObsName> (32 char or less, no blanks) *
// *       <Index>3</Index>        (opt nonneg integer: default 0)    *
// *       <Description>....</Description>    (optional)              *
// *       <Simple/>      (if simple observable)                      *
// *       <Arg>RealPart</Arg> or <Arg>Re</Arg>                       *
// *           or <Arg>ImaginaryPart</Arg> or <Arg>Im</Arg>           *
// *     </MCObservable>                                              *
// *                                                                  *
// *                                                                  * 
// ********************************************************************


// ********************************************************************
// *                                                                  *
// *   Implementation notes:                                          *
// *                                                                  *
// *   All observables are encoded in a std::vector<unsigned int>.    *
// *   The first unsigned integer icode[0] contains the following     *
// *   information:                                                   *
// *                                                                  *
// *   Content of icode[0]:                                           *
// *        rightmost bit:  0 --> real part,  1 --> imaginary part    *
// *   next rightmost bit:  0 --> simple,     1 --> nonsimple         *
// *   next rightmost bit:  0 --> standard,   1 --> nonstandard       *
// *    remaining 29 bits:                                            *
// *        if standard:                                              *
// *            0 = Vacuum, 1 = VEV,   2 = CorrelatorAtTimeInfo       *
// *        if nonstandard:                                           *
// *            rightmost 16 bits:  id number in static map           *
// *             leftmost 13 bits:  unsigned integer index            *
// *                                                                  *
// *   For standard observables, the remaining elements of the        *
// *   icode[] vector match those of OperatorInfo and                 *
// *   CorrelatorAtTimeInfo which speeds up encoding keys for         *
// *   file access.                                                   *
// *                                                                  *
// *   For nonstandard observables, a static map "m_encodings" is     *
// *   used for the encoding, with the key given by the content of    *
// *   the <Name> tag, which must be 32 characters or less with no    *
// *   blanks, tabs, or newlines.   "m_encodings" associates an       *
// *   unsigned integer with a name string.  The static vector        *
// *   "m_decodings" contains the name string and the description     *
// *   strings associated with each integer code.  The optional       *
// *   <Description> tag must be valid XML is used only for           *
// *   informational and outputting purposes.  Once an MCObsInfo      *
// *   is created for a particular <Name>, the <Description> cannot   *
// *   be changed and need not be present in subsequent constructors. *
// *                                                                  *
// *                                                                  *
// ********************************************************************


class MCObsInfo
{

   std::vector<unsigned int> icode;

 public:

   MCObsInfo();

   MCObsInfo(XMLHandler& xml_in);

   MCObsInfo(const OperatorInfo& opinfo, ComplexArg arg=RealPart);   // a VEV

   MCObsInfo(const OperatorInfo& sinkop, const OperatorInfo& sourceop, 
             int timeval, bool hermitianmatrix=true, ComplexArg arg=RealPart,
             bool subtractvev=false);  // correlator

   MCObsInfo(const CorrelatorAtTimeInfo& corrinfo, 
             ComplexArg arg=RealPart);

   MCObsInfo(const CorrelatorInfo& corrinfo, int timeval, 
             bool hermitianmatrix=true, ComplexArg arg=RealPart,
             bool subtractvev=false);

   MCObsInfo(const std::string& obsname, uint index=0, bool simple=false,
             ComplexArg arg=RealPart);

   MCObsInfo(const MCObsInfo& B) : icode(B.icode) {}

   MCObsInfo& operator=(const MCObsInfo& B)
    {icode=B.icode; return *this;}

   ~MCObsInfo(){}

   void setToRealPart();

   void setToImaginaryPart();

   void resetObsIndex(uint ind);  // for nonstandard observables


    // output functions
    
   bool isVacuum() const;

   bool isVEV() const;

   bool isCorrelatorAtTime() const;

   bool isHermitianCorrelatorAtTime() const;

   bool isRealPart() const;

   bool isImaginaryPart() const;

   bool isSimple() const;

   bool isNonSimple() const;

   bool isStandard() const;

   bool isNonStandard() const;


   OperatorInfo getVEVInfo() const;

   void getVEVInfo(OperatorInfo& vop) const;

   CorrelatorAtTimeInfo getCorrelatorAtTimeInfo() const;

   void getCorrelatorAtTimeInfo(CorrelatorAtTimeInfo& ctinfo ) const;

   OperatorInfo getCorrelatorSourceInfo() const;

   OperatorInfo getCorrelatorSinkInfo() const;

   unsigned int getCorrelatorTimeIndex() const;

   CorrelatorInfo getCorrelatorInfo() const;

   void getCorrelatorInfo(CorrelatorInfo& cinfo) const;

   std::string getObsName() const;

   uint getObsIndex() const;


   std::string output(bool longform=true, int indent=0) const;  // XML output 

   std::string str() const;  // XML output 

   void output(XMLHandler& xmlout, bool longform=true) const;  // XML output


   bool operator==(const MCObsInfo& rhs) const;

   bool operator!=(const MCObsInfo& rhs) const;

   bool operator<(const MCObsInfo& rhs) const;


 private:

   static std::map<std::string,uint> m_encodings;

   static std::vector<std::pair<std::string,std::string> > m_decodings;


   void encode(const std::vector<uint>& precode, unsigned int optype, 
               bool simple, ComplexArg arg);

   void encode(const std::string& name, uint index, const std::string& description,
               bool simple, ComplexArg arg);

   void set_real_part();

   void set_imag_part();

   void set_index(uint ind);

   void set_arg(ComplexArg arg);

   bool read_arg_type(XMLHandler& xmlin, ComplexArg& arg);

   void assert_corrtype(const std::string& msg="") const;

   void assert_vevtype(const std::string& msg="") const;

   uint get_obs_index() const;

   uint get_obs_code() const;

   friend class CorrelatorAtTimeInfo;

};


// **************************************************
#endif  

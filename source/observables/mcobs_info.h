#ifndef MCOBS_INFO_H
#define MCOBS_INFO_H

#include "xml_handler.h"
#include "operator_info.h"
#include "correlator_info.h"
#include "scalar_defs.h"
#include "io_handler.h"


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
// *                for one time separation (no vev subtraction)      *
// *    (2) the real or imaginary part of a vacuum expectation value  *
// *                of a single operator                              *
// *   Other observables are referred to as "nonsimple".              *
// *                                                                  *
// *   The class "MCObsInfo" is meant to encompass all observables    *
// *   of interest.  Observables can be classified as "primary" or    *
// *   "secondary":  "primary" refers to operator VEVs and            *
// *   correlators of field operators, which can be of type           *
// *   "BasicLapH" or "GenIrrep"; "secondary" refers to fit           *
// *   parameters, and other user-defined observables.                *
// *                                                                  *
// *   "primary" observables:                                         *
// *                                                                  *
// *     For "primary" observables, the constructor can take          *
// *     an XMLHandler as its single argument, or there are           *
// *     other constructor versions that use other classes, such as   *
// *     "OperatorInfo" and "CorrelatorAtTimeInfo", as arguments.     *
// *     Construction by XML content requires XML in the following    *
// *     format:                                                      *
// *                                                                  *
// *     <MCObservable>                                               *
// *       <VEV>                                                      *
// *          <Operator>..</Operator> (or other operator tags)        *
// *          <Reweight/>      (optional)                             *
// *       </VEV>                                                     *
// *       <Arg>RealPart</Arg> or <Arg>Re</Arg>                       *
// *           or <Arg>ImaginaryPart</Arg> or <Arg>Im</Arg>           *
// *     </MCObservable>                                              *
// *                                                                  *
// *     <MCObservable>                                               *
// *       <Correlator>                                               *
// *         <Source>                                                 *
// *            <Operator>..</Operator>  (or other operator tags)     *
// *         </Source>                                                *
// *         <Sink>                                                   *
// *            <Operator>..</Operator>  (or other operator tags)     *
// *         </Sink>                                                  *
// *         <TimeIndex>..</TimeIndex>                                *
// *         <HermitianMatrix/>    (optional)                         *
// *         <SubtractVEV/>     (optional)                            *
// *         <Reweight/>      (optional)                              *
// *       </Correlator>                                              *
// *       <Arg>RealPart</Arg> or <Arg>Re</Arg>                       *
// *           or <Arg>ImaginaryPart</Arg> or <Arg>Im</Arg>           *
// *     </MCObservable>                                              *
// *                                                                  *
// *     If the <HermitianMatrix/> tag is given, this facilitates     *
// *     input which automatically averages using the complex         *
// *     conjugate elements, if available.  If <Arg> is omitted, then *
// *     the real part is assumed.                                    *
// *                                                                  *
// *     See the file "operator_info.h" for more information about    *
// *     the XML format needed inside the <Operator> or               *
// *     <OperatorString> tags.                                       *
// *                                                                  * 
// *                                                                  * 
// *   "secondary" observables:                                       *
// *                                                                  *
// *     "secondary" observables are specified using a <ObsName>      *
// *     tag, unsigned integer <Index>, as well as an optional <Arg>  *
// *     tag (assumed RealPart if absent) and an optional <Simple/>   *
// *     tag (if absent, the observable is assumed to be nonsimple).  *
// *                                                                  *
// *     For these observables, there is a constructor which takes an *
// *     XMLHandler as its single argument, and another version which *
// *     takes the <ObsName> string, <Index> integer, and so on,      *
// *     Construction by XML content requires XML in the              *
// *     following format:                                            *
// *                                                                  *
// *     <MCObservable>                                               *
// *       <ObsName>T1up_Energy</ObsName> (32 char or less, no blanks)*
// *       <Index>3</Index>        (opt nonneg integer: default 0)    *
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
// *   next rightmost bit:  0 --> primary,    1 --> secondary         *
// *    remaining 29 bits:                                            *
// *        if primary:                                               *
// *            0 = Vacuum, 1 = VEV,   2 = CorrelatorAtTimeInfo,      *
// *            3 = ReweightingFactor                                 *
// *        if secondary:                                             *
// *            the unsigned integer index                            *
// *                                                                  *
// *   For "primary" observables, the remaining elements of the       *
// *   icode[] vector match those of OperatorInfo and                 *
// *   CorrelatorAtTimeInfo which speeds up encoding keys for         *
// *   file access.                                                   *
// *                                                                  *
// *   For "secondary" observables, the remaining elements of the     *
// *   icode[] vector contain the maximum 32-character "ObsName"      *
// *   string converted byte-by-byte to integers by ASCII code.       *
// *                                                                  *
// *   If a particular observable is simple, then it can be           *
// *   computed on a single field configuration and stored in bins;   *
// *   a non-simple observable cannot have a single bin file.  If a   *
// *   particular observable is a "primary", then it can be           *
// *   constructed from one or more bin files in a simple known       *
// *   way; if "secondary", then only jackknife and bootstrap         *
// *   files can be used since its construction in terms of several   *
// *   bin files is unknown and perhaps complicated.  Whether or not  *
// *   a "primary" observable is "BasicLapH" or "GenIrrep"            *
// *   determines how the quantity is read from file(s).              *
// *                                                                  *
// *                                                                  *
// ********************************************************************


class MCObsInfo
{

   std::vector<unsigned int> icode;

 public:

   MCObsInfo(bool reweighting_factor=false);

   MCObsInfo(XMLHandler& xml_in);

   MCObsInfo(const OperatorInfo& opinfo, ComplexArg arg=RealPart,
             bool reweight=false);   // an Op VEV

   MCObsInfo(const OperatorInfo& sinkop, const OperatorInfo& sourceop, 
             int timeval, bool hermitianmatrix=true, ComplexArg arg=RealPart,
             bool subtractvev=false, bool reweight=false);  // correlator

   MCObsInfo(const CorrelatorAtTimeInfo& corrinfo, 
             ComplexArg arg=RealPart);

   MCObsInfo(const CorrelatorInfo& corrinfo, int timeval, 
             bool hermitianmatrix=true, ComplexArg arg=RealPart,
             bool subtractvev=false, bool reweight=false);

   MCObsInfo(const std::string& obsname, uint index=0, bool simple=false,
             ComplexArg arg=RealPart);

   MCObsInfo(const MCObsInfo& B) : icode(B.icode) {}

   MCObsInfo& operator=(const MCObsInfo& B)
    {icode=B.icode; return *this;}

   ~MCObsInfo(){}

   void setToRealPart();

   void setToImaginaryPart();

   void resetObsIndex(uint ind);  // for secondary observables


    // output functions
    
   bool isVacuum() const;

   bool isVEV() const;

   bool isCorrelatorAtTime() const;

   bool isHermitianCorrelatorAtTime() const;

   bool isReweightingFactor() const;

   bool isRealPart() const;

   bool isImaginaryPart() const;

   bool isSimple() const;

   bool isNonSimple() const;

   bool isPrimary() const;

   bool isSecondary() const;

   bool isBasicLapH() const;

   bool isGenIrrep() const;

   bool isVEVsubtractedCorrelatorAtTime() const;

   bool isReweightedCorrelatorAtTime() const;

   bool isReweightedVEV() const;


     // only allowed if is a vev or a non-vev subtracted correlator at time
   void setSimple();

   void setNotSimple();



     // routines below throw exception if inappropriate

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


   std::string output(bool longform=false, int indent=0) const;  // XML output 

   std::string str() const;  // XML output 

   void output(XMLHandler& xmlout, bool longform=false) const;  // XML output


   bool operator==(const MCObsInfo& rhs) const;

   bool operator!=(const MCObsInfo& rhs) const;

   bool operator<(const MCObsInfo& rhs) const;


   //  Routines below are used when MCObsInfo is a record key in
   //  an IOMap.  The IOMap class requires that every record key
   //  must occupy the same number of bytes.  

   static int numints() { return max_ints; }

   void copyTo(unsigned int *buf) const;

   explicit MCObsInfo(const unsigned int *buf);

   size_t numbytes() const {return max_ints*sizeof(unsigned int);}



 private:

   static const unsigned int max_ints = 24;

   void encode(const std::vector<uint>& precode, unsigned int optype, 
               bool simple, ComplexArg arg);

   void encode(const std::string& name, uint index,
               bool simple, ComplexArg arg);

   void set_real_part();

   void set_imag_part();

   void set_index(uint ind);

   void set_arg(ComplexArg arg);

   bool read_arg_type(XMLHandler& xmlin, ComplexArg& arg);

   void assert_corrtype(const std::string& msg="") const;

   void assert_vevtype(const std::string& msg="") const;

   std::string get_obs_name() const;

   uint get_obs_index() const;

   friend class CorrelatorAtTimeInfo;

};


// ***************************************************************

inline size_t numbytes(IOHandler& ioh, const MCObsInfo& rkey)
{ 
 return rkey.numbytes(); 
}

// ***************************************************************
#endif  

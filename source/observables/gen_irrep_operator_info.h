#ifndef GEN_IRREP_OPERATOR_INFO_H
#define GEN_IRREP_OPERATOR_INFO_H

#include <map>
#include <vector>
#include "momenta.h"
#include "args_handler.h"


// *******************************************************************
// *                                                                 *
// *   Objects of class "GenIrrepOperatorInfo" store identifying     *
// *   information about one particular QCD stationary-state         *
// *   operator of general but irreducible form.  For example,       *
// *   a "rotated" operator or other linear superpositions of the    *
// *   basic operators are general irreducible operators.  The XML   *
// *   format for specifying such an operator must be of the form:   *
// *                                                                 *
// *       <GIOperator>                                              *
// *          <Momentum>  0 1 1  </Momentum>                         *
// *          <ReferenceMomentum/>         (optional)                *
// *          <LGIrrep> T1gm </LGIrrep>      (optional)              *
// *          <LGIrrepRow> 3 </LGIrrepRow>   (optional)              *
// *          <Flavor> 1h -1 </Flavor>                               *
// *          <IDName>a_string_no_whitespace</IDName> (24 char max)  *
// *          <IDIndex> 2 </IDIndex> (0 if absent)                   *
// *       </GIOperator>                                             *
// *                                                                 *
// *   or of the form (see below)                                    *
// *                                                                 *
// *       <GIOperatorString> ... </GIOperatorString>                *
// *                                                                 *
// *   The "Flavor" tag must be 1 to 2 numbers separated by a space. *
// *   Two possibilities                                             *
// *     1) SU(3) - single integer specifying the SU(3)-flavor irrep *
// *     2) SU(2) - two numbers: Isospin and Strangeness             *
// *                                                                 *
// *   Construction can also be done by a short string.              *
// *   Example:                                                      *
// *     "Pref=(0,1,1) T1gm_3 flavor=1h,-1 a_string_no_whitespace 2" *
// *                                                                 *
// *   To specify that the operator has definite momentum, use P     *
// *   instead of Pref.                                              *
// *                                                                 *
// *   Note that the Pref provided will be stored in a canonical     *
// *   ordering such that |p_x| <= |p_y| <= |p_z|, so as to map all  *
// *   equivalent frames to one unique momentum.                     *
// *                                                                 *
// *   If the IDIndex token is absent, a value 0 is assumed.         *
// *                                                                 *
// *                                                                 *
// *******************************************************************


class GenIrrepOperatorInfo
{

   std::vector<unsigned int> icode;

#ifndef NO_CXX11
    GenIrrepOperatorInfo() = delete;
#else
    GenIrrepOperatorInfo();
#endif

 public:

   GenIrrepOperatorInfo(XMLHandler& xml_in);

   GenIrrepOperatorInfo(const std::string& opstring);

   GenIrrepOperatorInfo(const GenIrrepOperatorInfo& B) : icode(B.icode) {}

   GenIrrepOperatorInfo& operator=(const GenIrrepOperatorInfo& B)
    {icode=B.icode; return *this;}

   ~GenIrrepOperatorInfo(){}

   GenIrrepOperatorInfo& resetIDIndex(uint level);


    // output functions
    
   Momentum getMomentum() const;

   int getXMomentum() const;

   int getYMomentum() const;

   int getZMomentum() const;

   bool isReferenceMomentum() const;

   std::string getLGIrrep() const;

   unsigned int getLGIrrepRow() const;

   bool isSU3flavor() const;

   std::vector<std::string> getFlavor() const;

   std::string getIDName() const;

   unsigned int getIDIndex() const;



   std::string output(bool longform=false, int indent=0) const;  // XML output 

   std::string str() const;  // XML output 

   void output(XMLHandler& xmlout, bool longform=false) const;  // XML output

   std::string short_output() const;

   bool operator==(const GenIrrepOperatorInfo& rhs) const;

   bool operator!=(const GenIrrepOperatorInfo& rhs) const;

   bool operator<(const GenIrrepOperatorInfo& rhs) const;


 private:

   class Encoder
    {
     std::map<std::string, unsigned int> m_code;
     std::map<unsigned int, std::string> m_string;
     std::string m_codetype;
     Encoder(int ctype);
     Encoder(const Encoder& in);
     Encoder& operator=(const Encoder& in);
     ~Encoder(){}
     unsigned int encode(const std::string& description) const;
     std::string decode(unsigned int code) const;
     friend class GenIrrepOperatorInfo;
     void set_irreps();
     };

   static Encoder m_irreps;

   void assign(ArgsHandler& xt);
   void assign_from_string(const std::string& opstring);
   void momentum_from_string(const std::string& momstr, std::vector<int>& p);
   void encode(std::vector<int> mom, bool ref_mom, const std::string& irrep,
               uint irrepRow, const std::vector<std::string>& flavor,
               const std::string& name, uint index);

   static const unsigned int momt_bits = 27;
   static const unsigned int momj_bits = 8;
   static const unsigned int girr_bits = 3;
   static const unsigned int irrp_bits = 6;
   static const unsigned int irrw_bits = 3;
   static const unsigned int flav_bits = 12;
   static const unsigned int isop_bits = 6;
   static const unsigned int strange_bits = 5;
   static const unsigned int irrprwfl_bits = irrp_bits+irrw_bits+flav_bits;

   static const unsigned int momt_mask = 0x7FFFFFFu;
   static const unsigned int momj_mask = 0xFFu;
   static const unsigned int girr_mask = 0x7u;
   static const unsigned int irrp_mask = 0x3Fu;
   static const unsigned int irrw_mask = 0x7u;
   static const unsigned int flav_mask = 0xFFFu;
   static const unsigned int isop_mask = 0x3Fu;
   static const unsigned int strange_mask = 0x1Fu;
   static const unsigned int irrprwfl_mask = 0x7FFu;

   GenIrrepOperatorInfo(const std::vector<unsigned int>& incode)
     : icode(incode) {}

   static void resetIDIndex(uint level, uint& code);

   friend class MCObsInfo;
   friend class OperatorInfo;
   friend class BasicLapHOperatorInfo;
   friend class CorrelatorInfo;
   friend class CorrelatorAtTimeInfo;


};


// **************************************************
#endif  

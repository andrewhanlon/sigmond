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
// *          <Isospin> triplet </Isospin>                           *
// *          <Strangeness> -2 </Strangeness>  (default 0)           *
// *          <Momentum>  0 0 0  </Momentum>                         *
// *            or <MomentumSquared> 1 </MomentumSquared>            *
// *          <LGIrrep> T1gm </LGIrrep>                              *
// *          <LGIrrepRow> 3 </LGIrrepRow>   (optional)              *
// *          <IDName>a_string_no_whitespace</IDName> (24 char max)  *
// *          <IDIndex> 0 </IDIndex> (0 if absent)                   *
// *       </GIOperator>                                             *
// *                                                                 *
// *   or of the form (see below)                                    *
// *                                                                 *
// *       <GIOperatorString> ... </GIOperatorString>                *
// *                                                                 *
// *   The "Isospin" tag must take a value such as "singlet",        *
// *   "doublet", "triplet", etc.  The ID name cannot contain any    *
// *   white space.                                                  *
// *                                                                 *
// *   Construction can also be done by a short string.              *
// *   Example:                                                      *
// *     "isotriplet S=-1 P=(0,0,0) A1um_1 IDname 2"                 *
// *                                                                 *
// *   You can also specify the momentum squared and/or leave out    *
// *   the irrep row (e.g. if it was averaged over).                 *
// *   Example:                                                      *
// *     "isotriplet S=-1 PSQ=1 A2m IDname 2"                        *
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

   unsigned int getMomentumSquared() const;

   bool hasDefiniteMomentum() const;

   int getXMomentum() const;

   int getYMomentum() const;

   int getZMomentum() const;

   std::string getLGIrrep() const;

   unsigned int getLGIrrepRow() const;

   std::string getIsospin() const;

   int getStrangeness() const;

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
     void set_isospin();
     };

   static Encoder m_irreps, m_isospin;

   void assign(ArgsHandler& xt);
   void assign_from_string(const std::string& opstring);
   void momentum_from_string(const std::string& momstr, std::vector<int>& p);
   void encode(const std::string& isostr, int strangeness, const std::string& irrep, 
               uint irrepRow, const std::vector<int>& mom, const std::string& name,
               uint index);
   void encode(const std::string& isostr, int strangeness, const std::string& irrep, 
               uint irrepRow, uint mom_sqr, const std::string& name, uint index);

   static const unsigned int momt_bits = 24;
   static const unsigned int momj_bits = 7;
   static const unsigned int girr_bits = 3;
   static const unsigned int irrp_bits = 7;
   static const unsigned int isop_bits = 6;
   static const unsigned int strange_bits = 3;
   static const unsigned int irrw_bits = 4;
   static const unsigned int isirrw_bits = irrp_bits+isop_bits+irrw_bits;

   static const unsigned int momt_mask = 0xFFFFFFu;
   static const unsigned int momj_mask = 0x7Fu;
   static const unsigned int girr_mask = 0x7u;
   static const unsigned int irrp_mask = 0x7Fu;
   static const unsigned int isop_mask = 0x3Fu;
   static const unsigned int strange_mask = 0x7u;
   static const unsigned int irrw_mask = 0xFu;
   static const unsigned int isirrw_mask = 0x1FFFFu;

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

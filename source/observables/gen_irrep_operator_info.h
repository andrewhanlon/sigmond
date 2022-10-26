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
// *            or <Flavor> 1 -2 </Flavor>                           *
// *          <Momentum>  0 0 0  </Momentum>                         *
// *            or <MomentumSquared> 1 </MomentumSquared>            *
// *            or <ReferenceMomentum> 1 2 2 </ReferenceMomentum>    *
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
// *   Either the "Flavor" tag or the "Isospin/Strangeness" tag must *
// *   be present. The "Isospin" tag must take a value such as       *
// *   "singlet", "doublet", "triplet", etc. The "Flavor" tag must   *
// *   be 1 to 2 strings separated by a space.                       *
// *   Two possibilities                                             *
// *     1) SU(3) - single integer specifying the SU(3)-flavor irrep *
// *                if followed by '*', then it is the complex conj  *
// *                representation                                   *
// *     2) SU(2) - two numbers: Isospin and Strangeness             *
// *   The ID name cannot contain any white space.                   *
// *                                                                 *
// *   Construction can also be done by a short string.              *
// *   Examples:                                                     *
// *     "isotriplet S=-1 P=(0,0,0) A1um_1 IDname 2"                 *
// *          or                                                     *
// *        "Flavor=1,-1 P=(0,0,0) A1um_1 IDname 2"                  *
// *     "Flavor=27 P=(0,01) A2m particle1 3"                        *
// *                                                                 *
// *   You can also specify the momentum squared and/or leave out    *
// *   the irrep row (e.g. if it was averaged over).                 *
// *   Example:                                                      *
// *     "isotriplet S=-1 PSQ=1 A2m IDname 2"                        *
// *                                                                 *
// *   Or, you can even specify a reference momentum, which is       *
// *   userful when PSQ becomes ambiguous but you still want to      *
// *   indicate that equivalent momenta have been averaged over.     *
// *   Example:                                                      *
// *     "isotriplet S=-1 Pref=(1,2,2) A2m IDname 2"                 *
// *     "isotriplet S=-1 Pref=(0,0,3) A2m IDname 2"                 *
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
    
   bool hasDefiniteMomentum() const;

   bool hasMomentumSquared() const;

   bool hasReferenceMomentum() const;

   Momentum getMomentum() const;

   unsigned int getMomentumSquared() const;

   int getXMomentum() const;

   int getYMomentum() const;

   int getZMomentum() const;

   std::string getLGIrrep() const;

   unsigned int getLGIrrepRow() const;

   bool hasSU3Flavor() const;

   bool hasSU2Flavor() const;

   std::string getIsospin() const;

   int getStrangeness() const;

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
     void set_isospin();
     };

   static Encoder m_irreps, m_isospin;

   void assign(ArgsHandler& xt);
   void assign_from_string(const std::string& opstring);
   void momentum_from_string(const std::string& momstr, std::vector<int>& p);
   void encode(const std::vector<std::string>& flavor, const std::string& irrep, 
               uint irrepRow, std::vector<int>& mom, bool reference,
               const std::string& name, uint index);
   void encode(const std::vector<std::string>& flavor, const std::string& irrep,
               uint irrepRow, uint mom_sqr, const std::string& name, uint index);

   static const unsigned int flav_bits = 12;
   static const unsigned int su3flav_bits = 11;
   static const unsigned int strange_bits = 5;
   static const unsigned int momt_bits = 24;
   static const unsigned int momj_bits = 7;
   static const unsigned int girr_bits = 3;
   static const unsigned int irrp_bits = 6;
   static const unsigned int irrw_bits = 3;
   static const unsigned int irflav_bits = 19;

   static const unsigned int flav_mask = 0xFFFu;
   static const unsigned int su3flav_mask = 0x7FFu;
   static const unsigned int isop_mask = 0x3Fu;
   static const unsigned int strange_mask = 0x1Fu;
   static const unsigned int momt_mask = 0xFFFFFFu;
   static const unsigned int momj_mask = 0x7Fu;
   static const unsigned int girr_mask = 0x7u;
   static const unsigned int irrp_mask = 0x3Fu;
   static const unsigned int irrw_mask = 0x7u;
   static const unsigned int id_mask = 0x1FFFu;
   static const unsigned int irflav_mask = 0x7FFFFu;

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

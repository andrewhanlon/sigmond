#ifndef LAPH_OPERATOR_INFO_H
#define LAPH_OPERATOR_INFO_H

#include <map>
#include <vector>
#include "momenta.h"
#include "args_handler.h"


// *******************************************************************
// *                                                                 *
// *   Objects of class "OperatorInfo" store identifying information *
// *   about one particular QCD stationary-state operator.  The XML  *
// *   format for specifying such an operator must be of the form:   *
// *                                                                 *
// *       <Operator>                                                *
// *          <NumberOfHadrons> 2 </NumberOfHadrons>                 *
// *               ....                                              *
// *       </Operator>                                               *
// *                                                                 *
// *   or of the form (see below)                                    *
// *                                                                 *
// *       <OperatorString> ... </OperatorString>                    *
// *                                                                 *
// *   or of the form (see below)                                    *
// *                                                                 *
// *       <RotatedOperator> ... </RotatedOperator>                  *
// *                                                                 *
// *   or of the form (see below)                                    *
// *                                                                 *
// *       <RotatedOpString>name 4</RotatedOpString>                 *
// *                                                                 *
// *   The "Operator" tag must contain a "NumberOfHadrons" tag.      *
// *   The rest of the XML depends on the number of hadrons.         *
// *                                                                 *
// *   A "hadron" here means a baryon, a meson, or a glueball.       *
// *   It could also mean a tetraquark system, a pentaquark system,  *
// *   etc.  It is a gauge-invariance localized quantity.            *
// *                                                                 *
// *   Number of hadrons = 0                                         *
// *                                                                 *
// *     No further tags are required                                *
// *        (this is the "default" constructor).                     *
// *                                                                 *
// *   Number of hadrons = 1                                         *
// *                                                                 *
// *      <Hadron>                                                   *
// *         ....described below....                                 *
// *      </Hadron>                                                  *
// *      <LGIrrepRow> 3 </LGIrrepRow>                               *
// *                                                                 *
// *   Number of hadrons = 2,3,4,...                                 *
// *                                                                 *
// *      <Total>                                                    *
// *          <Isospin> triplet </Isospin>                           *
// *          <IsoCGId> 0 </IsoCGId> (assume zero if absent)         *
// *          <Momentum>  0 0 0  </Momentum>                         *
// *          <LGIrrep> T1gm </LGIrrep>                              *
// *          <LGCGId> 0 </LGCGId> (assume zero if absent)           *
// *          <LGIrrepRow> 3 </LGIrrepRow>                           *
// *      </Total>                                                   *
// *      <Hadron1>                                                  *
// *         ....                                                    *
// *      </Hadron1>                                                 *
// *      <Hadron2>                                                  *
// *         ....                                                    *
// *      </Hadron2>                                                 *
// *      ....                                                       *
// *                                                                 *
// *   In the "Total" tag, the "Isospin" tag must take a value such  *
// *   as "singlet", "doublet", "triplet", etc.  When there are      *
// *   three or more hadrons, an isospin Clebsch-Gordan occurrence   *
// *   identifying number "IsoCGId" must be specified: value 0 to    *
// *   one less than the number of times that LGIrrep occurs in      *
// *   the direct product of the single-hadron isospins.  A value    *
// *   of zero for this tag is assumed if the tag is absent.         *
// *   If the irrep occurs more than once in the Clebsch-Gordan      *
// *   series of the direct product of single-hadron irreps, then    *
// *   a little group Clebsch-Gordan identifying number "LGCGId"     *
// *   must be specified (value 0 to one less than the number of     *
// *   occurrences).  A value of 0 is assumed if absent.             *
// *                                                                 *
// *   Each of the constituent single-hadrons must be described in   *
// *   a separate tag "Hadron1", "Hadron2", etc.  If a single hadron *
// *   is a meson or a baryon, then the "Hadron" tag has the form    *
// *                                                                 *
// *         <Hadron1>                                               *
// *            <Flavor> eta </Flavor>                               *
// *            <Momentum>  0 0 0  </Momentum>                       *
// *            <LGIrrep> A1p </LGIrrep>                             *
// *            <SpatialType> DDL </SpatialType>                     *
// *            <SpatialIdNum> 4 </SpatialIdNum>                     *
// *            <DispLength> 3 </DispLength>                         *
// *         </Hadron1>                                              *
// *                                                                 *
// *   Allowed meson "Flavor" tag values are "pion", "eta", "phi",   *
// *   "kaon", and "kbar". Note that "pion" does NOT mean an actual  *
// *   pion, rather, it just means an isovector consisting of u,d    *
// *   quarks.  "pion" is just a shorter name than "isovector_du".   *
// *   Similarly, an "eta" means an isoscalar meson consisting of    *
// *   u,d quarks.  A "phi" is an isoscalar meson that is s-sbar.    *
// *   "kaon" refers to a strangeness S=1 meson, and "kbar" refers   *
// *   to a strangeness S=-1 meson. Allowed baryon "Flavor" tag      *
// *   values are "nucleon", "delta", "lambda", "sigma", "xi", and   *
// *   "omega".  "Momentum" is the momentum of the particle          *
// *   in a chosen "reference" term of the total operator.  The      *
// *   total momentum can be obtained by adding all of the single    *
// *   hadron momenta, and the total momentum is fixed.  The total   *
// *   operator is a superposition of terms that are rotations,      *
// *   parity-transformations, and G-parity transformations of the   *
// *   reference term. Each "Momentum" tag must contain three        *
// *   integers (in units of 2*Pi/L) which describe the momentum     *
// *   of the hadron in the reference term. "LGIrrep" specifies the  *
// *   irreducible representation of the little group corresponding  *
// *   to the hadron momentum.                                       *
// *                                                                 *
// *   If the single constituent is a glueball, then its "Hadron"    *
// *   tag must be of the form                                       *
// *                                                                 *
// *         <Hadron1>                                               *
// *            <Flavor> glueball </Flavor>                          *
// *            <Momentum>  0 0 0  </Momentum>                       *
// *            <LGIrrep> A1gp </LGIrrep>                            *
// *            <SpatialType> TrEig </SpatialType>                   *
// *         </Hadron1>                                              *
// *                                                                 *
// *   Construction can also be done by a short string:              *
// *   Examples:                                                     *
// *     "glueball P=(0,0,0) A1gp_1 TrEig"                           *
// *     "pion P=(0,0,0) A1um_1 SD_5"                                *
// *     "isotriplet_pion_pion A1um_1 CG_1 [P=(0,0,1) A1p LSD_1] [P=(0,0,-1) A2m TSD_2]"
// *                                                                 *
// *   If the "CG_1" token is absent, a value 0 is assumed.          *
// *                                                                 *
// *   These strings can also be used inside an <OperatorString>     *
// *   tag in XML format:                                            *
// *                                                                 *
// *       <OperatorString> ... </OperatorString>                    *
// *                                                                 *
// *   where the tag should contain a string in the format above.    *
// *                                                                 *
// *   Often, a so-called "rotated" operator is needed.  These are   *
// *   formed using linear superpositions of the usual operators     *
// *   described above.  The XML for specifying such an operator     *
// *   is                                                            *
// *       <RotatedOperator>                                         *
// *         <ObsName>a_string_no_blanks</ObsName>                   *
// *         <Level>4</Level>    (default 0)                         *
// *         <Description>....</Description>    (optional)           *
// *       </RotatedOperator>                                        *
// *   or                                                            *
// *       <RotatedOpString>name 4</RotatedOpString>                 *
// *                                                                 *
// *   A call to "getNumberOfHadrons()" returns 7, but this is       *
// *   notational only since such operators may not have a definite  *
// *   number of hadrons.  The optional <Description> tag must be    *
// *   valid XML is used only for informational and outputting       *
// *   purposes.                                                     *
// *                                                                 *
// *                                                                 *
// *******************************************************************


class OperatorInfo
{

   std::vector<unsigned int> icode;

 public:

   enum OpKind { LaphSpec, Rotated };

   OperatorInfo();

   OperatorInfo(XMLHandler& xml_in);

   OperatorInfo(const std::string& opstring, OpKind opkind = LaphSpec);

   OperatorInfo(const OperatorInfo& B) : icode(B.icode) {}

   OperatorInfo& operator=(const OperatorInfo& B)
    {icode=B.icode; return *this;}

   ~OperatorInfo(){}


    // output functions
    
   unsigned int getNumberOfHadrons() const;   // returns 7 for rotated operator

   unsigned int get_NumberOfHadrons() const;   // exception for rotated operator

   bool isVacuum() const;

   bool isGlueball() const;

   bool isMeson() const;

   bool isBaryon() const;

   bool isMesonMeson() const;

   bool isMesonBaryon() const; 

   bool isRotated() const;
   

   Momentum getMomentum() const;

   int getXMomentum() const;

   int getYMomentum() const;

   int getZMomentum() const;

   std::string getLGIrrep() const;

   unsigned int getLGClebschGordonIdNum() const;

   unsigned int getLGIrrepRow() const;

   std::string getIsospin() const;

   unsigned int getIsospinClebschGordonIdNum() const;

   std::string getFlavor() const;

   std::string getFlavorCode() const;

   std::string getRotatedName() const;

   unsigned int getRotatedLevel() const;
   
   OperatorInfo& resetRotatedLevel(uint level);


             // note: hadron_index = 1, 2, ...
             
   std::string getFlavor(unsigned int hadron_index) const;

   int getStrangeness(unsigned int hadron_index) const;
 
   bool isGlueball(unsigned int hadron_index) const;

   bool isMeson(unsigned int hadron_index) const;

   bool isBaryon(unsigned int hadron_index) const;

   bool isFermion(unsigned int hadron_index) const;
   
   bool isBoson(unsigned int hadron_index) const;

   std::string getLGIrrep(unsigned int hadron_index) const;

   std::string getSpatialType(unsigned int hadron_index) const;

   unsigned int getSpatialIdNumber(unsigned int hadron_index) const;

   unsigned int getDisplacementLength(unsigned int hadron_index) const;

   Momentum getMomentum(unsigned int hadron_index) const;

   int getXMomentum(unsigned int hadron_index) const;

   int getYMomentum(unsigned int hadron_index) const;

   int getZMomentum(unsigned int hadron_index) const;


   std::string output(bool longform=true, int indent=0) const;  // XML output 

   std::string str() const;  // XML output 

   void output(XMLHandler& xmlout, bool longform=true) const;  // XML output

   std::string short_output() const;

   bool operator==(const OperatorInfo& rhs) const;

   bool operator!=(const OperatorInfo& rhs) const;

   bool operator<(const OperatorInfo& rhs) const;


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
     friend class OperatorInfo;
     void set_flavor();
     void set_flavorcode();
     void set_spatial();
     void set_irreps();
     void set_isospin();
     void set_isospin_char();
     };

   static Encoder m_flavor, m_irreps, m_isospin, m_spatial, m_flavorcode, m_isochar;
   
   void check_hadron_index(unsigned int hadron_index, const char *msg) const;
   void xmlread_momentum(ArgsHandler& xmlh, const std::string& tag,
                         Momentum& mom, unsigned int& momcode);
   void xmlwrite_momentum(XMLHandler& xmlout, const std::string& tag,
                          unsigned int momcode) const;
   void xmlread_hadron(ArgsHandler& xml_in, const std::string& toptag,
                       unsigned int& momcode, unsigned int& hadroncode);
   void xmlwrite_hadron(XMLHandler& xml_in, const std::string& toptag,
                        const unsigned int& momcode, const unsigned int& hadroncode) const;
   void xmlread_total(ArgsHandler& xml_in, unsigned int& code);
   void xmlwrite_total(XMLHandler& xml_in, const unsigned int& code) const;

   void encode_momentum(const std::string& momstr, unsigned int& momcode);
   void encode_hadron(const std::string& flav, const std::string& irrep, 
                      const std::string& sptype, const std::string& spid, 
                      unsigned int& hadroncode);
   void encode_total(const std::string& totalisospin, const std::string& isocgid,
                     const std::string& totalirrep, const std::string& irreprow, 
                     const std::string& lgcgid, unsigned int& code);

   void shortwrite_hadron(unsigned int momcode, unsigned int hadroncode,
                          const std::string& irreprow, std::string& flavor, 
                          std::string& hadstring) const;
   void shortwrite_total(unsigned int code, std::string& totaliso, 
                         std::string& totirrep) const;

   void assign(const std::string& opstring);
   Momentum get_momentum(unsigned int momcode) const;
   int get_xmomentum(unsigned int momcode) const;
   int get_ymomentum(unsigned int momcode) const;
   int get_zmomentum(unsigned int momcode) const;
   std::string get_flavor(unsigned int hadroncode) const;
   std::string get_flavorcode(unsigned int hadroncode) const;
   bool is_glueball(unsigned int hadroncode) const;
   bool is_meson(unsigned int hadroncode) const;
   bool is_baryon(unsigned int hadroncode) const;
   bool is_fermion(unsigned int hadroncode) const;
   bool is_boson(unsigned int hadroncode) const;
   std::string get_lgirrep(unsigned int hadroncode) const;
   std::string get_spatial_type(unsigned int hadroncode) const;
   unsigned int get_spatial_id(unsigned int hadroncode) const;
   unsigned int get_disp_length(unsigned int hadroncode) const;

   int count(const std::string& astr, char delimiter) const;
   std::vector<std::string> split(const std::string& astr, char delimiter) const;
   std::string extract(const std::string& astr, char left, char right) const;

   static const unsigned int momt_bits = 24;
   static const unsigned int momj_bits = 7;
   static const unsigned int nhad_bits = 3;
   static const unsigned int dlen_bits = 3;
   static const unsigned int spid_bits = 13;
   static const unsigned int sptp_bits = 4;
   static const unsigned int irrp_bits = 7;
   static const unsigned int flav_bits = 5;
   static const unsigned int isop_bits = 6;
   static const unsigned int iscg_bits = 6;
   static const unsigned int lgcg_bits = 6;
   static const unsigned int irrw_bits = 4;

   static const unsigned int momt_mask = 0xFFFFFFu;
   static const unsigned int momj_mask = 0x7Fu;
   static const unsigned int nhad_mask = 0x7u;
   static const unsigned int dlen_mask = 0x7u;
   static const unsigned int spid_mask = 0x1FFFu;
   static const unsigned int sptp_mask = 0xFu;
   static const unsigned int irrp_mask = 0x7Fu;
   static const unsigned int flav_mask = 0x1Fu;
   static const unsigned int isop_mask = 0x3Fu;
   static const unsigned int iscg_mask = 0x3Fu;
   static const unsigned int lgcg_mask = 0x3Fu;
   static const unsigned int irrw_mask = 0xFu;

   static std::map<std::string,uint> m_rotated_encodings;
   static std::vector<std::pair<std::string,std::string> > m_rotated_decodings;

   void xmlread_rotated(ArgsHandler& xml_in);
   void assign_rotated(const std::string& opstring);
   void encode_rotated(const std::string& name, uint level, 
                       const std::string& description);
   unsigned int get_rotated_obs_code() const;
   std::string short_rotated_output() const;

   OperatorInfo(std::vector<unsigned int>::const_iterator inbegin,
                std::vector<unsigned int>::const_iterator inend) 
          : icode(inbegin,inend) {}

   static unsigned int codesize(unsigned int icode0)
    {unsigned int i0=(icode0 & nhad_mask);
     return (i0!=7)?((i0>1)?(2*i0+1):(i0+1)):1;}

   friend class MCObsInfo;
   friend class CorrelatorInfo;
   friend class CorrelatorAtTimeInfo;


};


// **************************************************
#endif  

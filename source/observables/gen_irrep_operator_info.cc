#include "gen_irrep_operator_info.h"
#include "multi_compare.h"
#include "args_handler.h"
#include "log_helper.h"
#include "encoder.h"
#include <stdexcept>

using namespace std;

  //  Encode operator information into 3 to 8 unsigned 32-bit integers
  //  Operator encoding:
  //       icode[0] (left to right) =
  //             5 bits empty, momentum in 24 bits (8 bits each momx, 
  //             momy, momz; first bit is sign, 7 bits remaining), 
  //             rightmost 3 bits set to 111 = 7 to indicate GIOperatorInfo
  //       icode[1] =
  //             15 bits ID index, 6 bits <Isospin>, 7 bits <LGIrrep>, 
  //             4 bits <LGIrrepRow>,  
  //       icode[2,...] = integer representation of ID string (max 6 ints)
  //



GenIrrepOperatorInfo::GenIrrepOperatorInfo(XMLHandler& xml_in)
{
 try{
 set<string> tags;
 tags.insert("GIOperator");
 tags.insert("GIOperatorString");
 ArgsHandler xin(xml_in,tags);
 string rtag=xin.getInputRootTag();
 if (rtag=="GIOperator"){
    assign(xin);
    return;}
 else if (rtag=="GIOperatorString"){
    assign_from_string(xin.getString("GIOperatorString"));
    return;}
 else
    throw(std::invalid_argument("Invalid XML input"));
 }
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("GenIrrepOperatorInfo construction failed: \n")
      +string(errmsg.what())+string("\nInput XML:")+xml_in.output()));}
}


    //   Examples:
    //     "isotriplet P=(0,0,0) A1um_1 IDname 2"

GenIrrepOperatorInfo::GenIrrepOperatorInfo(const std::string& opstring)
{
 assign_from_string(opstring);
}


 // ********************************************************************


void GenIrrepOperatorInfo::assign(ArgsHandler& xt)
{
 string isostr(xt.getString("Isospin"));
 string irrep(xt.getString("LGIrrep"));
 int irrepRow;
 xt.getInt("LGIrrepRow",irrepRow);
 vector<int> mom(xt.getIntVector("Momentum"));
 string name(xt.getName("IDName"));
 uint index=0;
 xt.getOptionalUInt("IDIndex",index);
 encode(isostr,irrep,irrepRow,mom,name,index);
}


void GenIrrepOperatorInfo::assign_from_string(const string& opstring)
{
 try{
 string opstr(tidyString(opstring));
 vector<string> tokens=split(opstr,' ');
 if ((tokens.size()!=4)&&(tokens.size()!=5)) throw(std::runtime_error(""));
 string isostr(tokens[0]); 
 size_t pos=isostr.find("iso");
 if (pos!=string::npos) isostr.erase(pos,3);
 vector<int> mom;
 momentum_from_string(tokens[1],mom);
 string name(tokens[3]);
 string irrep(tokens[2]);  
 uint index=0;
 if (tokens.size()==5){
    extract_from_string(tokens[4],index);}
 tokens=split(irrep,'_');
 if (tokens.size()!=2) throw(std::runtime_error(""));
 irrep=tokens[0];
 uint irrepRow;
 extract_from_string(tokens[1],irrepRow);
 encode(isostr,irrep,irrepRow,mom,name,index);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("Invalid GenIrrepOperatorInfo string: ")+opstring));}
}



void GenIrrepOperatorInfo::encode(const string& isostr, const string& irrep, 
                  uint irrepRow, const vector<int>& mom, const string& name,
                  uint index)
{
 unsigned int isocode=m_isospin.encode(isostr);
 unsigned int irrep_code=m_irreps.encode(irrep);
 if ((irrepRow<=0)||(irrepRow>int(irrw_mask))){
    throw(std::invalid_argument("Irrep row not currently supported"));}
 if (mom.size()!=3){
    throw(std::invalid_argument("Bad momentum"));}
 if (((unsigned int)abs(mom[0])>momj_mask)||((unsigned int)abs(mom[1])>momj_mask)
    ||((unsigned int)abs(mom[2])>momj_mask)){
    throw(std::invalid_argument("momentum component magnitude not currently supported"));}
 uint momcode=(mom[0]<0)?1:0; 
 momcode<<=momj_bits; momcode|=abs(mom[0]);
 momcode<<=1; momcode|=(mom[1]<0)?1:0;
 momcode<<=momj_bits; momcode|=abs(mom[1]);
 momcode<<=1; momcode|=(mom[2]<0)?1:0;
 momcode<<=momj_bits; momcode|=abs(mom[2]);
 const uint maxlength=24;
 if (name.length()>maxlength){
    throw(std::invalid_argument("GIOperator name too long"));}
 if (index>=32768){
    throw(std::invalid_argument("GIOperator index too large"));}
 vector<uint> namecode;
 encode_string_to_uints(name,maxlength,namecode);
 icode.resize(namecode.size()+2);
 icode[0]=momcode;
 icode[0]<<=girr_bits;
 icode[0]|=girr_mask;
 uint tcode=index;   
 tcode<<=isop_bits;  tcode|=isocode;
 tcode<<=irrp_bits;  tcode|=irrep_code;
 tcode<<=irrw_bits;  tcode|=irrepRow;
 icode[1]=tcode;
 for (uint k=0;k<namecode.size();k++)
    icode[k+2]=namecode[k];
}



string GenIrrepOperatorInfo::output(bool longform, int indent) const
{
 XMLHandler xmlout;
 output(xmlout,longform);
 return xmlout.output(indent);
}

string GenIrrepOperatorInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}

void GenIrrepOperatorInfo::output(XMLHandler& xmlout, bool longform) const
{
 string irrep(getLGIrrep());
 uint irrepRow=getLGIrrepRow();
 string isospin(getIsospin());
 string idname(getIDName());
 uint index=getIDIndex();
 Momentum P(getMomentum());
 if (!longform){
     string opstr="iso"+isospin+" P=("+make_string(P.x)+","+make_string(P.y)+","+make_string(P.z)+") "
                 +irrep+"_"+make_string(irrepRow)+" "+idname+" "+make_string(index);
    xmlout.set_root("GIOperatorString",opstr);
    return;}
 xmlout.set_root("GIOperator");
 xmlout.put_child("Isospin",isospin);
 xmlout.put_child("Momentum",make_string(P.x)+" "+make_string(P.y)+" "+make_string(P.z));
 xmlout.put_child("LGIrrep",irrep);
 xmlout.put_child("LGIrrepRow",make_string(irrepRow));
 xmlout.put_child("IDName",idname);
 xmlout.put_child("IDIndex",make_string(index));
}


string GenIrrepOperatorInfo::short_output() const
{
 string irrep(getLGIrrep());
 uint irrepRow=getLGIrrepRow();
 string isospin(getIsospin());
 string idname(getIDName());
 uint index=getIDIndex();
 Momentum P(getMomentum());
 string opstr=isospin+" P=("+make_string(P.x)+","+make_string(P.y)+","+make_string(P.z)+") "
              +irrep+"_"+make_string(irrepRow)+" "+idname+" "+make_string(index);
 return opstr;
}




bool GenIrrepOperatorInfo::operator==(const GenIrrepOperatorInfo& rhs) const
{
 return multiEqual(icode,rhs.icode);   
}

bool GenIrrepOperatorInfo::operator!=(const GenIrrepOperatorInfo& rhs) const
{
 return multiNotEqual(icode,rhs.icode);   
}

bool GenIrrepOperatorInfo::operator<(const GenIrrepOperatorInfo& rhs) const
{
 return multiLessThan(icode,rhs.icode);   
}


 // *******************************************************************


std::string GenIrrepOperatorInfo::getLGIrrep() const
{
 return m_irreps.decode((icode[1]>>irrw_bits) & irrp_mask);
}

unsigned int GenIrrepOperatorInfo::getLGIrrepRow() const
{
 return (icode[1]&irrw_mask);
}

std::string GenIrrepOperatorInfo::getIsospin() const
{
 return m_isospin.decode((icode[1]>>(irrp_bits+irrw_bits)) & isop_mask);
}


std::string GenIrrepOperatorInfo::getIDName() const
{
 vector<uint>::const_iterator it=icode.begin()+2;
 vector<uint> namecode(it,icode.end());
 return decode_uints_to_string(namecode);
}


unsigned int GenIrrepOperatorInfo::getIDIndex() const
{
 return icode[1]>>isirrw_bits; 
}


GenIrrepOperatorInfo& GenIrrepOperatorInfo::resetIDIndex(uint level)
{
 resetIDIndex(level,icode[1]);
 return *this;
}


Momentum GenIrrepOperatorInfo::getMomentum() const
{   
 unsigned int tmp=(icode[0]>>girr_bits);
 int pz=tmp & momj_mask;
 tmp>>=momj_bits;
 if ((tmp&0x1u)==1) pz=-pz;
 tmp>>=1;
 int py=tmp & momj_mask;
 tmp>>=momj_bits;
 if ((tmp&0x1u)==1) py=-py;
 tmp>>=1;
 int px=tmp & momj_mask;
 tmp>>=momj_bits;
 if ((tmp&0x1u)==1) px=-px;
 return Momentum(px,py,pz);
}

int GenIrrepOperatorInfo::getXMomentum() const
{
 unsigned int tmp=(icode[0]>>(girr_bits+2*momj_bits+2));
 int res=tmp & momj_mask;
 if (((tmp>>momj_bits)&0x1u)==1) return -res;
 return res;
}
  
int GenIrrepOperatorInfo::getYMomentum() const
{
 unsigned int tmp=(icode[0]>>(girr_bits+momj_bits+1));
 int res=tmp & momj_mask;
 if (((tmp>>momj_bits)&0x1u)==1) return -res;
 return res;
}

int GenIrrepOperatorInfo::getZMomentum() const
{
 unsigned int tmp=(icode[0]>>girr_bits);
 int res=tmp & momj_mask;
 if (((tmp>>momj_bits)&0x1u)==1) return -res;
 return res;
}



// **************************************************

       // static data needed
       
GenIrrepOperatorInfo::Encoder GenIrrepOperatorInfo::m_irreps(0);
GenIrrepOperatorInfo::Encoder GenIrrepOperatorInfo::m_isospin(1);


GenIrrepOperatorInfo::Encoder::Encoder(int ctype)
{
 if (ctype==0) set_irreps();
 else if (ctype==1) set_isospin();
}

void GenIrrepOperatorInfo::Encoder::set_irreps()
{
 m_codetype="irrep name";
 m_code["A1gp"]= 0;   m_string[ 0]="A1gp";
 m_code["A1"  ]= 1;   m_string[ 1]="A1"  ;
 m_code["A1g" ]= 2;   m_string[ 2]="A1g" ;
 m_code["A1gm"]= 3;   m_string[ 3]="A1gm";
 m_code["A1m" ]= 4;   m_string[ 4]="A1m" ;
 m_code["A1p" ]= 5;   m_string[ 5]="A1p" ;
 m_code["A1u" ]= 6;   m_string[ 6]="A1u" ;
 m_code["A1um"]= 7;   m_string[ 7]="A1um";
 m_code["A1up"]= 8;   m_string[ 8]="A1up";
 m_code["A2"  ]= 9;   m_string[ 9]="A2"  ;
 m_code["A2g" ]=10;   m_string[10]="A2g" ;
 m_code["A2gm"]=11;   m_string[11]="A2gm";
 m_code["A2gp"]=12;   m_string[12]="A2gp";
 m_code["A2m" ]=13;   m_string[13]="A2m" ;
 m_code["A2p" ]=14;   m_string[14]="A2p" ;
 m_code["A2u" ]=15;   m_string[15]="A2u" ;
 m_code["A2um"]=16;   m_string[16]="A2um";
 m_code["A2up"]=17;   m_string[17]="A2up";
 m_code["B1"  ]=18;   m_string[18]="B1"  ;
 m_code["B1m" ]=19;   m_string[19]="B1m" ;
 m_code["B1p" ]=20;   m_string[20]="B1p" ;
 m_code["B2"  ]=21;   m_string[21]="B2"  ;
 m_code["B2m" ]=22;   m_string[22]="B2m" ;
 m_code["B2p" ]=23;   m_string[23]="B2p" ;
 m_code["E"   ]=24;   m_string[24]="E"   ;
 m_code["Eg"  ]=25;   m_string[25]="Eg"  ;
 m_code["Egm" ]=26;   m_string[26]="Egm" ;
 m_code["Egp" ]=27;   m_string[27]="Egp" ;
 m_code["Em"  ]=28;   m_string[28]="Em"  ;
 m_code["Ep"  ]=29;   m_string[29]="Ep"  ;
 m_code["Eu"  ]=30;   m_string[30]="Eu"  ;
 m_code["Eum" ]=31;   m_string[31]="Eum" ;
 m_code["Eup" ]=32;   m_string[32]="Eup" ;
 m_code["F1"  ]=33;   m_string[33]="F1"  ;
 m_code["F2"  ]=34;   m_string[34]="F2"  ;
 m_code["G"   ]=35;   m_string[35]="G"   ;
 m_code["G1"  ]=36;   m_string[36]="G1"  ;
 m_code["G1g" ]=37;   m_string[37]="G1g" ;
 m_code["G1u" ]=38;   m_string[38]="G1u" ;
 m_code["G2"  ]=39;   m_string[39]="G2"  ;
 m_code["G2g" ]=40;   m_string[40]="G2g" ;
 m_code["G2u" ]=41;   m_string[41]="G2u" ;
 m_code["Hg"  ]=42;   m_string[42]="Hg"  ;
 m_code["Hu"  ]=43;   m_string[43]="Hu"  ;
 m_code["T1g" ]=44;   m_string[44]="T1g" ;
 m_code["T1gm"]=45;   m_string[45]="T1gm";
 m_code["T1gp"]=46;   m_string[46]="T1gp";
 m_code["T1u" ]=47;   m_string[47]="T1u" ;
 m_code["T1um"]=48;   m_string[48]="T1um";
 m_code["T1up"]=49;   m_string[49]="T1up";
 m_code["T2g" ]=50;   m_string[50]="T2g" ;
 m_code["T2gm"]=51;   m_string[51]="T2gm";
 m_code["T2gp"]=52;   m_string[52]="T2gp";
 m_code["T2u" ]=53;   m_string[53]="T2u" ;
 m_code["T2um"]=54;   m_string[54]="T2um";
 m_code["T2up"]=55;   m_string[55]="T2up";
}

void GenIrrepOperatorInfo::Encoder::set_isospin()
{
 m_codetype="total isospin";
     // Total isospin
 m_code["singlet"]=0;   m_string[0]="singlet";
 m_code["doublet"]=1;   m_string[1]="doublet";
 m_code["triplet"]=2;   m_string[2]="triplet";
 m_code["quartet"]=3;   m_string[3]="quartet";
 m_code["quintet"]=4;   m_string[4]="quintet";
 m_code["sextet" ]=5;   m_string[5]="sextet" ;
}


unsigned int GenIrrepOperatorInfo::Encoder::encode(const std::string& description) const
{
 map<string,unsigned int>::const_iterator it=m_code.find(description);
 if (it==m_code.end()) 
    throw(std::invalid_argument(string("Invalid ")+m_codetype
                  +string(" string in GenIrrepOperatorInfo")));
 return it->second;
}
 
std::string GenIrrepOperatorInfo::Encoder::decode(unsigned int code) const
{
 map<unsigned int,string>::const_iterator it=m_string.find(code);
 if (it==m_string.end()) 
    throw(std::invalid_argument(string("Invalid ")+m_codetype
             +string(" code in GenIrrepOperatorInfo")));
 return it->second;
}


vector<string> GenIrrepOperatorInfo::split(const string& astr, char delimiter) const
{
 vector<string> tokens;
 size_t lastpos=astr.find_first_not_of(delimiter);
 size_t pos=(lastpos==string::npos)?string::npos:astr.find_first_of(delimiter,lastpos+1);
 while (lastpos!=string::npos){
    if (pos==string::npos) pos=astr.length();
    tokens.push_back(astr.substr(lastpos,pos-lastpos));
    lastpos=astr.find_first_not_of(delimiter,pos+1);
    pos=(lastpos==string::npos)?string::npos:astr.find_first_of(delimiter,lastpos+1);}
 return tokens;
}


void GenIrrepOperatorInfo::momentum_from_string(const std::string& momstr, std::vector<int>& pmom)
{
 if (momstr.length()<9) throw(std::invalid_argument("Invalid momentum string"));
 if ((momstr[0]!='P')||(momstr[1]!='=')||(momstr[2]!='(')||(momstr[momstr.length()-1]!=')'))
    throw(std::invalid_argument("Invalid momentum string"));
 string mmm(momstr.substr(3,momstr.length()-4));
 vector<string> p=split(mmm,',');
 if (p.size()!=3){
    throw(std::invalid_argument("Invalid momentum string"));}
 pmom.resize(3);
 try{
    extract_from_string(p[0],pmom[0]);
    extract_from_string(p[1],pmom[1]);
    extract_from_string(p[2],pmom[2]);}
 catch(const std::exception& msg){
    throw(std::invalid_argument("Invalid momentum string"));}
}


void GenIrrepOperatorInfo::resetIDIndex(uint level, uint& kcode)
{
 uint keep=kcode & isirrw_mask;  // keep right 17 bits  
 kcode=(level<<isirrw_bits);
 kcode|=keep;
}

// ******************************************************************************

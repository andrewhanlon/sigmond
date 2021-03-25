#include "gen_irrep_operator_info.h"
#include "multi_compare.h"
#include "args_handler.h"
#include "log_helper.h"
#include "encoder.h"
#include <algorithm>
#include <stdexcept>

using namespace std;

  //  Encode operator information into 3 to 8 unsigned 32-bit integers
  //  Operator encoding:
  //       icode[0] (left to right) =
  //             1 bit for SU(3)-flavor (1) or Isospin/Strangeness (0),
  //             1 bit for momentum reference (1) or definite momentum (0),
  //             momentum in 27 bits (9 bits each momx, momy, momz;
  //             first bit is sign, 8 bits remaining),
  //             rightmost 3 bits set to 111 = 7 to indicate GIOperatorInfo
  //       icode[1] =
  //             11 bits ID index, 6 bits <LGIrrep>, 3 bits <LGIrrepRow>,  
  //             12 bits for flavor:
  //               * if SU(3) flavor, then the leftmost bit is 1 if a
  //                 complex conjugate representation, otherwise it
  //                 is 0.
  //               * if isopspin/strangeness, then the left 6 bits
  //                 are the isospin and the right 6 bits are the
  //                 strangeness (with the leftmost strangeness
  //                 bit being the sign).
  //
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
    //     "P=(0,1,1) T1gm_3 flavor=1h,-1 IDname 2"
    //     "Pref=(0,1,1) T1gm flavor=27 IDname 2"
    //     "Pref=(0,1,1) flavor=27 IDname 2"
    //     "Pref=(0,1,1) flavor=27 IDname"

GenIrrepOperatorInfo::GenIrrepOperatorInfo(const std::string& opstring)
{
 assign_from_string(opstring);
}


 // ********************************************************************


void GenIrrepOperatorInfo::assign(ArgsHandler& xt)
{
 vector<int> mom(xt.getIntVector("Momentum"));
 bool ref_mom = xt.queryTag("ReferenceMomentum");
 string irrep="NONE";
 xt.getOptionalString("LGIrrep",irrep);
 uint irrepRow=0;
 xt.getOptionalUInt("LGIrrepRow",irrepRow);
 vector<string> flavor(ArgsHandler::split(xt.getString("Flavor"),' '));
 string name(xt.getName("IDName"));
 uint index=0;
 xt.getOptionalUInt("IDIndex",index);
 encode(mom,ref_mom,irrep,irrepRow,flavor,name,index);
}


void GenIrrepOperatorInfo::assign_from_string(const string& opstring)
{
 try{
 string opstr(tidyString(opstring));
 vector<string> tokens=ArgsHandler::split(opstr,' ');
 if ((tokens.size()<3)||(tokens.size()>5)) throw(std::runtime_error(""));
 uint count=0;
 string mom_str = tokens.at(count++);
 size_t pos = mom_str.find("Pref=");
 bool ref_mom = (pos!=string::npos);
 if (ref_mom) mom_str.erase(pos+1,3);
 vector<int> mom;
 momentum_from_string(mom_str,mom);
 string irrep="NONE";
 uint irrepRow=0;
 string flavor_str = tokens.at(count++);
 if (flavor_str.find("flavor=")==string::npos){  // found irrep instead of flavor
    vector<string> irrep_tokens=ArgsHandler::split(flavor_str,'_');
    irrep=irrep_tokens[0];
    if (irrep_tokens.size()==2) extract_from_string(irrep_tokens[1],irrepRow);
    flavor_str = tokens.at(count++);
 }
 pos = flavor_str.find("flavor=");
 flavor_str.erase(pos,7);
 vector<string> flavor=ArgsHandler::split(flavor_str,',');
 string name(tokens.at(count++));
 uint index=0;
 if (tokens.size()==(count+1)){
    extract_from_string(tokens.at(count),index);}
 encode(mom,ref_mom,irrep,irrepRow,flavor,name,index);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("Invalid GenIrrepOperatorInfo string: ")+opstring));}
}



void GenIrrepOperatorInfo::encode(vector<int> mom, bool ref_mom, const string& irrep,
                                  uint irrepRow, const vector<string>& flavor,
                                  const string& name, uint index)
{
 if (mom.size()!=3){
    throw(std::invalid_argument("Bad momentum"));}
 if (((unsigned int)abs(mom[0])>momj_mask)||((unsigned int)abs(mom[1])>momj_mask)
    ||((unsigned int)abs(mom[2])>momj_mask)){
    throw(std::invalid_argument("momentum component magnitude not currently supported"));}
 if (ref_mom){
    mom[0] = abs(mom[0]);
    mom[1] = abs(mom[1]);
    mom[2] = abs(mom[2]);
    sort(mom.begin(), mom.end());}
 uint momcode=(mom[0]<0)?1:0; 
 momcode<<=momj_bits; momcode|=abs(mom[0]);
 momcode<<=1; momcode|=(mom[1]<0)?1:0;
 momcode<<=momj_bits; momcode|=abs(mom[1]);
 momcode<<=1; momcode|=(mom[2]<0)?1:0;
 momcode<<=momj_bits; momcode|=abs(mom[2]);

 uint irrep_code=m_irreps.encode(irrep);
 if (irrep=="NONE") irrepRow=0;

 if (irrepRow>int(irrw_mask)){
    throw(std::invalid_argument("Irrep row not currently supported"));}

 bool is_su3_flavor;
 uint flavcode=0u;
 if (flavor.size()==1){
    is_su3_flavor=true;
    string flav_str = flavor.at(0);
    size_t pos = flav_str.find("*");
    if (pos==(flav_str.size()-1)){
      flavcode=1u;
      flavcode<<=su3flav_bits;
      flav_str.erase(pos,1);}
    int su3_flavor=stoi(flav_str, &pos);
    if (pos!=flav_str.size()){
      throw(std::invalid_argument("Invalid SU(3) flavor irrep"));}
    if ((su3_flavor<=0)||(su3_flavor>int(su3flav_mask))){
      throw(std::invalid_argument("Invalid SU(3) flavor irrep"));}
    flavcode|=(uint)su3_flavor;}
 else if (flavor.size()==2){
   is_su3_flavor=false;
   string iso_str = flavor.at(0);
   size_t pos = iso_str.find("h");
   if (pos!=string::npos) iso_str.erase(pos,1);
   int isospin=stoi(iso_str);
   if (isospin<0) throw(std::invalid_argument("Isospin must be greater than zero"));
   if (pos==string::npos) isospin *= 2;
   if (isospin>int(isop_mask)) throw(std::invalid_argument("Isospin not currently supported"));
   flavcode=(uint)isospin;
   int strangeness=stoi(flavor.at(1));
   if (abs(strangeness)>int(strange_mask)) throw(std::invalid_argument("Strangeness not supported"));
   flavcode<<=1; flavcode|=(strangeness<0)?1:0;
   flavcode<<=strange_bits; flavcode|=abs(strangeness);}
 else {
    throw(std::invalid_argument("Invalid flavor"));}

 const uint maxlength=24;
 if (name.length()>maxlength){
    throw(std::invalid_argument("GIOperator name too long"));}
 if (index>int(irrprwfl_mask)){
    throw(std::invalid_argument("GIOperator index too large"));}
 vector<uint> namecode;
 encode_string_to_uints(name,maxlength,namecode);
 icode.resize(namecode.size()+2);
 icode[0]<<=1; if (is_su3_flavor) icode[0]|=1u;
 icode[0]<<=1; if (ref_mom) icode[0]|=1u;
 icode[0]<<=momt_bits; icode[0]|=momcode;
 icode[0]<<=girr_bits; icode[0]|=girr_mask;
 uint tcode=index;   
 tcode<<=irrp_bits; tcode|=irrep_code;
 tcode<<=irrw_bits; tcode|=irrepRow;
 tcode<<=flav_bits; tcode|=flavcode;

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
 if (!longform){
    xmlout.set_root("GIOperatorString",short_output());
    return;}

 xmlout.set_root("GIOperator");
 string irrep(getLGIrrep());
 uint irrepRow=getLGIrrepRow();
 vector<string> flavor(getFlavor());
 string idname(getIDName());
 uint index=getIDIndex();
 string mom_str;
 Momentum P(getMomentum());
 xmlout.put_child("Momentum",make_string(P.x)+" "+make_string(P.y)+" "+make_string(P.z));
 if (isReferenceMomentum()) xmlout.put_child("ReferenceMomentum");
 if (irrep!="None"){
    xmlout.put_child("LGIrrep",irrep);
    if (irrepRow!=0) xmlout.put_child("LGIrrepRow",make_string(irrepRow));}
 string flavor_str="";
 for (vector<string>::iterator flav_it=flavor.begin(); flav_it!=flavor.end(); ++flav_it){
    flavor_str+=*flav_it+" ";}
 xmlout.put_child("Flavor",flavor_str);
 xmlout.put_child("IDName",idname);
 xmlout.put_child("IDIndex",make_string(index));
}


string GenIrrepOperatorInfo::short_output() const
{
 string irrep(getLGIrrep());
 uint irrepRow=getLGIrrepRow();
 vector<string> flavor(getFlavor());
 string idname(getIDName());
 uint index=getIDIndex();
 string mom_str = "P";
 if (isReferenceMomentum()) mom_str += "ref";
 Momentum P(getMomentum());
 mom_str += "=("+make_string(P.x)+","+make_string(P.y)+","+make_string(P.z)+")";
 string irrep_str="";
 if (irrep!="NONE") {
    irrep_str = " "+irrep;
    if (irrepRow!=0) irrep_str += "_"+make_string(irrepRow);}

 string flavor_str="flavor=";
 for (vector<string>::iterator flav_it=flavor.begin(); flav_it!=flavor.end(); ++flav_it){
    flavor_str+=*flav_it+",";}
 flavor_str = flavor_str.substr(0, flavor_str.size()-1);

 string opstr=mom_str+irrep_str+" "+flavor_str+" "+idname+" "+make_string(index);
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

bool GenIrrepOperatorInfo::isReferenceMomentum() const
{
 uint tcode = icode[0]>>(momt_bits+girr_bits);
 return ((tcode&1u)==1);
}

std::string GenIrrepOperatorInfo::getLGIrrep() const
{
 return m_irreps.decode((icode[1]>>(irrw_bits+flav_bits)) & irrp_mask);
}

unsigned int GenIrrepOperatorInfo::getLGIrrepRow() const
{
 return ((icode[1]>>flav_bits) & irrw_mask);
}

bool GenIrrepOperatorInfo::isSU3flavor() const
{
 uint tcode = icode[0]>>(momt_bits+girr_bits+1);
 return ((tcode&1u)==1);
}

std::vector<std::string> GenIrrepOperatorInfo::getFlavor() const
{
 uint flavcode=(icode[1]&flav_mask);
 vector<string> flavor;
 if (isSU3flavor()){
   string flav_str=to_string((flavcode&su3flav_mask));
   if ((flavcode>>su3flav_bits)==1u) flav_str+="*";
   flavor.push_back(flav_str);}
 else{
   int strangeness=(flavcode & strange_mask);
   flavcode>>=strange_bits;
   if ((flavcode&1u)==1) strangeness=-strangeness;
   flavcode>>=1;
   uint isospin=(flavcode&isop_mask);
   string iso_str;
   if ((isospin%2)==0) iso_str = to_string(isospin/2);
   else iso_str = to_string(isospin)+"h";
   flavor.push_back(iso_str);
   flavor.push_back(to_string(strangeness));
 }
 return flavor;
}

std::string GenIrrepOperatorInfo::getIDName() const
{
 vector<uint>::const_iterator it=icode.begin()+2;
 vector<uint> namecode(it,icode.end());
 return decode_uints_to_string(namecode);
}


unsigned int GenIrrepOperatorInfo::getIDIndex() const
{
 return icode[1]>>irrprwfl_bits; 
}


GenIrrepOperatorInfo& GenIrrepOperatorInfo::resetIDIndex(uint level)
{
 resetIDIndex(level,icode[1]);
 return *this;
}



// **************************************************

       // static data needed
       
GenIrrepOperatorInfo::Encoder GenIrrepOperatorInfo::m_irreps(0);


GenIrrepOperatorInfo::Encoder::Encoder(int ctype)
{
 if (ctype==0) set_irreps();
}

void GenIrrepOperatorInfo::Encoder::set_irreps()
{
 m_codetype="irrep name";
 m_code["NONE"]= 0;   m_string[ 0]="NONE";
 m_code["A1gp"]= 1;   m_string[ 1]="A1gp";
 m_code["A1"  ]= 2;   m_string[ 2]="A1"  ;
 m_code["A1g" ]= 3;   m_string[ 3]="A1g" ;
 m_code["A1gm"]= 4;   m_string[ 4]="A1gm";
 m_code["A1m" ]= 5;   m_string[ 5]="A1m" ;
 m_code["A1p" ]= 6;   m_string[ 6]="A1p" ;
 m_code["A1u" ]= 7;   m_string[ 7]="A1u" ;
 m_code["A1um"]= 8;   m_string[ 8]="A1um";
 m_code["A1up"]= 9;   m_string[ 9]="A1up";
 m_code["A2"  ]=10;   m_string[10]="A2"  ;
 m_code["A2g" ]=11;   m_string[11]="A2g" ;
 m_code["A2gm"]=12;   m_string[12]="A2gm";
 m_code["A2gp"]=13;   m_string[13]="A2gp";
 m_code["A2m" ]=14;   m_string[14]="A2m" ;
 m_code["A2p" ]=15;   m_string[15]="A2p" ;
 m_code["A2u" ]=16;   m_string[16]="A2u" ;
 m_code["A2um"]=17;   m_string[17]="A2um";
 m_code["A2up"]=18;   m_string[18]="A2up";
 m_code["B1"  ]=19;   m_string[19]="B1"  ;
 m_code["B1m" ]=20;   m_string[20]="B1m" ;
 m_code["B1p" ]=21;   m_string[21]="B1p" ;
 m_code["B2"  ]=22;   m_string[22]="B2"  ;
 m_code["B2m" ]=23;   m_string[23]="B2m" ;
 m_code["B2p" ]=24;   m_string[24]="B2p" ;
 m_code["E"   ]=25;   m_string[25]="E"   ;
 m_code["Eg"  ]=26;   m_string[26]="Eg"  ;
 m_code["Egm" ]=27;   m_string[27]="Egm" ;
 m_code["Egp" ]=28;   m_string[28]="Egp" ;
 m_code["Em"  ]=29;   m_string[29]="Em"  ;
 m_code["Ep"  ]=30;   m_string[30]="Ep"  ;
 m_code["Eu"  ]=31;   m_string[31]="Eu"  ;
 m_code["Eum" ]=32;   m_string[32]="Eum" ;
 m_code["Eup" ]=33;   m_string[33]="Eup" ;
 m_code["F1"  ]=34;   m_string[34]="F1"  ;
 m_code["F2"  ]=35;   m_string[35]="F2"  ;
 m_code["G"   ]=36;   m_string[36]="G"   ;
 m_code["G1"  ]=37;   m_string[37]="G1"  ;
 m_code["G1g" ]=38;   m_string[38]="G1g" ;
 m_code["G1u" ]=39;   m_string[39]="G1u" ;
 m_code["G2"  ]=40;   m_string[40]="G2"  ;
 m_code["G2g" ]=41;   m_string[41]="G2g" ;
 m_code["G2u" ]=42;   m_string[42]="G2u" ;
 m_code["Hg"  ]=43;   m_string[43]="Hg"  ;
 m_code["Hu"  ]=44;   m_string[44]="Hu"  ;
 m_code["T1g" ]=45;   m_string[45]="T1g" ;
 m_code["T1gm"]=46;   m_string[46]="T1gm";
 m_code["T1gp"]=47;   m_string[47]="T1gp";
 m_code["T1u" ]=48;   m_string[48]="T1u" ;
 m_code["T1um"]=49;   m_string[49]="T1um";
 m_code["T1up"]=50;   m_string[50]="T1up";
 m_code["T2g" ]=51;   m_string[51]="T2g" ;
 m_code["T2gm"]=52;   m_string[52]="T2gm";
 m_code["T2gp"]=53;   m_string[53]="T2gp";
 m_code["T2u" ]=54;   m_string[54]="T2u" ;
 m_code["T2um"]=55;   m_string[55]="T2um";
 m_code["T2up"]=56;   m_string[56]="T2up";
}

unsigned int GenIrrepOperatorInfo::Encoder::encode(const std::string& description) const
{
 map<string,unsigned int>::const_iterator it=m_code.find(description);
 if (it==m_code.end()) {
    throw(std::invalid_argument(string("Invalid ")+m_codetype
                  +string(" string in GenIrrepOperatorInfo")));}
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


void GenIrrepOperatorInfo::momentum_from_string(const std::string& momstr, std::vector<int>& pmom)
{
 if (momstr.length()<9) throw(std::invalid_argument("Invalid momentum string"));
 if ((momstr[0]!='P')||(momstr[1]!='=')||(momstr[2]!='(')||(momstr[momstr.length()-1]!=')'))
    throw(std::invalid_argument("Invalid momentum string"));
 string mmm(momstr.substr(3,momstr.length()-4));
 vector<string> p=ArgsHandler::split(mmm,',');
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
 if (level>int(irrprwfl_mask)){
    throw(std::invalid_argument("GIOperator index too large"));}
 uint keep=kcode & irrprwfl_mask;  // keep right 21 bits  
 kcode=(level<<irrprwfl_bits);
 kcode|=keep;
}

// ******************************************************************************

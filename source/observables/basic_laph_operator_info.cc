#include "basic_laph_operator_info.h"
#include "multi_compare.h"
#include "args_handler.h"
#include "log_helper.h"
#include <stdexcept>


  //  Encode operator information into several unsigned 32-bit integers
  //
  //   Each single-hadron stored in two 32-bit unsigned integers:
  //      icode1 =  5 bits empty, momentum in 24 bits, 3 bits empty
  //                8 bits each momx, momy, momz (first bit is sign, 7 bits remaining), 
  //      icode2 = 3 bits displacement length, 13 bits spatial id num,
  //               4 bits spatial type, 7 bits irrep, 5 bits flavor
  //
  //  Operator encoding:
  //
  //      icode[2*k-2] contains icode1 of particle k=1,2,3,....
  //      icode[2*k-1] contains icode2 of particle k=1,2,3,....
  //
  //   number_of_hadrons = 0   
  //         icode[0]=0
  //
  //   number_of_hadrons = 1  (flavor code identifies as meson, baryon, tetraquark....)
  //         icode[0]=icode1 with rightmost 3 bits number_of_hadrons
  //                    and leftmost 5 bits containing 1 bit (color bit for tetraquarks, empty otherwise) 
  //                    and 4 bits LGIrrepRow
  //         icode[1]=icode2
  //
  //   number_of_hadrons = 2,3,4,5,6  
  //         icode[0]=icode1 of hadron 1 with rightmost 3 bits number_of_hadrons
  //         icode[1]=icode2 of hadron 1
  //         icode[2]=icode1 of hadron 2
  //         icode[3]=icode2 of hadron 2
  //           ....
  //         icode[2*number_of_hadrons]=total_code
  //
  //    the total_code has the format (left to right)
  //          3 bits empty,  6 bits <Isospin>, 6 bits <IsoCGId>,
  //          7 bits <LGIrrep>, 6 bits <LGCGId>, 4 bits <LGIrrepRow>
  //


using namespace std;



BasicLapHOperatorInfo::BasicLapHOperatorInfo() : icode(1)
{
 icode[0]=0;     // default is zero particles
}


BasicLapHOperatorInfo::BasicLapHOperatorInfo(XMLHandler& xml_in)
{

 try{
 set<string> tags;
 tags.insert("Operator");
 tags.insert("BLOperator");
 tags.insert("OperatorString");
 tags.insert("BLOperatorString");
 ArgsHandler xin(xml_in,tags);
 string rtag=xin.getInputRootTag();

 if ((rtag=="Operator")||(rtag=="BLOperator")){
       // first, read the number of hadrons (particles)
    unsigned int nhadrons;
    xin.getUInt("NumberOfHadrons",nhadrons);
    if (nhadrons>nhad_mask){
       throw(std::invalid_argument("Unsupported value of NumberOfHadrons"));}
    if (nhadrons==0){
       icode.resize(1); 
       icode[0]=0;} 
    else if (nhadrons==1){
       icode.resize(2);
       xmlread_hadron(xin,"Hadron",icode[0],icode[1]);
       icode[0]|=0x1u;
       int LGIrrepRow;
       xin.getInt("LGIrrepRow",LGIrrepRow);
       if ((LGIrrepRow<1)||(LGIrrepRow>int(irrw_mask))){
          throw(std::invalid_argument("Unsupported value of LGIrrepRow"));}
       icode[0]|=LGIrrepRow<<(momt_bits+nhad_bits);}
    else if (nhadrons>=2){
       icode.resize(2*nhadrons+1);
       for (unsigned int k=0;k<nhadrons;++k){
          ostringstream oss; oss << "Hadron" << k+1;
          xmlread_hadron(xin,oss.str(),icode[2*k],icode[2*k+1]);}
       icode[0]|=nhadrons;
       xmlread_total(xin,icode[2*nhadrons]);}
    else
       throw(std::invalid_argument("Unsupported number of hadrons"));
    return;}
 else if (rtag=="OperatorString"){
    assign(xin.getString("OperatorString"));
    return;}
 else if (rtag=="BLOperatorString"){
    assign(xin.getString("BLOperatorString"));
    return;}
 else
    throw(std::invalid_argument("Invalid XML input"));
 }
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("BasicLapHOperatorInfo construction failed: \n")
      +string(errmsg.what())+string("\nInput XML:")+xml_in.output()));}
}


    // converts specification of operator by a string into an BasicLapHOperatorInfo object
    // CG_0 is default, put in CG_1, etc for other LGCG id's

    //   Examples:
    //     "glueball P=(0,0,0) A1gp_1 TrEig"
    //     "pion P=(0,0,0) A1um_1 SD_5"
    //     "isotriplet_pion_pion A1um_1 CG_1 [P=(0,0,1) A1p LSD_1] [P=(0,0,-1) A2m TSD_2] 
    //     "tquuuu1p P=(0,0,0) A1um_1 QDX_1"

BasicLapHOperatorInfo::BasicLapHOperatorInfo(const std::string& opstring)
{
 assign(opstring);
}


 // ********************************************************************


unsigned int BasicLapHOperatorInfo::getNumberOfHadrons() const
{
 return (icode[0] & nhad_mask);   // rightmost 3 bits
}

unsigned int BasicLapHOperatorInfo::get_NumberOfHadrons() const 
{
 return (icode[0] & nhad_mask);   // rightmost 3 bits
}

bool BasicLapHOperatorInfo::isVacuum() const
{
 return (getNumberOfHadrons()==0);
}

bool BasicLapHOperatorInfo::isGlueball() const
{
 return ((getNumberOfHadrons()==1)&&(isGlueball(1)));
}

bool BasicLapHOperatorInfo::isMeson() const
{
 return ((getNumberOfHadrons()==1)&&(isMeson(1)));
}

bool BasicLapHOperatorInfo::isBaryon() const
{
 return ((getNumberOfHadrons()==1)&&(isBaryon(1)));
}

bool BasicLapHOperatorInfo::isTetraquark() const
{
 return ((getNumberOfHadrons()==1)&&(isTetraquark(1)));
}

bool BasicLapHOperatorInfo::isMesonMeson() const
{
 return ((getNumberOfHadrons()==2)&&(isMeson(1))&&(isMeson(2)));
}

bool BasicLapHOperatorInfo::isMesonBaryon() const
{
 return ((getNumberOfHadrons()==2)&&(isMeson(1))&&(isBaryon(2)));
}

bool BasicLapHOperatorInfo::isBaryonBaryon() const
{
 return ((getNumberOfHadrons()==2)&&(isBaryon(1))&&(isBaryon(2)));
}

Momentum BasicLapHOperatorInfo::getMomentum() const
{
 Momentum p(0,0,0),padd;
 for (unsigned int k=1;k<=get_NumberOfHadrons();k++){
    padd=getMomentum(k);
    p.x+=padd.x;
    p.y+=padd.y;
    p.z+=padd.z;}
 return p;
}

int BasicLapHOperatorInfo::getXMomentum() const
{
 int p=0;
 for (unsigned int k=1;k<=get_NumberOfHadrons();k++) 
    p+=getXMomentum(k);
 return p;
}

int BasicLapHOperatorInfo::getYMomentum() const
{
 int p=0;
 for (unsigned int k=1;k<=get_NumberOfHadrons();k++) 
    p+=getYMomentum(k);
 return p;
}

int BasicLapHOperatorInfo::getZMomentum() const
{
 int p=0;
 for (unsigned int k=1;k<=get_NumberOfHadrons();k++) 
    p+=getZMomentum(k);
 return p;
}

std::string BasicLapHOperatorInfo::getLGIrrep() const
{
 unsigned int nhadrons=get_NumberOfHadrons();
 if (nhadrons==0) return "A1gp";
 else if (nhadrons==1) return getLGIrrep(1);
 return m_irreps.decode((icode[2*nhadrons]>>(lgcg_bits+irrw_bits)) & irrp_mask);
}

unsigned int BasicLapHOperatorInfo::getLGClebschGordonIdNum() const
{
 unsigned int nhadrons=get_NumberOfHadrons();
 if (nhadrons<2) return 0;
 return (icode[2*nhadrons]>>irrw_bits) & lgcg_mask;
}

unsigned int BasicLapHOperatorInfo::getLGIrrepRow() const
{
 unsigned int nhadrons=get_NumberOfHadrons();
 if (nhadrons==0) return 0;
 else if (nhadrons==1) return ((icode[0]>>(nhad_bits+momt_bits))&irrw_mask);
 return (icode[2*nhadrons] & irrw_mask);
}

std::string BasicLapHOperatorInfo::getIsospin() const
{
 unsigned int nhadrons=get_NumberOfHadrons();
 if (nhadrons==0) return "singlet";
 else if (nhadrons==1) return get_flavor(icode[1]);
 return m_isospin.decode((icode[2*nhadrons]>>
       (iscg_bits+irrp_bits+lgcg_bits+irrw_bits)) & isop_mask);
}

unsigned int BasicLapHOperatorInfo::getIsospinClebschGordonIdNum() const
{
 unsigned int nhadrons=get_NumberOfHadrons();
 if (nhadrons<2) return 0;
 return ((icode[2*nhadrons]>>(irrp_bits+lgcg_bits+irrw_bits)) & iscg_mask);
}


std::string BasicLapHOperatorInfo::getFlavor() const
{
 unsigned int nhadrons=get_NumberOfHadrons();
 if (nhadrons==0) return "singlet";
 else if (nhadrons==1) return get_flavor(icode[1]);
 string flav(getIsospin());
 for (uint k=1;k<=nhadrons;++k) flav+="_"+getFlavor(k);
 return flav;
}

std::string BasicLapHOperatorInfo::getFlavorCode() const
{
 unsigned int nhadrons=get_NumberOfHadrons();
 if (nhadrons==0) return "0";
 else if (nhadrons==1) return get_flavorcode(icode[1]);
 string flav;
 for (uint k=1;k<=nhadrons;++k) flav+=get_flavorcode(icode[2*k-1]);
 flav+=m_isochar.decode((icode[2*nhadrons]>>
       (iscg_bits+irrp_bits+lgcg_bits+irrw_bits)) & isop_mask);
 return flav;
}

int BasicLapHOperatorInfo::getStrangeness() const
{
 unsigned int nhadrons=get_NumberOfHadrons();
 if (nhadrons==0) return 0;
 int strangeness = 0;
 for (uint k=1;k<=nhadrons;++k) strangeness+=getStrangeness(k);
 return strangeness;
}

int BasicLapHOperatorInfo::getTetraquarkColorType() const
{
 if (!isTetraquark()){
    cerr << "Invalid call to getTetraquarkColorType in BasicLapHOperatorInfo::"<<endl;
    throw(std::runtime_error("not Tetraquark in getTetraquarkColorType"));}
 return getTetraquarkColorType(1);
}


          // note: hadron_index = 1, 2, ...
         
std::string BasicLapHOperatorInfo::getFlavor(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"getFlavor");
 return get_flavor(icode[2*hadron_index-1]);
}

int BasicLapHOperatorInfo::getStrangeness(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"getStrangeness");
 string flav(get_flavor(icode[2*hadron_index-1]));
 if (flav==string("kaon")) return 1;
 else if ((flav==string("kbar"))||(flav==string("lambda"))||(flav==string("sigma"))) return -1;
 else if (flav==string("xi")) return -2;
 else if (flav==string("omega")) return -3;
 else return 0;
}

bool BasicLapHOperatorInfo::isGlueball(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"isGlueball");
 return is_glueball(icode[2*hadron_index-1]);
}

bool BasicLapHOperatorInfo::isMeson(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"isMeson");
 return is_meson(icode[2*hadron_index-1]);
}

bool BasicLapHOperatorInfo::isBaryon(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"isBaryon");
 return is_baryon(icode[2*hadron_index-1]);
}

bool BasicLapHOperatorInfo::isTetraquark(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"isTetraquark");
 return is_tetraquark(icode[2*hadron_index-1]);
}

bool BasicLapHOperatorInfo::isFermion(unsigned int hadron_index) const
 {
 check_hadron_index(hadron_index,"isFermion");
 return is_fermion(icode[2*hadron_index-1]);
}
  
bool BasicLapHOperatorInfo::isBoson(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"isBoson");
 return is_boson(icode[2*hadron_index-1]);
}

std::string BasicLapHOperatorInfo::getLGIrrep(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"getLGIrrep");
 return get_lgirrep(icode[2*hadron_index-1]);
}

std::string BasicLapHOperatorInfo::getSpatialType(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"getSpatialType");
 return get_spatial_type(icode[2*hadron_index-1]);
}

unsigned int BasicLapHOperatorInfo::getSpatialIdNumber(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"getSpatialIdNumber");
 return get_spatial_id(icode[2*hadron_index-1]);
}

unsigned int BasicLapHOperatorInfo::getDisplacementLength(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"getDisplacementLength");
 return get_disp_length(icode[2*hadron_index-1]);
}


Momentum BasicLapHOperatorInfo::getMomentum(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"getMomentum");
 return get_momentum(icode[2*hadron_index-2]);
}

int BasicLapHOperatorInfo::getXMomentum(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"getXMomentum");
 return get_xmomentum(icode[2*hadron_index-2]);
}

int BasicLapHOperatorInfo::getYMomentum(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"getYMomentum");
 return get_ymomentum(icode[2*hadron_index-2]);
}

int BasicLapHOperatorInfo::getZMomentum(unsigned int hadron_index) const
{
 check_hadron_index(hadron_index,"getZMomentum");
 return get_zmomentum(icode[2*hadron_index-2]);
}

int BasicLapHOperatorInfo::getTetraquarkColorType(unsigned int hadron_index) const
{
 if (!isTetraquark(hadron_index)){
    cerr << "Invalid call to getTetraquarkColorType in BasicLapHOperatorInfo::"<<endl;
    throw(std::runtime_error("not Tetraquark in getTetraquarkColorType"));}
 return 1-2*(icode[2*hadron_index-2]>>31);
}


string BasicLapHOperatorInfo::output(bool longform, int indent) const
{
 XMLHandler xmlout;
 output(xmlout,longform);
 return xmlout.output(indent);
}

string BasicLapHOperatorInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}

void BasicLapHOperatorInfo::output(XMLHandler& xmlout, bool longform) const
{
 if (!longform){
    xmlout.set_root("BLOperatorString",short_output());
    return;}
 xmlout.set_root("BLOperator");
 unsigned int nhadrons=getNumberOfHadrons();
 xmlout.put_child("NumberOfHadrons",make_string(nhadrons));
 if (nhadrons==0) return;
 else if (nhadrons==1){
    xmlwrite_hadron(xmlout,"Hadron",icode[0],icode[1]);
    xmlout.put_child("LGIrrepRow",make_string(getLGIrrepRow()));}
 else{
     xmlwrite_total(xmlout,icode[2*nhadrons]);
    for (unsigned int k=0;k<nhadrons;++k){
       ostringstream oss; oss << "Hadron" << k+1;
       xmlwrite_hadron(xmlout,oss.str(),icode[2*k],icode[2*k+1]);}}
}

bool BasicLapHOperatorInfo::operator==(const BasicLapHOperatorInfo& rhs) const
{
 return multiEqual(icode,rhs.icode);   
}

bool BasicLapHOperatorInfo::operator!=(const BasicLapHOperatorInfo& rhs) const
{
 return multiNotEqual(icode,rhs.icode);   
}

bool BasicLapHOperatorInfo::operator<(const BasicLapHOperatorInfo& rhs) const
{
 return multiLessThan(icode,rhs.icode);   
}


 // *******************************************************************



void BasicLapHOperatorInfo::check_hadron_index(unsigned int hadron_index, 
                                      const char *msg) const
{
 uint nhad=get_NumberOfHadrons();
 if ((hadron_index<1)||(hadron_index>nhad)){
    //cerr << "Invalid hadron index in BasicLapHOperatorInfo::"<<msg<<endl;
    throw(std::invalid_argument("Invalid hadron index"));}
}

   //  reads the momentum tag using tag name "tag", putting
   //  result in "mom" and encodes it into 24 bits in "momcode", then
   //  bitshift to left by 3

void BasicLapHOperatorInfo::xmlread_momentum(ArgsHandler& xmlh, const string& tag,
                                             Momentum& mom, unsigned int& momcode)
{
 vector<int> p(xmlh.getIntVector(tag));
 if (p.size()!=3){
    throw(std::invalid_argument("Bad momentum"));}
 mom.x=p[0];
 mom.y=p[1];
 mom.z=p[2];
 if (((unsigned int)abs(mom.x)>momj_mask)||((unsigned int)abs(mom.y)>momj_mask)
    ||((unsigned int)abs(mom.z)>momj_mask)){
    throw(std::invalid_argument("momentum component magnitude not currently supported"));}
 momcode=(mom.x<0)?1:0; 
 momcode<<=momj_bits; momcode|=abs(mom.x);
 momcode<<=1; momcode|=(mom.y<0)?1:0;
 momcode<<=momj_bits; momcode|=abs(mom.y);
 momcode<<=1; momcode|=(mom.z<0)?1:0;
 momcode<<=momj_bits; momcode|=abs(mom.z);
 momcode<<=nhad_bits;
}
    

void BasicLapHOperatorInfo::xmlwrite_momentum(XMLHandler& xmlout, const string& tag,
                                              unsigned int momcode) const
{
 Momentum mom=get_momentum(momcode);
 vector<int> p(3);  
 p[0]=mom.x; p[1]=mom.y; p[2]=mom.z;
 xmlout.put_child(tag,make_string(p));
}


     
void BasicLapHOperatorInfo::xmlread_hadron(ArgsHandler& xin, const string& toptag,
                                           unsigned int& momcode, unsigned int& hadroncode)
{
 ArgsHandler xh(xin,toptag);

 string flav(xh.getString("Flavor"));
 unsigned int fcode=m_flavor.encode(flav);

 uint colorcode=0;
 if ((fcode>=13)&&(fcode<=24)){  // is tetraquark
    int ctype=1;
    xh.getInt("ColorType",ctype);
    if (ctype==-1) colorcode=1;
    else if (ctype==1) colorcode=0;
    else throw(std::invalid_argument("Bad hadron operator input xml data: tetraquark with no ColorType"));}

 Momentum mom;
 xmlread_momentum(xh,"Momentum",mom,momcode);
 if (colorcode==1) momcode|=(colorcode<<31);

 string irrep(xh.getString("LGIrrep"));
 unsigned int irrepcode=m_irreps.encode(irrep);

 string spatialType(xh.getString("SpatialType"));
 unsigned int spcode=m_spatial.encode(spatialType);

 int spatialIdNum, dispLength; 
 if (is_glueball(fcode)){  // glueball
    spatialIdNum=0; dispLength=0;}
 else{
    xh.getInt("SpatialIdNum",spatialIdNum);
    xh.getInt("DispLength",dispLength);
    if ((dispLength<0)||(dispLength>int(dlen_mask))
       ||(spatialIdNum<0)||(spatialIdNum>int(spid_mask))){
       throw(std::invalid_argument("Bad hadron operator input xml data"));}
    if ((spatialType!="SS")&&(spatialType!="VI")&&(dispLength==0)){
       throw(std::invalid_argument("Bad hadron operator input xml data"));}
       // check if a single-site operator, then make the displacement length zero
    if ((spatialType=="SS")||(spatialType=="VI")) dispLength=0;}

    //  now do the encoding
 hadroncode=dispLength; 
 hadroncode<<=spid_bits;  hadroncode|=spatialIdNum; 
 hadroncode<<=sptp_bits;  hadroncode|=spcode; 
 hadroncode<<=irrp_bits;  hadroncode|=irrepcode;
 hadroncode<<=flav_bits;  hadroncode|=fcode;

    // do some rudimentary checks (not exhaustive)
 if (is_boson(hadroncode)){
    if (!((irrep[0]=='A')||(irrep[0]=='B')||(irrep[0]=='E')||(irrep[0]=='T'))){
       throw(std::invalid_argument("Bad bosonic irrep in input xml data"));}}
 else if (is_fermion(hadroncode)){
    if (!((irrep[0]=='G')||(irrep[0]=='F')||(irrep[0]=='H'))){
       throw(std::invalid_argument("Bad fermionic irrep in input xml data"));}}
 bool spcheck=true;
 if (spcode!=12){
    if (is_glueball(hadroncode))
       spcheck=((spcode>=9)&&(spcode<=11))||(spcode==14);
    else if (is_meson(hadroncode))
       spcheck=((spcode>=0)&&(spcode<=3))||(spcode==5)||(spcode==7)||(spcode==8);
    else if (is_baryon(hadroncode))
       spcheck=((spcode>=0)&&(spcode<=2))||((spcode>=4)&&(spcode<=6))||(spcode==8);
    else if (is_tetraquark(hadroncode))
       spcheck=(spcode==0)||(spcode==4)||(spcode==13);}
 if (!spcheck){
    throw(std::invalid_argument("Bad spatial type in input xml data"));}
}

void BasicLapHOperatorInfo::xmlwrite_hadron(XMLHandler& xmlout, const string& toptag,
                                            const unsigned int& momcode, 
                                            const unsigned int& hadroncode) const
{
 XMLHandler xmlhad;
 xmlhad.set_root(toptag);
 Momentum mom=get_momentum(momcode);

 unsigned int fcode=hadroncode;
 unsigned int ircode=fcode>>flav_bits;
 unsigned int spcode=ircode>>irrp_bits;
 unsigned int idcode=spcode>>sptp_bits;
 unsigned int dlcode=idcode>>spid_bits;
 fcode&=flav_mask; 
 ircode&=irrp_mask;
 spcode&=sptp_mask;
 idcode&=spid_mask;
 dlcode&=dlen_mask;

 xmlhad.put_child("Flavor",m_flavor.decode(fcode));
 vector<int> p(3);  
 p[0]=mom.x; p[1]=mom.y; p[2]=mom.z;
 xmlhad.put_child("Momentum",make_string(p));
 xmlhad.put_child("LGIrrep",m_irreps.decode(ircode));
 xmlhad.put_child("SpatialType",m_spatial.decode(spcode));
 if (!is_glueball(fcode)){
    xmlhad.put_child("SpatialIdNum",make_string(idcode));
    xmlhad.put_child("DispLength",make_string(dlcode));}
 if (is_tetraquark(fcode)){
    int ctype=1-2*int(momcode>>31);
    xmlhad.put_child("ColorType",make_string(ctype));}
 xmlout.put_child(xmlhad);
}

  //    the total_code has the format (left to right)
  //          4 bits empty,  6 bits <Isospin>, 6 bits <IsoCGId>,
  //          6 bits <LGIrrep>, 6 bits <LGCGId>, 4 bits <LGIrrepRow>

void BasicLapHOperatorInfo::xmlread_total(ArgsHandler& xin, unsigned int& code)
{
 ArgsHandler xt(xin,"Total");

 string isostr(xt.getString("Isospin"));
 unsigned int isocode=m_isospin.encode(isostr);

 int isoCGid=0;
 xt.getOptionalInt("IsoCGId",isoCGid);
 if ((isoCGid<0)||(isoCGid>int(iscg_mask))){
    throw(std::invalid_argument("Requested value of isospin CG id not supported"));}
 
 string irrep(xt.getString("LGIrrep"));
 unsigned int irrep_code=m_irreps.encode(irrep);
 
 int lgCGid=0;
 xt.getOptionalInt("LGCGId",lgCGid);
 if ((lgCGid<0)||(lgCGid>int(lgcg_mask))){
    throw(std::invalid_argument("Requested value of little group CG id not supported"));}

 int irrepRow;
 xt.getInt("LGIrrepRow",irrepRow);
 if ((irrepRow<=0)||(irrepRow>int(irrw_mask))){
    throw(std::invalid_argument("Irrep row not currently supported"));}

 Momentum mom;
 unsigned int momcode;
 xmlread_momentum(xt,"Momentum",mom,momcode);
 Momentum check=getMomentum();
 if ((check.x!=mom.x)||(check.y!=mom.y)||(check.z!=mom.z)){
    throw(std::invalid_argument("Total momentum does not equal sum of hadron momenta"));}

 code=isocode;
 code<<=iscg_bits;  code|=isoCGid;
 code<<=irrp_bits;  code|=irrep_code;
 code<<=lgcg_bits;  code|=lgCGid;
 code<<=irrw_bits;  code|=irrepRow;
} 


void BasicLapHOperatorInfo::xmlwrite_total(XMLHandler& xmlout, const unsigned int& code) const
{
 XMLHandler xmlt;
 xmlt.set_root("Total");
 unsigned int irreprow=code;
 unsigned int lgcgid=irreprow>>irrw_bits;
 unsigned int ircode=lgcgid>>lgcg_bits;
 unsigned int isocgid=ircode>>irrp_bits;
 unsigned int isocode=isocgid>>iscg_bits;
 irreprow&=irrw_mask;
 lgcgid&=lgcg_mask;
 ircode&=irrp_mask;
 isocgid&=iscg_mask;
 isocode&=isop_mask;

 xmlt.put_child("Isospin",m_isospin.decode(isocode));
 if (isocgid>0) xmlt.put_child("IsoCGId",make_string(isocgid));

 Momentum mom=getMomentum();
 vector<int> p(3);  
 p[0]=mom.x; p[1]=mom.y; p[2]=mom.z;
 xmlt.put_child("Momentum",make_string(p));

 xmlt.put_child("LGIrrep",m_irreps.decode(ircode));
 if (lgcgid>0) xmlt.put_child("LGCGId",make_string(lgcgid));
 xmlt.put_child("LGIrrepRow",make_string(irreprow));
 xmlout.put_child(xmlt);
}

  // ***************************************************************

Momentum BasicLapHOperatorInfo::get_momentum(unsigned int momcode) const
{   
 unsigned int tmp=(momcode>>nhad_bits);
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

int BasicLapHOperatorInfo::get_xmomentum(unsigned int momcode) const
{
 unsigned int tmp=(momcode>>(nhad_bits+2*momj_bits+2));
 int res=tmp & momj_mask;
 if (((tmp>>momj_bits)&0x1u)==1) return -res;
 return res;
}
  
int BasicLapHOperatorInfo::get_ymomentum(unsigned int momcode) const
{
 unsigned int tmp=(momcode>>(nhad_bits+momj_bits+1));
 int res=tmp & momj_mask;
 if (((tmp>>momj_bits)&0x1u)==1) return -res;
 return res;
}

int BasicLapHOperatorInfo::get_zmomentum(unsigned int momcode) const
{
 unsigned int tmp=(momcode>>nhad_bits);
 int res=tmp & momj_mask;
 if (((tmp>>momj_bits)&0x1u)==1) return -res;
 return res;
}


std::string BasicLapHOperatorInfo::get_flavor(unsigned int hadroncode) const
{
 unsigned int fcode=hadroncode & flav_mask;
 return m_flavor.decode(fcode);
}

std::string BasicLapHOperatorInfo::get_flavorcode(unsigned int hadroncode) const
{
 unsigned int fcode=hadroncode & flav_mask;
 return m_flavorcode.decode(fcode);
}

bool BasicLapHOperatorInfo::is_glueball(unsigned int hadroncode) const
{
 unsigned int fcode=hadroncode & flav_mask;
 return (fcode==1);
}

bool BasicLapHOperatorInfo::is_meson(unsigned int hadroncode) const
{
 unsigned int fcode=hadroncode & flav_mask;
 return ((fcode>=2)&&(fcode<=6));
}

bool BasicLapHOperatorInfo::is_baryon(unsigned int hadroncode) const
{
 unsigned int fcode=hadroncode & flav_mask;
 return ((fcode>=7)&&(fcode<=12));
}

bool BasicLapHOperatorInfo::is_tetraquark(unsigned int hadroncode) const
{
 unsigned int fcode=hadroncode & flav_mask;
 return ((fcode>=13)&&(fcode<=24));
}

bool BasicLapHOperatorInfo::is_fermion(unsigned int hadroncode) const
{
 unsigned int fcode=hadroncode & flav_mask;
 return ((fcode>=7)&&(fcode<=12));
}

bool BasicLapHOperatorInfo::is_boson(unsigned int hadroncode) const
{
 unsigned int fcode=hadroncode & flav_mask;
 return (fcode==1) || ((fcode>=2)&&(fcode<=6));
}

std::string BasicLapHOperatorInfo::get_lgirrep(unsigned int hadroncode) const
{
 unsigned int irrepcode=(hadroncode>>flav_bits) & irrp_mask;
 return m_irreps.decode(irrepcode);
}

std::string BasicLapHOperatorInfo::get_spatial_type(unsigned int hadroncode) const
{
 unsigned int spcode=(hadroncode>>(irrp_bits+flav_bits)) & sptp_mask;
 return m_spatial.decode(spcode);
}

unsigned int BasicLapHOperatorInfo::get_spatial_id(unsigned int hadroncode) const
{
 return (hadroncode>>(sptp_bits+irrp_bits+flav_bits)) & spid_mask;
}

unsigned int BasicLapHOperatorInfo::get_disp_length(unsigned int hadroncode) const
{
 return (hadroncode>>(spid_bits+sptp_bits+irrp_bits+flav_bits)) & dlen_mask;
}

// **************************************************

       // static data needed
       
BasicLapHOperatorInfo::Encoder BasicLapHOperatorInfo::m_flavor(0);
BasicLapHOperatorInfo::Encoder BasicLapHOperatorInfo::m_irreps(1);
BasicLapHOperatorInfo::Encoder BasicLapHOperatorInfo::m_isospin(2);
BasicLapHOperatorInfo::Encoder BasicLapHOperatorInfo::m_spatial(3);
BasicLapHOperatorInfo::Encoder BasicLapHOperatorInfo::m_flavorcode(4);
BasicLapHOperatorInfo::Encoder BasicLapHOperatorInfo::m_isochar(5);


BasicLapHOperatorInfo::Encoder::Encoder(int ctype)
{
 if      (ctype==0) set_flavor();
 else if (ctype==1) set_irreps();
 else if (ctype==2) set_isospin();
 else if (ctype==3) set_spatial();
 else if (ctype==4) set_flavorcode();
 else if (ctype==5) set_isospin_char();
}

void BasicLapHOperatorInfo::Encoder::set_flavor()
{
 m_codetype="flavor";
     // zero quark operators
 m_code["glueball"]= 1;  m_string[1]="glueball";
     // quark-antiquark operators
 m_code["pion"]= 2;  m_string[2]="pion";   // this means isovector du           
 m_code["eta" ]= 3;  m_string[3]="eta" ;   // this means isoscalar (uu+dd)
 m_code["phi" ]= 4;  m_string[4]="phi" ;   // this means isoscalar ss
 m_code["kaon"]= 5;  m_string[5]="kaon";   // strangeness = 1                   
 m_code["kbar"]= 6;  m_string[6]="kbar";   // strangeness = -1                  
     // three-quark operators
 m_code["nucleon"]= 7;  m_string[ 7]="nucleon";
 m_code["delta"  ]= 8;  m_string[ 8]="delta"  ;
 m_code["sigma"  ]= 9;  m_string[ 9]="sigma"  ;
 m_code["lambda" ]=10;  m_string[10]="lambda" ;
 m_code["xi"     ]=11;  m_string[11]="xi"     ;
 m_code["omega"  ]=12;  m_string[12]="omega"  ;
     // four-quark operators
 m_code["isosinglet_eta_eta"  ]=13;   m_string[13]="isosinglet_eta_eta"  ;
 m_code["isotriplet_eta_pion" ]=14;   m_string[14]="isotriplet_eta_pion" ;
 m_code["isosinglet_pion_pion"]=15;   m_string[15]="isosinglet_pion_pion";
 m_code["isotriplet_pion_pion"]=16;   m_string[16]="isotriplet_pion_pion";
 m_code["isoquintet_pion_pion"]=17;   m_string[17]="isoquintet_pion_pion";
 m_code["isodoublet_kaon_eta" ]=18;   m_string[18]="isodoublet_kaon_eta" ;
 m_code["isodoublet_kaon_pion"]=19;   m_string[19]="isodoublet_kaon_pion";
 m_code["isoquartet_kaon_pion"]=20;   m_string[20]="isoquartet_kaon_pion";
 m_code["isotriplet_phi_pion" ]=21;   m_string[21]="isotriplet_phi_pion" ;
 m_code["isosinglet_eta_phi"  ]=22;   m_string[22]="isosinglet_eta_phi"  ;
 m_code["isodoublet_kaon_phi" ]=23;   m_string[23]="isodoublet_kaon_phi" ;
 m_code["isosinglet_phi_phi"  ]=24;   m_string[24]="isosinglet_phi_phi"  ;
    // short forms
 m_code["tquuuu1p"]=13;  m_code["tquuuu1m"]=13;   // last integer is isospin
 m_code["tquudu3p"]=14;  m_code["tquudu3m"]=14;   //  1=singlet, 2=doublet,...
 m_code["tqdudu1p"]=15;  m_code["tqdudu1m"]=15; 
 m_code["tqdudu3p"]=16;  m_code["tqdudu3m"]=16; 
 m_code["tqdudu5p"]=17;  m_code["tqdudu5m"]=17; 
 m_code["tqsuuu2p"]=18;  m_code["tqsuuu2m"]=18; 
 m_code["tqsudu2p"]=19;  m_code["tqsudu2m"]=19; 
 m_code["tqsudu4p"]=20;  m_code["tqsudu4m"]=20; 
 m_code["tqssdu3p"]=21;  m_code["tqssdu3m"]=21; 
 m_code["tquuss1p"]=22;  m_code["tquuss1m"]=22; 
 m_code["tqsuss2p"]=23;  m_code["tqsuss2m"]=23; 
 m_code["tqssss1p"]=24;  m_code["tqssss1m"]=24; 

 m_string[25]="tquuuu1p";  m_string[37]="tquuuu1m"; 
 m_string[26]="tquudu3p";  m_string[38]="tquudu3m"; 
 m_string[27]="tqdudu1p";  m_string[39]="tqdudu1m"; 
 m_string[28]="tqdudu3p";  m_string[40]="tqdudu3m"; 
 m_string[29]="tqdudu5p";  m_string[41]="tqdudu5m"; 
 m_string[30]="tqsuuu2p";  m_string[42]="tqsuuu2m"; 
 m_string[31]="tqsudu2p";  m_string[43]="tqsudu2m"; 
 m_string[32]="tqsudu4p";  m_string[44]="tqsudu4m"; 
 m_string[33]="tqssdu3p";  m_string[45]="tqssdu3m"; 
 m_string[34]="tquuss1p";  m_string[46]="tquuss1m"; 
 m_string[35]="tqsuss2p";  m_string[47]="tqsuss2m"; 
 m_string[36]="tqssss1p";  m_string[48]="tqssss1m"; 
}

void BasicLapHOperatorInfo::Encoder::set_flavorcode()
{
 m_codetype="FC";
     // zero quark operators
 m_code["G"]= 1;  m_string[1]="G";   // glueball
     // quark-antiquark operators
 m_code["P"]= 2;  m_string[2]="P";   // pion - this means isovector du           
 m_code["E"]= 3;  m_string[3]="E";   // eta -- this means isoscalar (uu+dd)
 m_code["F"]= 4;  m_string[4]="F";   // phi -- this means isoscalar ss
 m_code["K"]= 5;  m_string[5]="K";   // kaon -- strangeness = 1                   
 m_code["k"]= 6;  m_string[6]="k";   // kbar -- strangeness = -1                  
     // three-quark operators
 m_code["N"]= 7;  m_string[ 7]="N";    // nucleon
 m_code["D"]= 8;  m_string[ 8]="D";    // delta
 m_code["S"]= 9;  m_string[ 9]="S";    // sigma
 m_code["L"]=10;  m_string[10]="L";    // lambda
 m_code["X"]=11;  m_string[11]="X";    // xi
 m_code["W"]=12;  m_string[12]="W";    // omega
     // four-quark operators
 m_code["EE0"]=13;   m_string[13]="EE0";
 m_code["EP2"]=14;   m_string[14]="EP2";
 m_code["PP0"]=15;   m_string[15]="PP0";
 m_code["PP2"]=16;   m_string[16]="PP2";
 m_code["PP4"]=17;   m_string[17]="PP4";
 m_code["KE1"]=18;   m_string[18]="KE1";
 m_code["KP1"]=19;   m_string[19]="KP1";
 m_code["KP3"]=20;   m_string[20]="KP3";
 m_code["FP2"]=21;   m_string[21]="FP2";
 m_code["EF0"]=22;   m_string[22]="EF0";
 m_code["KF1"]=23;   m_string[23]="KF1";
 m_code["FF0"]=24;   m_string[24]="FF0";
}

void BasicLapHOperatorInfo::Encoder::set_spatial()
{
 m_codetype="spatial type";
 m_code["SS"]=0;       m_string[0]="SS";
 m_code["SD"]=1;       m_string[1]="SD";
 m_code["LSD"]=2;      m_string[2]="LSD";
 m_code["TSD"]=3;      m_string[3]="TSD";
 m_code["DDI"]=4;      m_string[4]="DDI";
 m_code["DDL"]=5;      m_string[5]="DDL";
 m_code["TDT"]=6;      m_string[6]="TDT";
 m_code["TDU"]=7;      m_string[7]="TDU";
 m_code["TDO"]=8;      m_string[8]="TDO";
 m_code["TrEig"]=9;    m_string[9]="TrEig";
 m_code["TrW1Eig"]=10; m_string[10]="TrW1Eig";
 m_code["TrW2Eig"]=11; m_string[11]="TrW2Eig";
 m_code["VI"]=12;      m_string[12]="VI";    // variationally improved
 m_code["QDX"]=13;     m_string[13]="QDX";
 m_code["Plaq"]=14;    m_string[14]="Plaq";
}




void BasicLapHOperatorInfo::Encoder::set_irreps()
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

void BasicLapHOperatorInfo::Encoder::set_isospin()
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

void BasicLapHOperatorInfo::Encoder::set_isospin_char()
{
 m_codetype="isospin char";
     // Total isospin
 m_code["0"]=0;   m_string[0]="0";
 m_code["1"]=1;   m_string[1]="1";
 m_code["2"]=2;   m_string[2]="2";
 m_code["3"]=3;   m_string[3]="3";
 m_code["4"]=4;   m_string[4]="4";
 m_code["5"]=5;   m_string[5]="5";
}


unsigned int BasicLapHOperatorInfo::Encoder::encode(const std::string& description) const
{
 map<string,unsigned int>::const_iterator it=m_code.find(description);
 if (it==m_code.end()) 
    throw(std::invalid_argument(string("Invalid ")+m_codetype
                  +string(" string in BasicLapHOperatorInfo")));
 return it->second;
}
 
std::string BasicLapHOperatorInfo::Encoder::decode(unsigned int code) const
{
 map<unsigned int,string>::const_iterator it=m_string.find(code);
 if (it==m_string.end()) 
    throw(std::invalid_argument(string("Invalid ")+m_codetype
             +string(" code in BasicLapHOperatorInfo")));
 return it->second;
}


int BasicLapHOperatorInfo::count(const string& astr, char delimiter) const
{
 int cnt=0;
 size_t pos=astr.find(delimiter);
 while (pos!=string::npos){
    cnt++;
    pos=astr.find(delimiter,pos+1);}
 return cnt;
}


vector<string> BasicLapHOperatorInfo::split(const string& astr, char delimiter) const
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

string BasicLapHOperatorInfo::extract(const string& astr, char left, char right) const
{
 size_t lpos=astr.find_first_of(left);
 size_t rpos=(lpos==string::npos)?string::npos:astr.find_first_of(right,lpos+1);
 if (rpos!=string::npos){
    return astr.substr(lpos+1,rpos-lpos-1);}
 return string("");
}


    // converts specification of operator by a string into an BasicLapHOperatorInfo object
    // CG_0 is default, put in CG_1, etc for other LGCG id's

    //   Examples:
    //     "glueball P=(0,0,0) A1gp_1 TrEig"
    //     "pion P=(0,0,0) A1um_1 SD_5"
    //     "isotriplet_pion_pion A1um_1 CG_1 [P=(0,0,1) A1p LSD_1] [P=(0,0,-1) A2m TSD_2] 
    //     "tquuuu1p P=(0,0,0) A1um_1 QDX_1"


void BasicLapHOperatorInfo::assign(const std::string& opstring)
{
 try{
 string opstr(tidyString(opstring));
 size_t pos=opstr.find_first_of(" ");
 if (pos==string::npos)
    throw(std::invalid_argument("Invalid operator string"));
 string optype=opstr.substr(0,pos);
 int nus=count(optype,'_');

 if (optype=="glueball"){

    vector<string> tokens=split(opstr,' ');
    if (tokens.size()!=4) throw(std::invalid_argument("Invalid glueball operator string"));
    string momstr(tokens[1]);
    vector<string> tk=split(tokens[2],'_');
    if (tk.size()!=2) throw(std::invalid_argument("Invalid glueball operator string"));
    string irrep(tk[0]),irreprow(tk[1]),sptype(tokens[3]);
    icode.resize(2);
    encode_momentum(momstr,icode[0]);
    encode_hadron(optype,irrep,sptype,"0",icode[1]);
    icode[0]|=0x1u;
    int LGIrrepRow;
    extract_from_string(irreprow,LGIrrepRow);
    if ((LGIrrepRow<1)||(LGIrrepRow>int(irrw_mask))){
       throw(std::invalid_argument("Unsupported value of LGIrrepRow"));}
    icode[0]|=LGIrrepRow<<(momt_bits+nhad_bits);}

 else if (nus==0){

    vector<string> tokens=split(opstr,' ');
    if (tokens.size()!=4) throw(std::invalid_argument("Invalid single hadron operator string"));
    string flav(tokens[0]),momstr(tokens[1]);
    vector<string> tk=split(tokens[2],'_');
    if (tk.size()!=2) throw(std::invalid_argument("Invalid single hadron operator string"));
    string irrep(tk[0]),irreprow(tk[1]);
    tk=split(tokens[3],'_');
    if (tk.size()!=2) throw(std::invalid_argument("Invalid single hadron operator string"));
    string sptype(tk[0]),spid(tk[1]);
    icode.resize(2);
    encode_momentum(momstr,icode[0]);
    encode_hadron(flav,irrep,sptype,spid,icode[1]);
    if ((flav[0]=='t')&&(flav[1]=='q')){   // tetraquark adjustment for colortype
       if (flav[7]=='m') icode[0]|=(0x1u << 31);}
    icode[0]|=0x1u;
    int LGIrrepRow;
    extract_from_string(irreprow,LGIrrepRow);
    if ((LGIrrepRow<1)||(LGIrrepRow>int(irrw_mask))){
       throw(std::invalid_argument("Unsupported value of LGIrrepRow"));}
    icode[0]|=LGIrrepRow<<(momt_bits+nhad_bits);}
      
 else if (nus==2){

    vector<string> majortags=split(opstr,'[');
    if (majortags.size()!=3) throw(std::invalid_argument("Invalid two hadron operator string"));
    vector<string> tokens=split(majortags[0],' ');
    if ((tokens.size()!=2)&&(tokens.size()!=3))
       throw(std::invalid_argument("Invalid two hadron operator string"));
    string totaliso(tokens[0]);
    vector<string> tk=split(tokens[1],'_');
    if (tk.size()!=2) throw(std::invalid_argument("Invalid two hadron operator string"));
    string totalirrep(tk[0]),irreprow(tk[1]);
    string lgcgid="0";
    if (tokens.size()==3){
       vector<string> lg=split(tokens[2],'_');
       if (lg[0]!="CG") throw(std::invalid_argument("Invalid two hadron operator string"));
       lgcgid=lg[1];}
    tokens=split(totaliso,'_');
    if (tokens.size()!=3) throw(std::invalid_argument("Invalid two hadron operator string"));
    totaliso=tokens[0];
    size_t pos=totaliso.find("iso");
    if (pos!=string::npos) totaliso.erase(pos,3);
    string flav1(tokens[1]),flav2(tokens[2]);
    if ((flav1[0]=='t')||(flav2[0]=='t')) throw(std::invalid_argument("Invalid two hadron operator string"));
    pos=majortags[1].find("]");
    if (pos!=string::npos) majortags[1].erase(pos,1);
    vector<string> hadron1=split(majortags[1],' ');
    pos=majortags[2].find("]");
    if (pos!=string::npos) majortags[2].erase(pos,1);
    vector<string> hadron2=split(majortags[2],' ');
    string mom1str(hadron1[0]),irrep1(hadron1[1]);
    tk=split(hadron1[2],'_');
    string spid1(tk[1]),sptype1(tk[0]);
    string mom2str(hadron2[0]),irrep2(hadron2[1]);
    tk=split(hadron2[2],'_');
    string spid2(tk[1]),sptype2(tk[0]);
    icode.resize(2*nus+1);
    encode_momentum(mom1str,icode[0]);
    encode_hadron(flav1,irrep1,sptype1,spid1,icode[1]);
    encode_momentum(mom2str,icode[2]);
    encode_hadron(flav2,irrep2,sptype2,spid2,icode[3]);
    icode[0]|=nus;
    encode_total(totaliso,"0",totalirrep,irreprow,lgcgid,icode[2*nus]);}

 else{
    throw(std::invalid_argument("Unsupported operator by string"));}}
 catch(const std::exception& errmsg){
    //cerr << "Invalid operator string: "<<errmsg<<endl;
    throw(std::invalid_argument(string("Invalid operator string:")
                   +string(errmsg.what())));}
}




void BasicLapHOperatorInfo::encode_momentum(const std::string& momstr, unsigned int& momcode)
{
 if (momstr.length()<9) throw(std::invalid_argument("Invalid momentum string"));
 if ((momstr[0]!='P')||(momstr[1]!='=')||(momstr[2]!='(')||(momstr[momstr.length()-1]!=')'))
    throw(std::invalid_argument("Invalid momentum string"));
 string mmm(momstr.substr(3,momstr.length()-4));
 vector<string> p=split(mmm,',');
 if (p.size()!=3){
    throw(std::invalid_argument("Invalid momentum string"));}
 int px,py,pz;
 try{
    extract_from_string(p[0],px);
    extract_from_string(p[1],py);
    extract_from_string(p[2],pz);}
 catch(const std::exception& msg){
    throw(std::invalid_argument("Invalid momentum string"));}

 if (((unsigned int)abs(px)>momj_mask)||((unsigned int)abs(py)>momj_mask)
    ||((unsigned int)abs(pz)>momj_mask)){
    //cerr << "momentum component magnitude not currently supported"<<endl;
    throw(std::invalid_argument("momentum component magnitude not currently supported"));}
 momcode=(px<0)?1:0; 
 momcode<<=momj_bits; momcode|=abs(px);
 momcode<<=1; momcode|=(py<0)?1:0;
 momcode<<=momj_bits; momcode|=abs(py);
 momcode<<=1; momcode|=(pz<0)?1:0;
 momcode<<=momj_bits; momcode|=abs(pz);
 momcode<<=nhad_bits;
}


void BasicLapHOperatorInfo::encode_hadron(const std::string& flav, const std::string& irrep,
                                 const std::string& sptype, const std::string& spid, 
                                 unsigned int& hadroncode)
{
 try{
    unsigned int fcode=m_flavor.encode(tidyString(flav));
    unsigned int irrepcode=m_irreps.encode(tidyString(irrep));
    unsigned int spcode=m_spatial.encode(tidyString(sptype));
    unsigned int spidnum;
    extract_from_string(spid,spidnum);
    unsigned int dlen;
    if ((sptype=="SS")||(sptype=="VI")||(is_glueball(fcode))) dlen=0;
    else if ((flav=="pion")||(flav=="eta")||(flav=="phi")||(flav=="kaon")||(flav=="kbar")) dlen=3;
    else dlen=2;

    //  now do the encoding
    hadroncode=dlen; 
    hadroncode<<=spid_bits;  hadroncode|=spidnum; 
    hadroncode<<=sptp_bits;  hadroncode|=spcode; 
    hadroncode<<=irrp_bits;  hadroncode|=irrepcode;
    hadroncode<<=flav_bits;  hadroncode|=fcode;

    // do some rudimentary checks (not exhaustive)
    if (is_boson(hadroncode)){
       if (!((irrep[0]=='A')||(irrep[0]=='B')||(irrep[0]=='E')||(irrep[0]=='T'))){
          throw(std::invalid_argument("Bad bosonic irrep in input xml data"));}}
    else if (is_fermion(hadroncode)){
       if (!((irrep[0]=='G')||(irrep[0]=='F')||(irrep[0]=='H'))){
          throw(std::invalid_argument("Bad fermionic irrep in input xml data"));}}}
 catch(const std::exception& errmsg){
    //cerr << "Invalid operator string: "<<errmsg<<endl;
    throw(std::invalid_argument(string("Invalid operator string: ")
           +string(errmsg.what())));}
}



void BasicLapHOperatorInfo::encode_total(const std::string& totalisospin, const std::string& isocgid,
                                const std::string& totalirrep, const std::string& irreprow, 
                                const std::string& lgcgid, unsigned int& code)
{
 try{
    unsigned int isocode=m_isospin.encode(tidyString(totalisospin));
    unsigned int irrep_code=m_irreps.encode(tidyString(totalirrep));
    unsigned int isoCGid;
    extract_from_string(isocgid,isoCGid);
    unsigned int lgCGid;
    extract_from_string(lgcgid,lgCGid);
    unsigned int irrepRow;
    extract_from_string(irreprow,irrepRow);
    if ((isoCGid<0)||(isoCGid>int(iscg_mask))){
       throw(std::invalid_argument("Requested value of isospin CG id not supported"));}
    if ((lgCGid<0)||(lgCGid>int(lgcg_mask))){
       throw(std::invalid_argument("Requested value of little group CG id not supported"));}
    if ((irrepRow<=0)||(irrepRow>int(irrw_mask))){
       throw(std::invalid_argument("Irrep row not currently supported"));}
    code=isocode;
    code<<=iscg_bits;  code|=isoCGid;
    code<<=irrp_bits;  code|=irrep_code;
    code<<=lgcg_bits;  code|=lgCGid;
    code<<=irrw_bits;  code|=irrepRow;}
 catch(const std::exception& errmsg){
    //cerr << "Invalid operator string: "<<errmsg<<endl;
    throw(std::invalid_argument(string("Invalid operator total: ")
          +string(errmsg.what())));}
}



string BasicLapHOperatorInfo::short_output() const
{
 unsigned int nhadrons=getNumberOfHadrons();
 if (nhadrons==0) return string("");
 else if (nhadrons==1){
    string flavor,hadstring;
    shortwrite_hadron(icode[0],icode[1],make_string(getLGIrrepRow()),flavor,hadstring);
    return flavor+" "+hadstring;}
 else if (nhadrons==2){
    string flav,had,totiso,totirrep;
    shortwrite_total(icode[2*nhadrons],totiso,totirrep);
    shortwrite_hadron(icode[0],icode[1],"",flav,had);
    totiso+="_"+flav;
    totirrep+=" ["+had+"]";
    shortwrite_hadron(icode[2],icode[3],"",flav,had);
    totiso+="_"+flav;
    totirrep+=" ["+had+"]";
    return totiso+" "+totirrep;}
 else{
    throw(std::invalid_argument("Unsupported operator for short_output"));}
}





void BasicLapHOperatorInfo::shortwrite_hadron(unsigned int momcode, unsigned int hadroncode,
                                              const string& irreprow, string& flavor, string& hadstring) const
{
 Momentum mom=get_momentum(momcode);
 hadstring=string("P=(");
 hadstring+=make_string(mom.x)+","+make_string(mom.y)+","+make_string(mom.z)+string(") ");
 unsigned int fcode=hadroncode;
 unsigned int ircode=fcode>>flav_bits;
 unsigned int spcode=ircode>>irrp_bits;
 unsigned int idcode=spcode>>sptp_bits;
 fcode&=flav_mask;
 bool notglueball=!is_glueball(fcode);
 if (fcode>12){ // tetraquark
    if ((momcode&(0x1u<<31))==0) fcode+=12;
    else fcode+=24;}
 ircode&=irrp_mask;
 spcode&=sptp_mask;
 idcode&=spid_mask;
 flavor=m_flavor.decode(fcode);
 hadstring+=m_irreps.decode(ircode);
 if (!irreprow.empty()) hadstring+="_"+irreprow;
 hadstring+=" "+m_spatial.decode(spcode);
 if (notglueball){
    hadstring+="_"+make_string(idcode);}
}



void BasicLapHOperatorInfo::shortwrite_total(unsigned int code, string& totaliso, string& totirrep) const
{
 unsigned int irreprow=code;
 unsigned int lgcgid=irreprow>>irrw_bits;
 unsigned int ircode=lgcgid>>lgcg_bits;
 unsigned int isocgid=ircode>>irrp_bits;
 unsigned int isocode=isocgid>>iscg_bits;
 irreprow&=irrw_mask;
 lgcgid&=lgcg_mask;
 ircode&=irrp_mask;
 isocgid&=iscg_mask;
 isocode&=isop_mask;
 totaliso=string("iso")+m_isospin.decode(isocode);
// if (isocgid>0) 
 totirrep=m_irreps.decode(ircode)+"_"+make_string(irreprow);
 if (lgcgid>0) totirrep+=" CG_"+make_string(lgcgid);
}


// ******************************************************************************

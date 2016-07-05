#include "mcobs_info.h"
#include "multi_compare.h"

using namespace std;


// ***************************************************************


MCObsInfo::MCObsInfo() : icode(1)
{
 icode[0]=0;     // default is zero particles
}


MCObsInfo::MCObsInfo(XMLHandler& xml_in)
{
 try{
    XMLHandler xmlr(xml_in);
    ComplexArg arg=RealPart;
    xmlr.set_exceptions_on();
    xml_tag_assert(xmlr,"MCObservable","MCObsInfo");
    XMLHandler xmlb(xmlr,"MCObservable");
    if (xmlb.count("VEV")==1){
       XMLHandler xmlv(xmlb,"VEV");
       OperatorInfo vop(xmlv);
       read_arg_type(xmlb,arg);
       encode(vop.icode,1,true,arg);}
    else if (xmlb.count("Correlator")==1){
       XMLHandler xmlc(xmlb,"Correlator");
       CorrelatorAtTimeInfo corr(xmlc);
       read_arg_type(xmlb,arg);
       bool subvev=corr.isVEVsubtracted();
       encode(corr.icode,2,!subvev,arg);}
    else if (xmlb.count("ObsName")==1){
       string name;
       xmlreadchild(xmlb,"ObsName",name);
       uint index=0;
       xmlreadifchild(xmlb,"Index",index);
       bool simple=(xmlb.count("Simple")==1)?true:false;
       string description;
       if (xmlb.count("Description")==1){
          XMLHandler xmld(xmlb,"Description");
          description=xmld.str();}
       read_arg_type(xmlb,arg);
       encode(name,index,description,simple,arg);}
    else{ throw(std::invalid_argument("Error"));}}
 catch(const std::exception& msg){
    throw(std::invalid_argument((string("Invalid XML for MCObsInfo constructor")
        +string(msg.what())).c_str()));}
}


MCObsInfo::MCObsInfo(const OperatorInfo& opinfo, ComplexArg arg) 
{
 encode(opinfo.icode,1,true,arg);
}


MCObsInfo::MCObsInfo(const OperatorInfo& sinkop, const OperatorInfo& sourceop, 
                     int timeval, bool hermitianmatrix, ComplexArg arg,
                     bool subvev)
{
 CorrelatorAtTimeInfo corr(sinkop,sourceop,timeval,hermitianmatrix,subvev);
 encode(corr.icode,2,!subvev,arg);
}


MCObsInfo::MCObsInfo(const CorrelatorAtTimeInfo& corrinfo, 
                     ComplexArg arg) 
{
 encode(corrinfo.icode,2,!(corrinfo.isVEVsubtracted()),arg);
}


MCObsInfo::MCObsInfo(const CorrelatorInfo& corrinfo, int timeval, 
                     bool hermitianmatrix, ComplexArg arg, bool subvev)
{
 CorrelatorAtTimeInfo corrt(corrinfo,timeval,hermitianmatrix,subvev);
 encode(corrt.icode,2,!subvev,arg);
}


MCObsInfo::MCObsInfo(const string& obsname, uint index, bool simple,
                     ComplexArg arg)
{
 encode(obsname,index,"",simple,arg);
}


void MCObsInfo::setToRealPart()
{
 set_real_part();
}


void MCObsInfo::setToImaginaryPart()
{
 if (isVacuum()) throw(std::invalid_argument("Cannot set vacuum to imaginary part"));
 set_imag_part();
}


void MCObsInfo::resetObsIndex(uint ind)
{
 if (isStandard())
    throw(std::invalid_argument("Cannot reset index for standard observable"));
 set_index(ind);
}


bool MCObsInfo::isVacuum() const
{
 return (icode[0]==0u);
}

bool MCObsInfo::isVEV() const
{
 return ((icode[0]>>2)==2u);
}

bool MCObsInfo::isCorrelatorAtTime() const
{
 return ((icode[0]>>2)==4u);
}

bool MCObsInfo::isHermitianCorrelatorAtTime() const
{
 if ((icode[0]>>2)!=4u) return false;
 return CorrelatorAtTimeInfo::isHermitian(icode);
}

bool MCObsInfo::isRealPart() const
{
 return ((icode[0]&1u)==0);
}

bool MCObsInfo::isImaginaryPart() const
{
 return ((icode[0]&1u)==1);
}



bool MCObsInfo::isSimple() const
{
 return (((icode[0]>>1)&1u)==0);
}

bool MCObsInfo::isNonSimple() const
{
 return (((icode[0]>>1)&1u)!=0);
}


bool MCObsInfo::isStandard() const
{
 return (((icode[0]>>2)&1u)==0);
}

bool MCObsInfo::isNonStandard() const
{
 return (((icode[0]>>2)&1u)!=0);
}


OperatorInfo MCObsInfo::getVEVInfo() const
{
 assert_vevtype("getVEVInfo");
 return OperatorInfo(icode.begin()+1,icode.end());
}


void MCObsInfo::getVEVInfo(OperatorInfo& vop) const
{
 assert_vevtype("getVEVInfo");
 vop.icode.resize(icode.size()-1);
 std::copy(icode.begin()+1,icode.end(),vop.icode.begin());
}


CorrelatorAtTimeInfo MCObsInfo::getCorrelatorAtTimeInfo() const
{
 assert_corrtype("getCorrelatorAtTimeInfo");
 return CorrelatorAtTimeInfo(icode.begin()+1,icode.end());
}


void MCObsInfo::getCorrelatorAtTimeInfo(CorrelatorAtTimeInfo& ctinfo) const
{
 assert_corrtype("getCorrelatorAtTimeInfo");
 ctinfo.icode.resize(icode.size()-1);
 std::copy(icode.begin()+1,icode.end(),ctinfo.icode.begin());
}




OperatorInfo MCObsInfo::getCorrelatorSourceInfo() const
{
 assert_corrtype("getCorrelatorSourceInfo");
 CorrelatorAtTimeInfo corr(icode.begin()+1,icode.end());
 return corr.getSource();
}


OperatorInfo MCObsInfo::getCorrelatorSinkInfo() const
{
 assert_corrtype("getCorrelatorSinkInfo");
 CorrelatorAtTimeInfo corr(icode.begin()+1,icode.end());
 return corr.getSink();
}


unsigned int MCObsInfo::getCorrelatorTimeIndex() const
{
 assert_corrtype("getCorrelatorTimeIndex");
 CorrelatorAtTimeInfo corr(icode.begin()+1,icode.end());
 return corr.getTimeSeparation();
}


CorrelatorInfo MCObsInfo::getCorrelatorInfo() const
{
 assert_corrtype("getCorrelatorInfo");
 CorrelatorAtTimeInfo corr(icode.begin()+1,icode.end());
 return corr.getCorrelator();
}


void MCObsInfo::getCorrelatorInfo(CorrelatorInfo& cinfo) const
{
 assert_corrtype("getCorrelatorInfo");
 CorrelatorAtTimeInfo corr(icode.begin()+1,icode.end());
 cinfo=corr.getCorrelator();
}


string MCObsInfo::getObsName() const
{
 if (isStandard()){
    throw(std::invalid_argument("getObsName called for standard observable"));}
 string res;
 try{
    res=m_decodings[get_obs_code()].first;}
 catch(const std::exception& msg){
    throw(std::invalid_argument((string("Error in getObsName: ")+string(msg.what())).c_str()));}
 return res;
}


uint MCObsInfo::getObsIndex() const
{
 if (isStandard()){
    throw(std::invalid_argument("getObsName called for standard observable"));}
 return get_obs_index();
}


string MCObsInfo::output(bool longform, int indent) const
{
 XMLHandler xmlout;
 output(xmlout,longform);
 return xmlout.output(indent);
}

string MCObsInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}

void MCObsInfo::output(XMLHandler& xmlout, bool longform) const
{
 xmlout.set_root("MCObservable");
 if (isVEV()){
    xmlout.put_child("VEV");
    XMLHandler xmlop;
    OperatorInfo vop(getVEVInfo());
    vop.output(xmlop,longform);
    xmlout.seek_first_child();
    xmlout.put_child(xmlop);
    if (isRealPart()) xmlout.put_sibling("Arg","RealPart");
    else xmlout.put_sibling("Arg","ImaginaryPart");}
 else if (isCorrelatorAtTime()){
    XMLHandler xmlop;
    CorrelatorAtTimeInfo corop(getCorrelatorAtTimeInfo());
    corop.output(xmlop,longform);
    xmlout.put_child(xmlop);
    if (isRealPart()) xmlout.put_child("Arg","RealPart");
    else xmlout.put_child("Arg","ImaginaryPart");}
 else if (isVacuum()){
    xmlout.put_child("Vacuum");}
 else{
    uint namecode=get_obs_code();
    uint index=get_obs_index();
    xmlout.put_child("ObsName",m_decodings[namecode].first);
    xmlout.put_child("Index",make_string(index));
    if (isSimple()) xmlout.put_child("Simple");
    if (isRealPart()) xmlout.put_child("Arg","RealPart");
    else xmlout.put_child("Arg","ImaginaryPart");
    if (m_decodings[namecode].second.length()>0){
       XMLHandler xmld;
       xmld.set_from_string(m_decodings[namecode].second);
       xmlout.put_child(xmld);}}
}


bool MCObsInfo::operator==(const MCObsInfo& rhs) const
{
 return multiEqual(icode,rhs.icode);   
}

bool MCObsInfo::operator!=(const MCObsInfo& rhs) const
{
 return multiNotEqual(icode,rhs.icode);   
}

bool MCObsInfo::operator<(const MCObsInfo& rhs) const
{
 return multiLessThan(icode,rhs.icode);
}

                     //  private routines

void MCObsInfo::encode(const vector<uint>& precode, unsigned int optype, 
                       bool simple, ComplexArg arg)
{
 icode.resize(precode.size()+1);
 std::copy(precode.begin(),precode.end(),icode.begin()+1);
 uint tcode=optype; 
 tcode<<=2; 
 if (!simple) tcode|=1u; 
 tcode<<=1;
 if (arg==ImaginaryPart) tcode|=1u;
 icode[0]=tcode;
}


void MCObsInfo::encode(const string& name, uint index, const string& description,
                       bool simple, ComplexArg arg)
{
 if (name.length()>32){
    throw(std::invalid_argument("MCObsInfo name too long"));}
 if (name.find_first_of("\t\n ")!=string::npos){
    throw(std::invalid_argument("MCObsInfo name cannot contain space, tabs, newlines"));}
 if (index>=8192){
    throw(std::invalid_argument("MCObsInfo index too large"));}
 uint namecode=0;
 map<string,uint>::iterator mt=m_encodings.find(name);
 if (mt!=m_encodings.end()) namecode=mt->second;
 else{
    namecode=m_encodings.size();
    if (namecode>=65536){
       throw(std::invalid_argument("MCObsInfo too many names"));}
    m_encodings.insert(make_pair(name,namecode));
    m_decodings.push_back(make_pair(name,description));}
 uint tcode=index;
 tcode<<=16;
 tcode|=namecode;
 tcode<<=1;
 tcode|=1u;  // nonstandard
 tcode<<=1; 
 if (!simple) tcode|=1u; 
 tcode<<=1;
 if (arg==ImaginaryPart) tcode|=1u;
 icode.resize(1);
 icode[0]=tcode;
}

 

void MCObsInfo::set_real_part()
{
 icode[0] &= ~1u;            // clear the bit
}


void MCObsInfo::set_imag_part()
{
 icode[0] |= 1u;          // set the bit
}


void MCObsInfo::set_index(uint index)
{
 if (index>=8192){
    throw(std::invalid_argument("MCObsInfo index too large"));}
 uint tcode=index;
 tcode<<=19;
 tcode|=(icode[0]&524287u);
 icode[0]=tcode;
}




void MCObsInfo::set_arg(ComplexArg arg)
{ 
 if (arg==RealPart) set_real_part();
 else set_imag_part();
}


bool MCObsInfo::read_arg_type(XMLHandler& xmlin, ComplexArg& arg)
{
 arg=RealPart;
 XMLHandler xmlt(xmlin);
 int na=xml_child_tag_count(xmlt,"Arg");
 if (na>1) throw(std::invalid_argument("Invalid Arg tag"));
 else if (na==0) return false;
 string reply;
 if (xmlreadifchild(xmlt,"Arg",reply)){
    if ((reply=="ImaginaryPart")||(reply=="Im")) arg=ImaginaryPart;
    else if ((reply=="RealPart")||(reply=="Re")) arg=RealPart;
    else throw(std::invalid_argument("Invalid Arg tag"));
    return true;}
 return false;
}


void MCObsInfo::assert_corrtype(const std::string& msg) const
{
 if (!isCorrelatorAtTime()){
    throw(std::invalid_argument((string("Assertion failed fpr Correlator type in MCObsInfo ")
          +msg).c_str()));}
}


void MCObsInfo::assert_vevtype(const std::string& msg) const
{
 if (!isVEV()){
    throw(std::invalid_argument((string("Assertion failed for VEV type in MCObsInfo ")
            +msg).c_str()));}
}


uint MCObsInfo::get_obs_code() const
{
 uint namecode=(icode[0]>>3)&65535u;
 if (namecode>=m_decodings.size()){
    throw(std::invalid_argument("getObsName failed: invalid name code"));}
 return namecode;
}


uint MCObsInfo::get_obs_index() const
{
 return (icode[0]>>19);
}



   //   static maps

std::map<std::string,uint> MCObsInfo::m_encodings;

std::vector<std::pair<std::string,std::string> > MCObsInfo::m_decodings;


// ******************************************************************************

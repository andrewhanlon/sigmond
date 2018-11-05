#include "mcobs_info.h"
#include "multi_compare.h"
#include "encoder.h"

using namespace std;


// ***************************************************************


MCObsInfo::MCObsInfo(bool reweighting_factor) : icode(1)
{
 if (reweighting_factor)
    icode[0]=8u;
 else
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
       bool reweight=(xmlv.count("Reweight")>0)?true:false;
       read_arg_type(xmlb,arg);
       encode(vop.icode,2+reweight,!reweight,arg);}
    else if (xmlb.count("Correlator")==1){
       XMLHandler xmlc(xmlb,"Correlator");
       CorrelatorAtTimeInfo corr(xmlc);
       read_arg_type(xmlb,arg);
       bool subvev=corr.subtractVEV();
       bool reweight=corr.reweight();
       encode(corr.icode,4,!(subvev||reweight),arg);}
    else if (xmlb.count("ObsName")==1){
       string name;
       xmlreadchild(xmlb,"ObsName",name);
       uint index=0;
       xmlreadifchild(xmlb,"Index",index);
       bool simple=(xmlb.count("Simple")>0)?true:false;
       read_arg_type(xmlb,arg);
       encode(name,index,simple,arg);}
    else if (xmlb.count("ReweightingFactor")==1){
       icode.resize(1); icode[0]=8u;}
    else if (xmlb.count("Vacuum")==1){
       icode.resize(1); icode[0]=0;}
    else{ throw(std::invalid_argument("Error"));}}
 catch(const std::exception& msg){
    throw(std::invalid_argument(string("Invalid XML for MCObsInfo constructor: ")
        +string(msg.what())));}
}


MCObsInfo::MCObsInfo(const OperatorInfo& opinfo, ComplexArg arg,
                     bool reweight) 
{
 encode(opinfo.icode,2+reweight,!reweight,arg);
}


MCObsInfo::MCObsInfo(const OperatorInfo& sinkop, const OperatorInfo& sourceop, 
                     int timeval, bool hermitianmatrix, ComplexArg arg,
                     bool subvev, bool reweight)
{
 CorrelatorAtTimeInfo corr(sinkop,sourceop,timeval,hermitianmatrix,subvev,reweight);
 encode(corr.icode,4,!(subvev||reweight),arg);
}


MCObsInfo::MCObsInfo(const CorrelatorAtTimeInfo& corrinfo, 
                     ComplexArg arg) 
{
 encode(corrinfo.icode,4,!(corrinfo.subtractVEV()||corrinfo.reweight()),arg);
}


MCObsInfo::MCObsInfo(const CorrelatorInfo& corrinfo, int timeval, 
                     bool hermitianmatrix, ComplexArg arg, bool subvev,
                     bool reweight)
{
 CorrelatorAtTimeInfo corrt(corrinfo,timeval,hermitianmatrix,subvev,reweight);
 encode(corrt.icode,4,!(subvev||reweight),arg);
}


MCObsInfo::MCObsInfo(const string& obsname, uint index, bool simple,
                     ComplexArg arg)
{
 encode(obsname,index,simple,arg);
}


void MCObsInfo::setToRealPart()
{
 set_real_part();
}


void MCObsInfo::setToImaginaryPart()
{
 if (isVacuum()) throw(std::invalid_argument("Cannot set vacuum to imaginary part"));
 if (isReweightingFactor()) throw(std::invalid_argument("Cannot set reweighting factor to imaginary part"));
 set_imag_part();
}


void MCObsInfo::resetObsIndex(uint ind)
{
 if (isPrimary())
    throw(std::invalid_argument("Cannot reset index for primary observable"));
 set_index(ind);
}


bool MCObsInfo::isVacuum() const
{
 return (icode[0]==0u);
}

bool MCObsInfo::isVEV() const
{
 unsigned int shifted = icode[0] >> 2;
 return ((shifted==4u)||(shifted==6u));
}

bool MCObsInfo::isReweightedVEV() const
{
 return ((icode[0]>>2)==6u);
}

bool MCObsInfo::isCorrelatorAtTime() const
{
 return ((icode[0]>>2)==8u);
}

bool MCObsInfo::isHermitianCorrelatorAtTime() const
{
 if ((icode[0]>>2)!=8u) return false;
 return CorrelatorAtTimeInfo::isHermitian(icode);
}

bool MCObsInfo::isReweightingFactor() const
{
 return (icode[0]==8u);
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


bool MCObsInfo::isPrimary() const
{
 return (((icode[0]>>2)&1u)==0);
}

bool MCObsInfo::isSecondary() const
{
 return (((icode[0]>>2)&1u)!=0);
}


bool MCObsInfo::isBasicLapH() const
{
 if (isSecondary()) 
    return false;
 if (isVEV())
    return getVEVInfo().isBasicLapH();
 if ((isCorrelatorAtTime())||(isHermitianCorrelatorAtTime()))
    return (getCorrelatorSourceInfo().isBasicLapH())
         &&(getCorrelatorSinkInfo().isBasicLapH());
 return false;
}

bool MCObsInfo::isGenIrrep() const
{
 if (isSecondary()) 
    return false;
 if (isVEV())
    return getVEVInfo().isGenIrrep();
 if ((isCorrelatorAtTime())||(isHermitianCorrelatorAtTime()))
    return (getCorrelatorSourceInfo().isGenIrrep())
         &&(getCorrelatorSinkInfo().isGenIrrep());
 return false;
}

bool MCObsInfo::isVEVsubtractedCorrelatorAtTime() const
{
 if (!isCorrelatorAtTime())
    return false;

 return getCorrelatorAtTimeInfo().subtractVEV();
}

bool MCObsInfo::isReweightedCorrelatorAtTime() const
{
 if (!isCorrelatorAtTime())
    return false;

 return getCorrelatorAtTimeInfo().reweight();
}

void MCObsInfo::setSimple()
{
 if (isVEV() || (isCorrelatorAtTime() && !isVEVsubtractedCorrelatorAtTime()))
    icode[0] &= ~2u;
}

void MCObsInfo::setNotSimple()
{
 if (isVEV() || (isCorrelatorAtTime() && !isVEVsubtractedCorrelatorAtTime()))
    icode[0] |= 2u;
}

void MCObsInfo::setReweight()
{
  if ((isVEV()) && (!isReweightedVEV())){
    icode[0] |= 8u;
  }
  else if ((isCorrelatorAtTime()) && (!isReweightedCorrelatorAtTime())){
    CorrelatorAtTimeInfo corr = getCorrelatorAtTimeInfo();
    corr.resetReweight(true);
    std::copy(corr.icode.begin(),corr.icode.end(),icode.begin()+1);
  }
}

void MCObsInfo::setNotReweight()
{
  if (isReweightedVEV()){
    icode[0] &= ~8u;
  }
  else if (isReweightedCorrelatorAtTime()){
    CorrelatorAtTimeInfo corr = getCorrelatorAtTimeInfo();
    corr.resetReweight(false);
    std::copy(corr.icode.begin(),corr.icode.end(),icode.begin()+1);
  }
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
 if (isPrimary()){
    throw(std::invalid_argument("getObsName called for primary observable"));}
 return get_obs_name();
}


uint MCObsInfo::getObsIndex() const
{
 if (isPrimary()){
    throw(std::invalid_argument("getObsName called for primary observable"));}
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
    XMLHandler xmlop;
    OperatorInfo vop(getVEVInfo());
    vop.output(xmlop,longform);
    if (longform){
       xmlout.put_child("VEV");
       xmlout.seek_first_child();
       xmlout.put_child(xmlop);
       if (isReweightedVEV()) xmlout.put_child("Reweight");
       if (isRealPart()) xmlout.put_sibling("Arg","RealPart");
       else xmlout.put_sibling("Arg","ImaginaryPart");}
    else{
       xmlop.rename_tag("VEV");
       xmlout.put_child(xmlop);
       string infostr=isReweightedVEV()?"RW ":"";
       infostr+=isRealPart()?" Re":" Im"; 
       xmlout.put_child("Info",infostr);}}
 else if (isCorrelatorAtTime()){
    XMLHandler xmlop;
    CorrelatorAtTimeInfo corop(getCorrelatorAtTimeInfo());
    corop.output(xmlop,longform);
    if (longform){
       xmlout.put_child(xmlop);
       if (isRealPart()) xmlout.put_child("Arg","RealPart");
       else xmlout.put_child("Arg","ImaginaryPart");}
    else{
       xmlop.rename_tag("MCObservable");
       xmlop.seek_unique("Src"); xmlop.rename_tag("CorrSrc");
       xmlop.seek_unique("Snk"); xmlop.rename_tag("CorrSnk");
       xmlop.seek_unique("Info");   
       string infostr=xmlop.get_text_content();
       infostr+= isRealPart() ? " Re" : " Im"; 
       xmlop.seek_unique("Info"); xmlop.seek_next_node();   
       xmlop.set_text_content(infostr);
       xmlout=xmlop;}}
 else if (isVacuum()){
    xmlout.put_child("Vacuum");}
 else if (isReweightingFactor()){
    xmlout.put_child("ReweightingFactor");}
 else{
    string obsname=get_obs_name();
    uint index=get_obs_index();
    if (longform){
       xmlout.put_child("ObsName",obsname);
       xmlout.put_child("Index",make_string(index));
       if (isSimple()) xmlout.put_child("Simple");
       if (isRealPart()) xmlout.put_child("Arg","RealPart");
       else xmlout.put_child("Arg","ImaginaryPart");}
    else{
       string infostr(obsname);
       infostr+=" Index="+make_string(index);
       if (isSimple()) infostr+=" Simple";
       infostr+=(isRealPart() ? " Re" : " Im");
       xmlout.put_child("Info",infostr);}}
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


void MCObsInfo::encode(const string& obsname, uint index,
                       bool simple, ComplexArg arg)
{
 uint nchar=obsname.length();
 if (nchar>32){
    throw(std::invalid_argument("MCObsInfo name cannot be longer than 32 characters"));}
 vector<uint> namecode;
 encode_string_to_uints(obsname,32,namecode);
 icode.resize(namecode.size()+1);
 std::copy(namecode.begin(),namecode.end(),icode.begin()+1);
 uint tcode=index;
 tcode<<=1;
 tcode|=1u;  // secondary
 tcode<<=1; 
 if (!simple) tcode|=1u; 
 tcode<<=1;
 if (arg==ImaginaryPart) tcode|=1u;
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
 uint tcode=index;
 tcode<<=3;
 tcode|=(icode[0]&7u);
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
    throw(std::invalid_argument(string("Assertion failed for Correlator type in MCObsInfo ")
          +msg));}
}


void MCObsInfo::assert_vevtype(const std::string& msg) const
{
 if (!isVEV()){
    throw(std::invalid_argument(string("Assertion failed for VEV type in MCObsInfo ")
            +msg));}
}


string MCObsInfo::get_obs_name() const
{
 vector<uint> namecode(icode.begin()+1,icode.end());
 return decode_uints_to_string(namecode);
}


uint MCObsInfo::get_obs_index() const
{
 return (icode[0]>>3);
}


void MCObsInfo::copyTo(unsigned int *buf) const
{
 buf[0]=icode.size();
 if (buf[0]>=max_ints)
    throw(std::runtime_error("icode in MCObsInfo too large for record key output"));
 for (uint i=0;i<buf[0];i++) buf[i+1]=icode[i];
 for (uint i=buf[0]+1;i<max_ints;i++) buf[i]=0;
}


MCObsInfo::MCObsInfo(const unsigned int *buf)
{
 if (buf[0]>=max_ints) 
    throw(std::runtime_error("icode in MCObsInfo too large for record key output"));
 icode.resize(buf[0]);
 std::copy(buf+1,buf+1+buf[0],icode.begin());
}

// ******************************************************************************

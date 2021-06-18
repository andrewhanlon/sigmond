#include "mcobs_info.h"
#include "multi_compare.h"
#include "encoder.h"

using namespace std;


// ***************************************************************


MCObsInfo::MCObsInfo() : icode(1)
{
 icode[0]=0;     // default is zero particles
}


MCObsInfo::MCObsInfo(XMLHandler& xml_in)
{
 try{
    set<string> tags;
    tags.insert("MCObservable");
    tags.insert("MCObsInfo");
    tags.insert("MCObs");
    ArgsHandler xin(xml_in,tags);
    string argstr="Re";
    xin.getOptionalString("Arg",argstr);
    ComplexArg arg=RealPart;
    if ((argstr=="Re")||(argstr=="RealPart")) arg=RealPart;
    else if ((argstr=="Im")||(argstr=="ImaginaryPart")) arg=ImaginaryPart;
    else throw(std::invalid_argument("Invalid Arg tag in MCObsInfo"));
    tags.clear();
    tags.insert("VEV");
    tags.insert("Correlator");
    tags.insert("CorrelatorInfo");
    tags.insert("CorrT");
    tags.insert("ObsName");
    tags.insert("Info");
    tags.insert("Vacuum");
    ArgsHandler xim(xin,tags);
    if (xim.getInputRootTag()=="VEV"){
       OperatorInfo vop(xim.getItem<OperatorInfo>("VEV"));
       encode(vop.icode,1,true,arg);}
    else if (xim.getInputRootTag()=="Correlator"){
       CorrelatorAtTimeInfo corr(xin.getItem<CorrelatorAtTimeInfo>("Correlator"));
       bool subvev=corr.subtractVEV();
       encode(corr.icode,2,!subvev,arg);}
    else if (xim.getInputRootTag()=="CorrelatorInfo"){
       CorrelatorAtTimeInfo corr(xin.getItem<CorrelatorAtTimeInfo>("CorrelatorInfo"));
       bool subvev=corr.subtractVEV();
       encode(corr.icode,2,!subvev,arg);}
    else if (xim.getInputRootTag()=="CorrT"){
       CorrelatorAtTimeInfo corr(xin.getItem<CorrelatorAtTimeInfo>("CorrT"));
       bool subvev=corr.subtractVEV();
       encode(corr.icode,2,!subvev,arg);}
    else if (xim.getInputRootTag()=="ObsName"){
       string name=xin.getString("ObsName");
       uint index=0;
       xin.getOptionalUInt("Index",index);
       bool simple= xin.queryTag("Simple") ? true : false;
       encode(name,index,simple,arg);}
    else if (xim.getInputRootTag()=="Info"){
       string info=xin.getString("Info");
       assign_from_string(info);}
    else{ // vacuum
       icode.resize(1); icode[0]=0;}}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("MCObsInfo construction failed: \n")
      +string(errmsg.what())+string("\nInput XML:")+xml_in.output()));}
}


MCObsInfo::MCObsInfo(const OperatorInfo& opinfo, ComplexArg arg) 
{
 encode(opinfo.icode,1,true,arg);
}


MCObsInfo::MCObsInfo(const OperatorInfo& sinkop, const OperatorInfo& sourceop, 
                     int timesep, bool hermitianmatrix, ComplexArg arg,
                     bool subvev)
{
 CorrelatorAtTimeInfo corr(sinkop,sourceop,timesep,hermitianmatrix,subvev);
 encode(corr.icode,2,!subvev,arg);
}

MCObsInfo::MCObsInfo(const OperatorInfo& sinkop, const OperatorInfo& insertop, const OperatorInfo& sourceop, 
                     int timesep, int timeins, bool hermitianmatrix, ComplexArg arg,
                     bool subvev)
{
 CorrelatorAtTimeInfo corr(sinkop, insertop, sourceop, timesep, timeins, hermitianmatrix, subvev);
 encode(corr.icode, 2, !subvev, arg);
}

MCObsInfo::MCObsInfo(const CorrelatorAtTimeInfo& corrinfo, 
                     ComplexArg arg) 
{
 encode(corrinfo.icode,2,!(corrinfo.subtractVEV()),arg);
}


MCObsInfo::MCObsInfo(const CorrelatorInfo& corrinfo, int timesep, 
                     bool hermitianmatrix, ComplexArg arg, bool subvev)
{
 CorrelatorAtTimeInfo corrt(corrinfo,timesep,hermitianmatrix,subvev);
 encode(corrt.icode,2,!subvev,arg);
}

MCObsInfo::MCObsInfo(const CorrelatorInfo& corrinfo, int timesep, int timeins,
                     bool hermitianmatrix, ComplexArg arg, bool subvev)
{
 CorrelatorAtTimeInfo corrt(corrinfo, timesep, timeins, hermitianmatrix, subvev);
 encode(corrt.icode, 2, !subvev, arg);
}


void MCObsInfo::assign_from_string(const string& opstring)
{
 string opstr(tidyString(opstring));
 vector<string> tokens=ArgsHandler::split(opstr,' ');
 if ((tokens.size()>4)&&(tokens.size()<1)) throw(std::runtime_error(""));
 string name(tokens[0]);
 uint index=0;
 if (tokens.size()>1){
    extract_from_string(tokens[1],index);}
 bool simple;
 if (tokens.size()<3) simple=false;
 else if (tokens[2]=="s") simple=true;
 else if (tokens[2]=="n") simple=false;
 else throw(std::runtime_error(""));
 ComplexArg arg;
 if (tokens.size()<4) arg=RealPart;
 else if (tokens[3]=="re") arg=RealPart;
 else if (tokens[3]=="im") arg=ImaginaryPart;
 else throw(std::runtime_error(""));
 encode(name,index,simple,arg);
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

bool MCObsInfo::hasOperatorInsertion() const
{
 if ((icode[0]>>2)!=4u) return false;
 return CorrelatorAtTimeInfo::hasOperatorInsertion(icode);
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
 if (hasOperatorInsertion())
    return (getCorrelatorSourceInfo().isBasicLapH())
         &&(getCorrelatorInsertionInfo().isBasicLapH())
         &&(getCorrelatorSinkInfo().isBasicLapH());
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

bool MCObsInfo::isImagDiagOfHermCorr() const
{
 return isImaginaryPart() && isHermitianCorrelatorAtTime() && getCorrelatorAtTimeInfo().isSinkSourceSame();
}

bool MCObsInfo::hasNoRelatedFlip() const
{
 return (!isHermitianCorrelatorAtTime()) || (getCorrelatorAtTimeInfo().isSinkSourceSame());
}

MCObsInfo MCObsInfo::getTimeFlipped() const
{
 return isHermitianCorrelatorAtTime() ?
      MCObsInfo(getCorrelatorAtTimeInfo().getTimeFlipped(),isRealPart() ? RealPart : ImaginaryPart)
     : *this;
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


OperatorInfo MCObsInfo::getCorrelatorInsertionInfo() const
{
 assert_corrtype("getCorrelatorInsertionInfo");
 CorrelatorAtTimeInfo corr(icode.begin()+1,icode.end());
 return corr.getInsertion();
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

unsigned int MCObsInfo::getCorrelatorInsertionTimeIndex() const
{
 assert_corrtype("getCorrelatorInsertionTimeIndex");
 CorrelatorAtTimeInfo corr(icode.begin()+1,icode.end());
 return corr.getTimeInsertion();
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
       if (isRealPart()) xmlout.put_sibling("Arg","RealPart");
       else xmlout.put_sibling("Arg","ImaginaryPart");}
    else{
       xmlout.put_child("VEV");
       xmlout.seek_first_child();
       xmlout.put_child(xmlop);
       xmlout.put_sibling("Arg",isRealPart()?"Re":"Im");}}
 else if (isCorrelatorAtTime()){
    XMLHandler xmlop;
    CorrelatorAtTimeInfo corop(getCorrelatorAtTimeInfo());
    corop.output(xmlop,longform);
    xmlout.put_child(xmlop);
    if (longform){
       if (isRealPart()) xmlout.put_child("Arg","RealPart");
       else xmlout.put_child("Arg","ImaginaryPart");}
    else{
       if (isRealPart()) xmlout.put_child("Arg","Re");
       else xmlout.put_child("Arg","Im");}}
 else if (isVacuum()){
    xmlout.put_child("Vacuum");}
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
       infostr+=" "+make_string(index);
       if (isSimple()) infostr+=" s"; else infostr+=" n";
       infostr+=(isRealPart() ? " re" : " im");
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
 if ((nchar<3)||(nchar>64)){
    throw(std::invalid_argument("MCObsInfo name cannot be longer than 64 characters or less than 3"));}
 if ((obsname=="RealPart")||(obsname=="ImaginaryPart")
    ||(obsname=="HermMat")||(obsname=="SubtractVEV")||(obsname=="HermitianMatrix")
    ||(obsname=="Arg")||(obsname=="Info")||(obsname=="Time")||(obsname=="TimeIndex")
    ||(obsname=="Corr")||(obsname=="CorrT")||(obsname=="Correlator")||(obsname=="CorrelatorInfo")
    ||(obsname=="MCObsInfo")||(obsname=="MCObservable")||(obsname=="MCObs")||(obsname=="Operator")
    ||(obsname=="OperatorInfo")||(obsname=="SubVEV"))
    throw(std::invalid_argument("MCObsInfo disallowed name"));
 vector<uint> namecode;
 encode_string_to_uints(obsname,64,namecode);
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

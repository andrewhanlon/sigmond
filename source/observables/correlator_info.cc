#include "correlator_info.h"
#include "multi_compare.h"
#include <algorithm>
#include "mcobs_info.h"
#include <stdexcept>

  //  See "operator_info.h" for a discussion of how each operator is
  //  encoded into a vector "icode" of unsigned integers.  

  //  In a CorrelatorInfo object, the icode vectors of the source and sink
  //  operators are appended, source first, sink last.  An extra integer
  //  is added at the end which contains the size of the source icode.  The 
  //  vacuum operator is NOT allowed.   

  //  For a CorrelatorAtTimeInfo object, the extra integer is added at 
  //  the end which contains (left to right) 24 bits for the time, 
  //  6 bits for the size of the source icode,
  //  the VEV subtraction bit and the Hermiticity bit.


using namespace std;


// ****************************************************************


CorrelatorInfo::CorrelatorInfo(XMLHandler& xml_in)
{
 try{
    set<string> tags;
    tags.insert("Correlator");
    tags.insert("CorrelatorInfo");
    tags.insert("Corr");
    ArgsHandler xin(xml_in,tags);
    if (xin.getInputRootTag()!="Corr"){
       ArgsHandler xin1(xin,"Source");
       OperatorInfo source(xin1.getItem<OperatorInfo>("Operator"));
       ArgsHandler xin2(xin,"Sink");
       OperatorInfo sink(xin2.getItem<OperatorInfo>("Operator"));
       assign(sink,source);}
    else{
        XMLHandler xmls(xml_in,"Corr");
       string corrstr=xmls.get_text_content();
       assign_from_string(corrstr);}}
 catch(const std::exception& msg){
    throw(std::invalid_argument(string("Invalid XML for CorrelatorInfo constructor")
          +string(msg.what())));}
}


CorrelatorInfo::CorrelatorInfo(const OperatorInfo& sink, const OperatorInfo& source) 
{
 assign(sink,source);
}


OperatorInfo CorrelatorInfo::getSource() const
{
 return OperatorInfo(icode.begin(),icode.begin()+icode.back());
}


OperatorInfo CorrelatorInfo::getSink() const
{
 return OperatorInfo(icode.begin()+icode.back(),icode.begin()+icode.size()-1);
}


CorrelatorInfo CorrelatorInfo::getTimeFlipped() const
{
 vector<unsigned int> tmpvec(icode.size());
 interchange_ends(tmpvec,icode);
 return CorrelatorInfo(tmpvec);
}


bool CorrelatorInfo::isSinkSourceSame() const
{
 return (getSink()==getSource());
}


string CorrelatorInfo::output(bool longform, int indent) const
{
 XMLHandler xmlout;
 output(xmlout,longform);
 return xmlout.output(indent);
}

string CorrelatorInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}

void CorrelatorInfo::output(XMLHandler& xmlout, bool longform) const
{
 if (longform){
    xmlout.set_root("Correlator");
    xmlout.put_child("Source");
    xmlout.seek_first_child();
    XMLHandler xmlop;
    getSource().output(xmlop,longform);
    xmlout.put_child(xmlop);
    xmlout.put_sibling("Sink");
    getSink().output(xmlop,longform);
    xmlout.put_child(xmlop);}
 else{
    XMLHandler xmlop;
    OperatorInfo opinfo(getSource());
    opinfo.output(xmlop,longform);
    string srcstr( opinfo.isBasicLapH() ? " BL{" : " GI{"); 
    srcstr+=xmlop.get_text_content()+"}";
    opinfo=getSink();
    opinfo.output(xmlop,longform);
    string snkstr( opinfo.isBasicLapH() ? " BL{" : " GI{"); 
    snkstr+=xmlop.get_text_content()+"}";
    xmlout.set_root("Corr",snkstr+srcstr);}
}


bool CorrelatorInfo::operator==(const CorrelatorInfo& rhs) const
{
 return multiEqual(icode,rhs.icode);   
}

bool CorrelatorInfo::operator!=(const CorrelatorInfo& rhs) const
{
 return multiNotEqual(icode,rhs.icode);   
}

bool CorrelatorInfo::operator<(const CorrelatorInfo& rhs) const
{
 return multiLessThan(icode,rhs.icode);   
}



               //  private routines


void CorrelatorInfo::assign(const OperatorInfo& sink, const OperatorInfo& source) 
{
 if ((source.icode[0]==0)||(sink.icode[0]==0)){
    throw(std::invalid_argument("Cannot have vacuum operator in CorrelatorInfo"));}
 uint sourcesize=source.icode.size();
 uint sinksize=sink.icode.size();
 icode.resize(sourcesize+sinksize+1);
 std::copy(source.icode.begin(),source.icode.end(),icode.begin());
 std::copy(sink.icode.begin(),sink.icode.end(),icode.begin()+sourcesize);
 icode.back()=sourcesize;
}

void CorrelatorInfo::assign_from_string(const std::string& corrstr)
{
 vector<string> tokens=ArgsHandler::split(tidyString(corrstr),'}');
 if (tokens.size()!=2)
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
 size_t pos1=tokens[0].find("BL{");
 size_t pos2=tokens[0].find("GI{");
 if ((pos1==string::npos)&&(pos2==string::npos))
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
 size_t pos = (pos1!=string::npos) ? pos1 : pos2;
 OperatorInfo::OpKind opkind = (pos1!=string::npos) ? OperatorInfo::BasicLapH : OperatorInfo::GenIrrep;
 OperatorInfo snkOp(tokens[0].substr(pos+3),opkind);
 pos1=tokens[1].find("BL{");
 pos2=tokens[1].find("GI{");
 if ((pos1==string::npos)&&(pos2==string::npos))
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
 pos = (pos1!=string::npos) ? pos1 : pos2;
 opkind = (pos1!=string::npos) ? OperatorInfo::BasicLapH : OperatorInfo::GenIrrep;
 OperatorInfo srcOp(tokens[1].substr(pos+3),opkind);
 assign(snkOp,srcOp);
}

void CorrelatorInfo::interchange_ends(std::vector<unsigned int>& outcode,
                                      const std::vector<unsigned int>& incode) const
{
 uint isrcsize=incode.back();
 uint isnksize=incode.size()-1-isrcsize;
 outcode.back()=isnksize;
 std::copy(incode.begin(),incode.begin()+isrcsize,outcode.begin()+isnksize);
 std::copy(incode.begin()+isrcsize,incode.end()-1,outcode.begin());
}






// ******************************************************************************






CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(XMLHandler& xml_in)
{
 try{
    set<string> tags;
    tags.insert("Correlator");
    tags.insert("CorrelatorInfo");
    tags.insert("CorrT");
    ArgsHandler xin(xml_in,tags);
    if (xin.getInputRootTag()!="CorrT"){
       ArgsHandler xin1(xin,"Source");
       OperatorInfo source(xin1.getItem<OperatorInfo>("Operator"));
       ArgsHandler xin2(xin,"Sink");
       OperatorInfo sink(xin2.getItem<OperatorInfo>("Operator"));
       unsigned int timeval;
       xin.getUInt("TimeIndex",timeval);
       bool hermitian=xin.getBool("HermitianMatrix");
       bool subvev=xin.getBool("SubtractVEV");
       assign(sink,source,timeval,hermitian,subvev);}
    else{
       XMLHandler xmls(xml_in,"CorrT");
       string corrstr=xmls.get_text_content();
       assign_from_string(corrstr);}}
 catch(const std::exception& msg){
    throw(std::invalid_argument(string("Invalid XML for CorrelatorAtTimeInfo constructor")
            +string(msg.what())));}
}



CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(const OperatorInfo& sink, 
                             const OperatorInfo& source,
                             int timeval, bool hermitianmatrix,
                             bool subvev) 
{
 assign(sink,source,timeval,hermitianmatrix,subvev);
}



CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(const CorrelatorInfo& corr, 
                        int timeval, bool hermitianmatrix, bool subvev)
{
 assign(corr,timeval,hermitianmatrix,subvev);
}



CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(const MCObsInfo& obsinfo)
               : icode(obsinfo.icode.begin()+1,obsinfo.icode.end())
{
 obsinfo.assert_corrtype("Constructing CorrelatorAtTimeInfo from MCObsInfo");
}


CorrelatorAtTimeInfo& CorrelatorAtTimeInfo::resetTimeSeparation(int timeval)
{
 if (timeval<0){
   throw(std::invalid_argument("Nonnegative time separation required in CorrelatorAtTimeInfo"));}
 uint tcode=timeval; tcode<<=8;
 tcode|= (icode.back()&255u);
 icode.back()=tcode;
 return *this;
}


CorrelatorAtTimeInfo& CorrelatorAtTimeInfo::resetSubtractVEV(bool subvev)
{
 if (subvev) icode.back()|=(1u<<1);
 else icode.back()&=~(1u<<1);
 return *this;
}



CorrelatorInfo CorrelatorAtTimeInfo::getCorrelator() const
{
 return CorrelatorInfo(icode,get_source_size());
}

OperatorInfo CorrelatorAtTimeInfo::getSource() const
{
 return OperatorInfo(icode.begin(),icode.begin()+get_source_size());
}

OperatorInfo CorrelatorAtTimeInfo::getSink() const
{
 return OperatorInfo(icode.begin()+get_source_size(),icode.begin()+icode.size()-1);
}

bool CorrelatorAtTimeInfo::isSinkSourceSame() const
{
 return (getSink()==getSource());
}

CorrelatorAtTimeInfo CorrelatorAtTimeInfo::getTimeFlipped() const
{
 vector<unsigned int> tmpcode(icode.size());
 uint isrcsize=(icode.back()>>2)&0x3Fu;
 uint isnksize=icode.size()-1-isrcsize;
 tmpcode.back()=(icode.back() & 0xFFFFFF03u) | (isnksize<<2);
 std::copy(icode.begin(),icode.begin()+isrcsize,tmpcode.begin()+isnksize);
 std::copy(icode.begin()+isrcsize,icode.end()-1,tmpcode.begin());
 return CorrelatorAtTimeInfo(tmpcode.begin(),tmpcode.end());
}


string CorrelatorAtTimeInfo::output(bool longform, int indent) const
{
 XMLHandler xmlout;
 output(xmlout,longform);
 return xmlout.output(indent);
}

string CorrelatorAtTimeInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}

void CorrelatorAtTimeInfo::output(XMLHandler& xmlout, bool longform) const
{
 if (longform){
    xmlout.set_root("Correlator");
    xmlout.put_child("Source");
    xmlout.seek_first_child();
    XMLHandler xmlop;
    getSource().output(xmlop,longform);
    xmlout.put_child(xmlop);
    xmlout.put_sibling("Sink");
    getSink().output(xmlop,longform);
    xmlout.put_child(xmlop);
    xmlout.put_sibling("TimeIndex",make_string(getTimeSeparation()));
    if (isHermitianMatrix()) 
       xmlout.put_sibling("HermitianMatrix");
    if (subtractVEV()) 
       xmlout.put_sibling("SubtractVEV");}
 else{
    XMLHandler xmlop;
    OperatorInfo opinfo(getSource());
    opinfo.output(xmlop,longform);
    string srcstr( opinfo.isBasicLapH() ? " BL{" : " GI{"); 
    srcstr+=xmlop.get_text_content()+"}";
    opinfo=getSink();
    opinfo.output(xmlop,longform);
    string snkstr( opinfo.isBasicLapH() ? " BL{" : " GI{"); 
    snkstr+=xmlop.get_text_content()+"}";
    string infostr(" time=");
    infostr+=make_string(getTimeSeparation());
    if (isHermitianMatrix()) 
       infostr+=" HermMat";
    if (subtractVEV()) 
       infostr+=" SubVEV";
    xmlout.set_root("CorrT",snkstr+srcstr+infostr);}
}


bool CorrelatorAtTimeInfo::operator==(const CorrelatorAtTimeInfo& rhs) const
{
 return multiEqual(icode,rhs.icode);   
}

bool CorrelatorAtTimeInfo::operator!=(const CorrelatorAtTimeInfo& rhs) const
{
 return multiNotEqual(icode,rhs.icode);   
}

bool CorrelatorAtTimeInfo::operator<(const CorrelatorAtTimeInfo& rhs) const
{
 return multiLessThan(icode,rhs.icode);   
}



               // private routines

void CorrelatorAtTimeInfo::assign(const OperatorInfo& sink, const OperatorInfo& source,
                                  int timeval, bool hermitianmatrix,
                                  bool subvev) 
{
 if ((source.icode[0]==0)||(sink.icode[0]==0)){
    throw(std::invalid_argument("Cannot have vacuum operator in CorrelatorAtTimeInfo"));}
 uint sourcesize=source.icode.size();
 uint sinksize=sink.icode.size();
 icode.resize(sourcesize+sinksize+1);
 std::copy(source.icode.begin(),source.icode.end(),icode.begin());
 std::copy(sink.icode.begin(),sink.icode.end(),icode.begin()+sourcesize);
 set_time_herm_vev(sourcesize,timeval,hermitianmatrix,subvev);
}


void CorrelatorAtTimeInfo::assign(const CorrelatorInfo& corr, int timeval,
                                  bool hermitianmatrix, bool subvev) 
{
 icode.resize(corr.icode.size());
 std::copy(corr.icode.begin(),corr.icode.end()-1,icode.begin());
 set_time_herm_vev(corr.icode.back(),timeval,hermitianmatrix,subvev);
}


void CorrelatorAtTimeInfo::assign_from_string(const std::string& corrstr)
{
 vector<string> tokens=ArgsHandler::split(tidyString(corrstr),'}');
 if (tokens.size()!=3)
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
 size_t pos1=tokens[0].find("BL{");
 size_t pos2=tokens[0].find("GI{");
 if ((pos1==string::npos)&&(pos2==string::npos))
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
 size_t pos = (pos1!=string::npos) ? pos1 : pos2;
 OperatorInfo::OpKind opkind = (pos1!=string::npos) ? OperatorInfo::BasicLapH : OperatorInfo::GenIrrep;
 OperatorInfo snkOp(tokens[0].substr(pos+3),opkind);
 pos1=tokens[1].find("BL{");
 pos2=tokens[1].find("GI{");
 if ((pos1==string::npos)&&(pos2==string::npos))
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
 pos = (pos1!=string::npos) ? pos1 : pos2;
 opkind = (pos1!=string::npos) ? OperatorInfo::BasicLapH : OperatorInfo::GenIrrep;
 OperatorInfo srcOp(tokens[1].substr(pos+3),opkind);
    // now parse the time, herm, subvev
 tokens=ArgsHandler::split(tidyString(tokens[2]),' ');
 uint timeval=0;
 vector<string> tt=ArgsHandler::split(tokens[0],'=');
 if (tt[0]!="time")
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
 extract_from_string(tt[1],timeval);
 bool herm=false;
 uint count=1;
 if ((tokens.size()>count)&&(tokens[count]=="HermMat")){ herm=true; ++count;}
 bool subvev=false;
 if ((tokens.size()>count)&&(tokens[count]=="SubVEV")){ subvev=true; ++count;}
 if (count!=tokens.size())
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
 assign(snkOp,srcOp,timeval,herm,subvev); 
}


void CorrelatorAtTimeInfo::set_time_herm_vev(uint srcsize, int timeval,
                                             bool hermitianmatrix, bool subvev)
{
 if (timeval<0){
    throw(std::invalid_argument("Nonnegative time separation required in CorrelatorAtTimeInfo"));}
 if (srcsize>=64){
    throw(std::invalid_argument("Too large source size in CorrelatorAtTimeInfo"));}
 uint tcode=timeval; tcode<<=6; 
 tcode|=srcsize; tcode<<=1;
 if (subvev) tcode|=1u; 
 tcode<<=1;
 if (hermitianmatrix) tcode|=1u;
 icode.back()=tcode;
}


// ******************************************************************************

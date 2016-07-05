#include "correlator_info.h"
#include "multi_compare.h"
#include <algorithm>
#include "mcobs_info.h"
#include <stdexcept>

  //  See "operator_info.h" for a discussion of how each operator is
  //  encoded into a vector "icode" of unsigned integers.  For a
  //  zero hadron operator (vacuum), icode has length 1; for a single
  //  hadron operator, icode has length 2; for n>=2 hadrons in the
  //  operator, icode has length 2*n+1.

  //  In a CorrelatorInfo object, the icode vectors of the source and sink
  //  operators are appended, source first, sink last.  The vacuum
  //  operator is NOT allowed.   For a CorrelatorAtTimeInfo object,
  //  another integer is added at the end which contains the time,
  //  the Hermiticity bit, and the VEV subtraction bit.


using namespace std;


// ****************************************************************


CorrelatorInfo::CorrelatorInfo(XMLHandler& xml_in)
{
 try{
    set<string> tags;
    tags.insert("Correlator");
    tags.insert("CorrelatorInfo");
    ArgsHandler xin(xml_in,tags);
    ArgsHandler xin1(xin,"Source");
    OperatorInfo source(xin1.getItem<OperatorInfo>("Operator"));
    ArgsHandler xin2(xin,"Sink");
    OperatorInfo sink(xin2.getItem<OperatorInfo>("Operator"));
    icode.resize(sink.icode.size()+source.icode.size());
    assign(icode,sink,source);}
 catch(const std::exception& msg){
    throw(std::invalid_argument((string("Invalid XML for CorrelatorInfo constructor")
          +string(msg.what())).c_str()));}
}


CorrelatorInfo::CorrelatorInfo(const OperatorInfo& sink, const OperatorInfo& source) 
               : icode(sink.icode.size()+source.icode.size())
{
 assign(icode,sink,source);
}


OperatorInfo CorrelatorInfo::getSource() const
{
 return OperatorInfo(sourcebegin(icode),sourceend(icode));
}


OperatorInfo CorrelatorInfo::getSink() const
{
 return OperatorInfo(sinkbegin(icode),sinkend(icode));
}


CorrelatorInfo CorrelatorInfo::getTimeFlipped() const
{
 vector<unsigned int> tmpvec(icode.size());
 interchange_ends(tmpvec,icode);
 return CorrelatorInfo(tmpvec.begin(),tmpvec.end());
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
 xmlout.set_root("Correlator");
 xmlout.put_child("Source");
 xmlout.seek_first_child();
 XMLHandler xmlop;
 getSource().output(xmlop,longform);
 xmlout.put_child(xmlop);
 xmlout.put_sibling("Sink");
 getSink().output(xmlop,longform);
 xmlout.put_child(xmlop);
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


void CorrelatorInfo::assign(std::vector<unsigned int>& outcode,
                            const OperatorInfo& sink, const OperatorInfo& source) 
{
 if ((source.icode[0]==0)||(sink.icode[0]==0)){
    throw(std::invalid_argument("Cannot have vacuum operator in CorrelatorInfo"));}
 std::copy(source.icode.begin(),source.icode.end(),sourcebegin(outcode));
 std::copy(sink.icode.begin(),sink.icode.end(),sinkbegin(outcode));
}


void CorrelatorInfo::interchange_ends(std::vector<unsigned int>& outcode,
                                      const std::vector<unsigned int>& incode)
{
 unsigned int n=incode.size();
 unsigned int sk=OperatorInfo::codesize(incode[0]);
 std::copy(incode.begin(),incode.begin()+sk,outcode.begin()+(n-sk));
 std::copy(incode.begin()+sk,incode.end(),outcode.begin());
}






// ******************************************************************************






CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(XMLHandler& xml_in)
{
 try{
    set<string> tags;
    tags.insert("Correlator");
    tags.insert("CorrelatorInfo");
    ArgsHandler xin(xml_in,tags);
    ArgsHandler xin1(xin,"Source");
    OperatorInfo source(xin1.getItem<OperatorInfo>("Operator"));
    ArgsHandler xin2(xin,"Sink");
    OperatorInfo sink(xin2.getItem<OperatorInfo>("Operator"));
    unsigned int timeval;
    xin.getUInt("TimeIndex",timeval);
    bool hermitian=xin.getBool("HermitianMatrix");
    bool subvev=xin.getBool("SubtractVEV");
    icode.resize(sink.icode.size()+source.icode.size()+1);
    assign(icode,sink,source,timeval,hermitian,subvev);}
 catch(const std::exception& msg){
    throw(std::invalid_argument((string("Invalid XML for CorrelatorAtTimeInfo constructor")
            +string(msg.what())).c_str()));}
}



CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(const OperatorInfo& sink, 
                             const OperatorInfo& source,
                             int timeval, bool hermitianmatrix,
                             bool subvev) 
               : icode(sink.icode.size()+source.icode.size()+1)
{
 assign(icode,sink,source,timeval,hermitianmatrix,subvev);
}



CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(const CorrelatorInfo& corr, 
                        int timeval, bool hermitianmatrix, bool subvev)
               : icode(corr.icode.size()+1)
{
 assign(icode,corr,timeval,hermitianmatrix,subvev);
}



CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(const MCObsInfo& obsinfo)
               : icode(obsinfo.icode.begin()+1,obsinfo.icode.end())
{
 obsinfo.assert_corrtype("Constructing CorrelatorAtTimeInfo from MCObsInfo");
}


CorrelatorAtTimeInfo& CorrelatorAtTimeInfo::resetTimeSeparation(int timeval)
{
 reset_time(icode,timeval);
 return *this;
}


CorrelatorAtTimeInfo& CorrelatorAtTimeInfo::resetVEVSubtracted(bool subvev)
{
 reset_vev(icode,subvev);
 return *this;
}

CorrelatorInfo CorrelatorAtTimeInfo::getCorrelator() const
{
 return CorrelatorInfo(sourcebegin(icode),sinkend(icode));
}

OperatorInfo CorrelatorAtTimeInfo::getSource() const
{
 return OperatorInfo(sourcebegin(icode),sourceend(icode));
}

OperatorInfo CorrelatorAtTimeInfo::getSink() const
{
 return OperatorInfo(sinkbegin(icode),sinkend(icode));
}

unsigned int CorrelatorAtTimeInfo::getTimeSeparation() const
{
 return extract_time(icode);
}

bool CorrelatorAtTimeInfo::isHermitianMatrix() const
{
 return isHermitian(icode);
}

bool CorrelatorAtTimeInfo::isVEVsubtracted() const
{
 return isVEVsubtracted(icode);
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
 if (isVEVsubtracted()) 
    xmlout.put_sibling("SubtractVEV");
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

void CorrelatorAtTimeInfo::assign(std::vector<unsigned int>& outcode,
                                  const OperatorInfo& sink, const OperatorInfo& source,
                                  int timeval, bool hermitianmatrix,
                                  bool subvev) 
{
 if ((source.icode[0]==0)||(sink.icode[0]==0)){
    throw(std::invalid_argument("Cannot have vacuum operator in CorrelatorAtTimeInfo"));}
 std::copy(source.icode.begin(),source.icode.end(),sourcebegin(outcode));
 std::copy(sink.icode.begin(),sink.icode.end(),sinkbegin(outcode));
 set_time_herm_vev(outcode,timeval,hermitianmatrix,subvev);
}


void CorrelatorAtTimeInfo::assign(std::vector<unsigned int>& outcode,
                                  const CorrelatorInfo& corr, int timeval,
                                  bool hermitianmatrix, bool subvev) 
{
 std::copy(corr.icode.begin(),corr.icode.end(),outcode.begin());
 set_time_herm_vev(outcode,timeval,hermitianmatrix,subvev);
}


// ******************************************************************************

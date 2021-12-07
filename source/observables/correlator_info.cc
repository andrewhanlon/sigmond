#include "correlator_info.h"
#include "multi_compare.h"
#include <algorithm>
#include "mcobs_info.h"
#include <stdexcept>

  //  See "operator_info.h" for a discussion of how each operator is
  //  encoded into a vector "icode" of unsigned integers.  

  //  OLD:
  //  In a CorrelatorInfo object, the icode vectors of the source and sink
  //  operators are appended, source first, sink last.  An extra integer
  //  is added at the end which contains the size of the source icode.  The 
  //  vacuum operator is NOT allowed.   

  //  For a CorrelatorAtTimeInfo object, the extra integer is added at 
  //  the end which contains (left to right) 24 bits for the time, 
  //  6 bits for the size of the source icode,
  //  the VEV subtraction bit and the Hermiticity bit.

  //  NEW:
  //  In a CorrelatorInfo object, the icode vectors of the source, insertion, and sink
  //  operators are appended, source first, insertion in the middle, and sink last.
  //  An extra integer is added at the end which contains the sizes of the insertion icode (16 leftmost bits)
  //  and source icode (16 rightmost bits).  The vacuum operator is NOT allowed.   

  //  For a CorrelatorAtTimeInfo object, the extra integer is added at 
  //  the end which contains (left to right) 10 bits for insertion time,
  //  10 bits for time separation, 5 bits for the size of the insertion icode, 
  //  5 bits for the size of the source icdoe,
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
    ArgsHandler xin(xml_in, tags);
    if (xin.getInputRootTag() != "Corr") {
      ArgsHandler xin1(xin, "Source");
      OperatorInfo source(xin1.getItem<OperatorInfo>("Operator"));
      ArgsHandler xin2(xin, "Sink");
      OperatorInfo sink(xin2.getItem<OperatorInfo>("Operator"));
      if (xin.queryTag("Insertion")) {
        ArgsHandler xin3(xin, "Insertion");
        OperatorInfo ins(xin3.getItem<OperatorInfo>("Operator"));
        assign(sink, ins, source);
      }
      else {
        assign(sink, source);
      }
    }
    else {
      XMLHandler xmls(xml_in, "Corr");
      string corrstr = xmls.get_text_content();
      assign_from_string(corrstr);
    }
  }
  catch(const std::exception& msg) {
    throw(std::invalid_argument(string("Invalid XML for CorrelatorInfo constructor")
          +string(msg.what())));
  }
}


CorrelatorInfo::CorrelatorInfo(const OperatorInfo& sink, const OperatorInfo& source) 
{
  assign(sink, source);
}

CorrelatorInfo::CorrelatorInfo(const OperatorInfo& sink, const OperatorInfo& insertion, const OperatorInfo& source) 
{
  assign(sink, insertion, source);
}

OperatorInfo CorrelatorInfo::getSource() const
{
  uint src_size = icode.back() & 0xFFFFu;
  return OperatorInfo(icode.begin(), icode.begin() + src_size);
}

OperatorInfo CorrelatorInfo::getInsertion() const
{
  if (!hasInsertion()) {
    throw(std::runtime_error(string("Correlator has no insertion time")));
  }
  uint src_size = icode.back() & 0xFFFFu;
  uint ins_size = icode.back() >> 16;
  return OperatorInfo(icode.begin() + src_size, icode.begin() + src_size + ins_size);
}

OperatorInfo CorrelatorInfo::getSink() const
{
  uint src_size = icode.back() & 0xFFFFu;
  uint ins_size = icode.back() >> 16;
  return OperatorInfo(icode.begin() + src_size + ins_size, icode.begin() + icode.size() - 1);
}

bool CorrelatorInfo::hasInsertion() const
{
  return (icode.back() >> 16);
}

CorrelatorInfo CorrelatorInfo::getTimeFlipped() const
{
  vector<unsigned int> tmpvec(icode.size());
  interchange_ends(tmpvec,icode);
  return CorrelatorInfo(tmpvec);
}


bool CorrelatorInfo::isSinkSourceSame() const
{
  return (getSink() == getSource());
}

bool CorrelatorInfo::isBackwards() const
{
 return getSink().isBackwards();
}

void CorrelatorInfo::setBackwards()
{
  OperatorInfo op_snk = getSink();
  OperatorInfo op_src = getSource();
  op_snk.setBackwards();
  op_src.setBackwards();
  if (hasInsertion()) {
    OperatorInfo op_ins = getInsertion();
    op_ins.setBackwards();
    CorrelatorInfo temp(op_snk, op_ins, op_src);
    icode = temp.icode;
  }
  else {
    CorrelatorInfo temp(op_snk, op_src);
    icode = temp.icode;
  }
}

void CorrelatorInfo::setForwards()
{
  OperatorInfo op_snk = getSink();
  OperatorInfo op_src = getSource();
  op_snk.setForwards();
  op_src.setForwards();
  if (hasInsertion()) {
    OperatorInfo op_ins = getInsertion();
    op_ins.setForwards();
    CorrelatorInfo temp(op_snk, op_ins, op_src);
    icode = temp.icode;
  }
  else {
    CorrelatorInfo temp(op_snk, op_src);
    icode = temp.icode;
  }
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
  if (longform) {
    xmlout.set_root("Correlator");
    xmlout.put_child("Source");
    xmlout.seek_first_child();
    XMLHandler xmlop;
    getSource().output(xmlop, longform);
    xmlout.put_child(xmlop);
    if (hasInsertion()) {
      xmlout.put_sibling("Insertion");
      getInsertion().output(xmlop, longform);
      xmlout.put_child(xmlop);
    }
    xmlout.put_sibling("Sink");
    getSink().output(xmlop, longform);
    xmlout.put_child(xmlop);
  }
  else {
    XMLHandler xmlop;
    OperatorInfo opinfo(getSource());
    opinfo.output(xmlop, longform);
    string srcstr(opinfo.isBasicLapH() ? " BL{" : " GI{"); 
    srcstr += xmlop.get_text_content() + "}";
    string insstr = "";
    if (hasInsertion()) {
      opinfo = getInsertion();
      opinfo.output(xmlop, longform);
      insstr = opinfo.isBasicLapH() ? " BL{" : " GI{";
      insstr += xmlop.get_text_content() + "}";
    }
    opinfo = getSink();
    opinfo.output(xmlop, longform);
    string snkstr(opinfo.isBasicLapH() ? " BL{" : " GI{"); 
    snkstr += xmlop.get_text_content() + "}";
    xmlout.set_root("Corr", snkstr + insstr + srcstr);}
}


bool CorrelatorInfo::operator==(const CorrelatorInfo& rhs) const
{
 return multiEqual(icode, rhs.icode);   
}

bool CorrelatorInfo::operator!=(const CorrelatorInfo& rhs) const
{
 return multiNotEqual(icode, rhs.icode);   
}

bool CorrelatorInfo::operator<(const CorrelatorInfo& rhs) const
{
 return multiLessThan(icode, rhs.icode);   
}



               //  private routines


void CorrelatorInfo::assign(const OperatorInfo& sink, const OperatorInfo& source) 
{
  if ((source.icode[0]==0) || (sink.icode[0]==0)) {
    throw(std::invalid_argument("Cannot have vacuum operator in CorrelatorInfo"));
  }
  /*
  if (source.isBackwards() != sink.isBackwards()) {
     throw(std::invalid_argument("Cannot have sink and source have different 'isBackwards'"));
  }
  */
  uint sourcesize = source.icode.size();
  uint sinksize = sink.icode.size();
  icode.resize(sourcesize + sinksize + 1);
  std::copy(source.icode.begin(), source.icode.end(), icode.begin());
  std::copy(sink.icode.begin(), sink.icode.end(), icode.begin() + sourcesize);
  icode.back() = sourcesize;
}

void CorrelatorInfo::assign(const OperatorInfo& sink, const OperatorInfo& insert, const OperatorInfo& source) 
{
  if ((source.icode[0] == 0) || (insert.icode[0] == 0) || (sink.icode[0] == 0)) {
    throw(std::invalid_argument("Cannot have vacuum operator in CorrelatorInfo"));
  }
  /*
  if (source.isBackwards() != sink.isBackwards()) {
    throw(std::invalid_argument("Cannot have sink and source have different 'isBackwards'"));
  }
  */
  uint sourcesize = source.icode.size();
  uint insertsize = insert.icode.size();
  uint sinksize = sink.icode.size();
  icode.resize(sourcesize + insertsize + sinksize + 1);
  std::copy(source.icode.begin(), source.icode.end(), icode.begin());
  std::copy(insert.icode.begin(), insert.icode.end(), icode.begin() + sourcesize);
  std::copy(sink.icode.begin(), sink.icode.end(), icode.begin() + sourcesize + insertsize);
  icode.back() = insertsize; icode.back() <<= 16;
  icode.back() |= sourcesize;
}

void CorrelatorInfo::assign_from_string(const std::string& corrstr)
{
  vector<string> tokens=ArgsHandler::split(tidyString(corrstr),'}');
  if ((tokens.size() != 2) || (tokens.size() != 3)) {
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
  }
  size_t pos_brace = tokens[0].find("{");
  size_t pos_gi = tokens[0].find("GI{");
  if (pos_brace == string::npos) {
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
  }
  OperatorInfo::OpKind opkind = (pos_gi != string::npos) ? OperatorInfo::GenIrrep : OperatorInfo::BasicLapH;
  OperatorInfo snkOp(tokens[0].substr(pos_brace+1), opkind);
  pos_brace = tokens[1].find("{");
  pos_gi = tokens[1].find("GI{");
  if (pos_brace == string::npos) {
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
  }
  opkind = (pos_gi!=string::npos) ? OperatorInfo::GenIrrep : OperatorInfo::BasicLapH;
  if (tokens.size() == 2) {
    OperatorInfo srcOp(tokens[1].substr(pos_brace+1), opkind);
    assign(snkOp, srcOp);
  }
  else {
    OperatorInfo insOp(tokens[1].substr(pos_brace+1), opkind);
    pos_brace = tokens[2].find("{");
    pos_gi = tokens[2].find("GI{");
    if (pos_brace == string::npos) {
      throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
    }
    opkind = (pos_gi!=string::npos) ? OperatorInfo::GenIrrep : OperatorInfo::BasicLapH;
    OperatorInfo srcOp(tokens[2].substr(pos_brace+1), opkind);
    assign(snkOp, insOp, srcOp);
  }
}

void CorrelatorInfo::interchange_ends(std::vector<unsigned int>& outcode,
                                      const std::vector<unsigned int>& incode) const
{
 uint isrcsize = incode.back() & 0xFFFFu;
 uint iinssize = incode.back() >> 16;
 uint isnksize = incode.size() - iinssize - isrcsize - 1;
 outcode.back() = iinssize; outcode.back() <<= 16;
 outcode.back() |= isnksize;
 std::copy(incode.begin(), incode.begin()+isrcsize, outcode.begin()+isnksize+iinssize);
 std::copy(incode.begin()+isrcsize, incode.begin()+isrcsize+iinssize, outcode.begin()+isnksize);
 std::copy(incode.begin()+isrcsize+iinssize, incode.end()-1, outcode.begin());
}






// ******************************************************************************






CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(XMLHandler& xml_in)
{
  try {
    set<string> tags;
    tags.insert("Correlator");
    tags.insert("CorrelatorInfo");
    tags.insert("CorrT");
    ArgsHandler xin(xml_in, tags);
    if (xin.getInputRootTag() != "CorrT") {
      unsigned int timesep;
      xin.getUInt("TimeIndex", timesep);
      bool hermitian = xin.getBool("HermitianMatrix");
      bool subvev = xin.getBool("SubtractVEV");
      ArgsHandler xin1(xin, "Source");
      OperatorInfo source(xin1.getItem<OperatorInfo>("Operator"));
      ArgsHandler xin2(xin, "Sink");
      OperatorInfo sink(xin2.getItem<OperatorInfo>("Operator"));
      if (xin.queryTag("Insertion")) {
        ArgsHandler xin3(xin, "Insertion");
        OperatorInfo ins(xin3.getItem<OperatorInfo>("Operator"));
        unsigned int timeins;
        xin.getUInt("TimeInsertion", timeins);
        assign(sink, ins, source, timesep, timeins, hermitian, subvev);
      }
      else {
        assign(sink, source, timesep, hermitian, subvev);
      }
    }
    else {
      XMLHandler xmls(xml_in,"CorrT");
      string corrstr=xmls.get_text_content();
      assign_from_string(corrstr);
    }
  }
  catch(const std::exception& msg) {
    throw(std::invalid_argument(string("Invalid XML for CorrelatorAtTimeInfo constructor")
            +string(msg.what())));
  }
}



CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(const OperatorInfo& sink, 
                             const OperatorInfo& source,
                             int timesep, bool hermitianmatrix,
                             bool subvev) 
{
 assign(sink, source, timesep, hermitianmatrix, subvev);
}

CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(const OperatorInfo& sink, 
                             const OperatorInfo& insertion,
                             const OperatorInfo& source,
                             int timesep, int timeins,
                             bool hermitianmatrix, bool subvev) 
{
 assign(sink, insertion, source, timesep, timeins, hermitianmatrix, subvev);
}


CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(const CorrelatorInfo& corr, 
                        int timesep, bool hermitianmatrix, bool subvev)
{
 assign(corr, timesep, hermitianmatrix, subvev);
}

CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(const CorrelatorInfo& corr, 
                        int timesep, int timeins, bool hermitianmatrix, bool subvev)
{
 assign(corr, timesep, timeins, hermitianmatrix, subvev);
}


CorrelatorAtTimeInfo::CorrelatorAtTimeInfo(const MCObsInfo& obsinfo)
               : icode(obsinfo.icode.begin() + 1, obsinfo.icode.end())
{
 obsinfo.assert_corrtype("Constructing CorrelatorAtTimeInfo from MCObsInfo");
}


CorrelatorAtTimeInfo& CorrelatorAtTimeInfo::resetTimeSeparation(int timesep)
{
  if (timesep < 0) {
    throw(std::invalid_argument("Nonnegative time separation required in CorrelatorAtTimeInfo"));
  }
  uint tcode = timesep; tcode <<= 12;
  tcode |= (icode.back() & 0xFFC00FFFu);
  icode.back() = tcode;
  return *this;
}

CorrelatorAtTimeInfo& CorrelatorAtTimeInfo::resetTimeInsertion(int timeins)
{
  if (timeins < 0) {
    throw(std::invalid_argument("Nonnegative time insertion required in CorrelatorAtTimeInfo"));
  }
  uint tcode = timeins; tcode <<= 22;
  tcode |= (icode.back() & 0x3FFFFFu);
  icode.back() = tcode;
  return *this;
}


CorrelatorAtTimeInfo& CorrelatorAtTimeInfo::resetSubtractVEV(bool subvev)
{
  if (subvev) icode.back() |= (1u << 1);
  else icode.back() &= ~(1u << 1);
  return *this;
}



CorrelatorInfo CorrelatorAtTimeInfo::getCorrelator() const
{
  return CorrelatorInfo(icode, get_source_size(), get_insertion_size());
}

OperatorInfo CorrelatorAtTimeInfo::getSource() const
{
  return OperatorInfo(icode.begin(), icode.begin() + get_source_size());
}

OperatorInfo CorrelatorAtTimeInfo::getInsertion() const
{
  if (!hasInsertion()) {
    throw(std::runtime_error(string("Correlator has no insertion time")));
  }
  return OperatorInfo(icode.begin() + get_source_size(), icode.begin() + get_source_size() + get_insertion_size());
}

OperatorInfo CorrelatorAtTimeInfo::getSink() const
{
  return OperatorInfo(icode.begin() + get_source_size() + get_insertion_size(), icode.begin() + icode.size() - 1);
}

bool CorrelatorAtTimeInfo::isSinkSourceSame() const
{
  return (getSink() == getSource());
}

bool CorrelatorAtTimeInfo::isBackwards() const
{
  return getSink().isBackwards();
}

void CorrelatorAtTimeInfo::setBackwards()
{
  CorrelatorInfo corr = getCorrelator();
  corr.setBackwards();
  if (hasInsertion()) {
    CorrelatorAtTimeInfo temp(corr, getTimeSeparation(), getTimeInsertion(), isHermitianMatrix(), subtractVEV());
    icode = temp.icode;
  }
  else {
    CorrelatorAtTimeInfo temp(corr, getTimeSeparation(), isHermitianMatrix(), subtractVEV());
    icode = temp.icode;
  }
}

void CorrelatorAtTimeInfo::setForwards()
{
  CorrelatorInfo corr = getCorrelator();
  corr.setForwards();
  if (hasInsertion()) {
    CorrelatorAtTimeInfo temp(corr, getTimeSeparation(), getTimeInsertion(), isHermitianMatrix(), subtractVEV());
    icode = temp.icode;
  }
  else {
    CorrelatorAtTimeInfo temp(corr, getTimeSeparation(), isHermitianMatrix(), subtractVEV());
    icode = temp.icode;
  }
}

CorrelatorAtTimeInfo CorrelatorAtTimeInfo::getTimeFlipped() const
{
 vector<unsigned int> tmpcode(icode.size());
 uint isrcsize = (icode.back() >> 2) & 0x1Fu;
 uint iinssize = (icode.back() >> 7) & 0x1Fu;
 uint isnksize = icode.size() - isrcsize - iinssize - 1;
 tmpcode.back() = (icode.back() & 0xFFFFFF83u) | (isnksize << 2);
 std::copy(icode.begin(), icode.begin()+isrcsize, tmpcode.begin()+isnksize+iinssize);
 std::copy(icode.begin()+isrcsize, icode.begin()+isrcsize+iinssize, tmpcode.begin()+isnksize);
 std::copy(icode.begin()+isrcsize+iinssize, icode.end()-1, tmpcode.begin());
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
  if (longform) {
    xmlout.set_root("Correlator");
    xmlout.put_child("Source");
    xmlout.seek_first_child();
    XMLHandler xmlop;
    getSource().output(xmlop, longform);
    xmlout.put_child(xmlop);
    if (hasInsertion()) {
      xmlout.put_sibling("Insertion");
      getInsertion().output(xmlop, longform);
      xmlout.put_child(xmlop);
    }
    xmlout.put_sibling("Sink");
    getSink().output(xmlop, longform);
    xmlout.put_child(xmlop);
    xmlout.put_sibling("TimeIndex", make_string(getTimeSeparation()));
    if (hasInsertion()) {
      xmlout.put_sibling("TimeInsIndex", make_string(getTimeInsertion()));
    }
    if (isHermitianMatrix()) {
      xmlout.put_sibling("HermitianMatrix");
    }
    if (subtractVEV()) {
      xmlout.put_sibling("SubtractVEV");
    }
  }
  else {
    XMLHandler xmlop;
    OperatorInfo opinfo(getSource());
    opinfo.output(xmlop,longform);
    string srcstr(opinfo.isBasicLapH() ? " BL{" : " GI{"); 
    srcstr += xmlop.get_text_content()+"}";
    string insstr = "";
    if (hasInsertion()) {
      opinfo = getInsertion();
      opinfo.output(xmlop, longform);
      insstr = opinfo.isBasicLapH() ? " BL{" : " GI{";
      insstr += xmlop.get_text_content() + "}";
    }
    opinfo = getSink();
    opinfo.output(xmlop,longform);
    string snkstr(opinfo.isBasicLapH() ? " BL{" : " GI{"); 
    snkstr += xmlop.get_text_content()+"}";
    string infostr(" time=");
    infostr += make_string(getTimeSeparation());
    if (hasInsertion()) {
      infostr += " timeins=";
      infostr += make_string(getTimeInsertion());
    }
    if (isHermitianMatrix()) {
      infostr+=" HermMat";
    }
    if (subtractVEV()) {
      infostr+=" SubVEV";
    }
    xmlout.set_root("CorrT", snkstr + insstr + srcstr + infostr);}
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
                                  int timesep, bool hermitianmatrix,
                                  bool subvev) 
{
  if ((source.icode[0] == 0) || (sink.icode[0] == 0)) {
    throw(std::invalid_argument("Cannot have vacuum operator in CorrelatorAtTimeInfo"));
  }
  /*
  if (source.isBackwards() != sink.isBackwards()) {
     throw(std::invalid_argument("Cannot have sink and source have different 'isBackwards'"));
  }
  */
  uint sourcesize = source.icode.size();
  uint sinksize = sink.icode.size();
  icode.resize(sourcesize + sinksize + 1);
  std::copy(source.icode.begin(),source.icode.end(),icode.begin());
  std::copy(sink.icode.begin(),sink.icode.end(),icode.begin()+sourcesize);
  set_time_herm_vev(sourcesize, timesep, hermitianmatrix, subvev);
}

void CorrelatorAtTimeInfo::assign(const OperatorInfo& sink, const OperatorInfo& insertion,
                                  const OperatorInfo& source,
                                  int timesep, int timeins,
                                  bool hermitianmatrix, bool subvev) 
{
  if ((source.icode[0] == 0) || (insertion.icode[0] == 0) || (sink.icode[0] == 0)) {
    throw(std::invalid_argument("Cannot have vacuum operator in CorrelatorAtTimeInfo"));
  }
  /*
  if (source.isBackwards() != sink.isBackwards()) {
     throw(std::invalid_argument("Cannot have sink and source have different 'isBackwards'"));
  }
  */
  uint sourcesize = source.icode.size();
  uint insertsize = insertion.icode.size();
  uint sinksize = sink.icode.size();
  icode.resize(sourcesize + insertsize + sinksize + 1);
  std::copy(source.icode.begin(), source.icode.end(), icode.begin());
  std::copy(insertion.icode.begin(), insertion.icode.end(), icode.begin() + sourcesize);
  std::copy(sink.icode.begin(), sink.icode.end(), icode.begin() + sourcesize + insertsize);
  set_time_herm_vev(sourcesize, insertsize, timesep, timeins, hermitianmatrix, subvev);
}


void CorrelatorAtTimeInfo::assign(const CorrelatorInfo& corr, int timesep,
                                  bool hermitianmatrix, bool subvev) 
{
  if (corr.hasInsertion()) {
    throw(std::runtime_error(string("Correlator needs insertion time")));
  }
  icode.resize(corr.icode.size());
  std::copy(corr.icode.begin(), corr.icode.end() - 1, icode.begin());
  set_time_herm_vev(corr.icode.back(), timesep, hermitianmatrix, subvev);
}

void CorrelatorAtTimeInfo::assign(const CorrelatorInfo& corr, int timesep, int timeins,
                                  bool hermitianmatrix, bool subvev) 
{
  if (!corr.hasInsertion()) {
    throw(std::runtime_error(string("Correlator has no insertion time")));
  }
  icode.resize(corr.icode.size());
  std::copy(corr.icode.begin(), corr.icode.end() - 1, icode.begin());
  uint srcsize = corr.icode.back() & 0xFFFFu;
  uint inssize = (corr.icode.back() >> 16) & 0xFFFFu;
  set_time_herm_vev(srcsize, inssize, timesep, timeins, hermitianmatrix, subvev);
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
 uint timesep=0;
 vector<string> tt=ArgsHandler::split(tokens[0],'=');
 if (tt[0]!="time")
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
 extract_from_string(tt[1],timesep);
 bool herm=false;
 uint count=1;
 if ((tokens.size()>count)&&(tokens[count]=="HermMat")){ herm=true; ++count;}
 bool subvev=false;
 if ((tokens.size()>count)&&(tokens[count]=="SubVEV")){ subvev=true; ++count;}
 if (count!=tokens.size())
    throw(std::invalid_argument("Invalid string to construct CorrelatorInfo"));
 assign(snkOp,srcOp,timesep,herm,subvev); 
}


void CorrelatorAtTimeInfo::set_time_herm_vev(uint srcsize, int timesep,
                                             bool hermitianmatrix, bool subvev)
{
  if (timesep < 0) {
    throw(std::invalid_argument("Nonnegative time separation required in CorrelatorAtTimeInfo"));
  }
  if (srcsize >= 32) {
    throw(std::invalid_argument("Too large source size in CorrelatorAtTimeInfo"));
  }
  uint tcode = timesep; tcode <<= 10; 
  tcode |= srcsize; tcode <<= 1;
  if (subvev) tcode |= 1u; 
  tcode <<= 1;
  if (hermitianmatrix) tcode |= 1u;
  icode.back() = tcode;
}

void CorrelatorAtTimeInfo::set_time_herm_vev(uint srcsize, uint inssize,
                                             int timesep, int timeins,
                                             bool hermitianmatrix, bool subvev)
{
  if (timesep < 0) {
    throw(std::invalid_argument("Nonnegative time separation required in CorrelatorAtTimeInfo"));
  }
  if (timeins < 0) {
    throw(std::invalid_argument("Nonnegative time insertion required in CorrelatorAtTimeInfo"));
  }
  if (srcsize >= 32) {
    throw(std::invalid_argument("Too large source size in CorrelatorAtTimeInfo"));
  }
  if (inssize >= 32) {
    throw(std::invalid_argument("Too large insertion size in CorrelatorAtTimeInfo"));
  }
  uint tcode = timeins; tcode <<= 10; 
  tcode |= timesep; tcode <<= 5; 
  tcode |= inssize; tcode <<= 5;
  tcode |= srcsize; tcode <<= 1;
  if (subvev) tcode |= 1u; 
  tcode <<= 1;
  if (hermitianmatrix) tcode |= 1u;
  icode.back() = tcode;
}

// ******************************************************************************

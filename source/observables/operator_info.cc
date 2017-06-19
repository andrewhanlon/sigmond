#include "operator_info.h"
#include "args_handler.h"
#include "log_helper.h"
#include <stdexcept>

using namespace std;



OperatorInfo::OperatorInfo() : icode(1)
{
 icode[0]=0;     // default is zero particles
}


OperatorInfo::OperatorInfo(XMLHandler& xml_in)
{
 try{
 set<string> tags;
 tags.insert("Operator");
 tags.insert("OperatorString");
 tags.insert("BLOperator");
 tags.insert("BLOperatorString");
 tags.insert("GIOperator");
 tags.insert("GIOperatorString");
 ArgsHandler xin(xml_in,tags);
 string rtag=xin.getInputRootTag();
 if ((rtag=="Operator")||(rtag=="OperatorString")
   ||(rtag=="BLOperator")||(rtag=="BLOperatorString")){
    BasicLapHOperatorInfo buf(xml_in);
    icode=buf.icode;}
 else if ((rtag=="GIOperator")||(rtag=="GIOperatorString")){
    GenIrrepOperatorInfo buf(xml_in);
    icode=buf.icode;}
 else
    throw(std::invalid_argument("Invalid XML input"));
 }
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("OperatorInfo construction failed: \n")
      +string(errmsg.what())+string("\nInput XML:")+xml_in.output()));}
}



OperatorInfo::OperatorInfo(const std::string& opstring, OperatorInfo::OpKind opkind)
{
 if (opkind==BasicLapH){
    BasicLapHOperatorInfo temp(opstring);
    icode=temp.icode;}
 else{
    GenIrrepOperatorInfo temp(opstring);
    icode=temp.icode;}
}


 // ********************************************************************


bool OperatorInfo::isBasicLapH() const
{
 return ((icode[0] & 0x7u)!=7);
}
   
bool OperatorInfo::isGenIrrep() const
{
 return ((icode[0] & 0x7u)==7);
}

BasicLapHOperatorInfo OperatorInfo::getBasicLapH() const
{
 if (isBasicLapH()) 
    return BasicLapHOperatorInfo(icode);
 throw(std::runtime_error("OperatorInfo cannot call getBasicLapH"));
}

GenIrrepOperatorInfo OperatorInfo::getGenIrrep() const
{
 if (isGenIrrep()) 
    return GenIrrepOperatorInfo(icode);
 throw(std::runtime_error("OperatorInfo cannot call getGenIrrep"));
}


OperatorInfo& OperatorInfo::resetGenIrrepIDIndex(uint index)
{
 if (isGenIrrep()){
    GenIrrepOperatorInfo::resetIDIndex(index,icode[1]);}
 else
   throw(std::runtime_error("OperatorInfo cannot call resetGenIrrepIDIndex"));
 return *this;
}


string OperatorInfo::output(bool longform, int indent) const
{
 XMLHandler xmlout;
 output(xmlout,longform);
 return xmlout.output(indent);
}

string OperatorInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}

void OperatorInfo::output(XMLHandler& xmlout, bool longform) const
{
 if (isBasicLapH()){ 
    BasicLapHOperatorInfo temp(icode);
    temp.output(xmlout,longform);}   
 else if (isGenIrrep()){
    GenIrrepOperatorInfo temp(icode);
    temp.output(xmlout,longform);}   
}

string OperatorInfo::short_output() const
{
 if (isBasicLapH()){
    BasicLapHOperatorInfo temp(icode);
    return temp.short_output();}   
 else if (isGenIrrep()){
    GenIrrepOperatorInfo temp(icode);
    return temp.short_output();}
 return string("");
}

bool OperatorInfo::operator==(const OperatorInfo& rhs) const
{
 return multiEqual(icode,rhs.icode);   
}

bool OperatorInfo::operator!=(const OperatorInfo& rhs) const
{
 return multiNotEqual(icode,rhs.icode);   
}

bool OperatorInfo::operator<(const OperatorInfo& rhs) const
{
 return multiLessThan(icode,rhs.icode);   
}


 // *******************************************************************

#include "correlator_matrix_info.h"
#include "multi_compare.h"
#include <algorithm>
#include "mcobs_info.h"

  //  See "operator_info.h" for a discussion of how each operator is
  //  encoded into a vector "icode" of unsigned integers.  For a
  //  zero hadron operator (vacuum), icode has length 1; for a single
  //  hadron operator, icode has length 2; for n>=2 hadrons in the
  //  operator, icode has length 2*n+1.

  //  In a CorrelatorMatrixInfo object, the icode vectors of the source and sink
  //  operators are appended, source first, sink last.  The vacuum
  //  operator is NOT allowed.   For a CorrelatorAtTimeInfo object,
  //  another integer is added at the end which contains the time,
  //  the Hermiticity bit, and the VEV subtraction bit.


using namespace std;


map<string,CorrelatorMatrixInfo> CorrelatorMatrixInfo::m_cormat_namemap;


// ****************************************************************


CorrelatorMatrixInfo::CorrelatorMatrixInfo(XMLHandler& xml_in)
{
 m_current=0;
 try{
    XMLHandler xmlf(xml_in,"CorrelatorMatrixInfo");
    int namecount=xmlf.count("Name");
    if ((namecount>1)||((namecount==1)&&(xmlf.count_children()>1)))
       throw(std::invalid_argument("Multiple <Name> tags or single <Name> with other tags"));
    if (namecount==1){
       string name; xmlreadchild(xmlf,"Name",name);
       name=tidyString(name);
       map<string,CorrelatorMatrixInfo>::const_iterator it=m_cormat_namemap.find(name);
       if (it==m_cormat_namemap.end())
          throw(std::invalid_argument("<Name> tag not associated with any object"));
       m_opinfos=it->second.m_opinfos;
       m_hermitian=it->second.m_hermitian;
       m_vevsubt=it->second.m_vevsubt;}
    else{
       namecount=xmlf.count("AssignName");
       if (namecount>1) throw(std::invalid_argument("Multiple <AssignName> tags"));
       list<string> tagnames;
       tagnames.push_back("Operator");
       tagnames.push_back("OperatorString");
       tagnames.push_back("BLOperator");
       tagnames.push_back("BLOperatorString");
       tagnames.push_back("GIOperator");
       tagnames.push_back("GIOperatorString");
       list<XMLHandler> opxml=xmlf.find_among_children(tagnames);
       for (list<XMLHandler>::iterator
          ot=opxml.begin();ot!=opxml.end();++ot)
             m_opinfos.insert(OperatorInfo(*ot));
       m_hermitian=(xml_tag_count(xmlf,"HermitianMatrix")>0);
       m_vevsubt=(xml_tag_count(xmlf,"SubtractVEV")>0);
       if (namecount==1){
          string name; xmlreadchild(xmlf,"AssignName",name);
          setName(name);}}}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument(string("Invalid XML for CorrelatorMatrixInfo constructor: ")
         +string(errmsg.what())));}
}


void CorrelatorMatrixInfo::setName(const string& corrmatname)
{
 string name=tidyString(corrmatname);
 if (name.length()>256){
    throw(std::invalid_argument("CorrelatorMatrixInfo assign name too long"));}
 if (name.find_first_of("\t\n ")!=string::npos){
    throw(std::invalid_argument("CorrelatorMatrixInfo assign name cannot contain space, tabs, newlines"));}
 if (name.empty()){
    throw(std::invalid_argument("CorrelatorMatrixInfo assign name is empty"));}
 map<string,CorrelatorMatrixInfo>::const_iterator it=m_cormat_namemap.find(name);
 if (it!=m_cormat_namemap.end())
    throw(std::invalid_argument("CorrelatorMatrixInfo <AssignName> tag already associated with an object"));
 m_cormat_namemap.insert(make_pair(name,*this));
}


string CorrelatorMatrixInfo::output(bool longform, int indent) const
{
 XMLHandler xmlout;
 output(xmlout,longform);
 return xmlout.output(indent);
}

string CorrelatorMatrixInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}

void CorrelatorMatrixInfo::output(XMLHandler& xmlout, bool longform) const
{
 xmlout.set_root("CorrelatorMatrixInfo");
 if (m_hermitian){
    xmlout.put_child("HermitianMatrix");}
 if (m_vevsubt){
    xmlout.put_child("SubtractVEV");}
 for (set<OperatorInfo>::const_iterator ot=m_opinfos.begin();
      ot!=m_opinfos.end();ot++){
    XMLHandler xmlop;
    ot->output(xmlop,longform);
    xmlout.put_child(xmlop);}
}


bool CorrelatorMatrixInfo::operator==(const CorrelatorMatrixInfo& rhs) const
{
 return multiEqual(m_hermitian, rhs.m_hermitian,  m_vevsubt, rhs.m_vevsubt)
     && multiEqual(m_opinfos,rhs.m_opinfos);   
}

bool CorrelatorMatrixInfo::operator!=(const CorrelatorMatrixInfo& rhs) const
{
 return multiNotEqual(m_hermitian, rhs.m_hermitian,  m_vevsubt, rhs.m_vevsubt)
     || multiNotEqual(m_opinfos,rhs.m_opinfos);   
}

bool CorrelatorMatrixInfo::operator<(const CorrelatorMatrixInfo& rhs) const
{
 return multiLessThan(m_hermitian, rhs.m_hermitian,  m_vevsubt, rhs.m_vevsubt)
   || ( multiEqual(m_hermitian, rhs.m_hermitian,  m_vevsubt, rhs.m_vevsubt)
       && multiLessThan(m_opinfos,rhs.m_opinfos) );   
}


// ******************************************************************************

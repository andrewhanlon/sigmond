#include "diag_corr_set.h"
#include "args_handler.h"
#include "correlator_matrix_info.h"

using namespace std;

 // ******************************************************************

DiagonalCorrelatorSet::DiagonalCorrelatorSet(XMLHandler& xmlin)
{
 try{
 XMLHandler xin(xmlin,"DiagonalCorrelatorSet");
 ArgsHandler gin(xin);    
 m_subvev=false;
 gin.getOptionalBool("SubtractVEV",m_subvev);
 if (gin.queryTag("Sequential")){
    ArgsHandler gseq(gin,"Sequential");
    XMLHandler xmlf;
    gseq.getInput(xmlf);
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
          m_opset.push_back(OperatorInfo(*ot));}
 else if (gin.queryTag("RotatedSequential")){
    ArgsHandler grot(gin,"RotatedSequential");
    GenIrrepOperatorInfo rop(grot.getItem<GenIrrepOperatorInfo>("RotatedSequential"));
    uint nlevels=grot.getUInt("NumberOfLevels");
    for (uint level=0;level<nlevels;level++){
       m_opset.push_back(OperatorInfo(rop.resetIDIndex(level)));}}
 else if (gin.queryTag("CorrelatorMatrixInfo")){
    CorrelatorMatrixInfo cm(gin.getItem<CorrelatorMatrixInfo>("CorrelatorMatrixInfo"));
    const set<OperatorInfo>& opsref=cm.getOperators();
    m_opset.assign(opsref.begin(),opsref.end());   
    m_subvev=cm.subtractVEV();}
 else{
    uint opindex;
    std::list<XMLHandler> oplist=xin.find_among_children("DiagonalCorrelator");
    uint count=0;
    for (list<XMLHandler>::iterator it=oplist.begin();it!=oplist.end();it++){
       m_opset.push_back(OperatorInfo(*it));
       xmlreadchild(*it,"OperatorIndex",opindex);
       if (opindex!=count) throw(std::runtime_error("Error creating DiagonalCorrelatorSet from XML"));
       count++;}}}
 catch(std::exception& xp){
    throw;}
}


void DiagonalCorrelatorSet::outputOperatorInfos(XMLHandler& xmlout)
{
 xmlout.set_root("DiagonalCorrelatorSet");
 XMLHandler xmlop;
 for (uint k=0;k<m_opset.size();k++){
    XMLHandler xmlc("DiagonalCorrelator");
    m_opset[k].output(xmlop);
    xmlc.put_child("OperatorIndex",make_string(k));
    xmlc.put_child(xmlop);
    xmlout.put_child(xmlc);}
}


void DiagonalCorrelatorSet::output(XMLHandler& xmlout)
{
 outputOperatorInfos(xmlout);
 if (m_subvev) xmlout.put_child("SubtractVEV");
}


void DiagonalCorrelatorSet::addCorrelator(const OperatorInfo& anop)
{
 for (uint k=0;k<m_opset.size();k++)
    if (m_opset[k]==anop) 
       throw(std::runtime_error("Cannot add OperatorInfo already in DiagonalCorrelatorSet"));
 m_opset.push_back(anop);
}

CorrelatorInfo DiagonalCorrelatorSet::getCorrelatorInfo(uint opnum) const
{
 return CorrelatorInfo(m_opset.at(opnum),m_opset.at(opnum));
}


 // ******************************************************************

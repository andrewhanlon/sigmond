#include "obs_get_handler.h"
#include "filelist_info.h"
#include "correlator_matrix_info.h"

using namespace std;

namespace LaphEnv {

// *************************************************************************


MCObsGetHandler::MCObsGetHandler(XMLHandler& xmlin)
                     : m_corrdh(0), m_vevdh(0), m_ensembleptr(0)
{

 XMLHandler xmlr(xmlin);
 list<FileListInfo> corrinputfiles;
 list<FileListInfo> vevinputfiles;
 xml_tag_assert(xmlr,"MCObservables");

    //  read the input file lists to find the data

 if (xmlr.query_unique_to_among_children("CorrelatorData")){
    XMLHandler xmlp(xmlr,"CorrelatorData");
    {list<XMLHandler> infxml;
    infxml=xmlp.find("FileListInfo");  
    for (list<XMLHandler>::iterator 
        it=infxml.begin();it!=infxml.end();++it)
       corrinputfiles.push_back(FileListInfo(*it));}}

 if (xmlr.query_unique_to_among_children("VEVData")){
    XMLHandler xmlp(xmlr,"VEVData");
    {list<XMLHandler> infxml;
    infxml=xmlp.find("FileListInfo");  
    for (list<XMLHandler>::iterator 
        it=infxml.begin();it!=infxml.end();++it)
       vevinputfiles.push_back(FileListInfo(*it));}}
 bool nodata=(corrinputfiles.empty())&&(vevinputfiles.empty());

    //  read requested ensemble info (optionally given)

 if (xml_tag_count(xmlr,"MCEnsembleInfo")>1){
    throw(std::invalid_argument("Multiple <MCEnsembleInfo> tags not allowed"));}
 try{
    MCEnsembleInfo ensemble(xmlr);
    m_ensembleptr=new MCEnsembleInfo(ensemble);}
 catch(const std::exception &xp){
    m_ensembleptr=0;
    if (nodata) throw(std::invalid_argument("No data and missing <MCEnsembleInfo> tag"));}

 bool use_checksums=false;
 if (xml_child_tag_count(xmlr,"UseCheckSums")>1){
    use_checksums=true;}

 if (nodata) return;

    //  read requested observables (optionally given)

 try{

 set<CorrelatorInfo> corrSetNoSym;
 set<CorrelatorInfo> corrSetSym;
 set<OperatorInfo> vevSet;

 int speccount=xml_tag_count(xmlr,"Specifications");
 if (speccount>1){
    throw(std::invalid_argument("Multiple <Specifications> tags not allowed"));}
 else if (speccount==1){

    bool allhermitian=false;
    XMLHandler xmlspec(xmlr,"Specifications");
    if (xml_tag_count(xmlspec,"AllHermitian")>=1){
       allhermitian=true;}
    set<CorrelatorInfo>& cset=(allhermitian)? corrSetSym : corrSetNoSym;

          // individual correlators
    read_correlators(xmlspec,"Correlator",false,cset,vevSet);
          // individual correlators with VEV subtraction
    read_correlators(xmlspec,"CorrelatorWithVEV",true,cset,vevSet);
          // individual VEVs
    read_vevs(xmlspec,"VEV",vevSet);
          // Hermitian correlation matrices
    read_correlator_matrices(xmlspec,"HermitianCorrelationMatrix",
                             true,false,corrSetSym,vevSet);
          // Hermitian correlation matrices with VEVs
    read_correlator_matrices(xmlspec,"HermitianCorrelationMatrixWithVEVs",
                             true,true,corrSetSym,vevSet);
          // non-Hermitian correlation matrices
    read_correlator_matrices(xmlspec,"CorrelationMatrix",
                             allhermitian,false,cset,vevSet);
          // non-Hermitian correlation matrices with VEVs
    read_correlator_matrices(xmlspec,"CorrelationMatrixWithVEVs",
                             allhermitian,true,cset,vevSet);
    }

 bool cspec=(!corrSetNoSym.empty())||(!corrSetSym.empty());
 bool vspec=!vevSet.empty();
 if ((cspec)&&(corrinputfiles.empty()))
    throw(std::invalid_argument("No correlator input files but correlators requested"));
 if ((vspec)&&(vevinputfiles.empty()))
    throw(std::invalid_argument("No VEV input files but VEVs requested"));

 if ((!vevinputfiles.empty()) && (vspec || !cspec)){
    m_vevdh=new VEVDataHandler(
                         vevinputfiles,vevSet,m_ensembleptr,use_checksums);}
 if ((!corrinputfiles.empty()) && (cspec || !vspec)){
    m_corrdh=new CorrelatorDataHandler(
                         corrinputfiles,corrSetNoSym,
                         corrSetSym,m_ensembleptr,use_checksums);}

 if ((m_corrdh==0)&&(m_vevdh==0))
    throw(std::invalid_argument("No data: nothing to do"));

 if (m_ensembleptr==0){
    if (m_corrdh!=0) m_ensembleptr=new MCEnsembleInfo(m_corrdh->getEnsemble());
    else m_ensembleptr=new MCEnsembleInfo(m_vevdh->getEnsemble());}

 if ((m_vevdh!=0)&&(m_corrdh!=0)){
    if (m_vevdh->getEnsemble()!=m_corrdh->getEnsemble()){
       throw(std::invalid_argument("VEV and Correlator Ensembles do not match: fatal"));}}
 }
 catch(const std::exception& msg){
    throw;}
}



MCObsGetHandler::~MCObsGetHandler()
{
 delete m_corrdh;
 delete m_vevdh;
 delete m_ensembleptr;
}


string MCObsGetHandler::getEnsembleId() const
{
 return m_ensembleptr->getId();
}


MCEnsembleInfo MCObsGetHandler::getEnsembleInfo() const
{
 return *m_ensembleptr;
}


unsigned int MCObsGetHandler::getNumberOfMeasurements()
{
 return m_ensembleptr->getNumberOfConfigs();
}


void MCObsGetHandler::getData(const MCObsInfo& obsinfo, 
                              int serial_index, Scalar& data)
{
 if (obsinfo.isHermitianCorrelatorAtTime()){
    if (m_corrdh==0) throw(std::invalid_argument((string("getData failed for ")+obsinfo.str()).c_str()));
    m_corrdh->getSymData(obsinfo.getCorrelatorAtTimeInfo(),serial_index,data);}
 else if (obsinfo.isCorrelatorAtTime()){
    if (m_corrdh==0) throw(std::invalid_argument((string("getData failed for ")+obsinfo.str()).c_str()));
    m_corrdh->getData(obsinfo.getCorrelatorAtTimeInfo(),serial_index,data);}
 else if (obsinfo.isVEV()){
    if (m_vevdh==0) throw(std::invalid_argument((string("getData failed for ")+obsinfo.str()).c_str()));
    m_vevdh->getData(obsinfo.getVEVInfo(),serial_index,data);}
}


bool MCObsGetHandler::queryData(const MCObsInfo& obsinfo, int serial_index)
{
 if (obsinfo.isHermitianCorrelatorAtTime()){
    if (m_corrdh==0) return false;
    return m_corrdh->querySymData(obsinfo.getCorrelatorAtTimeInfo(),serial_index);}
 else if (obsinfo.isCorrelatorAtTime()){
    if (m_corrdh==0) return false;
    return m_corrdh->queryData(obsinfo.getCorrelatorAtTimeInfo(),serial_index);}
 else if (obsinfo.isVEV()){
    if (m_vevdh==0) return false;
    return m_vevdh->queryData(obsinfo.getVEVInfo(),serial_index);}
 return false;
}



void MCObsGetHandler::close()
{
 if (m_corrdh) m_corrdh->close();
 if (m_vevdh) m_vevdh->close();
}


void MCObsGetHandler::getFileMap(XMLHandler& xmlout) const
{
 if (m_corrdh){
    m_corrdh->getFileMap(xmlout);}
 else{
    xmlout.set_root("FileMap");}
 if (m_vevdh){
    XMLHandler xmlv;
    m_vevdh->getFileMap(xmlv);
    list<XMLHandler> ventries=xmlv.find("Entry");
    for (list<XMLHandler>::const_iterator 
         it=ventries.begin();it!=ventries.end();it++)
       xmlout.put_child(*it);
    }
}


set<OperatorInfo> MCObsGetHandler::getVEVInfos() const
{
 if (m_vevdh){
    return m_vevdh->getOperatorSet();}
 return set<OperatorInfo>();
}


set<CorrelatorInfo> MCObsGetHandler::getCorrelatorInfos() const
{
 if (m_corrdh){
    return m_corrdh->getCorrelatorSet();}
 return set<CorrelatorInfo>();
}




   //   private members

           // read individual correlators with/without VEV subtraction

void MCObsGetHandler::read_correlators(XMLHandler& xmlin, 
                          const string& tagname, bool vevs, 
                          set<CorrelatorInfo>& corrSet,
                          set<OperatorInfo>& vevSet)
{
 try{
    list<XMLHandler> corrxml;
    corrxml=xmlin.find_among_children(tagname);
    for (list<XMLHandler>::iterator 
       it=corrxml.begin();it!=corrxml.end();++it){
       CorrelatorInfo corr(*it);
       corrSet.insert(corr);
       if (vevs){
          vevSet.insert(corr.getSource());
          vevSet.insert(corr.getSink());}}}
 catch(const std::exception& xp){}
}


void MCObsGetHandler::read_vevs(XMLHandler& xmls, 
                                const string& tagname,
                                set<OperatorInfo>& vevSet)
{
 try{
    list<XMLHandler> vevxml;
    vevxml=xmls.find_among_children(tagname);
    for (list<XMLHandler>::iterator 
       it=vevxml.begin();it!=vevxml.end();++it){
       OperatorInfo vev(*it);
       vevSet.insert(vev);}}
 catch(const std::exception& xp){}
}


          // read correlation matrices

void MCObsGetHandler::read_correlator_matrices(XMLHandler& xmlin,
                          const string& tagname,
                          bool hermitian, bool vevs, 
                          set<CorrelatorInfo>& corrSet,
                          set<OperatorInfo>& vevSet)
{
 try{
    list<XMLHandler> corrxml;
    corrxml=xmlin.find_among_children(tagname);
    for (list<XMLHandler>::iterator 
       it=corrxml.begin();it!=corrxml.end();++it){
       list<XMLHandler> opxml,opxml2;
       opxml=it->find_among_children("Operator");
       opxml2=it->find_among_children("OperatorString");
       opxml.splice(opxml.end(),opxml2);
       set<OperatorInfo> opSet;
       for (list<XMLHandler>::iterator
          ot=opxml.begin();ot!=opxml.end();++ot)
             opSet.insert(OperatorInfo(*ot));
       for (set<OperatorInfo>::const_iterator
          ita=opSet.begin();ita!=opSet.end();ita++){
          if (vevs){ vevSet.insert(*ita);}
          for (set<OperatorInfo>::const_iterator
             itb=(hermitian)?ita:opSet.begin();itb!=opSet.end();itb++){
             CorrelatorInfo corr(*ita,*itb);
             corrSet.insert(corr);}}
       int namecount=it->count("AssignName");
       if (namecount>1) throw(std::invalid_argument("Multiple <AssignName> tags"));
       if (namecount==1){
          string name; xmlreadchild(*it,"AssignName",name);
          CorrelatorMatrixInfo corrkeep(opSet,hermitian,vevs);
          corrkeep.setName(name);}}}
 catch(const std::exception& xp){}
}


// ***************************************************************************************
}
 

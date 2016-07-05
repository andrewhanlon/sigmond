#include "vev_data_handler.h"

using namespace std;

namespace LaphEnv {

// *************************************************************************


VEVDataHandler::VEVDataHandler(
                const list<FileListInfo>& inputfiles,
                const set<OperatorInfo>& opSet,
                const MCEnsembleInfo *ensemble,
                bool use_checksums)
        : m_getter(0), m_ensembleptr(0)
{
 try{

 if (ensemble!=0) m_ensembleptr=new MCEnsembleInfo(*ensemble);

   //  determine the stubs vector and the fileMapper which 
   //  maps an Operator info -> stub + suffix)

 string fileId("Laph--VEVFile");   
 map<FileKey,pair<int,int> > fileMapper;
 vector<string> stubs(inputfiles.size());
 IOMap<RecordKey,DataType> iom;
 string header;
 int stubcount=0;
 for (list<FileListInfo>::const_iterator ft=inputfiles.begin();
  ft!=inputfiles.end();++ft,stubcount++){
  stubs[stubcount]=ft->getFileStub();
  for (int suffix=ft->getMinFileNumber();suffix<=ft->getMaxFileNumber();suffix++){
    string fname=ft->getFileName(suffix);
    if (!fileExists(fname)) continue;
    if (iom.peekHeader(header,fname,fileId)){
       XMLHandler xmlh;
       xmlh.set_from_string(header);
       string numtype;
       xmlread(xmlh,"NumberType",numtype,"CorrelatorDataHandler");
       check_number_type(numtype);
       MCEnsembleInfo currmc(xmlh);
       if (m_ensembleptr==0){
          m_ensembleptr=new MCEnsembleInfo(currmc);}
       else if (currmc!=(*m_ensembleptr)){
          throw(std::invalid_argument("An incorrect Monte Carlo ensemble was encountered"));}
       FileKey opinfo(xmlh);
       if ((opSet.empty())||(opSet.find(opinfo)!=opSet.end()))
          fileMapper.insert(make_pair(opinfo,make_pair(stubcount,suffix)));}
    }}
 if (m_ensembleptr==0) throw(std::invalid_argument("LapH ensemble not found"));

 bool cmissing=false;
 string errmsg;
 for (set<OperatorInfo>::const_iterator ct=opSet.begin();ct!=opSet.end();ct++){
    if (fileMapper.find(*ct)==fileMapper.end()){
       errmsg+="Following Operator not found in data files:\n";
       errmsg+=ct->output()+"\n";
       cmissing=true;}}
 if (cmissing) throw(std::invalid_argument(errmsg.c_str()));

 m_getter=new DataGetHandlerMF<FileKey,RecordKey,DataType>(
               stubs,fileMapper,fileId,maxgetopen,cleanfrac,use_checksums);
 }
 catch(const std::exception& errmsg){
    string out("Invalid VEVDataHandler initialization: ");
    out+=string(errmsg.what())+string("\n");
    delete m_ensembleptr;
    throw(std::invalid_argument(out.c_str()));}
}



VEVDataHandler::~VEVDataHandler()
{
 delete m_ensembleptr;
 delete m_getter;
}




set<VEVDataHandler::FileKey> VEVDataHandler::getOperatorSet() const
{
 return m_getter->getFileKeys();
}


MCEnsembleInfo VEVDataHandler::getEnsemble() const
{
 return (*m_ensembleptr);
}


string VEVDataHandler::getEnsembleId() const
{
 return m_ensembleptr->getId();
}


unsigned int VEVDataHandler::getNumberOfMeasurements()
{
 return m_ensembleptr->getNumberOfConfigs();
}


void VEVDataHandler::getData(const OperatorInfo& mckey, 
                             int serial_index, Scalar& data)
{
 InScalar buffer;
 m_getter->getData(mckey,serial_index,buffer);
 data=buffer;
}


bool VEVDataHandler::queryData(const OperatorInfo& mckey, 
                               int serial_index)
{
 return m_getter->queryData(mckey,serial_index);
}


void VEVDataHandler::close()
{
 m_getter->close();
}


void VEVDataHandler::getFileMap(XMLHandler& xmlout) const
{
 XMLHandler xmlv;
 m_getter->getFileMap(xmlv);
 xmlout.set_root("FileMap");
 list<XMLHandler> ventries=xmlv.find("Entry");
 for (list<XMLHandler>::iterator 
      it=ventries.begin();it!=ventries.end();it++){
    it->rename_tag("VEV");
    XMLHandler xmlt("Entry");
    xmlt.put_child(*it);
    xmlout.put_child(xmlt);}
}


std::set<OperatorInfo> VEVDataHandler::getFileKeys() const
{
 return m_getter->getFileKeys();
}


std::set<VEVDataHandler::RecordKey> VEVDataHandler::getKeys(const OperatorInfo& fkey)
{
 return m_getter->getKeys(fkey);
}


void VEVDataHandler::outputKeys(XMLHandler& xmlout)
{
 m_getter->outputKeys(xmlout);
}


// ***************************************************************************************
}

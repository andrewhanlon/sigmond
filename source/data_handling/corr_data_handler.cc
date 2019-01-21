#include "corr_data_handler.h"

using namespace std;

namespace LaphEnv {


// *************************************************************************


BLCorrelatorDataHandler::BLCorrelatorDataHandler(
                const list<FileListInfo>& inputfiles,
                const set<CorrelatorInfo>& corrSetNoSym,
                const set<CorrelatorInfo>& corrSetSym,
                const MCEnsembleInfo *ensemble,
                bool use_checksums)
        : m_getter(0), m_ensembleptr(0)
{
 try{

 if (ensemble!=0) m_ensembleptr=new MCEnsembleInfo(*ensemble);

   //  determine the stubs vector and the fileMapper which 
   //  maps a Correlator info -> stub + suffix)

 string fileId("Laph--CorrelatorFile");
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
       xmlread(xmlh,"NumberType",numtype,"BLCorrelatorDataHandler");
       check_number_type(numtype);
       MCEnsembleInfo currmc(xmlh);
       if (m_ensembleptr==0){
          m_ensembleptr=new MCEnsembleInfo(currmc);}
       else if (currmc!=(*m_ensembleptr)){
          throw(std::invalid_argument("An incorrect Monte Carlo ensemble was encountered"));}
       FileKey corrinfo(xmlh);
       if ( ((corrSetNoSym.empty())&&(corrSetSym.empty()))
          || (corrSetNoSym.find(corrinfo)!=corrSetNoSym.end())
          || (corrSetSym.find(corrinfo)!=corrSetSym.end())
          || (corrSetSym.find(corrinfo.getTimeFlipped())!=corrSetSym.end()) )
          fileMapper.insert(make_pair(corrinfo,make_pair(stubcount,suffix)));}
    }}
 if (m_ensembleptr==0) throw(std::invalid_argument("LapH ensemble not found"));

 bool cmissing=false;
 string errmsg;
 for (set<CorrelatorInfo>::const_iterator ct=corrSetNoSym.begin();ct!=corrSetNoSym.end();ct++){
    if (fileMapper.find(*ct)==fileMapper.end()){
       errmsg+="Following correlator not found in data files:\n";
       errmsg+=ct->output()+"\n";
       cmissing=true;}}
 for (set<CorrelatorInfo>::const_iterator ct=corrSetSym.begin();ct!=corrSetSym.end();ct++){
    if   ((fileMapper.find(*ct)==fileMapper.end())
       && (fileMapper.find(ct->getTimeFlipped())==fileMapper.end())){
       errmsg+="Following correlator and time-flipped not found in data files:\n";
       errmsg+=ct->output()+"\n";
       cmissing=true;}}
 if (cmissing) throw(std::invalid_argument(errmsg));

 m_getter=new LapHDataGetHandlerMF<FileKey,RecordKey,DataType>(
               stubs,fileMapper,fileId,maxgetopen,cleanfrac,use_checksums);

 }
 catch(const std::exception& errmsg){
    string out("Invalid BLCorrelatorDataHandler initialization: ");
    out+=string(errmsg.what())+string("\n");
    delete m_ensembleptr;
    throw(std::invalid_argument(out));}
}


BLCorrelatorDataHandler::~BLCorrelatorDataHandler()
{
 delete m_ensembleptr;
 delete m_getter;
}



set<CorrelatorInfo> BLCorrelatorDataHandler::getCorrelatorSet() const
{
 return m_getter->getFileKeys();
}


MCEnsembleInfo BLCorrelatorDataHandler::getEnsemble() const
{
 return (*m_ensembleptr);
}


string BLCorrelatorDataHandler::getEnsembleId() const
{
 return m_ensembleptr->getId();
}


unsigned int BLCorrelatorDataHandler::getNumberOfMeasurements()
{
 return m_ensembleptr->getNumberOfMeasurements();
}


void BLCorrelatorDataHandler::getData(const CorrelatorAtTimeInfo& mckey, 
                                    int serial_index, Scalar& data)
{
 FileKey fkey(mckey.getCorrelator());
 RecordKey rkey(mckey.getTimeSeparation(),serial_index);
 DataType buffer;
 m_getter->getData(fkey,rkey,buffer);
 data=buffer[0];
 for (unsigned int k=1;k<buffer.size();k++) data+=buffer[k];
}


void BLCorrelatorDataHandler::getSymData(const CorrelatorAtTimeInfo& mckey, 
                                       int serial_index, Scalar& data)
{
 CorrelatorInfo fkey(mckey.getCorrelator());
 RecordKey rkey(mckey.getTimeSeparation(),serial_index);
 DataType buffer;
 bool flag=false;
 if (m_getter->queryData(fkey,rkey)){
    m_getter->getData(fkey,rkey,buffer);
    data=buffer[0];
    for (unsigned int k=1;k<buffer.size();k++) data+=buffer[k];
    if (fkey.isSinkSourceSame()) return;
    flag=true;}
 CorrelatorInfo ffkey(fkey.getTimeFlipped());
 if (m_getter->queryData(ffkey,rkey)){
    m_getter->getData(ffkey,rkey,buffer);
    Scalar data2;
    data2=buffer[0];
    for (unsigned int k=1;k<buffer.size();k++) data2+=buffer[k];
    if (flag){
       data+=conjugate(data2); data*=0.5;}
    else{
       data=conjugate(data2);}
    flag=true;}
 if (!flag){
    throw(std::invalid_argument(string("getSymData failed for ")+mckey.str()));}
}


bool BLCorrelatorDataHandler::queryData(const CorrelatorAtTimeInfo& mckey, 
                                      int serial_index)
{
 FileKey fkey(mckey.getCorrelator());
 RecordKey rkey(mckey.getTimeSeparation(),serial_index);
 return m_getter->queryData(fkey,rkey);
}


bool BLCorrelatorDataHandler::querySymData(const CorrelatorAtTimeInfo& mckey, 
                                         int serial_index)
{
 CorrelatorInfo fkey(mckey.getCorrelator());
 RecordKey rkey(mckey.getTimeSeparation(),serial_index);
 if (m_getter->queryData(fkey,rkey)) return true;
 CorrelatorInfo ffkey(fkey.getTimeFlipped());
 return m_getter->queryData(ffkey,rkey);
}


void BLCorrelatorDataHandler::close()
{
 m_getter->close();
}


void BLCorrelatorDataHandler::getFileMap(XMLHandler& xmlout) const
{
 m_getter->getFileMap(xmlout);
}


std::set<CorrelatorInfo> BLCorrelatorDataHandler::getFileKeys() const
{
 return m_getter->getFileKeys();
}


std::set<BLCorrelatorDataHandler::RecordKey> 
    BLCorrelatorDataHandler::getKeys(const CorrelatorInfo& fkey)
{
 return m_getter->getKeys(fkey);
}


std::string BLCorrelatorDataHandler::getFileName(const CorrelatorInfo& fkey)
{
 return m_getter->getFileName(fkey);
}


void BLCorrelatorDataHandler::outputKeys(XMLHandler& xmlout)
{
 m_getter->outputKeys(xmlout);
}


// ***************************************************************************************
}
 

#include "ensemble_info.h"
#include <stdexcept>
using namespace std;

 // *************************************************************


std::string MCEnsembleInfo::m_known_ensembles_filename=DEFAULTENSFILE;


 // *************************************************************

MCEnsembleInfo::MCEnsembleInfo(XMLHandler& xml_in)
{
 XMLHandler xmlr(xml_in);
 xmlread(xmlr,"MCEnsembleInfo",m_id,"MCEnsembleInfo");
 initialize();
}


MCEnsembleInfo::MCEnsembleInfo(const string& in_id) : m_id(in_id)
{
 initialize();
}


void MCEnsembleInfo::initialize()
{
 m_is_weighted=false;
 try{
 XMLHandler xmlf;
 xmlf.set_from_file(m_known_ensembles_filename);
 XMLHandler xmlens(xmlf,"Infos");
 list<XMLHandler> xmlt(xmlens.find("EnsembleInfo"));
 for (list<XMLHandler>::iterator it=xmlt.begin();it!=xmlt.end();++it){
    string eid;
    xmlread(*it,"Id",eid,"MCEnsembleInfo");
    if (eid==m_id){
       xmlread(*it,"NStreams",n_streams,"MCEnsembleInfo");
       xmlread(*it,"NMeas",n_meas,"MCEnsembleInfo");
       xmlread(*it,"NSpace",n_x,"MCEnsembleInfo");
       xmlread(*it,"NTime",n_t,"MCEnsembleInfo");
       n_y=n_z=n_x;
       if (it->count("Weighted")==1){
          m_is_weighted=true;}
      return;}}}
 catch(const std::exception& xp){}
 if (!parse(m_id)){
   //cerr << "Invalid MCEnsembleInfo id string"<<endl;
   throw(std::invalid_argument("Invalid MCEnsembleInfo id string"));}
}


bool MCEnsembleInfo::parse(const string& idstr)
{
 m_id=tidyString(idstr);
 vector<size_t> pos;
 for (uint k=0;k<idstr.length();k++){
    if (idstr[k]=='|') pos.push_back(k);}
 if (pos.size()!=6) return false;
 try{
 extract_from_string(idstr.substr(pos[0]+1,pos[1]-pos[0]-1),n_meas);
 extract_from_string(idstr.substr(pos[1]+1,pos[2]-pos[1]-1),n_streams);
 extract_from_string(idstr.substr(pos[2]+1,pos[3]-pos[2]-1),n_x);
 extract_from_string(idstr.substr(pos[3]+1,pos[4]-pos[3]-1),n_y);
 extract_from_string(idstr.substr(pos[4]+1,pos[5]-pos[4]-1),n_z);
 extract_from_string(idstr.substr(pos[5]+1,string::npos),n_t);}
 catch(const std::exception& xp){
   return false;}
 return true;
}


MCEnsembleInfo::MCEnsembleInfo(const std::string& id, uint num_meas, uint num_streams,
                               uint nx, uint ny, uint nz, uint nt)
    : m_id(tidyString(id)), n_meas(num_meas), n_streams(num_streams),
      n_x(nx), n_y(ny), n_z(nz), n_t(nt)
{
 m_is_weighted=false;
 ostringstream oss;
 oss <<"|"<<n_meas<<"|"<<n_streams<<"|"<<n_x<<"|"<<n_y<<"|"<<n_z<<"|"<<n_t;
 m_id+=oss.str();
}


MCEnsembleInfo::MCEnsembleInfo(const MCEnsembleInfo& fin)
   : m_id(fin.m_id), n_meas(fin.n_meas), n_streams(fin.n_streams),
     n_x(fin.n_x), n_y(fin.n_y), n_z(fin.n_z), n_t(fin.n_t),
     m_is_weighted(fin.m_is_weighted)
 {}


MCEnsembleInfo& MCEnsembleInfo::operator=(const MCEnsembleInfo& fin)
{
 m_id=fin.m_id;
 n_meas=fin.n_meas;
 n_streams=fin.n_streams;
 n_t=fin.n_t;
 n_x=fin.n_x;
 n_y=fin.n_y;
 n_z=fin.n_z;
 m_is_weighted=fin.m_is_weighted;
 return *this;
}


void MCEnsembleInfo::checkEqual(const MCEnsembleInfo& rhs) const
{
 if (m_id!=rhs.m_id){
    cerr << "MCEnsembleInfo checkEqual failed"<<endl;
    cerr << "LHS:"<<endl<<output()<<endl<<"RHS:"<<endl<<rhs.output()<<endl;
    throw(std::invalid_argument("MCEnsembleInfo checkEqual failed"));}
}


void MCEnsembleInfo::output(XMLHandler& xmlout) const
{
 xmlout.set_root("MCEnsembleInfo",m_id);
}

string MCEnsembleInfo::output(int indent) const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.output(indent);
}

string MCEnsembleInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}

void MCEnsembleInfo::getWeights(std::vector<double>& weights) const
{
 XMLHandler xmlf;
 xmlf.set_from_file(m_known_ensembles_filename);
 XMLHandler xmlcls(xmlf,"CLSEnsembleWeights");
 list<XMLHandler> xmlw(xmlcls.find("Ensemble"));
 for (list<XMLHandler>::iterator it=xmlw.begin();it!=xmlw.end();++it){
    string eid;
    xmlread(*it,"Id",eid,"MCEnsembleInfo");
    if (eid==m_id){
       xmlread(*it,"Weights",weights,"MCEnsembleInfo");
       if (weights.size()<n_meas)
          throw(std::invalid_argument("Insufficient number of weights"));
       return;}}
 throw(std::invalid_argument("Weights not available"));
}

// ***************************************************************

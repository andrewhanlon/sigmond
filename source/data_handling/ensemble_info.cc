#include "ensemble_info.h"
#include <stdexcept>
using namespace std;

namespace LaphEnv {

 // *************************************************************


MCEnsembleInfo::MCEnsembleInfo(XMLHandler& xml_in)
{
 XMLHandler xmlr(xml_in);
 xmlread(xmlr,"MCEnsembleInfo",m_id,"MCEnsembleInfo");

 if (m_id=="clover_s24_t128_ud840_s743"){
    n_streams=4;  n_configs=551;}
 else if (m_id=="clover_s24_t128_ud860_s743"){
    n_streams=1;  n_configs=584;}
 else if (m_id=="clover_s32_t256_ud860_s743"){
    n_streams=2;  n_configs=412;}
 else if (m_id=="clover_s16_t128_ud840_s743"){
    n_streams=1;  n_configs=100;}
 else{
    cerr << "Invalid MCEnsembleInfo id string"<<endl;
    throw(std::invalid_argument("Invalid MCEnsembleInfo id string"));}
}


MCEnsembleInfo::MCEnsembleInfo(const MCEnsembleInfo& fin) 
   : m_id(fin.m_id), n_configs(fin.n_configs),n_streams(fin.n_streams) {}


MCEnsembleInfo& MCEnsembleInfo::operator=(const MCEnsembleInfo& fin)
{
 m_id=fin.m_id;
 n_configs=fin.n_configs;
 n_streams=fin.n_streams;
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


uint MCEnsembleInfo::getLatticeTimeExtent() const
{
 if (m_id=="clover_s24_t128_ud840_s743") return 128;
 else if (m_id=="clover_s24_t128_ud860_s743") return 128;
 else if (m_id=="clover_s32_t256_ud860_s743") return 256;
 else if (m_id=="clover_s16_t128_ud840_s743") return 128;
 else throw(std::invalid_argument("invalid"));
}

uint MCEnsembleInfo::getLatticeXExtent() const
{
 if (m_id=="clover_s24_t128_ud840_s743") return 24;
 else if (m_id=="clover_s24_t128_ud860_s743") return 24;
 else if (m_id=="clover_s32_t256_ud860_s743") return 32;
 else if (m_id=="clover_s16_t128_ud840_s743") return 16;
 else throw(std::invalid_argument("invalid"));
}

uint MCEnsembleInfo::getLatticeYExtent() const
{
 if (m_id=="clover_s24_t128_ud840_s743") return 24;
 else if (m_id=="clover_s24_t128_ud860_s743") return 24;
 else if (m_id=="clover_s32_t256_ud860_s743") return 32;
 else if (m_id=="clover_s16_t128_ud840_s743") return 16;
 else throw(std::invalid_argument("invalid"));
}

uint MCEnsembleInfo::getLatticeZExtent() const
{
 if (m_id=="clover_s24_t128_ud840_s743") return 24;
 else if (m_id=="clover_s24_t128_ud860_s743") return 24;
 else if (m_id=="clover_s32_t256_ud860_s743") return 32;
 else if (m_id=="clover_s16_t128_ud840_s743") return 16;
 else throw(std::invalid_argument("invalid"));
}


// ***************************************************************
}


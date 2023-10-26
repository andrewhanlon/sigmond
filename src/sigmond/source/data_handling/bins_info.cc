#include "bins_info.h"
#include "args_handler.h"
#include <stdexcept>
#include <algorithm>
using namespace std;

 // *************************************************************


MCBinsInfo::MCBinsInfo(XMLHandler& xml_in)
{
 uint rebin=1;
 vector<int> ovec;
 try{
 ArgsHandler xmli(xml_in);
 if (xmli.queryTag("MCBinsInfo")){
    ArgsHandler xmlr(xmli,"MCBinsInfo");
    ArgsHandler xmlen(xmlr,"MCEnsembleInfo");
    XMLHandler xmlenn; xmlen.getInput(xmlenn);
    m_ensemble=new MCEnsembleInfo(xmlenn);
    if (xmlr.queryTag("TweakEnsemble")){
       ArgsHandler xmlt(xmlr,"TweakEnsemble");
       xmlt.getOptionalUInt("Rebin",rebin);
       if (xmlt.queryTag("Omissions")){
          ovec=xmlt.getIntVector("Omissions");}
       uint keepfirst=0;
       xmlt.getOptionalUInt("KeepFirst",keepfirst);
       for (uint k=0;k<keepfirst;++k) ovec.push_back(k);
       uint keeplast=m_ensemble->getNumberOfMeasurements()-1;
       xmlt.getOptionalUInt("KeepLast",keeplast);
       for (uint k=keeplast+1;k<m_ensemble->getNumberOfMeasurements();k++) 
          ovec.push_back(k);}}
 else if (xmli.queryTag("MCEnsembleInfo")){
    XMLHandler xmlenn; xmli.getInput(xmlenn);
    m_ensemble=new MCEnsembleInfo(xmlenn);}
 else{
    throw(std::invalid_argument("Invalid XML "));}}
 catch(const std::exception& xp){
    cerr << "Invalid MCBinsInfo construction"<<endl;
    throw(std::invalid_argument(string("Invalid MCBinsInfo construction: ")+xp.what()));}

 m_nmeasures=m_ensemble->getNumberOfMeasurements();
 setRebin(rebin);
 if (!(ovec.empty())){
    set<int> omits(ovec.begin(),ovec.end());
    addOmissions(omits);} 
}


MCBinsInfo::MCBinsInfo(const MCBinsInfo& fin) 
   :  m_ensemble(0), m_nmeasures(fin.m_nmeasures), m_rebin(fin.m_rebin),
      m_omit(fin.m_omit), m_nbins(fin.m_nbins) 
{
 m_ensemble=new MCEnsembleInfo(*(fin.m_ensemble));
}


MCBinsInfo& MCBinsInfo::operator=(const MCBinsInfo& fin)
{
 m_ensemble=new MCEnsembleInfo(*(fin.m_ensemble));
 m_nmeasures=fin.m_nmeasures; 
 m_rebin=fin.m_rebin;     
 m_omit=fin.m_omit;      
 m_nbins=fin.m_nbins;     
 return *this;
}


MCBinsInfo::MCBinsInfo(const MCEnsembleInfo& ens)
{
 m_ensemble=new MCEnsembleInfo(ens);
 m_nmeasures=ens.getNumberOfMeasurements(); 
 m_nbins=m_nmeasures;
 m_rebin=1;  
}


void MCBinsInfo::setRebin(int rebin)
{
 if ((rebin>0)&&(rebin<int(m_nmeasures))){
    m_rebin=rebin;
    reset_nbins();}
 else{
    cout << "Invalid value for rebin"<<endl;
    throw(std::invalid_argument("Invalid rebin value"));}
}


void MCBinsInfo::addOmission(int index)
{
 if ((index>=0)&&(index<int(m_nmeasures))){
    m_omit.insert(index);
    reset_nbins();}
}


void MCBinsInfo::addOmissions(set<int> indices)
{
 for (set<int>::const_iterator it=indices.begin();it!=indices.end();it++){
    if ((*it>=0)&&(*it<int(m_nmeasures))) m_omit.insert(*it);}
 reset_nbins();
}


void MCBinsInfo::clearOmissions()
{
 m_omit.clear();
 reset_nbins();
}


void MCBinsInfo::reset_nbins()
{
 m_nbins=(m_nmeasures-m_omit.size())/m_rebin;
}



uint MCBinsInfo::getLatticeTimeExtent() const
{
 return m_ensemble->getLatticeTimeExtent();
}

uint MCBinsInfo::getLatticeXExtent() const
{
 return m_ensemble->getLatticeXExtent();
}

uint MCBinsInfo::getLatticeYExtent() const
{
 return m_ensemble->getLatticeYExtent();
}

uint MCBinsInfo::getLatticeZExtent() const
{
 return m_ensemble->getLatticeZExtent();
}



void MCBinsInfo::checkEqual(const MCBinsInfo& rhs) const
{
 if (*this!=rhs){
    cerr << "MCBinsInfo checkEqual failed"<<endl;
    cerr << "LHS:"<<endl<<output()<<endl<<"RHS:"<<endl<<rhs.output()<<endl;
    throw(std::invalid_argument("MCBinsInfo checkEqual failed"));}
}


bool MCBinsInfo::operator==(const MCBinsInfo& rhs) const 
{
 return ((*m_ensemble==*(rhs.m_ensemble))&&(m_rebin==rhs.m_rebin)
            &&(m_omit==rhs.m_omit));
}

bool MCBinsInfo::operator!=(const MCBinsInfo& rhs) const 
{
 return ((*m_ensemble!=*(rhs.m_ensemble))||(m_rebin!=rhs.m_rebin)
            ||(m_omit!=rhs.m_omit));
}

   //  rhs matches ensemble and omissions, but rhs rebin factor 
   //  can be multiple of rebin factor of this object

bool MCBinsInfo::isConsistentWith(const MCBinsInfo& rhs) const 
{
 if (*m_ensemble!=*(rhs.m_ensemble)) return false;
 if (m_omit!=rhs.m_omit) return false;
 if ((rhs.m_rebin%m_rebin)!=0) return false;
 return true;
}


void MCBinsInfo::output(XMLHandler& xmlout) const
{
 xmlout.set_root("MCBinsInfo");
 XMLHandler xmltmp;
 m_ensemble->output(xmltmp);
 xmlout.put_child(xmltmp);
 xmlout.put_child("NumberOfMeasurements",make_string(m_nmeasures));
 xmlout.put_child("NumberOfBins",make_string(m_nbins));
 if ((m_rebin>1)||(!m_omit.empty())){
    XMLHandler xmltw("TweakEnsemble");
    if (m_rebin>1) xmltw.put_child("Rebin",make_string(m_rebin));
    if (!(m_omit.empty())){
       string omstr;
       for (set<unsigned int>::const_iterator it=m_omit.begin();it!=m_omit.end();it++)
          omstr+=make_string(*it)+" ";
       omstr.erase(omstr.length()-1);
       xmltw.put_child("Omissions",omstr);}
    xmlout.put_child(xmltw);}
}

string MCBinsInfo::output(int indent) const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.output(indent);
}

string MCBinsInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}


// ***************************************************************


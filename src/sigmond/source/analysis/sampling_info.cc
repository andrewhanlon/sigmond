#include "sampling_info.h"
#include "args_handler.h"
#include <stdexcept>
using namespace std;

 // *************************************************************


MCSamplingInfo::MCSamplingInfo() 
   :     icode(0), rngseed(0)
{}


MCSamplingInfo::MCSamplingInfo(XMLHandler& xml_in)
{
 try{
    ArgsHandler xmlr(xml_in,"MCSamplingInfo");
    bool jack=xmlr.queryTag("Jackknife");
    bool boot=xmlr.queryTag("Bootstrapper");
    if (jack==boot){
       throw(std::invalid_argument("Invalid mode in MCSamplingInfo"));}
    if (jack){
       icode=0; rngseed=0; return;}
    ArgsHandler xmlb(xmlr,"Bootstrapper");
    uint num_resamplings=1024;
    uint bootseed=0;
    uint bootskip=64;
    xmlb.getOptionalUInt("NumberResamplings",num_resamplings);
    xmlb.getOptionalUInt("Seed",bootseed);
    xmlb.getOptionalUInt("BootSkip",bootskip);
    encode_bootstrap(num_resamplings,bootseed,bootskip);}
 catch(std::exception& xp){
    throw(std::invalid_argument(string("MCSamplingInfo creation failed: ")
          +xp.what()));}
}


MCSamplingInfo::MCSamplingInfo(const MCSamplingInfo& fin) 
   :     icode(fin.icode), rngseed(fin.rngseed)
{}


MCSamplingInfo& MCSamplingInfo::operator=(const MCSamplingInfo& fin)
{
 icode=fin.icode;
 rngseed=fin.rngseed;
 return *this;
}


MCSamplingInfo::MCSamplingInfo(uint nbootsamp, unsigned long bootseed, uint bootskip)
{
 encode_bootstrap(nbootsamp,bootseed,bootskip);
}



bool MCSamplingInfo::isJackknifeMode() const
{
 return (icode==0);
}


bool MCSamplingInfo::isBootstrapMode() const
{
 return (icode!=0);
}


SamplingMode MCSamplingInfo::getSamplingMode() const
{
 return (icode==0) ? Jackknife : Bootstrap;
}


void MCSamplingInfo::setToJackknifeMode()
{
 icode=0; rngseed=0;
}


void MCSamplingInfo::setToBootstrapMode(const Bootstrapper& boot)
{
 encode_bootstrap(boot.getNumberOfResamplings(),boot.getRNGSeed(),
                  boot.getSkipValue());
}


unsigned int MCSamplingInfo::getNumberOfReSamplings(const MCBinsInfo& binfo) const
{
 if (icode==0) return binfo.getNumberOfBins();
 return icode&0xFFFFFu;
}

unsigned long MCSamplingInfo::getRNGSeed() const   // returns zero for jackknife mode
{
 return rngseed;
}


unsigned int MCSamplingInfo::getSkipValue() const   // returns zero for jackknife mode
{
 return (icode>>20);
}


void MCSamplingInfo::checkEqual(const MCSamplingInfo& rhs) const
{
 if ((icode!=rhs.icode)||(rngseed!=rhs.rngseed)){
    cerr << "MCSamplingInfo checkEqual failed"<<endl;
    cerr << "LHS:"<<endl<<output()<<endl<<"RHS:"<<endl<<rhs.output()<<endl;
    throw(std::invalid_argument("MCSamplingInfo checkEqual failed"));}
}


bool MCSamplingInfo::operator==(const MCSamplingInfo& rhs) const 
{
 return ((icode==rhs.icode)&&(rngseed==rhs.rngseed));
}

bool MCSamplingInfo::operator!=(const MCSamplingInfo& rhs) const 
{
 return ((icode!=rhs.icode)||(rngseed!=rhs.rngseed));
}


void MCSamplingInfo::output(XMLHandler& xmlout) const
{
 xmlout.set_root("MCSamplingInfo");
 if (icode==0){
    xmlout.put_child("Jackknife");
    return;}
 XMLHandler xmlb("Bootstrapper");
 uint nsamp=icode&0xFFFFFu;
 uint skip=icode>>20;
 xmlb.put_child("NumberResamplings",make_string(nsamp));
 xmlb.put_child("Seed",make_string(rngseed));
 xmlb.put_child("BootSkip",make_string(skip));
 xmlout.put_child(xmlb);
}


string MCSamplingInfo::output(int indent) const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.output(indent);
}


string MCSamplingInfo::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}



void MCSamplingInfo::encode_bootstrap(unsigned int num_resamplings, 
                                      unsigned long seed, unsigned int skip_value)
{
 if (num_resamplings<1)
    throw(std::invalid_argument("Number of resamplings > 0 required in MCSamplingInfo"));
 if (num_resamplings>=1048576)
    throw(std::invalid_argument("Number of resamplings too large in MCSamplingInfo"));
 if (skip_value>=4096)
    throw(std::invalid_argument("Skip value too large in MCSamplingInfo"));
 icode=skip_value;
 icode<<=20;
 icode|=num_resamplings;
 rngseed=seed;
}




// ***************************************************************


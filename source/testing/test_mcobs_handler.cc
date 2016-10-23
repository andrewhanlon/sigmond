#include "mcobs_handler.h"
#include "test_obs_get_handler.h"
#include "test_mcobs_handler.h"
#include <cstdio>
#include <ctime>
#include <map>
#include <cstdlib>
#include <sstream>

using namespace std;
using namespace LaphEnv;

// subroutines from test_obs_get_handler.cc

bool get_binflag(int nconfig, const MCObsInfo& mcobs, const set<int>& omit,
                 const map<MCObsInfo,set<int> >& fmissing)
{
 bool binflag;
 map<MCObsInfo,set<int> >::const_iterator bt=fmissing.find(mcobs);
 if (bt==fmissing.end())
    binflag=false; 
 else{
    binflag=true;
    const set<int> *bptr=&(bt->second);
    for (int serind=0;serind<nconfig;++serind)
       if ((omit.find(serind)==omit.end())&&(bptr->find(serind)!=bptr->end())){
          binflag=false; break;}}
 return binflag;
}


 // **********************************


BinMapper::BinMapper(int nconfigs, int rebin, const set<int>& omit)
{
 int kcf=0;
 list<int> confs;
 for (int k=0;k<nconfigs;++k){
    if (omit.find(k)==omit.end()){
       if (kcf==0){
          confs.clear(); confs.push_back(k); ++kcf;}
       else{
          confs.push_back(k); ++kcf;}
       if (kcf==rebin){
          configs.push_back(confs);
          kcf=0;}}}
}



bool getVEVBin(VEVCorrect& VC, const BinMapper& BM, 
               const MCObsInfo& mcobs, int bin, double& result)
{
 try{
    const list<int>& configs=BM.getBinConfigs(bin);
    result=0.0;
    double temp;
    for (list<int>::const_iterator it=configs.begin();it!=configs.end();++it){
        if (!VC.getCorrect(mcobs,*it,temp)) return false;
        result+=temp;}
    result/=double(configs.size());}
 catch(const std::exception& xp){
    return false;}
 return true;
}


bool getCorBin(CorrCorrect& CC, const BinMapper& BM, 
               const MCObsInfo& mcobs, int bin, double& result)
{
 try{
    const list<int>& configs=BM.getBinConfigs(bin);
    result=0.0;
    double temp;
    for (list<int>::const_iterator it=configs.begin();it!=configs.end();++it){
        if (!CC.getCorrect(mcobs,*it,temp)) return false;
        result+=temp;}
    result/=double(configs.size());}
 catch(const std::exception& xp){
    return false;}
 return true;
}

void do_an_mcobs(const MCObsInfo& obstest, MCObsHandler& MC, SamplingMode mode,
                 uint sampling_loop)
{
 cout <<endl<<endl<<"DO AN MCOBS"<<endl<<endl;
 cout << "Observable: "<<obstest.output()<<endl;
 MC.setSamplingMode(mode);
 if (mode==Jackknife) cout << "mode is Jackknife"<<endl;
 else cout << "mode is Bootstrap"<<endl;
 try{
 cout <<"Query all samplings = "<<MC.queryFullAndSamplings(obstest)<<endl;
 cout <<"Full value = "<< MC.getFullSampleValue(obstest)<<endl;
 MC.begin(); cout << "loop = "<<sampling_loop<<endl;
 for (uint k=0;k<sampling_loop;k++) ++MC;
 cout <<"Current value = "<< MC.getCurrentSamplingValue(obstest)<<endl;}
 catch(std::exception& xp){
   cout << "Exception caught"<<endl;}
 cout << "size of bins in memory = "<< MC.getBinsInMemorySize() <<endl;
 cout << "size of jackknife samplings in memory = "<< MC.getJackknifeSamplingsInMemorySize()<<endl;
 cout << "size of bootstrap samplings in memory = "<< MC.getBootstrapSamplingsInMemorySize()<<endl;
 MC.setToJackknifeMode(); cout << "Jackknife mode"<<endl;
 cout << "size of current samplings in memory = "<< MC.getCurrentModeSamplingsInMemorySize()<<endl;
 MC.setToBootstrapMode(); cout << "Bootstrap mode"<<endl;
 cout << "size of current samplings in memory = "<< MC.getCurrentModeSamplingsInMemorySize()<<endl;
 cout << "size of jackknife samplings in memory = "<< MC.getSamplingsInMemorySize(Jackknife)<<endl;
 cout << "size of bootstrap samplings in memory = "<< MC.getSamplingsInMemorySize(Bootstrap)<<endl;
}

void check_sizes(MCObsHandler& MC, SamplingMode mode)
{
 cout <<endl<<endl<<"CHECK SIZES"<<endl<<endl;
 MC.setSamplingMode(mode);
 if (mode==Jackknife) cout << "mode is Jackknife"<<endl;
 else cout << "mode is Bootstrap"<<endl;
 cout << "size of bins in memory = "<< MC.getBinsInMemorySize() <<endl;
 cout << "size of jackknife samplings in memory = "<< MC.getJackknifeSamplingsInMemorySize()<<endl;
 cout << "size of bootstrap samplings in memory = "<< MC.getBootstrapSamplingsInMemorySize()<<endl;
 MC.setToJackknifeMode(); cout << "Jackknife mode"<<endl;
 cout << "size of current samplings in memory = "<< MC.getCurrentModeSamplingsInMemorySize()<<endl;
 MC.setToBootstrapMode(); cout << "Bootstrap mode"<<endl;
 cout << "size of current samplings in memory = "<< MC.getCurrentModeSamplingsInMemorySize()<<endl;
 cout << "size of jackknife samplings in memory = "<< MC.getSamplingsInMemorySize(Jackknife)<<endl;
 cout << "size of bootstrap samplings in memory = "<< MC.getSamplingsInMemorySize(Bootstrap)<<endl;
}

  // *************************************************************************

     //  This tests basic getting functions

void testMCObsHandler(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestMCObsHandler")==0)
 return;

 cout << endl << "Starting TestMCObsHandler"<<endl;
 

 try{
 XMLHandler xmlr(xml_in,"TestMCObsHandler");

   //   get the omissions and rebin factor

 set<int> omit;
 int rebin=1;
 XMLHandler xmlc(xmlr,"ConfigManip");
 xmlreadchild(xmlc,"Rebin",rebin);
 list<XMLHandler> xmlo=xmlc.find("Omit");
 for (list<XMLHandler>::iterator tt=xmlo.begin();tt!=xmlo.end();tt++){
    int anomission;
    xmlread(*tt,"Omit",anomission,"Tester");
    omit.insert(anomission);}

 VEVCorrect VC;
 bool hermitian=true;
 CorrCorrect CC;
 if (hermitian) CC.setHermitian();

    // make the fake data files with possibly missing data

 list<XMLHandler> xmlf=xmlr.find("MakeFakeCorrelatorFile");
 for (list<XMLHandler>::iterator it=xmlf.begin();it!=xmlf.end();++it){
    MCEnsembleInfo mcens(*it);  
    string fname;
    xmlreadchild(*it,"FileName",fname);  
    CorrelatorInfo corr(*it);
    int tmin,tmax;
    double Aval,Bval,Cval;
    xmlreadchild(*it,"MinTime",tmin);
    xmlreadchild(*it,"MaxTime",tmax);
    xmlreadchild(*it,"Aval",Aval);
    xmlreadchild(*it,"Bval",Bval);
    xmlreadchild(*it,"Cval",Cval);
    set<int> missing;
    list<XMLHandler> xmlo=it->find("Omit");
    for (list<XMLHandler>::iterator tt=xmlo.begin();tt!=xmlo.end();tt++){
       int anomission;
       xmlread(*tt,"Omit",anomission,"Tester");
       missing.insert(anomission);}
    cout << it->output()<<endl;
    make_fake_correlator_file(mcens,fname,corr,tmin,tmax,
                              Aval,Bval,Cval,missing,CC);}


    // make the fake vev files with possibly missing data

 list<XMLHandler> xmlv=xmlr.find("MakeFakeVEVFile");
 for (list<XMLHandler>::iterator it=xmlv.begin();it!=xmlv.end();++it){
    MCEnsembleInfo mcens(*it);  
    string fname;
    xmlreadchild(*it,"FileName",fname);  
    OperatorInfo opinfo(*it);
    double Aval,Cval;
    xmlreadchild(*it,"Aval",Aval);
    xmlreadchild(*it,"Cval",Cval);
    set<int> missing;
    list<XMLHandler> xmlo=it->find("Omit");
    for (list<XMLHandler>::iterator tt=xmlo.begin();tt!=xmlo.end();tt++){
       int anomission;
       xmlread(*tt,"Omit",anomission,"Tester");
       missing.insert(anomission);}
    cout << it->output()<<endl;
    make_fake_vev_file(mcens,fname,opinfo,Aval,Cval,missing,VC);}

 cout <<endl<<endl<<"Starting Read Tests:"<<endl<<endl;
 XMLHandler xmlrdtest(xmlr,"DoReadTests");
 MCBinsInfo bininfo(xmlrdtest);
 cout << bininfo.output()<<endl;
 MCSamplingInfo sampinfo(xmlrdtest);
 cout << sampinfo.output()<<endl;

 cout << "Firing up the MCObsGetHandler..."<<endl;
 MCObsGetHandler MCOH(xmlrdtest,bininfo,sampinfo); 

 cout << endl<<endl;
 cout << "Number of measurements in ensemble = "<<MCOH.getNumberOfMeasurements()<<endl;
 cout << "Ensemble ID = "<<MCOH.getEnsembleId()<<endl;
 cout << "FileMap:"<<endl;
 XMLHandler xmlout;
 MCOH.getFileMap(xmlout);
 cout << endl<<endl<<"***********************"<<endl<<endl;
 cout << xmlout.output()<<endl;
 cout << endl<<endl<<"***********************"<<endl<<endl;


 set<OperatorInfo> vevinfos=MCOH.getVEVInfos();
 for (set<OperatorInfo>::iterator it=vevinfos.begin();it!=vevinfos.end();++it)
    cout << endl<<"VEVINFO:"<<endl<<it->output()<<endl;

 set<CorrelatorInfo> corrinfos=MCOH.getCorrelatorInfos();
 for (set<CorrelatorInfo>::iterator it=corrinfos.begin();it!=corrinfos.end();++it)
    cout << endl<<"CORRINFO:"<<endl<<it->output()<<endl;

 XMLHandler xmlopl(xmlrdtest,"OpList");
 list<XMLHandler> xmlop=xmlopl.find("OperatorString");
 list<OperatorInfo> oplist;
 for (list<XMLHandler>::iterator it=xmlop.begin();it!=xmlop.end();++it)
    oplist.push_back(OperatorInfo(*it));
 cout << endl<<endl<<"OpList"<<endl;
 for (list<OperatorInfo>::iterator it=oplist.begin();it!=oplist.end();++it)
    cout << "operator: "<<it->output()<<endl;


 cout << endl<<endl<<"NOW firing up MCObsHandler!!"<<endl<<endl;

 MCObsHandler MC(MCOH,true);
 uint nconfig=MC.getNumberOfMeasurements();
 uint nbins=MC.getNumberOfBins();
 BinMapper BM(nconfig,rebin,omit);

 cout << "LatticeTExtent = "<<MC.getLatticeTimeExtent()<<endl;
 cout << "LatticeXExtent = "<<MC.getLatticeXExtent()<<endl;
 cout << "LatticeYExtent = "<<MC.getLatticeYExtent()<<endl;
 cout << "LatticeZExtent = "<<MC.getLatticeZExtent()<<endl;
 cout << "number of measurements = "<<MC.getNumberOfMeasurements()<<endl;
 cout << "number of bins = "<<MC.getNumberOfBins()<<endl;
 cout << "number of bootstrap resamplings = "<<MC.getNumberOfBootstrapResamplings()<<endl;
 cout << "rebinning factor = "<<MC.getRebinFactor()<<endl;
 const set<uint>& theomissions=MC.getOmissions();
 cout << "Omissions: ";
 for (set<uint>::const_iterator it=theomissions.begin();it!=theomissions.end();++it)
    cout << " "<<*it;
 cout << endl<<endl;
 SamplingMode mode=MC.getDefaultSamplingMode();
 if (mode==Jackknife) cout << "Default sampling mode is Jackknife"<<endl;
 else if (mode==Bootstrap) cout << "Default sampling mode is Bootstrap"<<endl;
 cout << "Sampling Info: "<<MC.getSamplingInfo().output()<<endl;
 const Bootstrapper& bptr=MC.getBootstrapper();
 cout << "Bootstrap is precompute mode? "<<bptr.isPrecomputeMode()<<endl;


 cout << "SOME FIRST TIME TESTING"<<endl<<endl;

 XMLHandler xmlfiles; MC.getFileMap(xmlfiles);
 cout << "Filemap from MCObsHandler: "<<xmlfiles.output()<<endl;
 check_sizes(MC,Jackknife);

 MCObsInfo obstest(oplist.front());
 do_an_mcobs(obstest,MC,Jackknife,5);
 do_an_mcobs(obstest,MC,Jackknife,16);
 do_an_mcobs(obstest,MC,Bootstrap,32);
 obstest.setToImaginaryPart();
 do_an_mcobs(obstest,MC,Jackknife,3);

 list<OperatorInfo>::const_iterator it=oplist.begin();
 OperatorInfo opA(*it); it++; OperatorInfo opB(*it);
 CorrelatorInfo cftest(opA,opB);
 CorrelatorAtTimeInfo ctest(cftest,6,true);
 obstest=MCObsInfo(ctest,RealPart);
 do_an_mcobs(obstest,MC,Jackknife,5);
 do_an_mcobs(obstest,MC,Jackknife,16);
 do_an_mcobs(obstest,MC,Bootstrap,32);
 obstest.setToImaginaryPart();
 do_an_mcobs(obstest,MC,Jackknife,3);

 it++; OperatorInfo opC(*it);
 CorrelatorInfo cftest2(opA,opC);
 CorrelatorAtTimeInfo ctest2(cftest2,4,true,true);
 obstest=MCObsInfo(ctest2,RealPart);
 do_an_mcobs(obstest,MC,Jackknife,5);
 do_an_mcobs(obstest,MC,Jackknife,16);
 do_an_mcobs(obstest,MC,Bootstrap,32);
 obstest.setToImaginaryPart();
 do_an_mcobs(obstest,MC,Jackknife,3);

 MC.clearSamplings();
 check_sizes(MC,Jackknife);
 cout << endl<<endl;



#ifdef COMPLEXNUMBERS
    int kkmax=2;
#else
    int kkmax=1;
#endif

 int goodq=0, badq=0;
 int goodc=0, badc=0;
 int goode=0, bade=0;
 cout << "Query VEV bins"<<endl;
 cout << "Number of VEV operators = "<<oplist.size()<<endl;

 for (list<OperatorInfo>::iterator it=oplist.begin();it!=oplist.end();++it){
    MCObsInfo mcobs(*it);
    for (int kk=0;kk<kkmax;++kk){
    if (kk>0) mcobs.setToImaginaryPart();
    //cout <<endl<< "VEV query test for "<<mcobs.output()<<endl<<endl;
    vector<double> vevbins;
    bool binflag=true;
    for (uint bin=0;bin<nbins;bin++){
       double val;
       if (getVEVBin(VC,BM,mcobs,bin,val)) vevbins.push_back(val);
       else break;}
    if (vevbins.size()!=nbins){ vevbins.clear(); binflag=false;}
    bool qflag=MC.queryBins(mcobs);
    if (qflag!=binflag){ cout << "VEV query bin failed"<<endl; badq++;}
    else goodq++;
    if (qflag){
       if (!binflag) badc+=nbins;
       else{
          const Vector<double>& res=MC.getBins(mcobs);
          for (uint bin=0;bin<nbins;++bin){
             double crt=vevbins[bin];  
          //cout << crt<<" "<<res[bin]<<" "<<MC.getBin(mcobs,bin)<<endl;
          if (!compare_floats(res[bin],crt,1e-5)) badc++; else goodc++;
          if (!compare_floats(MC.getBin(mcobs,bin),crt,1e-5)) badc++; else goodc++;}}}
    else{
       try{
          const Vector<double>& res=MC.getBins(mcobs); cout << res[0]<<endl;
          bade+=nbins;}
       catch(const std::exception& xp){
          //cout << "good exception caught"<<endl; 
          goode+=nbins;}
       for (uint bin=0;bin<nbins;++bin){
          try{
             double temp=MC.getBin(mcobs,bin); cout << temp<<endl;
             bade++;}
          catch(const std::exception& xp){ // cout << "good exception caught"<<endl; 
             goode++;}}}}}
 cout << "Number of successful vev bin queries = "<<goodq<<endl;
 cout << "Number of successful vev comparisons = "<<goodc<<endl;
 cout << "Number of correct vev exceptions caught = "<<goode<<endl;
 cout << "Number of incorrect vev bin queries = "<<badq<<endl;
 cout << "Number of incorrect vev comparisons = "<<badc<<endl;
 cout << "Number of incorrect vev exceptions caught = "<<bade<<endl;


 goodq=0; badq=0;
 goodc=0; badc=0;
 goode=0; bade=0;
 cout << endl<<"Query correlator bins"<<endl;
 for (list<OperatorInfo>::iterator it1=oplist.begin();it1!=oplist.end();++it1){
 for (list<OperatorInfo>::iterator it2=oplist.begin();it2!=oplist.end();++it2){
    CorrelatorInfo cf(*it1,*it2);
    for (int tind=0;tind<=32;tind++){
       CorrelatorAtTimeInfo ct(cf,tind,hermitian);
       MCObsInfo mcobs(ct);
#ifdef COMPLEXNUMBERS
       int kkmax=2;
#else
       int kkmax=1;
#endif
       for (int kk=0;kk<kkmax;++kk){
       if (kk>0) mcobs.setToImaginaryPart();
       vector<double> corbins;
       bool binflag=true;
       for (uint bin=0;bin<nbins;bin++){
          double val;
          if (getCorBin(CC,BM,mcobs,bin,val)) corbins.push_back(val);
          else break;}
       if (corbins.size()!=nbins){ corbins.clear(); binflag=false;}
       bool qflag=MC.queryBins(mcobs);
       if (qflag!=binflag){ cout << "corr query bin failed"<<endl; badq++;}
       else goodq++;
       if (qflag){
          if (!binflag) badc+=nbins;
          else{
             const Vector<double>& res=MC.getBins(mcobs);
             for (uint bin=0;bin<nbins;++bin){
                double crt=corbins[bin];  
               // cout << crt<<" "<<res[bin]<<" "<<MC.getBin(mcobs,bin)<<endl;
                if (!compare_floats(res[bin],crt,1e-5)){ badc++; cout << "diff = "<<res[bin]-crt<<endl;} else goodc++;
                if (!compare_floats(MC.getBin(mcobs,bin),crt,1e-5)) badc++; else goodc++;}}}
       else{
          try{
             const Vector<double>& res=MC.getBins(mcobs); cout << res[0]<<endl;
             bade+=nbins;}
          catch(const std::exception& xp){
             //cout << "good exception caught"<<endl; 
             goode+=nbins;}
          for (uint bin=0;bin<nbins;++bin){
             try{
                double temp=MC.getBin(mcobs,bin); cout << temp<<endl;
                bade++;}
             catch(const std::exception& xp){ // cout << "good exception caught"<<endl; 
                goode++;}}}}}}}
 cout << "Number of successful corr bin queries = "<<goodq<<endl;
 cout << "Number of successful corr comparisons = "<<goodc<<endl;
 cout << "Number of correct corr exceptions caught = "<<goode<<endl;
 cout << "Number of incorrect corr bin queries = "<<badq<<endl;
 cout << "Number of incorrect corr comparisons = "<<badc<<endl;
 cout << "Number of incorrect corr exceptions caught = "<<bade<<endl;


 Vector<double> newvalues(nbins);
 for (uint bin=0;bin<nbins;bin++)
    newvalues[bin]=1.4*bin-0.43;
 MCObsInfo obskey(*(corrinfos.begin()),5);
 MC.putBins(obskey,newvalues);
 cout << endl<<"Did put bins"<<endl;

 MC.setToJackknifeMode(); 
 double fullmean=MC.getFullSampleValue(obskey);
 cout << "fullmean = "<<fullmean<<endl;
 int count=0;
 for (MC.begin();!MC.end();++MC){
    cout << MC.getCurrentSamplingValue(obskey)<<endl; count++;}
 cout << "count = "<<count<<endl;


 cout <<endl<<endl;

// MC.setBootstrapper(888,32,54,false);
 {const Bootstrapper& bt=MC.getBootstrapper();
 cout << "Number of bootstrap objects = "<<bt.getNumberOfObjects()<<endl;
 cout << "Number of bootstrap resamplings = "<<bt.getNumberOfResamplings()<<endl;
 cout << "Bootstrap seed = "<<bt.getRNGSeed()<<endl;
 cout << "Bootstrap skip value = "<<bt.getSkipValue()<<endl;
 cout << "Bootstrap current resampling = "<<bt.getCurrentResamplingCount()<<endl;
 cout << "Boostrap precompute? "<<bt.isPrecomputeMode()<<endl;}

 MC.setToBootstrapMode(); 
 count=0;
 for (MC.begin();!MC.end();++MC){
    cout << MC.getCurrentSamplingValue(obskey)<<endl; count++;}
 cout << "count = "<<count<<endl;


 MCObsInfo obskeyA("MapleSyrup",7,true,RealPart);
 count=0;
 cout << "Doing putCurrentSamplingValue"<<endl;
 for (MC.begin();!MC.end();++MC){
    cout << "count = "<<count;
    MC.putCurrentSamplingValue(obskeyA,4.0,true);
    MC.putCurrentSamplingValue(obskeyA,2.0*double(count++),true);
    cout << "  value = "<<MC.getCurrentSamplingValue(obskeyA); 
    cout << "  query all = "<<MC.queryFullAndSamplings(obskeyA)<<endl;}


/*
 cout << endl<<"Rebin to 3"<<endl;
 MC.setRebin(3);
 MC.clearOmissions();
 {const Bootstrapper& bt=MC.getBootstrapper();
 cout << "Number of bootstrap objects = "<<bt.getNumberOfObjects()<<endl;
 cout << "Number of bootstrap resamplings = "<<bt.getNumberOfResamplings()<<endl;
 cout << "Bootstrap seed = "<<bt.getRNGSeed()<<endl;
 cout << "Bootstrap skip value = "<<bt.getSkipValue()<<endl;
 cout << "Bootstrap current resampling = "<<bt.getCurrentResamplingCount()<<endl;
 cout << "Boostrap precompute? "<<bt.isPrecomputeMode()<<endl;}

 MC.setRebin(1);
 MC.addOmission(54);
 MC.addOmission(524);
 MC.addOmission(354);
 MC.addOmission(222);
 MC.addOmission(584);

 {const Bootstrapper& bt=MC.getBootstrapper();
 cout << "Number of bootstrap objects = "<<bt.getNumberOfObjects()<<endl;
 cout << "Number of bootstrap resamplings = "<<bt.getNumberOfResamplings()<<endl;
 cout << "Bootstrap seed = "<<bt.getRNGSeed()<<endl;
 cout << "Bootstrap skip value = "<<bt.getSkipValue()<<endl;
 cout << "Bootstrap current resampling = "<<bt.getCurrentResamplingCount()<<endl;
 cout << "Boostrap precompute? "<<bt.isPrecomputeMode()<<endl;} 
*/
 }
 catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}

}


// ***********************************************************************


double jack_covariances(const Vector<double> ajack, const Vector<double> bjack)
{
 uint njack=ajack.size();
 double avga=0.0, avgb=0.0;
 for (uint k=0;k<njack;k++){
    avga+=ajack[k];
    avgb+=bjack[k];}
 avga/=njack;
 avgb/=njack;
 double cov=0.0;
 for (uint k=0;k<njack;k++)
    cov+=(ajack[k]-avga)*(bjack[k]-avgb);
 return (1.0-1.0/njack)*cov;
}


double calc_simple_covariance(const RVector& dvals1,
                              const RVector& dvals2)
{
 unsigned int n=dvals1.size();
 double r1=1.0/double(n);
 double r2=1.0/(double(n-1));
 double dm1=0.0, dm2=0.0, dm12=0.0;
 for (unsigned int k=0;k<n;k++){
    dm1+=dvals1[k];
    dm2+=dvals2[k];
    dm12+=dvals1[k]*dvals2[k];}
 dm1*=r1; dm2*=r1; dm12*=r1;
 return (dm12-dm1*dm2)*r2;
}



     //  This tests jackknifing, bootstrapping on simple and vev-subtracted
     //  data

void testMCObsHandler2(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestMCObsHandler2")==0)
 return;

 cout << endl << "Starting TestMCObsHandler2"<<endl;
 

 try{
 XMLHandler xmlr(xml_in,"TestMCObsHandler2");

 XMLHandler xmlrdtest(xmlr,"DoReadTests");
 MCBinsInfo bininfo(xmlrdtest);
 cout << bininfo.output()<<endl;
 MCSamplingInfo sampinfo(xmlrdtest);
 cout << sampinfo.output()<<endl;
 MCObsGetHandler MCOH(xmlrdtest,bininfo,sampinfo); 

 MCObsHandler MC(MCOH);
 cout << endl<<endl;
 cout << "Number of measurements in ensemble = "<<MC.getNumberOfMeasurements()<<endl;
 cout << "Ensemble ID = "<<MCOH.getEnsembleId()<<endl;
 cout << "Number of bins = "<<MC.getNumberOfBins()<<endl;
 uint nbins = MC.getNumberOfBins();

    //   make some fake data

 Vector<double> bins1(nbins);
 Vector<double> bins2(nbins);
 double mean1=23.546;
 double spread1=4.536;
 double mean2=32.812;
 double spread2=6.783;
 for (uint k=0;k<nbins;++k){
    double r=double((rand()%200000)-100000)/double(100000);
    double r2=double((rand()%200000)-100000)/double(6000000);
    bins1[k]=mean1+spread1*r;
    bins2[k]=mean2+spread2*(r+r2);}
 MCObsInfo obskey1("Cat",0,true,RealPart);
 MCObsInfo obskey2("Dog",0,true,RealPart);
 MC.putBins(obskey1,bins1);
 MC.putBins(obskey2,bins2);

 cout.precision(15);
 uint success=0, fail=0;
 double avg1=0.0;
 for (uint i=0;i<nbins;++i) avg1+=bins1[i];
 avg1/=nbins;
 cout <<obskey1.output()<<endl
      <<"Full sample mean = "<<MC.getFullSampleValue(obskey1)<<" should be "<<avg1<<endl;
 double avg2=0.0;
 for (uint i=0;i<nbins;++i) avg2+=bins2[i];
 avg2/=nbins;
 cout <<obskey2.output()<<endl
      <<"Full sample mean = "<<MC.getFullSampleValue(obskey2)<<" should be "<<avg2<<endl;

 double eps=1e-8;
 if (compare_floats(MC.getFullSampleValue(obskey1),avg1,eps)) success++; 
 else fail++; 
 if (compare_floats(MC.getFullSampleValue(obskey2),avg2,eps)) success++; 
 else fail++; 

 double cov11=calc_simple_covariance(bins1,bins1);
 double cov12=calc_simple_covariance(bins1,bins2);
 double cov22=calc_simple_covariance(bins2,bins2);

 cout << "1,1, covariance = "<<MC.getCovariance(obskey1,obskey1,Jackknife)<<" should be "<<cov11<<endl;
 cout << "1,2, covariance = "<<MC.getCovariance(obskey1,obskey2,Jackknife)<<" should be "<<cov12<<endl;
 cout << "2,1, covariance = "<<MC.getCovariance(obskey2,obskey1,Jackknife)<<" should be "<<cov12<<endl;
 cout << "2,2, covariance = "<<MC.getCovariance(obskey2,obskey2,Jackknife)<<" should be "<<cov22<<endl;
 if (compare_floats(MC.getCovariance(obskey1,obskey1,Jackknife),cov11,eps)) success++; else fail++; 
 if (compare_floats(MC.getCovariance(obskey1,obskey2,Jackknife),cov12,eps)) success++; else fail++; 
 if (compare_floats(MC.getCovariance(obskey2,obskey1,Jackknife),cov12,eps)) success++; else fail++; 
 if (compare_floats(MC.getCovariance(obskey2,obskey2,Jackknife),cov22,eps)) success++; else fail++; 

 MC.setToJackknifeMode();

 int jackcount=-1;
 for (MC.begin();!MC.end();++MC,jackcount++){

 uint nb=nbins-((jackcount>=0)?1:0); 
 avg1=0.0;
 for (uint i=0;i<nbins;++i) if (int(i)!=jackcount) avg1+=bins1[i];
 avg1/=nb;
 avg2=0.0;
 for (uint i=0;i<nbins;++i) if (int(i)!=jackcount) avg2+=bins2[i];
 avg2/=nb;

 if (compare_floats(MC.getCurrentSamplingValue(obskey1),avg1,eps)) success++; 
 else fail++; 
 if (compare_floats(MC.getCurrentSamplingValue(obskey2),avg2,eps)) success++; 
 else fail++; 
/*
 cov11=0.0; cov12=0.0; cov22=0.0;
 for (uint i=0;i<nbins;++i){ 
    if (int(i)!=jackcount) {
    cov11+=(bins1[i]-avg1)*(bins1[i]-avg1);
    cov12+=(bins1[i]-avg1)*(bins2[i]-avg2);
    cov22+=(bins2[i]-avg2)*(bins2[i]-avg2);}}
 cov11/=(nb*(nb-1));
 cov12/=(nb*(nb-1));
 cov22/=(nb*(nb-1));

 diff=std::abs(MC.getCurrentSamplingCovariance(obskey1,obskey1)-cov11);
 if (diff<1e-8)  success++; else fail++; 
 diff=std::abs(MC.getCurrentSamplingCovariance(obskey1,obskey2)-cov12);
 if (diff<1e-8)  success++; else fail++; 
 diff=std::abs(MC.getCurrentSamplingCovariance(obskey2,obskey1)-cov12);
 if (diff<1e-8)  success++; else fail++; 
 diff=std::abs(MC.getCurrentSamplingCovariance(obskey2,obskey2)-cov22);
 if (diff<1e-8)  success++; else fail++; 
*/
 }

// MC.setBootstrapper(84,0,24,true);
 MC.setToBootstrapMode();

 int count=0;
 const Bootstrapper& Bref =MC.getBootstrapper();
 Bootstrapper BB(Bref.getNumberOfObjects(),Bref.getNumberOfResamplings(),
                 Bref.getRNGSeed(),Bref.getSkipValue());
 for (MC.begin();!MC.end();++MC,++count){

 Vector<unsigned int> indexmapper;
 if (count>0){
    indexmapper=BB.getResampling(count-1);}
 else{
    indexmapper.resize(nbins);
    for (uint k=0;k<nbins;++k) indexmapper[k]=k;}

 avg1=0.0;
 for (uint i=0;i<nbins;++i) avg1+=bins1[indexmapper[i]];
 avg1/=nbins;
 avg2=0.0;
 for (uint i=0;i<nbins;++i) avg2+=bins2[indexmapper[i]];
 avg2/=nbins;
 if (compare_floats(MC.getCurrentSamplingValue(obskey1),avg1,eps)) success++; 
 else fail++; 
 if (compare_floats(MC.getCurrentSamplingValue(obskey2),avg2,eps)) success++; 
 else fail++; 
/*
 cov11=0.0; cov12=0.0; cov22=0.0;
 for (uint i=0;i<nbins;++i){
    uint ii=indexmapper[i];
    cov11+=(bins1[ii]-avg1)*(bins1[ii]-avg1);
    cov12+=(bins1[ii]-avg1)*(bins2[ii]-avg2);
    cov22+=(bins2[ii]-avg2)*(bins2[ii]-avg2);}
 cov11/=(nbins*(nbins-1));
 cov12/=(nbins*(nbins-1));
 cov22/=(nbins*(nbins-1));

 diff=std::abs(MC.getCurrentSamplingCovariance(obskey1,obskey1)-cov11);
 if (diff<1e-8)  success++; else fail++; 
 diff=std::abs(MC.getCurrentSamplingCovariance(obskey1,obskey2)-cov12);
 if (diff<1e-8)  success++; else fail++; 
 diff=std::abs(MC.getCurrentSamplingCovariance(obskey2,obskey1)-cov12);
 if (diff<1e-8)  success++; else fail++; 
 diff=std::abs(MC.getCurrentSamplingCovariance(obskey2,obskey2)-cov22);
 if (diff<1e-8)  success++; else fail++; 
*/
 }


     //   now with vev subtractions


    //   make some fake data

 Vector<double> corr1rebins(nbins);
 Vector<double> corr1imbins(nbins);
 Vector<double> vevsrc1rebins(nbins);
 Vector<double> vevsrc1imbins(nbins);
 Vector<double> vevsnk1rebins(nbins);
 Vector<double> vevsnk1imbins(nbins);

 Vector<double> corr2rebins(nbins);
 Vector<double> corr2imbins(nbins);
 Vector<double> vevsrc2rebins(nbins);
 Vector<double> vevsrc2imbins(nbins);
 Vector<double> vevsnk2rebins(nbins);
 Vector<double> vevsnk2imbins(nbins);

 double vevsrc1re=45.453;
 double vevsrc1im=-2.34;
 double vevsnk1re=22.173;
 double vevsnk1im=5.234;
 double corr1re=1048.56;
 double corr1im=295.21;

 double vevsrc2re=18.43;
 double vevsrc2im=-0.45;
 double vevsnk2re=8.43;
 double vevsnk2im=2.22;
 double corr2re=210.56;
 double corr2im=54.21;

 spread1=2.254;
 spread2=6.333;
 for (uint k=0;k<nbins;++k){
    double r1=double((rand()%200000)-100000)/double(100000);
    double r2=double((rand()%200000)-100000)/double(100000);
    double rr1=double((rand()%200000)-100000)/double(6000000);
    double rr2=double((rand()%200000)-100000)/double(6000000);
    double rr3=double((rand()%200000)-100000)/double(6000000);
    double rr4=double((rand()%200000)-100000)/double(6000000);

    corr1rebins[k]=corr1re+spread1*r1;
    corr1imbins[k]=corr1im+spread1*r2;
    vevsrc1rebins[k]=vevsrc1re+spread1*(r1+rr1);
    vevsrc1imbins[k]=vevsrc1im+spread1*(r2+rr2);
    vevsnk1rebins[k]=vevsnk1re+spread1*(r1+rr3);
    vevsnk1imbins[k]=vevsnk1im+spread1*(r2+rr4);

    corr2rebins[k]=corr2re+spread1*r1;
    corr2imbins[k]=corr2im+spread1*r2;
    vevsrc2rebins[k]=vevsrc2re+spread2*(r1+rr1);
    vevsrc2imbins[k]=vevsrc2im+spread2*(r2+rr2);
    vevsnk2rebins[k]=vevsnk2re+spread2*(r1+rr3);
    vevsnk2imbins[k]=vevsnk2im+spread2*(r2+rr4); }

#ifndef COMPLEXNUMBERS
 for (uint k=0;k<nbins;++k){
    corr1imbins[k]=0.0;
    vevsrc1imbins[k]=0.0;
    vevsnk1imbins[k]=0.0;

    corr2imbins[k]=0.0;
    vevsrc2imbins[k]=0.0;
    vevsnk2imbins[k]=0.0; }
#endif

 OperatorInfo op1("glueball P=(0,0,0) A1gp_1 TrEig");
 OperatorInfo op2("eta P=(0,0,0) A1gp_1 SD_2");
 OperatorInfo op3("phi P=(0,0,0) A1gp_1 TDO_5");
 OperatorInfo op4("isosinglet_pion_pion A1gp_1 [P=(0,0,0) A1um SS_0] [P=(0,0,0) A1um SS_0]");

 MCObsInfo corvevkeyre1(op1,op2,3,false,RealPart,true);
 MCObsInfo corvevkeyim1(op1,op2,3,false,ImaginaryPart,true);
 MCObsInfo corvevkeyre2(op3,op4,3,false,RealPart,true);
 MCObsInfo corvevkeyim2(op3,op4,3,false,ImaginaryPart,true);

 MCObsInfo corkeyre1(op1,op2,3,false,RealPart,false);
 MCObsInfo corkeyim1(op1,op2,3,false,ImaginaryPart,false);
 MCObsInfo corkeyre2(op3,op4,3,false,RealPart,false);
 MCObsInfo corkeyim2(op3,op4,3,false,ImaginaryPart,false);

 MCObsInfo vevsnkkeyre1(op1,RealPart);
 MCObsInfo vevsnkkeyim1(op1,ImaginaryPart);
 MCObsInfo vevsrckeyre1(op2,RealPart);
 MCObsInfo vevsrckeyim1(op2,ImaginaryPart);

 MCObsInfo vevsnkkeyre2(op3,RealPart);
 MCObsInfo vevsnkkeyim2(op3,ImaginaryPart);
 MCObsInfo vevsrckeyre2(op4,RealPart);
 MCObsInfo vevsrckeyim2(op4,ImaginaryPart);

 MC.putBins(corkeyre1,corr1rebins);
 MC.putBins(corkeyim1,corr1imbins);
 MC.putBins(vevsrckeyre1,vevsrc1rebins);
 MC.putBins(vevsrckeyim1,vevsrc1imbins);
 MC.putBins(vevsnkkeyre1,vevsnk1rebins);
 MC.putBins(vevsnkkeyim1,vevsnk1imbins);

 MC.putBins(corkeyre2,corr2rebins);
 MC.putBins(corkeyim2,corr2imbins);
 MC.putBins(vevsrckeyre2,vevsrc2rebins);
 MC.putBins(vevsrckeyim2,vevsrc2imbins);
 MC.putBins(vevsnkkeyre2,vevsnk2rebins);
 MC.putBins(vevsnkkeyim2,vevsnk2imbins);


 vevsrc1re=0.0;
 vevsrc1im=0.0;
 vevsnk1re=0.0;
 vevsnk1im=0.0;
 corr1re=0.0;
 corr1im=0.0;
 vevsrc2re=0.0;
 vevsrc2im=0.0;
 vevsnk2re=0.0;
 vevsnk2im=0.0;
 corr2re=0.0;
 corr2im=0.0;
 for (uint k=0;k<nbins;++k){
    vevsrc1re+=vevsrc1rebins[k];
    vevsrc1im+=vevsrc1imbins[k];
    vevsnk1re+=vevsnk1rebins[k];
    vevsnk1im+=vevsnk1imbins[k];
    corr1re+=corr1rebins[k];
    corr1im+=corr1imbins[k];
    vevsrc2re+=vevsrc2rebins[k];
    vevsrc2im+=vevsrc2imbins[k];
    vevsnk2re+=vevsnk2rebins[k];
    vevsnk2im+=vevsnk2imbins[k];
    corr2re+=corr2rebins[k];
    corr2im+=corr2imbins[k];}
 vevsrc1re/=nbins;
 vevsrc1im/=nbins;
 vevsnk1re/=nbins;
 vevsnk1im/=nbins;
 corr1re/=nbins;
 corr1im/=nbins;
 vevsrc2re/=nbins;
 vevsrc2im/=nbins;
 vevsnk2re/=nbins;
 vevsnk2im/=nbins;
 corr2re/=nbins;
 corr2im/=nbins;

 //cout<<endl<<endl;
 complex<double> corr1(corr1re,corr1im);
 complex<double> vevsrc1(vevsrc1re,vevsrc1im);
 complex<double> vevsnk1(vevsnk1re,vevsnk1im);
 corr1-=conj(vevsrc1)*vevsnk1;  //cout << corr1<<endl;

 complex<double> corr2(corr2re,corr2im);
 complex<double> vevsrc2(vevsrc2re,vevsrc2im);
 complex<double> vevsnk2(vevsnk2re,vevsnk2im);
 corr2-=conj(vevsrc2)*vevsnk2;  //cout << corr2<<endl;

 double tester=MC.getFullSampleValue(corvevkeyre1);
 tester=MC.getFullSampleValue(corvevkeyre1);
 if (compare_floats(tester,real(corr1),eps)) success++; 
 else fail++; 
 tester=MC.getFullSampleValue(corvevkeyre2);
 if (compare_floats(tester,real(corr2),eps)) success++; 
 else fail++; 
#ifdef COMPLEXNUMBERS
 tester=MC.getFullSampleValue(corvevkeyim1);
 if (compare_floats(tester,imag(corr1),eps)) success++; 
 else fail++; 
 tester=MC.getFullSampleValue(corvevkeyim2);
 if (compare_floats(tester,imag(corr2),eps)) success++; 
 else fail++; 
#endif 

 Vector<double> corr1rejack(nbins);
 Vector<double> corr1imjack(nbins);
 Vector<double> vevsrc1rejack(nbins);
 Vector<double> vevsrc1imjack(nbins);
 Vector<double> vevsnk1rejack(nbins);
 Vector<double> vevsnk1imjack(nbins);

 Vector<double> corr2rejack(nbins);
 Vector<double> corr2imjack(nbins);
 Vector<double> vevsrc2rejack(nbins);
 Vector<double> vevsrc2imjack(nbins);
 Vector<double> vevsnk2rejack(nbins);
 Vector<double> vevsnk2imjack(nbins);

 for (uint k=0;k<nbins;++k){
    vevsrc1rejack[k]=(nbins*vevsrc1re-vevsrc1rebins[k])/(nbins-1);
    vevsrc1imjack[k]=(nbins*vevsrc1im-vevsrc1imbins[k])/(nbins-1);
    vevsnk1rejack[k]=(nbins*vevsnk1re-vevsnk1rebins[k])/(nbins-1);
    vevsnk1imjack[k]=(nbins*vevsnk1im-vevsnk1imbins[k])/(nbins-1);
    corr1rejack[k]=(nbins*corr1re-corr1rebins[k])/(nbins-1);
    corr1imjack[k]=(nbins*corr1im-corr1imbins[k])/(nbins-1);
    vevsrc2rejack[k]=(nbins*vevsrc2re-vevsrc2rebins[k])/(nbins-1);
    vevsrc2imjack[k]=(nbins*vevsrc2im-vevsrc2imbins[k])/(nbins-1);
    vevsnk2rejack[k]=(nbins*vevsnk2re-vevsnk2rebins[k])/(nbins-1);
    vevsnk2imjack[k]=(nbins*vevsnk2im-vevsnk2imbins[k])/(nbins-1);
    corr2rejack[k]=(nbins*corr2re-corr2rebins[k])/(nbins-1);
    corr2imjack[k]=(nbins*corr2im-corr2imbins[k])/(nbins-1);}


 MC.setToJackknifeMode();

 Vector<double> corrvev1rejack(nbins);
 Vector<double> corrvev1imjack(nbins);
 Vector<double> corrvev2rejack(nbins);
 Vector<double> corrvev2imjack(nbins);


 for (uint k=0;k<nbins;++k){
 
    if (k==0) MC.begin();
    ++MC;

    //cout<<endl<<endl;
    corr1=complex<double>(corr1rejack[k],corr1imjack[k]);
    vevsrc1=complex<double>(vevsrc1rejack[k],vevsrc1imjack[k]);
    vevsnk1=complex<double>(vevsnk1rejack[k],vevsnk1imjack[k]);
    corr1-=conj(vevsrc1)*vevsnk1; 

    corr2=complex<double>(corr2rejack[k],corr2imjack[k]);
    vevsrc2=complex<double>(vevsrc2rejack[k],vevsrc2imjack[k]);
    vevsnk2=complex<double>(vevsnk2rejack[k],vevsnk2imjack[k]);
    corr2-=conj(vevsrc2)*vevsnk2; 

    tester=MC.getCurrentSamplingValue(corvevkeyre1);
    tester=MC.getCurrentSamplingValue(corvevkeyre1); 
    if (compare_floats(tester,real(corr1),eps)) success++; else fail++; 
    tester=MC.getCurrentSamplingValue(corvevkeyre2);
    if (compare_floats(tester,real(corr2),eps)) success++; else fail++; 
#ifdef COMPLEXNUMBERS
    tester=MC.getCurrentSamplingValue(corvevkeyim1);
    if (compare_floats(tester,imag(corr1),eps)) success++; else fail++; 
    tester=MC.getCurrentSamplingValue(corvevkeyim2);
    if (compare_floats(tester,imag(corr2),eps)) success++; else fail++; 
#endif
    }

// MC.setBootstrapper(128,0,16,false);
 MC.setToBootstrapMode();

 count=0;
 for (MC.begin();!MC.end();++MC,++count){

 Vector<unsigned int> indexmapper;
 if (count>0){
    indexmapper=BB.getResampling(count-1);}
 else{
    indexmapper.resize(nbins);
    for (uint k=0;k<nbins;++k) indexmapper[k]=k;}

 vevsrc1re=0.0;
 vevsrc1im=0.0;
 vevsnk1re=0.0;
 vevsnk1im=0.0;
 corr1re=0.0;
 corr1im=0.0;
 vevsrc2re=0.0;
 vevsrc2im=0.0;
 vevsnk2re=0.0;
 vevsnk2im=0.0;
 corr2re=0.0;
 corr2im=0.0;
 for (uint kk=0;kk<nbins;++kk){
    uint k=indexmapper[kk];
    vevsrc1re+=vevsrc1rebins[k];
    vevsrc1im+=vevsrc1imbins[k];
    vevsnk1re+=vevsnk1rebins[k];
    vevsnk1im+=vevsnk1imbins[k];
    corr1re+=corr1rebins[k];
    corr1im+=corr1imbins[k];
    vevsrc2re+=vevsrc2rebins[k];
    vevsrc2im+=vevsrc2imbins[k];
    vevsnk2re+=vevsnk2rebins[k];
    vevsnk2im+=vevsnk2imbins[k];
    corr2re+=corr2rebins[k];
    corr2im+=corr2imbins[k];}
 vevsrc1re/=nbins;
 vevsrc1im/=nbins;
 vevsnk1re/=nbins;
 vevsnk1im/=nbins;
 corr1re/=nbins;
 corr1im/=nbins;
 vevsrc2re/=nbins;
 vevsrc2im/=nbins;
 vevsnk2re/=nbins;
 vevsnk2im/=nbins;
 corr2re/=nbins;
 corr2im/=nbins;

 //cout<<endl<<endl;
 complex<double> corr1(corr1re,corr1im);
 complex<double> vevsrc1(vevsrc1re,vevsrc1im);
 complex<double> vevsnk1(vevsnk1re,vevsnk1im);
 corr1-=conj(vevsrc1)*vevsnk1;  //cout << corr1<<endl;

 complex<double> corr2(corr2re,corr2im);
 complex<double> vevsrc2(vevsrc2re,vevsrc2im);
 complex<double> vevsnk2(vevsnk2re,vevsnk2im);
 corr2-=conj(vevsrc2)*vevsnk2;  //cout << corr2<<endl;

 double tester=MC.getCurrentSamplingValue(corvevkeyre1);
 tester=MC.getCurrentSamplingValue(corvevkeyre1);
 if (compare_floats(tester,real(corr1),eps)) success++; else fail++; 
 tester=MC.getCurrentSamplingValue(corvevkeyre2);
 if (compare_floats(tester,real(corr2),eps)) success++; else fail++; 
#ifdef COMPLEXNUMBERS
 tester=MC.getCurrentSamplingValue(corvevkeyim1);
 if (compare_floats(tester,imag(corr1),eps)) success++; else fail++; 
 tester=MC.getCurrentSamplingValue(corvevkeyim2);
 if (compare_floats(tester,imag(corr2),eps)) success++; else fail++; 
#endif

 for (uint kk=0;kk<nbins;++kk){
    uint k=indexmapper[kk];
    vevsrc1rejack[kk]=(nbins*vevsrc1re-vevsrc1rebins[k])/(nbins-1);
    vevsrc1imjack[kk]=(nbins*vevsrc1im-vevsrc1imbins[k])/(nbins-1);
    vevsnk1rejack[kk]=(nbins*vevsnk1re-vevsnk1rebins[k])/(nbins-1);
    vevsnk1imjack[kk]=(nbins*vevsnk1im-vevsnk1imbins[k])/(nbins-1);
    corr1rejack[kk]=(nbins*corr1re-corr1rebins[k])/(nbins-1);
    corr1imjack[kk]=(nbins*corr1im-corr1imbins[k])/(nbins-1);
    vevsrc2rejack[kk]=(nbins*vevsrc2re-vevsrc2rebins[k])/(nbins-1);
    vevsrc2imjack[kk]=(nbins*vevsrc2im-vevsrc2imbins[k])/(nbins-1);
    vevsnk2rejack[kk]=(nbins*vevsnk2re-vevsnk2rebins[k])/(nbins-1);
    vevsnk2imjack[kk]=(nbins*vevsnk2im-vevsnk2imbins[k])/(nbins-1);
    corr2rejack[kk]=(nbins*corr2re-corr2rebins[k])/(nbins-1);
    corr2imjack[kk]=(nbins*corr2im-corr2imbins[k])/(nbins-1);

    corr1=complex<double>(corr1rejack[kk],corr1imjack[kk]);
    vevsrc1=complex<double>(vevsrc1rejack[kk],vevsrc1imjack[kk]);
    vevsnk1=complex<double>(vevsnk1rejack[kk],vevsnk1imjack[kk]);
    corr1-=conj(vevsrc1)*vevsnk1; 

    corr2=complex<double>(corr2rejack[kk],corr2imjack[kk]);
    vevsrc2=complex<double>(vevsrc2rejack[kk],vevsrc2imjack[kk]);
    vevsnk2=complex<double>(vevsnk2rejack[kk],vevsnk2imjack[kk]);
    corr2-=conj(vevsrc2)*vevsnk2; 

    corrvev1rejack[kk]=real(corr1);
    corrvev1imjack[kk]=imag(corr1);
    corrvev2rejack[kk]=real(corr2);
    corrvev2imjack[kk]=imag(corr2);}

   // calculate covariances
/*
 double correct;
 correct=jack_covariances(corrvev1rejack,corrvev1rejack); tester=MC.getCurrentSamplingCovariance(corvevkeyre1,corvevkeyre1); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev1rejack,corrvev2rejack); tester=MC.getCurrentSamplingCovariance(corvevkeyre1,corvevkeyre2); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev1rejack,corrvev1rejack); tester=MC.getCurrentSamplingCovariance(corvevkeyre1,corvevkeyre1); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev2rejack,corrvev2rejack); tester=MC.getCurrentSamplingCovariance(corvevkeyre2,corvevkeyre2); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
#ifdef COMPLEXNUMBERS
 correct=jack_covariances(corrvev1rejack,corrvev1imjack); tester=MC.getCurrentSamplingCovariance(corvevkeyre1,corvevkeyim1); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev1imjack,corrvev1rejack); tester=MC.getCurrentSamplingCovariance(corvevkeyim1,corvevkeyre1); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev1imjack,corrvev1imjack); tester=MC.getCurrentSamplingCovariance(corvevkeyim1,corvevkeyim1); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev1rejack,corrvev2imjack); tester=MC.getCurrentSamplingCovariance(corvevkeyre1,corvevkeyim2); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev1imjack,corrvev2rejack); tester=MC.getCurrentSamplingCovariance(corvevkeyim1,corvevkeyre2); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev1imjack,corrvev2imjack); tester=MC.getCurrentSamplingCovariance(corvevkeyim1,corvevkeyim2); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev1rejack,corrvev1imjack); tester=MC.getCurrentSamplingCovariance(corvevkeyre1,corvevkeyim1); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev1imjack,corrvev1rejack); tester=MC.getCurrentSamplingCovariance(corvevkeyim1,corvevkeyre1); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev1imjack,corrvev1imjack); tester=MC.getCurrentSamplingCovariance(corvevkeyim1,corvevkeyim1); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev2rejack,corrvev2imjack); tester=MC.getCurrentSamplingCovariance(corvevkeyre2,corvevkeyim2); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev2imjack,corrvev2rejack); tester=MC.getCurrentSamplingCovariance(corvevkeyim2,corvevkeyre2); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
 correct=jack_covariances(corrvev2imjack,corrvev2imjack); tester=MC.getCurrentSamplingCovariance(corvevkeyim2,corvevkeyim2); diff=std::abs(tester-correct); if (diff<1e-8)  success++; else fail++; 
#endif
cout << "done"<<endl; */
 }

 cout <<endl<<endl;
 cout << "Number of successful tests = "<<success<<endl;
 cout << "Number of failures = "<<fail<<endl;

 }
 catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}

}

string output_status(const map<MCObsInfo,string>& statusmap)
{
 ostringstream fout;
 for (map<MCObsInfo,string>::const_iterator st=statusmap.begin();st!=statusmap.end();st++){
    string obs(st->first.str());
    size_t pos=obs.find("<Info>");
    obs.erase(0,pos+6);
    pos=obs.find("</Info>");
    obs.erase(pos,obs.length()-pos);
    fout << "  ("<<obs<<" "<<st->second<<") ";}
 return fout.str();
}

void make_samplings(uint nbins, uint nsamp, SamplingMode mode,
                    RVector& samps, double mean, double errfactor)
{
 double total=0.0;
 RVector bins(nbins);
 for (uint k=0;k<nbins;k++){
    bins[k]=mean+get_random_double()*errfactor;
    total+=bins[k];}
 double avg=total/double(nbins);
 if (mode==Jackknife){
    samps.resize(nbins+1);
    samps[0]=avg;
    for (uint k=0;k<nbins;k++)
       samps[k+1]=(total-bins[k])/double(nbins-1);}
 else{
    samps.resize(nsamp+1);
    samps[0]=avg;
    for (uint k=1;k<=nsamp;k++){
       double val=0.0;
       for (uint kk=0;kk<nbins;kk++)
          val+=bins[rand()%nbins];
       samps[k]=val/double(nbins);}}

}


// ***********************************************************************

     //  This tests input and output

void testMCObsHandlerIO(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestMCObsHandlerIO")==0)
 return;

 cout << endl << "Starting TestMCObsHandlerIO"<<endl;

 
 try{
 XMLHandler xmlq(xml_in,"TestMCObsHandlerIO");
 XMLHandler xmlr(xmlq,"Setup");
 MCBinsInfo bininfo(xmlr);
 cout << bininfo.output()<<endl;
 MCSamplingInfo sampinfo(xmlr);
 cout << sampinfo.output()<<endl;
 MCObsGetHandler MCOH(xmlr,bininfo,sampinfo); 
 MCObsHandler MC(MCOH);

 cout << endl<<endl;
 cout << "Number of measurements in ensemble = "<<MC.getNumberOfMeasurements()<<endl;
 cout << "Ensemble ID = "<<MCOH.getEnsembleId()<<endl;
 cout << "Number of bins = "<<MC.getNumberOfBins()<<endl;

 SamplingMode defaultmode=MC.getDefaultSamplingMode();
 if (defaultmode==Jackknife) cout << "sampling mode = Jackknife"<<endl;
 else cout << "sampling mode = Bootstrap"<<endl;

 uint nbins=MC.getNumberOfBins();
 uint nsamp=nbins;
 if (defaultmode==Bootstrap)
   nsamp=MC.getNumberOfBootstrapResamplings();



 cout << "Set Carrot 1"<<endl;
 MCObsInfo obs1("Carrot",1); 
 RVector samps;
 make_samplings(nbins,nsamp,defaultmode,samps,65.0,1.0);
 uint count=0;
 for (MC.begin();!MC.end();++MC,count++){
    MC.putCurrentSamplingValue(obs1,samps[count]);} 
 cout << MC.getEstimate(obs1).output()<<endl;

 cout << "Set Carrot 2"<<endl;
 MCObsInfo obs2("Carrot",2); 
 double mean=43.0;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs2,data);} 
 cout << MC.getEstimate(obs2).output()<<endl;

 cout << "Set Rice 0"<<endl;
 MCObsInfo obs3("Rice",0);
 mean=-12.0;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs3,data);} 
 cout << MC.getEstimate(obs3).output()<<endl;

// MC.setToBootstrapMode();        
 cout << "Set Potato 4 Bootstrap"<<endl;
 MCObsInfo obs4("Potato",4); 
 mean=44.4;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs4,data);} 
 cout << MC.getEstimate(obs4).output()<<endl;

 MC.begin(); for (uint k=0;k<16;k++) ++MC;
 cout << "Bootstrap sampling 16 = "<<MC.getCurrentSamplingValue(obs4)<<endl;
 MC.setToJackknifeMode();        
 cout << "Set Potato 4 Jackknife"<<endl;
 mean=86.7;
 for (MC.begin();!MC.end();++MC){
    double data=mean+3.0*get_random_double();
    MC.putCurrentSamplingValue(obs4,data);} 
 MC.begin(); for (uint k=0;k<16;k++) ++MC;
 cout << "Jackknife sampling 16 = "<<MC.getCurrentSamplingValue(obs4)<<endl;
 cout << MC.getEstimate(obs4).output()<<endl;

 MCEstimate est4;
/* MC.setToBootstrapMode();        
 est4=MC.getEstimate(obs4);
 cout << "Obs4 estimate: "<<est4.output()<<endl;
 est4=MC.getBootstrapEstimate(obs4);
 cout << "Obs4 estimate: "<<est4.output()<<endl;
 est4=MC.getEstimate(obs4,Bootstrap);
 cout << "Obs4 estimate: "<<est4.output()<<endl;
 est4=MC.getJackknifeEstimate(obs4);
 cout << "Obs4 estimate: "<<est4.output()<<endl;
 est4=MC.getEstimate(obs4,Jackknife);
 cout << "Obs4 estimate: "<<est4.output()<<endl;
 MC.begin(); for (uint k=0;k<16;k++) ++MC;
 cout << "Bootstrap sampling 16 = "<<MC.getCurrentSamplingValue(obs4)<<endl;
 MC.setToJackknifeMode();        
 MC.begin(); for (uint k=0;k<16;k++) ++MC;
 cout << "Jackknife sampling 16 = "<<MC.getCurrentSamplingValue(obs4)<<endl;
*/
 MC.setToDefaultSamplingMode();

 MCObsInfo obs5("Banana",3); 
 mean=100.0;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs5,data);} 
 cout << MC.getEstimate(obs5).output()<<endl;

 MCObsInfo obs6("Orange",2); 
 mean=330.0;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs6,data);} 
 cout << MC.getEstimate(obs6).output()<<endl;

 MCObsInfo obs9("Orange",8); 
 mean=888.0;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs9,data);} 
 cout << MC.getEstimate(obs9).output()<<endl;

 string filename("dummy_samplings.0");
 set<MCObsInfo> obskeys;
 MCObsInfo obs55("WhiteRabbit",10); 

 map<MCObsInfo,string> statusmap;
 obskeys.insert(obs55); statusmap.insert(make_pair(obs55,"not avail"));
 obskeys.insert(obs1); statusmap.insert(make_pair(obs1,"OK"));
 obskeys.insert(obs3); statusmap.insert(make_pair(obs3,"OK"));
 obskeys.insert(obs6); statusmap.insert(make_pair(obs6,"OK"));
 XMLHandler xmllog;
 bool overwrite=true;

 cout <<endl<< "Write Attempt 1"<<endl<<endl;
 cout <<"Correct behavior below:"<<endl;
 cout << output_status(statusmap)<<endl;
 MC.writeSamplingValuesToFile(obskeys,filename,xmllog,overwrite);
 cout << xmllog.output()<<endl<<endl;

 cout <<endl<< "Write Attempt 2"<<endl<<endl;
 cout <<"Correct behavior below:"<<endl;
 statusmap.clear();statusmap.insert(make_pair(obs55,"not avail"));
 statusmap.insert(make_pair(obs1,"no overwrite"));
 statusmap.insert(make_pair(obs3,"no overwrite"));
 statusmap.insert(make_pair(obs6,"no overwrite"));
 cout << output_status(statusmap)<<endl;
 try{
 MC.writeSamplingValuesToFile(obskeys,"dummy_samplings.0",xmllog);}
 catch(const std::exception& msg){
    cout << "error: "<<msg.what()<<endl;}
 cout << xmllog.output()<<endl<<endl;

 obskeys.insert(obs4);
 statusmap.insert(make_pair(obs4,"OK"));
 cout <<endl<< "Write Attempt 3"<<endl<<endl;
 cout <<"Correct behavior below:"<<endl;
 cout << output_status(statusmap)<<endl;
 MC.writeSamplingValuesToFile(obskeys,"dummy_samplings.0",xmllog);
 cout << xmllog.output()<<endl<<endl;

 MC.setToDefaultSamplingMode();

 mean=0.0;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs1,data);
    data+=0.3125;} 
 cout << MC.getEstimate(obs1).output()<<endl;

 mean=43.0;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs2,data);
    data+=2.23; data*=0.854;} 
 cout << MC.getEstimate(obs2).output()<<endl;

 MCObsInfo obs7("Rice",7);
 mean=-12.0;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs7,data);
    data+=1.823; data*=-1.12;} 
 cout << obs7.output()<<MC.getEstimate(obs7).output()<<endl;

 MCObsInfo obs8("Potato",2);
 mean=0.0;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs8,data);} 
 cout << MC.getEstimate(obs8).output()<<endl;
 
 mean=10.0;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs5,data);} 
 cout << MC.getEstimate(obs5).output()<<endl;

 mean=33.0;
 for (MC.begin();!MC.end();++MC){
    double data=mean+2.0*get_random_double();
    MC.putCurrentSamplingValue(obs6,data);} 
 cout << MC.getEstimate(obs6).output()<<endl;




 obskeys.insert(obs7);
 obskeys.insert(obs8);
 cout <<endl<< "Write Attempt 4"<<endl<<endl;
 if (fileExists("dummy_samplings.2")) cout << "file dummy_samplings.2 exists already"<<endl;
 statusmap.clear();statusmap.insert(make_pair(obs55,"not avail"));
 statusmap.insert(make_pair(obs1,"OK"));
 statusmap.insert(make_pair(obs3,"OK"));
 statusmap.insert(make_pair(obs4,"OK"));
 statusmap.insert(make_pair(obs6,"OK"));
 statusmap.insert(make_pair(obs7,"OK"));
 statusmap.insert(make_pair(obs8,"OK"));

 cout <<"Correct behavior below:"<<endl;
 cout << output_status(statusmap)<<endl;
 MC.writeSamplingValuesToFile(obskeys,"dummy_samplings.2",xmllog);
 cout << xmllog.output()<<endl<<endl;

 cout <<endl<< "Write Attempt 5"<<endl<<endl;
 cout <<"Correct behavior below:"<<endl;
 cout <<"(Carrot 1 OK) (Rice 0 not avail) (Rice 7 OK) (Potato 2 OK)  "<<endl
      << "(Potato 4 not avail) (Orange 2 OK) (WhiteRabbit 10 not avail)"<<endl<<endl;
 MC.writeSamplingValuesToFile(obskeys,"dummy_samplings.1",xmllog,true);
 cout << xmllog.output()<<endl<<endl;

 obskeys.erase(obs1);
 obskeys.erase(obs3);
 obskeys.erase(obs4);
 obskeys.insert(obs5);
 obskeys.insert(obs2);
 obskeys.insert(obs9);
 MC.eraseSamplings(obs5);
 MC.eraseData(MCObsInfo("WhiteRabbit",10));
 MC.eraseData(obs2);

 obskeys.insert(MCObsInfo("WhiteRabbit",10));
 cout <<endl<< "Write Attempt 6"<<endl<<endl;
 cout <<"Correct behavior below:"<<endl;
 cout << "OK: (Orange 8)"<<endl;
 cout <<"Already: (Rice 0) (Orange 2) (Potato 4)"<<endl;
 cout << "Not avail: (Rice 7) (Potato 2) (Banana 3) (Carrot 2) (WhiteRabbit 10)"<<endl<<endl;
 MC.writeSamplingValuesToFile(obskeys,filename,xmllog);
 cout << xmllog.output()<<endl<<endl;


 MC.clearData();

 
 cout <<endl<< "Read Attempt 1"<<endl<<endl;
 cout <<"Correct behavior:   open and read from "<<filename<<endl;
 MC.readSamplingValuesFromFile(filename,xmllog);
 cout << xmllog.output()<<endl;

 for (set<MCObsInfo>::const_iterator it=obskeys.begin();it!=obskeys.end();it++){
   cout << "Reading "<<it->output()<<endl;
   cout << "query = "<<MC.queryFullAndSamplings(*it)<<endl;
   try{ MCEstimate est=MC.getEstimate(*it); cout << est.output()<<endl;}
   catch(...){ cout << "no get:exception"<<endl;}}
 
// cout <<endl<< "Read Attempt 2"<<endl<<endl;
// cout <<"Correct behavior:   fail since already read"<<endl;
// try{
//    MC.readSamplingValuesFromFile(filename,xmllog);}
// catch(std::exception& xp){
//    cout << "exception caught: "<<xp.what()<<endl;}
// cout << xmllog.output()<<endl;

// cout <<endl<< "Read Attempt 3"<<endl<<endl;
// cout <<"Correct behavior: "<<endl;
// cout << "Read OK: (Carrot 1) (Rice 0) (Orange 2) (Potato 4) (Orange 8)"<<endl;
// try{
//    MC.readSamplingValuesFromFile(filename,xmllog);}
// catch(std::exception& xp){
//    cout << "exception caught: "<<xp.what()<<endl;}
// cout << xmllog.output()<<endl;

 MCEstimate est4B=MC.getEstimate(obs4,defaultmode);
 cout << "(Potato 4) estimate:"<<est4B.output()<<endl;
 cout << "(Potato 4) original estimate:"<<est4.output()<<endl; 

// obskeys.clear();
// obskeys.insert(MCObsInfo("WhiteRabbit",10));
// obskeys.insert(obs4);
// MC.readSamplingValuesFromFile(obskeys,filename,xmllog);
// cout << xmllog.output()<<endl;


 cout << endl<<" ************************************* "<<endl<<endl;
 cout << "Testing bins read from/write to file"<<endl<<endl;

 Vector<double> dummybins(nbins);

 for (uint k=0;k<nbins;k++){
    dummybins[k]=1.2317*(rand()%4096);}
 MCObsInfo binkey1("Monkey",8,true);
 cout << "bin key 1: "<<binkey1.output()<<endl;
 MC.putBins(binkey1,dummybins);
 est4=MC.getEstimate(binkey1,defaultmode);
 cout << "estimate: "<<est4.output()<<endl;

 for (uint k=0;k<nbins;k++){
    dummybins[k]=0.03184*(rand()%2048);}
 MCObsInfo binkey2("Elephant",3,true);
 cout << "bin key 2: "<<binkey2.output()<<endl;
 MC.putBins(binkey2,dummybins);
 est4=MC.getEstimate(binkey2,defaultmode);
 cout << "estimate: "<<est4.output()<<endl;

 for (uint k=0;k<nbins;k++){
    dummybins[k]=0.001237*(rand()%3333);}
 MCObsInfo binkey3("Giraffe",5,true);
 cout << "bin key 3: "<<binkey3.output()<<endl;
 MC.putBins(binkey3,dummybins);
 est4=MC.getEstimate(binkey3,defaultmode);
 cout << "estimate: "<<est4.output()<<endl;

 for (uint k=0;k<nbins;k++){
    dummybins[k]=0.07773*(rand()%2048);}
 MCObsInfo binkey4("Donkey",6,true);
 cout << "bin key 4: "<<binkey4.output()<<endl;
 MC.putBins(binkey4,dummybins);
 est4=MC.getEstimate(binkey4,defaultmode);
 cout << "estimate: "<<est4.output()<<endl;

 MCObsInfo binkey5("Donkey",11,true);
 cout << "bin key 5: "<<binkey5.output()<<endl;

 cout << "Write bins test 1:"<<endl;
 obskeys.clear();
 obskeys.insert(binkey1);
 obskeys.insert(binkey2);
 MC.writeBinsToFile(obskeys,"dummybinfile.0",xmllog,true);
 cout << "Correct behavior: (Monkey 8 OK) (Elephant 3 OK)"<<endl;
 cout << xmllog.output()<<endl;

 cout << "Write bins test 2:"<<endl;
 obskeys.insert(binkey3);
 MC.writeBinsToFile(obskeys,"dummybinfile.0",xmllog,true);
 cout << "Correct behavior: (Monkey 8 OK) (Elephant 3 OK) (Giraffe 5 OK)"<<endl;
 cout << xmllog.output()<<endl;

 cout << "Write bins test 3:"<<endl;
 MC.writeBinsToFile(obskeys,"dummybinfile.0",xmllog,false);
 cout << "Correct behavior: "<<endl;
 cout << "Already: (Monkey 8) (Elephant 3) (Giraffe 5)"<<endl;
 cout << xmllog.output()<<endl;

 cout << "Write bins test 4:"<<endl;
 obskeys.insert(binkey4);
 obskeys.insert(binkey5);
 MC.writeBinsToFile(obskeys,"dummybinfile.0",xmllog,false);
 cout << "Correct behavior: "<<endl;
 cout << "Already: (Monkey 8) (Elephant 3) (Giraffe 5)"<<endl;
 cout << "Not available: (Donkey 11)"<<endl;
 cout << "OK: (Donkey 6)"<<endl;
 cout << xmllog.output()<<endl;

 MC.clearData();
/*
 cout << "Read bins test 1:"<<endl;
 MC.readBinsFromFile("",xmllog);
 cout << xmllog.output()<<endl;

// cout << "Read bins test 2:"<<endl;
// MC.readBinsFromFile("farmer.1",xmllog);
// cout << xmllog.output()<<endl;

 cout << "Read bins test 3:"<<endl;
 MC.readBinsFromFile("dummybinfile.0",xmllog);
 cout << xmllog.output()<<endl;

 MC.clearData();
 obskeys.clear();
 obskeys.insert(binkey4);
 obskeys.insert(binkey5);
 cout << "Read bins test 4:"<<endl;
 MC.readBinsFromFile(obskeys,"dummybinfile.0",xmllog);
 cout << xmllog.output()<<endl;



 MC.clearSamplings();
 MC.setToBootstrapMode();        

 MCObsInfo tobs1("Carrot",1);
 for (MC.begin();!MC.end();++MC){
    data=43.0+double((rand()%4096)-2048)/3217.0;
    MC.putCurrentSamplingValue(tobs1,data);} 

 MCObsInfo tobs2("Carrot",2);
 for (MC.begin();!MC.end();++MC){
    data=21.0+double((rand()%4096)-2048)/3217.0;
    MC.putCurrentSamplingValue(tobs2,data);}

 MCObsInfo tobs3("Rice",0);
 for (MC.begin();!MC.end();++MC){
    data=87.0+double((rand()%4096)-2048)/1217.0;
    MC.putCurrentSamplingValue(tobs3,data);} 

 MCObsInfo tobs4("Potato",4);
 for (MC.begin();!MC.end();++MC){
    data=44.4+double((rand()%4096)-2048)/1217.0;
    MC.putCurrentSamplingValue(tobs4,data);} 
 
 MCObsInfo tobs5("Banana",3);
 for (MC.begin();!MC.end();++MC){
    data=10.0+double((rand()%4096)-2048)/8217.0;
    MC.putCurrentSamplingValue(tobs5,data);} 

 MCObsInfo tobs6("Orange",2);
 for (MC.begin();!MC.end();++MC){
    data=33.0+double((rand()%4096)-2048)/8217.0;
    MC.putCurrentSamplingValue(tobs6,data);} 

 MCObsInfo tobs9("Orange",8);
 for (MC.begin();!MC.end();++MC){
    data=8.0+double((rand()%4096)-2048)/217.0;
    MC.putCurrentSamplingValue(tobs9,data);}

 set<MCObsInfo> tobskeys;
 tobskeys.insert(tobs1);
 tobskeys.insert(tobs2);
 tobskeys.insert(tobs3);
 tobskeys.insert(tobs4);
 tobskeys.insert(tobs5);
 tobskeys.insert(tobs6);
 tobskeys.insert(tobs9);
 overwrite=true;

 MC.writeSamplingValuesToFile(tobskeys,"TestDataBoot",Bootstrap,xmllog,overwrite);
 cout << xmllog.output()<<endl<<endl;


 MC.setToJackknifeMode();        

 for (MC.begin();!MC.end();++MC){
    data=43.0+double((rand()%4096)-2048)/3217.0;
    MC.putCurrentSamplingValue(tobs1,data);} 

 for (MC.begin();!MC.end();++MC){
    data=21.0+double((rand()%4096)-2048)/3217.0;
    MC.putCurrentSamplingValue(tobs2,data);}

 for (MC.begin();!MC.end();++MC){
    data=87.0+double((rand()%4096)-2048)/1217.0;
    MC.putCurrentSamplingValue(tobs3,data);} 

 for (MC.begin();!MC.end();++MC){
    data=44.4+double((rand()%4096)-2048)/1217.0;
    MC.putCurrentSamplingValue(tobs4,data);} 
 
 for (MC.begin();!MC.end();++MC){
    data=10.0+double((rand()%4096)-2048)/8217.0;
    MC.putCurrentSamplingValue(tobs5,data);} 

 for (MC.begin();!MC.end();++MC){
    data=33.0+double((rand()%4096)-2048)/8217.0;
    MC.putCurrentSamplingValue(tobs6,data);} 

 for (MC.begin();!MC.end();++MC){
    data=8.0+double((rand()%4096)-2048)/217.0;
    MC.putCurrentSamplingValue(tobs9,data);}

 overwrite=true;

 MC.writeSamplingValuesToFile(tobskeys,"TestDataJack",Jackknife,xmllog,overwrite);
 cout << xmllog.output()<<endl<<endl;


 MCObsInfo bobs1("Carrot",1,true);
 MCObsInfo bobs2("Carrot",2,true);
 MCObsInfo bobs3("Rice",0,true);
 MCObsInfo bobs4("Potato",4,true);
 MCObsInfo bobs5("Banana",3,true);
 MCObsInfo bobs6("Orange",2,true);
 MCObsInfo bobs9("Orange",8,true);
 set<MCObsInfo> bobskeys;
 bobskeys.insert(bobs1);
 bobskeys.insert(bobs2);
 bobskeys.insert(bobs3);
 bobskeys.insert(bobs4);
 bobskeys.insert(bobs5);
 bobskeys.insert(bobs6);
 bobskeys.insert(bobs9);

 for (uint k=0;k<nbins;k++)
    dummybins[k]=43.0+double((rand()%4096)-2048)/3217.0;
 MC.putBins(bobs1,dummybins);

 for (uint k=0;k<nbins;k++)
    dummybins[k]=21.0+double((rand()%4096)-2048)/3217.0;
 MC.putBins(bobs2,dummybins);

 for (uint k=0;k<nbins;k++)
    dummybins[k]=87.0+double((rand()%4096)-2048)/1217.0;
 MC.putBins(bobs3,dummybins);

 for (uint k=0;k<nbins;k++)
    dummybins[k]=44.4+double((rand()%4096)-2048)/1217.0;
 MC.putBins(bobs4,dummybins);
 
 for (uint k=0;k<nbins;k++)
    dummybins[k]=10.0+double((rand()%4096)-2048)/8217.0;
 MC.putBins(bobs5,dummybins);

 for (uint k=0;k<nbins;k++)
    dummybins[k]=33.0+double((rand()%4096)-2048)/8217.0;
 MC.putBins(bobs6,dummybins);

 for (uint k=0;k<nbins;k++)
    dummybins[k]=8.0+double((rand()%4096)-2048)/217.0;
 MC.putBins(bobs9,dummybins);

 overwrite=true;

 MC.writeBinsToFile(bobskeys,"TestDataBins",xmllog,overwrite);
 cout << xmllog.output()<<endl<<endl;
*/
 }
 catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}

}

// ***********************************************************************

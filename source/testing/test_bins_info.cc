#include "bins_info.h"
#include <cstdio>
#include <ctime>
#include <map>
#include <cstdlib>
#include <iostream>

using namespace std;


  // *************************************************************************


void testMCBinsInfo(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestMCBinsInfo")==0)
 return;

 cout << endl << "Starting TestMCBinsInfo"<<endl;
 
 
 try{
 XMLHandler xmlr(xml_in,"TestMCBinsInfo");

 MCBinsInfo Mb(xmlr);
 cout << Mb.output()<<endl;

 const MCEnsembleInfo& mcens=Mb.getMCEnsembleInfo();
 cout << "ensemble Id = "<<mcens.getId()<<endl;
 cout << "number of streams = "<<mcens.getNumberOfStreams()<<endl;

 cout << "number of meas = "<<Mb.getNumberOfMeasurements()<<endl;
 cout << "number of bins = "<<Mb.getNumberOfBins()<<endl;
 cout << "rebin factor = "<<Mb.getRebinFactor()<<endl;
 const set<uint>& omits=Mb.getOmissions();
 cout << "omits: ";
 for (set<uint>::const_iterator it=omits.begin();it!=omits.end();it++)
    cout << " "<<*it;
 cout << endl;
 cout << "number of configs = "<<Mb.getNumberOfConfigs()<<endl;
 cout << "Nt = "<<Mb.getLatticeTimeExtent()<<endl;
 cout << "Nx = "<<Mb.getLatticeXExtent()<<endl;
 cout << "Ny = "<<Mb.getLatticeYExtent()<<endl;
 cout << "Nz = "<<Mb.getLatticeZExtent()<<endl;


 MCBinsInfo Mb2(Mb);
 cout << Mb2.output()<<endl;

 cout << "equal? "<<(Mb2==Mb)<<endl;
 cout << "not equal? "<<(Mb2!=Mb)<<endl;

 Mb.addOmission(218);
 Mb.setRebin(5);
 cout << Mb.output()<<endl;

 cout << "equal? "<<(Mb2==Mb)<<endl;
 cout << "not equal? "<<(Mb2!=Mb)<<endl;

 Mb2=Mb;
 cout << Mb2.output()<<endl;
 cout << "equal? "<<(Mb2==Mb)<<endl;
 cout << "not equal? "<<(Mb2!=Mb)<<endl;

    //  test isConsistentWith

 uint nmeas=2003; uint nstream=1;
 uint ns=24; uint nt=48;
 MCEnsembleInfo mctest("DummyG12",nmeas,nstream,ns,ns,ns,nt);

 MCBinsInfo b1test(mctest);
 uint rebin1=10;
 b1test.setRebin(rebin1);
 b1test.addOmission(5);
 b1test.addOmission(17);
 b1test.addOmission(547);
 cout << b1test.output()<<endl;

 cout << "Test 1:"<<endl;
 MCBinsInfo b2test(mctest);
 uint rebin2=20;
 b2test.setRebin(rebin2);
 b2test.addOmission(5);
 b2test.addOmission(17);
 b2test.addOmission(547);
// b2test.addOmission(333);
 cout << b2test.output()<<endl;
 cout << "b1test.isConsistentWith(b2test)? "<<b1test.isConsistentWith(b2test)<<endl;

   //   get the omissions and rebin factor
/*
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

 XMLHandler xmlrdtest(xmlr,"DoReadTests");
 MCObsGetHandler MCOH(xmlrdtest); 

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

 list<XMLHandler> xmlop=xmlrdtest.find("OperatorString");
 list<OperatorInfo> oplist;
 for (list<XMLHandler>::iterator it=xmlop.begin();it!=xmlop.end();++it)
    oplist.push_back(OperatorInfo(*it));


 cout << endl<<endl<<"NOW firing up MCObsHandler!!"<<endl<<endl;

 MCObsHandler MC(MCOH,rebin,omit);
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

#ifdef COMPLEXNUMBERS
    int kkmax=2;
#else
    int kkmax=1;
#endif

 int goodq=0, badq=0;
 int goodc=0, badc=0;
 int goode=0, bade=0;
 cout << "Query VEV bins"<<endl;
 for (list<OperatorInfo>::iterator it=oplist.begin();it!=oplist.end();++it){
    MCObsInfo mcobs(*it);
    for (int kk=0;kk<kkmax;++kk){
    if (kk>0) mcobs.setToImaginaryPart();
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
          if (std::abs(res[bin]-crt)>1e-6*std::abs(crt)) badc++; else goodc++;
          if (std::abs(MC.getBin(mcobs,bin)-crt)>1e-6*std::abs(crt)) badc++; else goodc++;}}}
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
                //cout << crt<<" "<<res[bin]<<" "<<MC.getBin(mcobs,bin)<<endl;
                if (std::abs(res[bin]-crt)>1e-5*std::abs(crt)){ badc++; cout << "diff = "<<res[bin]-crt<<endl;} else goodc++;
                if (std::abs(MC.getBin(mcobs,bin)-crt)>1e-5*std::abs(crt)) badc++; else goodc++;}}}
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

 double fullmean=MC.getFullSampleValue(obskey);
 cout << "fullmean = "<<fullmean<<endl;
 MC.setToJackknifeMode(); 
 int count=0;
 for (MC.begin();!MC.end();++MC){
    cout << MC.getCurrentSamplingValue(obskey)<<endl; count++;}
 cout << "count = "<<count<<endl;


 cout <<endl<<endl;

 MC.setBootstrapper(888,32,54,false);
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

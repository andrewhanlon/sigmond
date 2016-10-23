#include <cstdio>
#include <ctime>
#include <algorithm>
#include "test_obs_get_handler.h"
#include "bins_handler.h"
#include "samplings_handler.h"

using namespace std;
using namespace LaphEnv;


#ifdef SINGLEPRECISION
const double eps=1e-5;
const double ceps=5e-5;
#else
const double eps=1e-11;
const double ceps=5e-11;
#endif
const double deps=1e-11;

double get_random_double()
{
 return double(rand() % 1048576)/1048576.0 -0.5;
}


void testMCObsGetHandler(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestMCObsGetHandler")==0)
 return;

 cout << endl << "Starting TestMCObsGetHandler"<<endl;

 try{
 XMLHandler xmlr(xml_in,"TestMCObsGetHandler");
 MCBinsInfo bininfo(xmlr);
 cout << "Bins Info: "<<bininfo.output()<<endl;
 MCSamplingInfo sampinfo(xmlr);
 cout << "Sampling Info: "<<sampinfo.output()<<endl;

 MCObsGetHandler MCOH(xmlr,bininfo,sampinfo); 

 cout << endl<<endl;
 cout << "Number of measurements in ensemble = "<<MCOH.getNumberOfMeasurements()<<endl;
 cout << "Ensemble ID = "<<MCOH.getEnsembleId()<<endl;

 cout << "FileMap:"<<endl;
 XMLHandler xmlout;
 MCOH.getFileMap(xmlout);
 cout << endl<<endl<<"***********************"<<endl<<endl;
 cout << xmlout.output()<<endl;
 cout << endl<<endl<<"***********************"<<endl<<endl;

 XMLHandler xmli(xmlr,"GetSerialIndices");

 cout << xmli.output()<<endl;

 list<XMLHandler> getinds=xmli.find("SerialIndex");
 vector<unsigned int> serinds(getinds.size());
 cout << "serinds size = "<<serinds.size()<<endl;

 unsigned int count=0;
 for (list<XMLHandler>::iterator it=getinds.begin();it!=getinds.end();it++,count++){
    cout << it->output()<<endl;
    xmlread(*it,"SerialIndex",serinds[count],"testMCObsGetHandler");
    cout << "serinds["<<count<<"] = "<<serinds[count]<<endl;}
 xmli.clear();

 Scalar data;
 XMLHandler xmlg(xmlr,"GetKeys");

 cout << xmlg.output()<<endl;
 list<XMLHandler> getkeys=xmlg.find("MCObservable");
 cout << "getkeys size = "<<getkeys.size()<<endl;

 cout << endl<<endl<<"***********************************"<<endl<<endl;
 for (list<XMLHandler>::iterator it=getkeys.begin();it!=getkeys.end();it++){
    MCObsInfo obs(*it);
    cout << "Observable:"<<obs.output()<<endl;
    for (unsigned int k=0;k<serinds.size();k++){
       try{
          MCOH.getData(obs,serinds[k],data);
          cout << "result for serial index "<<serinds[k]<<" = "<<data<<endl;}
       catch(std::exception& xp){
          cout << "exception caught for "<<serinds[k]<<" "<<xp.what()<<endl;}
       cout << "query = "<<MCOH.queryData(obs,serinds[k])<<endl;}}



/*
 
 set<CorrelatorInfo> corrSet=CH.getCorrelatorSet();
 dcmplx result;
 for (set<CorrelatorInfo>::const_iterator 
      ct=corrSet.begin();ct!=corrSet.end();ct++){
    cout << "CorrelatorInfo: "<<endl;
    cout << ct->output()<<endl;
    CorrelatorAtTimeInfo ctinfo(*ct,0);
    for (int t=0;t<9;t++)
    for (unsigned int serial=0;serial<CH.getNumberOfMeasurements();serial++){
       ctinfo.resetTimeSeparation(t);
       try{CH.getData(ctinfo,serial,result);
       cout << "t="<<t<<" config="<<serial<<"  result: "<<result<<endl;}
       catch(const std::exception& xp){cout << "t="<<t<<" config="<<serial<<"  result not found"<<endl;}}}



 MCDataHandler<BLCorrelatorDataHandler,BLCorrelatorDataHandler::MCDataKey,dcmplx> mch(CH);
 cout << endl<<endl<<" MC Data Handler NOW"<<endl<<endl;
 cout << "number of MC measurements = "<<mch.getNumberOfMeasurements()<<endl;

 for (set<CorrelatorInfo>::const_iterator 
      ct=corrSet.begin();ct!=corrSet.end();ct++){
    cout << "CorrelatorInfo: "<<endl;
    cout << ct->output()<<endl;
    CorrelatorAtTimeInfo ctinfo(*ct,0);
    for (int t=3;t<9;t++){
       ctinfo.resetTimeSeparation(t);
       result=mch.getMean(ctinfo);
       cout << "t="<<t<<" mean: "<<result<<endl;}}
*/
 }
 catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}

}


// ***********************************************************************


void VEVCorrect::addValue(const OperatorInfo& opinfo, int serind, const InScalar& invalue)
{
 map<OperatorInfo,map<int,InScalar> >::iterator it=values.find(opinfo);
 if (it!=values.end()) 
    it->second.insert(make_pair(serind,invalue));
 else{ 
    map<int,InScalar> newmap; 
    newmap.insert(make_pair(serind,invalue)); 
    values.insert(make_pair(opinfo,newmap));}
}

bool VEVCorrect::getCorrect(const MCObsInfo& mcobs, int serind, double& result) const
{
 if (!mcobs.isVEV()) return false;
 map<OperatorInfo,map<int,InScalar> >::const_iterator it=values.find(mcobs.getVEVInfo());
 if (it==values.end()) return false;
 map<int,InScalar>::const_iterator st=it->second.find(serind);
 if (st==it->second.end()) return false;
 InScalar res=st->second;
 if (mcobs.isRealPart()) result=real(res);
 else result=imag(res);
 return true;
}

void increment_serial_index(uint& serind, set<uint>& omit)
{
 serind++;
 while (omit.find(serind)!=omit.end()) serind++;
}

void start_serial_index(uint& serind, set<uint>& omit)
{
 serind=0;
 while (omit.find(serind)!=omit.end()) serind++;
}

bool compare_vectors(const RVector& v1, const RVector& v2, double factor)
{
 if (v1.size()!=v2.size()) return false;
 for (uint k=0;k<v1.size();k++)
    if (std::abs(v1[k]-v2[k])>factor*std::max(std::abs(0.5*(v1[k]+v2[k])),1.0)) 
       return false;
 return true;
}


bool compare_floats(const double& v1, const double& v2, double factor)
{
 if (std::abs(v1-v2)>factor*std::max(std::abs(0.5*(v1+v2)),1.0)){ cout <<"MISMATCH: "<< v1<<" "<<v2<<endl; return false;}
 return true;
}


bool VEVCorrect::getCorrectBins(const MCObsInfo& mcobs, const MCBinsInfo& bininfo, 
                                RVector& result) const
{
 uint nbins=bininfo.getNumberOfBins();
 uint rebin=bininfo.getRebinFactor();
 set<uint> omit(bininfo.getOmissions());
 uint serind; start_serial_index(serind,omit);
 result.clear(); 
 double val;
 vector<double> buffer;
 for (uint bin=0;bin<nbins;bin++){
    double bval=0.0;
    for (uint j=0;j<rebin;j++){
       if (!(getCorrect(mcobs,serind,val))) {return false;}
       bval+=val; increment_serial_index(serind,omit);}
    bval/=double(rebin);
    buffer.push_back(bval);}
 result=RVector(buffer);
 return true;
}


void CorrCorrect::addValue(const CorrelatorAtTimeInfo& corrinfo, int serind, const InScalar& invalue)
{
 CorrelatorAtTimeInfo cf(corrinfo.getSink(),corrinfo.getSource(),
                         corrinfo.getTimeSeparation(),false,false);
 map<CorrelatorAtTimeInfo,map<int,InScalar> >::iterator it=values.find(cf);
 if (it!=values.end()) 
    it->second.insert(make_pair(serind,invalue));
 else{ 
    map<int,InScalar> newmap; 
    newmap.insert(make_pair(serind,invalue)); 
    values.insert(make_pair(cf,newmap));}
}

bool CorrCorrect::getCorrect(const MCObsInfo& mcobs, int serind, double& result) const
{  
 if ((!mcobs.isCorrelatorAtTime())&&(!mcobs.isHermitianCorrelatorAtTime())) return false;
 OperatorInfo sink(mcobs.getCorrelatorSinkInfo());
 OperatorInfo source(mcobs.getCorrelatorSourceInfo());
 int tind=mcobs.getCorrelatorTimeIndex();
 CorrelatorAtTimeInfo cf1(sink,source,tind,false,false);
 bool flag1=true;
 double res1=0.0;
 map<CorrelatorAtTimeInfo,map<int,InScalar> >::const_iterator it=values.find(cf1);
 if (it==values.end()) flag1=false;
 else{
    map<int,InScalar>::const_iterator st=it->second.find(serind);
    if (st==it->second.end()) flag1=false;
    else{
       InScalar res=st->second;
       if (mcobs.isRealPart()) res1=real(res);
       else res1=imag(res);
       flag1=true;}}
 if ((!hermitian)||(sink==source)){
    result=res1; 
    return flag1;}
 CorrelatorAtTimeInfo cf2(source,sink,tind,false,false);
 bool flag2=true;
 double res2=0.0;
 it=values.find(cf2);
 if (it==values.end()) flag2=false;
 else{
    map<int,InScalar>::const_iterator st=it->second.find(serind);
    if (st==it->second.end()) flag2=false;
    else{
       InScalar res=st->second;
       if (mcobs.isRealPart()) res2=real(res);
       else res2=-imag(res);
       flag2=true;}}
 if (flag1&&flag2) result=0.5*(res1+res2);
 else if (flag1) result=res1;
 else if (flag2) result=res2;
 return (flag1||flag2);
}


bool CorrCorrect::getCorrectBins(const MCObsInfo& mcobs, const MCBinsInfo& bininfo, 
                                 RVector& result) const
{
 uint nbins=bininfo.getNumberOfBins();
 uint rebin=bininfo.getRebinFactor();
 set<uint> omit(bininfo.getOmissions());
 uint serind; start_serial_index(serind,omit);
 result.clear(); 
 double val;
 vector<double> buffer;
 for (uint bin=0;bin<nbins;bin++){
    double bval=0.0;
    for (uint j=0;j<rebin;j++){
       if (!(getCorrect(mcobs,serind,val))) {return false;}
       bval+=val; increment_serial_index(serind,omit);}
    bval/=double(rebin);
    buffer.push_back(bval);}
 result=RVector(buffer);
 return true;
}

void corr_data_assign(Array<InScalar>& data, int tind, int serind, double Aval, 
                      double Bval, double Cval)
{
 double g=Aval*exp(-Bval*tind)*(1.0+Cval*serind)*get_random_double();
#ifdef COMPLEXNUMBERS
 data[0]=InScalar(g, 0.12*g);
 data[1]=InScalar(2.4*g,0.21);
 data[2]=InScalar(3.6*g,0.45);
#else
 data[0]=g;
 data[1]=2.4*g;
 data[2]=3.6*g;
#endif
}

string get_number_type()
{
#ifdef DOUBLEPRECISION
#ifdef COMPLEXNUMBERS
 return string("ComplexDoublePrecision");
#else
 return string("RealDoublePrecision");
#endif
#endif

#ifdef SINGLEPRECISION
#ifdef COMPLEXNUMBERS
 return string("ComplexSinglePrecision");
#else
 return string("RealSinglePrecision");
#endif
#endif
}




void make_fake_correlator_file(const MCEnsembleInfo& ens, const string& filename,
                               const CorrelatorInfo& corr,
                               int tmin, int tmax, double Aval, double Bval, 
                               double Cval, const set<int>& omit,
                               CorrCorrect& CC)
{
 XMLHandler xmlh;
 xmlh.set_root("CorrelatorHandlerDataFile");
 xmlh.put_child("NumberType",get_number_type());
 XMLHandler xmlt; ens.output(xmlt); xmlh.put_child(xmlt);
 corr.output(xmlt,true); xmlh.put_child(xmlt);
 cout << "Header:"<<xmlh.output()<<endl;
 IOMap<LaphEnv::BLCorrelatorDataHandler::RecordKey,Array<InScalar> > iom;
 iom.openNew(filename,"Laph--CorrelatorFile",xmlh.str(),false);
 Array<InScalar> data(3);
 for (int tind=tmin;tind<=tmax;++tind){
    CorrelatorAtTimeInfo ct(corr,tind,false);
    for (int serind=0;serind<int(ens.getNumberOfMeasurements());++serind){
       if (omit.find(serind)==omit.end()){
       LaphEnv::BLCorrelatorDataHandler::RecordKey ckey(tind,serind);
       corr_data_assign(data,tind,serind,Aval,Bval,Cval);
       iom.put(ckey,data);
       InScalar ddata=data[0]+data[1]+data[2];
       CC.addValue(ct,serind,ddata);}}}
 iom.close(); 
}


void vev_data_assign(InScalar& data, int serind, double Aval, double Cval)
{
 double g=Aval*(1.0+Cval*serind)*get_random_double();
#ifdef COMPLEXNUMBERS
 data=InScalar(7.0*g, 0.12*g+0.66);
#else
 data=7.0*g;
#endif
}


void make_fake_vev_file(const MCEnsembleInfo& ens, const string& filename,
                        const OperatorInfo& opinfo, double Aval, double Cval, 
                        const set<int>& omit, VEVCorrect& VC)
{
 XMLHandler xmlh;
 xmlh.set_root("VEVHandlerDataFile");
 xmlh.put_child("NumberType",get_number_type());
 XMLHandler xmlt; ens.output(xmlt); xmlh.put_child(xmlt);
 opinfo.output(xmlt,true); xmlh.put_child(xmlt);
 cout << "Header:"<<xmlh.output()<<endl;
 IOMap<LaphEnv::BLVEVDataHandler::RecordKey,InScalar > iom;
 iom.openNew(filename,"Laph--VEVFile",xmlh.str(),false);
 InScalar data;
 for (int serind=0;serind<int(ens.getNumberOfMeasurements());++serind){
    if (omit.find(serind)==omit.end()){
       LaphEnv::BLVEVDataHandler::RecordKey ckey(serind);
       vev_data_assign(data,serind,Aval,Cval);
       iom.put(ckey,data);
       VC.addValue(opinfo,serind,data);}}
 iom.close(); 
}


void bin_data_assign(double& data, int serind, double Aval, double Cval)
{
 double g=Aval*(1.0+Cval*serind)*get_random_double();
 data=3.3*g;
}


void make_fake_bin_file(const MCObsInfo& obskey, const string& filename,
                        const MCBinsInfo& bininfo, double Aval, double Cval, 
                        BinsCorrect& BC)
{
 BinsPutHandler SP(bininfo,filename,true,false);
 double data;
 uint nbins=bininfo.getNumberOfBins();
 uint rebin=bininfo.getRebinFactor();
 set<uint> omit(bininfo.getOmissions());
 vector<double> buffer;
 for (int serind=0;serind<int(bininfo.getNumberOfMeasurements());++serind){
    if (omit.find(serind)==omit.end()){
       bin_data_assign(data,serind,Aval,Cval);
       buffer.push_back(data);}}
 if ((buffer.size()/rebin)!=nbins)
    throw(std::runtime_error("Nbins mismatch in make_fake_bin_file"));
 if (rebin>1){
    vector<double> temp(buffer);
    buffer.clear();
    uint count=0;
    for (uint bin=0;bin<nbins;bin++){
       double binvalue=temp[count++];
       for (uint k=1;k<rebin;k++)
          binvalue+=temp[count++];
       binvalue/=double(rebin);
       buffer.push_back(binvalue);}}
 if (buffer.size()!=nbins)
    throw(std::runtime_error("Nbins mismatch in make_fake_bin_file"));
 SP.putData(obskey,RVector(buffer));
 BC.addValue(obskey,RVector(buffer));
}


void BinsCorrect::addValue(const MCObsInfo& mcobs, const RVector& inbins)
{
 map<MCObsInfo,RVector >::iterator it=values.find(mcobs);
 if (it!=values.end()) 
    it->second=inbins;
 else{ 
    values.insert(make_pair(mcobs,inbins));}
}

bool BinsCorrect::getCorrect(const MCObsInfo& mcobs, RVector& result) const
{
 map<MCObsInfo,RVector >::const_iterator it=values.find(mcobs);
 if (it==values.end()){ result.clear(); return false;}
 result=it->second;
 return true;
}



void samplings_data_assign(double& data, int serind, double Aval, double Cval)
{
 double g=Aval*(1.0+Cval*serind)*get_random_double();
 data=0.4172*g;
}


void make_fake_samplings_file(const MCObsInfo& obskey, const std::string& filename,
                              const MCBinsInfo& bininfo, const MCSamplingInfo& sampinfo, 
                              double Aval, double Cval, SamplingsCorrect& SC)
{
 SamplingsPutHandler SP(bininfo,sampinfo,filename,true,false);
 double data;
 uint nsamp=sampinfo.getNumberOfReSamplings(bininfo);
 vector<double> buffer;
   // first is full sampling, then do the samplings
 for (int serind=0;serind<=int(nsamp);++serind){
    samplings_data_assign(data,serind,Aval,Cval);
    buffer.push_back(data);}
 SP.putData(obskey,RVector(buffer));
 SC.addValue(obskey,RVector(buffer));
}


void SamplingsCorrect::addValue(const MCObsInfo& mcobs, const RVector& insamp)
{
 map<MCObsInfo,RVector >::iterator it=values.find(mcobs);
 if (it!=values.end()) 
    it->second=insamp;
 else{ 
    values.insert(make_pair(mcobs,insamp));}
}

bool SamplingsCorrect::getCorrect(const MCObsInfo& mcobs, RVector& result) const
{
 map<MCObsInfo,RVector >::const_iterator it=values.find(mcobs);
 if (it==values.end()){ result.clear(); return false;}
 result=it->second;
 return true;
}




// *******************************************************************


void testMCObsGetHandlerFake(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestMCObsGetHandlerFake")==0)
 return;

 cout << endl << "Starting TestMCObsGetHandlerFake"<<endl;
 srand (time(NULL));
 
 try{
 XMLHandler xmlr(xml_in,"TestMCObsGetHandlerFake");

 VEVCorrect VC;
 set<int> omit;
 bool hermitian=true;
 CorrCorrect CC;
 if (hermitian) CC.setHermitian();
 BinsCorrect BC;
 SamplingsCorrect SC;

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
    omit.clear();
    list<XMLHandler> xmlo=it->find("Omit");
    for (list<XMLHandler>::iterator tt=xmlo.begin();tt!=xmlo.end();tt++){
       int anomission;
       xmlread(*tt,"Omit",anomission,"Tester");
       omit.insert(anomission);}
    cout << it->output()<<endl;
    make_fake_correlator_file(mcens,fname,corr,tmin,tmax,
                              Aval,Bval,Cval,omit,CC);}


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
    omit.clear();
    list<XMLHandler> xmlo=it->find("Omit");
    for (list<XMLHandler>::iterator tt=xmlo.begin();tt!=xmlo.end();tt++){
       int anomission;
       xmlread(*tt,"Omit",anomission,"Tester");
       omit.insert(anomission);}
    cout << it->output()<<endl;
    make_fake_vev_file(mcens,fname,opinfo,Aval,Cval,omit,VC);}

    // make the fake bin files with possibly missing data

 list<XMLHandler> xmlb=xmlr.find("MakeFakeBinFile");
 for (list<XMLHandler>::iterator it=xmlb.begin();it!=xmlb.end();++it){
    MCBinsInfo bininfo(*it);
    string fname;
    xmlreadchild(*it,"FileName",fname);  
    double Aval,Cval;
    xmlreadchild(*it,"Aval",Aval);
    xmlreadchild(*it,"Cval",Cval);
    MCObsInfo obskey(*it);
    cout << it->output()<<endl;
    make_fake_bin_file(obskey,fname,bininfo,Aval,Cval,BC);}

    // make the fake sampling files

 list<XMLHandler> xmls=xmlr.find("MakeFakeSamplingFile");
 for (list<XMLHandler>::iterator it=xmls.begin();it!=xmls.end();++it){
    MCBinsInfo bininfo(*it);
    MCSamplingInfo sampinfo(*it);
    string fname;
    xmlreadchild(*it,"FileName",fname);  
    double Aval,Cval;
    xmlreadchild(*it,"Aval",Aval);
    xmlreadchild(*it,"Cval",Cval);
    MCObsInfo obskey(*it);
    cout << it->output()<<endl;
    make_fake_samplings_file(obskey,fname,bininfo,sampinfo,Aval,Cval,SC);}



 XMLHandler xmlrdtest(xmlr,"DoReadTests");
 MCBinsInfo bininfo(xmlrdtest);
 cout << "Bins Info: "<<bininfo.output()<<endl;
 MCSamplingInfo sampinfo(xmlrdtest);
 cout << "Sampling Info: "<<sampinfo.output()<<endl;

 MCObsGetHandler MCOH(xmlrdtest,bininfo,sampinfo);   //  <=====  the main object being tested!!!

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

 XMLHandler xmloplist(xmlrdtest,"OpList");
 list<XMLHandler> xmlop=xmloplist.find("OperatorString");
 list<OperatorInfo> oplist;
 for (list<XMLHandler>::iterator it=xmlop.begin();it!=xmlop.end();++it)
    oplist.push_back(OperatorInfo(*it));

 XMLHandler xmlblist(xmlrdtest,"BinObsList");
 list<XMLHandler> xmlbm=xmlblist.find("MCObservable");
 list<MCObsInfo> obslist;
 for (list<XMLHandler>::iterator it=xmlbm.begin();it!=xmlbm.end();++it)
    obslist.push_back(MCObsInfo(*it));

 XMLHandler xmlslist(xmlrdtest,"SamplingObsList");
 list<XMLHandler> xmlsm=xmlslist.find("MCObservable");
 list<MCObsInfo> obsslist;
 for (list<XMLHandler>::iterator it=xmlsm.begin();it!=xmlsm.end();++it)
    obsslist.push_back(MCObsInfo(*it));

#ifdef COMPLEXNUMBERS
 int kkmax=2;
#else
 int kkmax=1;
#endif

 cout << "VEV tests"<<endl;
 int goodq=0; int badq=0;
 int goodc=0; int badc=0;
 double corvalue;
 Scalar getvalue;
 for (list<OperatorInfo>::iterator it=oplist.begin();it!=oplist.end();++it){
    MCObsInfo mcobs(*it);
    for (int kk=0;kk<kkmax;++kk){
       if (kk>0) mcobs.setToImaginaryPart();
       for (int serind=0;serind<=800;serind++){
          bool cflag=VC.getCorrect(mcobs,serind,corvalue);
          if (MCOH.queryData(mcobs,serind)==cflag) goodq++;
          else{ cout << "ERROR: VEV real serind query result incorrect"<<endl; badq++;}
          if (cflag){
             MCOH.getData(mcobs,serind,getvalue);
             double gval=(kk==0)?real(getvalue):imag(getvalue);
             if (compare_floats(gval,corvalue,eps)) goodc++;
             else {cout << "vev Mismatch! correctvalue = "<<corvalue
                   <<"  getvalue  = "<<gval<<endl; badc++;}}}
       RVector vbincorrect;
       bool bflag=VC.getCorrectBins(mcobs,bininfo,vbincorrect);
       if (MCOH.queryBins(mcobs)==bflag) goodq++;
       else{ cout << "ERROR: VEV bins query result incorrect"<<endl; badq++;}
       if (bflag){
           RVector vbins;
           MCOH.getBins(mcobs,vbins);
           if (compare_vectors(vbins,vbincorrect,eps)) goodc++;
           else {cout << "VEV bins Mismatch! "<<endl; badc++;}}}
#ifdef COMPLEXNUMBERS
    if (MCOH.queryBins(mcobs)){
       mcobs.setToRealPart(); RVector vbinre1, vbinim1, vbinre2, vbinim2;
       MCOH.getBins(mcobs,vbinre1);
       mcobs.setToImaginaryPart();
       MCOH.getBins(mcobs,vbinim1);
       MCOH.getBinsComplex(mcobs,vbinre2,vbinim2);
       bool check=compare_vectors(vbinre1,vbinre2,eps) 
               && compare_vectors(vbinim1,vbinim2,eps);
       if (check) goodc++;
       else {cout << "VEV bins Mismatch! "<<endl; badc++;}}
#endif
    }
 cout << "Number of successful vev queries = "<<goodq<<endl;
 cout << "Number of successful vev comparisons = "<<goodc<<endl;
 cout << "Number of incorrect vev queries = "<<badq<<endl;
 cout << "Number of incorrect vev comparisons = "<<badc<<endl;


 cout << endl<<"Correlator tests"<<endl;
 goodq=0; badq=0;
 goodc=0; badc=0;
 for (list<OperatorInfo>::iterator it1=oplist.begin();it1!=oplist.end();++it1){
 for (list<OperatorInfo>::iterator it2=oplist.begin();it2!=oplist.end();++it2){
    CorrelatorInfo cf(*it1,*it2);
    for (int tind=0;tind<=32;tind++){
       CorrelatorAtTimeInfo ct(cf,tind,hermitian);
       MCObsInfo mcobs(ct);
       for (int kk=0;kk<kkmax;++kk){
          if (kk>0) mcobs.setToImaginaryPart();
          for (int serind=0;serind<=800;serind++){
             bool cflag=CC.getCorrect(mcobs,serind,corvalue);
             if (MCOH.queryData(mcobs,serind)==cflag) goodq++;
             else{ cout << "ERROR: corr serind query result incorrect"<<endl; badq++;}
             if (cflag){
                MCOH.getData(mcobs,serind,getvalue);
                double gval=(kk==0)?real(getvalue):imag(getvalue);
                if (compare_floats(gval,corvalue,ceps)) goodc++;
                else {cout << "corr Mismatch! correctvalue = "<<corvalue
                      <<"  getvalue  = "<<gval<<endl; badc++;}}}
       RVector cbincorrect;
       bool cflag=CC.getCorrectBins(mcobs,bininfo,cbincorrect);
       if (MCOH.queryBins(mcobs)==cflag) goodq++;
       else{ cout << "ERROR: CORR bins query result incorrect"<<endl; badq++;}
       if (cflag){
           RVector cbins;
           MCOH.getBins(mcobs,cbins);
           if (compare_vectors(cbins,cbincorrect,ceps)) goodc++;
           else {cout << "CORR bins Mismatch! "<<endl; badc++;}}}
#ifdef COMPLEXNUMBERS
    if (MCOH.queryBins(mcobs)){
       mcobs.setToRealPart(); RVector cbinre1, cbinim1, cbinre2, cbinim2;
       MCOH.getBins(mcobs,cbinre1);
       mcobs.setToImaginaryPart();
       MCOH.getBins(mcobs,cbinim1);
       MCOH.getBinsComplex(mcobs,cbinre2,cbinim2);
       bool check=compare_vectors(cbinre1,cbinre2,ceps) 
               && compare_vectors(cbinim1,cbinim2,ceps);
       if (check) goodc++;
       else {cout << "CORR bins Mismatch! "<<endl; badc++;}}
#endif
    }}}
 cout << "Number of successful corr queries = "<<goodq<<endl;
 cout << "Number of successful corr comparisons = "<<goodc<<endl;
 cout << "Number of incorrect corr queries = "<<badq<<endl;
 cout << "Number of incorrect corr comparisons = "<<badc<<endl;


 cout <<endl<< "BIN obs tests"<<endl;
 goodq=0; badq=0;
 goodc=0; badc=0;
 for (list<MCObsInfo>::iterator it=obslist.begin();it!=obslist.end();++it){
    MCObsInfo mcobs(*it);
    RVector bincorrect;
    bool cflag=BC.getCorrect(mcobs,bincorrect);
    if (MCOH.queryBins(mcobs)==cflag) goodq++;
    else{ cout << "ERROR: BIN obs query result incorrect"<<endl; badq++;}
    if (cflag){
       RVector bins;
       MCOH.getBins(mcobs,bins);
       if (compare_vectors(bins,bincorrect,deps)) goodc++;
       else {cout << "bin Mismatch! "<<endl; badc++;}}}
 cout << "Number of successful bin queries = "<<goodq<<endl;
 cout << "Number of successful bin comparisons = "<<goodc<<endl;
 cout << "Number of incorrect bin queries = "<<badq<<endl;
 cout << "Number of incorrect bin comparisons = "<<badc<<endl;


 cout <<endl<< "SAMPLING obs tests"<<endl;
 goodq=0; badq=0;
 goodc=0; badc=0;
 for (list<MCObsInfo>::iterator it=obsslist.begin();it!=obsslist.end();++it){
    MCObsInfo mcobs(*it);
    RVector sampcorrect;
    bool cflag=SC.getCorrect(mcobs,sampcorrect);
    if (MCOH.querySamplings(mcobs)==cflag) goodq++;
    else{ cout << "ERROR: SAMPLING obs query result incorrect"<<endl; badq++;}
    if (cflag){
       RVector samp;
       MCOH.getSamplings(mcobs,samp);
       if (compare_vectors(samp,sampcorrect,deps)) goodc++;
       else {cout << "bin Mismatch! "<<endl; badc++;}}}
 cout << "Number of successful bin queries = "<<goodq<<endl;
 cout << "Number of successful bin comparisons = "<<goodc<<endl;
 cout << "Number of incorrect bin queries = "<<badq<<endl;
 cout << "Number of incorrect bin comparisons = "<<badc<<endl;


 }
 catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}
}




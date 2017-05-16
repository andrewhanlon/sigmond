#include "mcobs_handler.h"
#include "test_obs_get_handler.h"
#include "test_mcobs_handler.h"
#include "correlator_matrix_info.h"
#include "task_utils.h"
#include <cstdio>
#include <ctime>
#include <map>
#include <cstdlib>

using namespace std;
using namespace LaphEnv;


  // *************************************************************************

     //  This tests reading a correlation matrix at one time value


void testCorrMatEstimates(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestCorrelatorMatrixEstimates")==0)
 return;

 cout << endl << "Starting TestCorrelatorMatrixEstimates"<<endl;
  /* initialize random seed: */
 srand (time(NULL));

 
 try{
 XMLHandler xmlr(xml_in,"TestCorrelatorMatrixEstimates");

   //  get the correlator matrix info

 CorrelatorMatrixInfo cormatinfo(xmlr);
 cout << cormatinfo.output(false)<<endl;
 bool herm=cormatinfo.isHermitian();
 if (!herm)
    throw(std::invalid_argument("Test required Hermitian Correlator Matrix"));
 bool vev=cormatinfo.subtractVEV();

   //   get the omissions and rebin factor

 set<int> omit;
 int rebin=1;
 XMLHandler xmlc(xmlr,"ConfigManip");
 xmlreadifchild(xmlc,"Rebin",rebin);
 list<XMLHandler> xmlo=xmlc.find("Omit");
 for (list<XMLHandler>::iterator tt=xmlo.begin();tt!=xmlo.end();tt++){
    int anomission;
    xmlread(*tt,"Omit",anomission,"Tester");
    omit.insert(anomission);}


    // make the fake data files with possibly missing data

 XMLHandler xmlf(xmlr,"MakeFakeCorrelatorFiles");
 MCBinsInfo binsinfo(xmlf);  
 MCEnsembleInfo mcens(binsinfo.getMCEnsembleInfo());
 string corrstub,vevstub("dummy");
 xmlreadchild(xmlf,"CorrFileStub",corrstub);
 int tmin,tmax;
 xmlreadchild(xmlf,"MinTime",tmin);
 xmlreadchild(xmlf,"MaxTime",tmax);
 double minenergy,stepenergy,Zavg,Zspread,vevavg,vevspread;
 xmlreadchild(xmlf,"MinEnergy",minenergy);
 xmlreadchild(xmlf,"StepEnergy",stepenergy);
 xmlreadchild(xmlf,"Zavg",Zavg);
 xmlreadchild(xmlf,"Zspread",Zspread);
 if (vev){
    xmlreadchild(xmlf,"VEVFileStub",vevstub);
    xmlreadchild(xmlf,"VEVavg",vevavg);
    xmlreadchild(xmlf,"VEVspread",vevspread);}

 cout << mcens.output()<<endl;
 cout << "Correlator File stub  = "<<corrstub<<endl;
 if (vev) cout << "VEV file stub = "<<vevstub<<endl;
 cout << "min time = "<<tmin<<endl;
 cout << "max time = "<<tmax<<endl;
 cout << "MinEnergy = "<<minenergy<<endl;
 cout << "StepEnergy = "<<stepenergy<<endl;
 cout << "Zavg = "<<Zavg<<endl;
 cout << "Zspread = "<<Zspread<<endl;
 if (vev){
    cout << "VEVavg = "<< vevavg<<endl;
    cout << "VEVspread = "<< vevspread<<endl;}

 int nops=cormatinfo.getNumberOfOperators();
 FileListInfo finfocor(corrstub,0,nops*nops-1);
 CorrCorrect CMcorrect;
 if (herm) CMcorrect.setHermitian();
 const set<OperatorInfo>& opset=cormatinfo.getOperators();
 const std::set<int> missing;

 cout <<endl<< "Making fake correlator files"<<endl;
 int suffix=0;
 int row=0;
 int col;
 for (set<OperatorInfo>::const_iterator itrow=opset.begin();itrow!=opset.end();itrow++,row++){
    col=0;
    for (set<OperatorInfo>::const_iterator itcol=opset.begin();itcol!=opset.end();itcol++,col++,suffix++){
       string fname(finfocor.getFileName(suffix));
       cout << "row  = "<<row<<" col = "<<col<<" file = "<<fname<<endl;
       CorrelatorInfo corr(*itrow,*itcol);
       cout << corr.output(false)<<endl;
       double Bval=minenergy+suffix*stepenergy;
       double Aval=Zavg+double((rand()%4096)-2048)/2048.0*Zspread;
       double Cval=double((rand()%4096)-2048)/4096.0*Zspread;
       cout << "Aval = "<<Aval<<endl;
       cout << "Bval = "<<Bval<<endl;
       cout << "Cval = "<<Cval<<endl;
       make_fake_correlator_file(mcens,fname,corr,tmin,tmax,
                                 Aval,Bval,Cval,missing,CMcorrect);}}

 VEVCorrect Vcorrect;
 FileListInfo finfovev(vevstub,0,nops-1);
 if (vev){

    cout <<endl<< "Making fake VEV files"<<endl;
    int suffix=0;
    for (set<OperatorInfo>::const_iterator it=opset.begin();it!=opset.end();it++,suffix++){
       string fname(finfovev.getFileName(suffix));
       cout << " vev file = "<<fname<<endl;
       cout << it->output(false)<<endl;
       double Aval=Zavg+double((rand()%4096)-2048)/2048.0*Zspread;
       double Cval=double((rand()%4096)-2048)/4096.0*Zspread;
       cout << "Aval = "<<Aval<<endl;
       cout << "Cval = "<<Cval<<endl;
       make_fake_vev_file(mcens,fname,*it,Aval,Cval,missing,Vcorrect);}}

    //  fake data made, now read the correlator matrix elements and check they are correct

 XMLHandler xmlrdtest("MCObservables");
 XMLHandler xmltmp,xmltmp2;
 mcens.output(xmltmp);
 xmlrdtest.put_child(xmltmp);
 finfocor.output(xmltmp);
 xmltmp2.set_root("BLCorrelatorData");
 xmltmp2.put_child(xmltmp);
 xmlrdtest.put_child(xmltmp2);
 if (vev){
    finfovev.output(xmltmp);
    xmltmp2.set_root("BLVEVData");
    xmltmp2.put_child(xmltmp);
    xmlrdtest.put_child(xmltmp2);}

 cout <<"*************************************"<<endl;
 cout << xmlrdtest.output()<<endl;
 MCSamplingInfo sampinfo(578,332,14);

 MCObsGetHandler MCOH(xmlrdtest,binsinfo,sampinfo); 

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

 MCObsHandler MC(MCOH);
 uint nconfig=MC.getNumberOfMeasurements();
 uint nbins=MC.getNumberOfBins();
 BinMapper BM(nconfig,rebin,omit);

 cout << "number of configs = "<<nconfig<<endl;
 cout << "number of bins = "<<nbins<<endl;

 MC.setToBootstrapMode();
 MC.setSamplingBegin();

#ifdef COMPLEXNUMBERS
 CVector vevvals(nops);
 for (int k=0;k<nops;k++) vevvals[k]=complex<double>(0.0,0.0);
#else
 RVector vevvals(nops);
 for (int k=0;k<nops;k++) vevvals[k]=0.0;
#endif
 if (vev){
    int index=0;
    for (set<OperatorInfo>::const_iterator it=opset.begin();it!=opset.end();it++,index++){
       MCObsInfo mcobs(*it);
#ifdef COMPLEXNUMBERS
       int kkmax=2;
#else
       int kkmax=1;
#endif
       vector<double> res(kkmax);
       for (int kk=0;kk<kkmax;++kk){
       if (kk>0) mcobs.setToImaginaryPart();
       Vector<double> vevbins(nbins);
       RVector vevbinread(MC.getBins(mcobs));
       for (uint bin=0;bin<nbins;bin++){
          double val;
          if (getVEVBin(Vcorrect,BM,mcobs,bin,val)){
             vevbins[bin]=val;
             if (std::abs(vevbinread[bin]-val)>1e-5)
                cout << "VEV bin is wrong: bin = "<<bin<<" read "
                     <<vevbinread[bin]<<"  correct = "<<val<<endl;}
          else break;}
       MCObsInfo obskeytemp("VEVValue",kk,true);
       MC.putBins(obskeytemp,vevbins);
       res[kk]=MC.getCurrentSamplingValue(obskeytemp);
       MC.eraseData(obskeytemp);}
#ifdef COMPLEXNUMBERS
       vevvals[index]=complex<double>(res[0],res[1]);
#else
       vevvals[index]=res[0];
#endif
       }}

#ifdef COMPLEXNUMBERS
 CMatrix cormat_estimates;
 ComplexHermitianMatrix corherm_estimates;
 CVector vev_estimates;
 complex<double> correct,cres;
 map<uint,double > all_correct_re,all_correct_im;
 map<uint,double>::const_iterator allit;
#else
 RMatrix cormat_estimates;
 RealSymmetricMatrix corherm_estimates;
 RVector vev_estimates;
 double correct,cres;
 map<uint,double> all_correct_re;
 map<uint,double>::const_iterator allit;
#endif

 bool erase_corr_data_after=true;

 if (vev) getHermCorrelatorMatrixVEVs_CurrentSampling(&MC,cormatinfo,vev_estimates);

 for (int timeval=tmin;timeval<=tmax;timeval++){
    getHermCorrelatorMatrixAtTime_CurrentSampling(&MC,cormatinfo,timeval,corherm_estimates);
    if (erase_corr_data_after) eraseHermCorrelatorMatrixAtTime(&MC,cormatinfo,timeval);

       // check if erasing worked

    cout << "Query first VEV real part = "<<MC.queryBinsInMemory(MCObsInfo(*(opset.begin()),RealPart))<<endl;
    cout << "Query first VEV imag part = "<<MC.queryBinsInMemory(MCObsInfo(*(opset.begin()),ImaginaryPart))<<endl;
    CorrelatorAtTimeInfo c00info(*(opset.begin()),*(opset.begin()),timeval,herm);
    cout << "Query C(0,0)[tval] real part = "<<MC.queryBinsInMemory(MCObsInfo(c00info,RealPart))<<endl;
    cout << "Query C(0,0)[tval] imag part = "<<MC.queryBinsInMemory(MCObsInfo(c00info,ImaginaryPart))<<endl;

    if (vev){
       cout << "vev_estimates size = "<<vev_estimates.size()<<endl;
       row=0;
       for (set<OperatorInfo>::const_iterator itrow=opset.begin();itrow!=opset.end();itrow++,row++){
          correct=vevvals[row];
          cres=vev_estimates[row];
          cout << "vev["<<row<<"]:  "<<cres;
          cout << "   correct = "<<correct;
          double diff=std::abs(cres-correct);
          if (std::abs(correct)>1.0) diff/=std::abs(correct);
          cout << "   diff = "<<diff;
          if (diff>1e-6) cout << "  MISMATCH";
          cout << endl;}}

    row=0;
    for (set<OperatorInfo>::const_iterator itrow=opset.begin();itrow!=opset.end();itrow++,row++){
       col=0;
       for (set<OperatorInfo>::const_iterator itcol=opset.begin();itcol!=opset.end();itcol++,col++){
       CorrelatorAtTimeInfo ctinfo(*itrow,*itcol,timeval,herm);
       MCObsInfo mcobs(ctinfo);
#ifdef COMPLEXNUMBERS
       int kkmax=2;
#else
       int kkmax=1;
#endif
       vector<double> res(kkmax);
       for (int kk=0;kk<kkmax;++kk){
       if (kk>0) mcobs.setToImaginaryPart();
       Vector<double> corbins(nbins);
       RVector corbinread(MC.getBins(mcobs));
       for (uint bin=0;bin<nbins;bin++){
          double val;
          if (getCorBin(CMcorrect,BM,mcobs,bin,val)){
             corbins[bin]=val;
             if (std::abs(corbinread[bin]-val)>1e-5)
                cout << "CORR bin is wrong: bin = "<<bin<<" read "
                     <<corbinread[bin]<<"  correct = "<<val<<endl;}
          else break;}
       MCObsInfo obskeytemp("CorrValue",kk,true);
       MC.putBins(obskeytemp,corbins);
       res[kk]=MC.getCurrentSamplingValue(obskeytemp);
       MC.eraseData(obskeytemp);}
#ifdef COMPLEXNUMBERS
       correct=complex<double>(res[0],res[1])-vevvals[row]*conjugate(vevvals[col]);
       all_correct_re.insert(make_pair(row*100000+col*1000+timeval,correct.real()));
       all_correct_im.insert(make_pair(row*100000+col*1000+timeval,correct.imag()));
       if ((herm)&&(row==col)) correct=complex<double>(correct.real(),0.0);
#else
       correct=res[0]-vevvals[row]*vevvals[col];
       all_correct_re.insert(make_pair(row*100000+col*1000+timeval,correct));
#endif
       if (herm) cres=corherm_estimates(row,col);
       else cres=cormat_estimates(row,col);
       cout << "("<<row<<","<<col<<")["<<timeval<<"]:  "<<cres;
       cout << "   correct = "<<correct;
       double diff=std::abs(cres-correct);
       if (std::abs(correct)>1.0) diff/=std::abs(correct);
       cout << "   diff = "<<diff;
       if (diff>1e-6) cout << "  MISMATCH";
       cout << endl;
    }}}

 double correctre;
 cout << endl<<endl;
 row=0;
 for (set<OperatorInfo>::const_iterator itrow=opset.begin();itrow!=opset.end();itrow++,row++){
    col=0;
    for (set<OperatorInfo>::const_iterator itcol=opset.begin();itcol!=opset.end();itcol++,col++){
       map<int,MCEstimate> results_re,results_im;
       getCorrelatorEstimates(&MC,CorrelatorInfo(*itrow,*itcol),herm,vev,RealPart, 
                              Bootstrap,results_re);
       for (map<int,MCEstimate>::const_iterator it=results_re.begin();it!=results_re.end();it++){
          double res=it->second.getFullEstimate();
          cout << "Re corr("<<row<<","<<col<<")[t="<<it->first<<"] = "<<res;
          allit=all_correct_re.find(row*100000+col*1000+it->first);
          if (allit!=all_correct_re.end()) correctre=allit->second;
          else throw(std::invalid_argument("Oops, something went wrong"));
          cout << "  correct = "<<correctre;
          double diff=std::abs(res-correctre);
          if (std::abs(correctre)>1.0) diff/=std::abs(correctre);
          cout << "   diff = "<<diff;
          if (diff>1e-6) cout << "  MISMATCH";
          cout << endl;}
#ifdef COMPLEXNUMBERS
       getCorrelatorEstimates(&MC,CorrelatorInfo(*itrow,*itcol),herm,vev,ImaginaryPart, 
                              Bootstrap,results_im);
       for (map<int,MCEstimate>::const_iterator it=results_im.begin();it!=results_im.end();it++){
          double res=it->second.getFullEstimate();
          cout << "Im corr("<<row<<","<<col<<")[t="<<it->first<<"] = "<<res;
          allit=all_correct_im.find(row*100000+col*1000+it->first);
          if (allit!=all_correct_im.end()) correctre=allit->second;
          else throw(std::invalid_argument("Oops, something went wrong"));
           cout << "  correct = "<<correctre;
          double diff=std::abs(res-correctre);
          if (std::abs(correctre)>1.0) diff/=std::abs(correctre);
          cout << "   diff = "<<diff;
          if (diff>1e-6) cout << "  MISMATCH";
          cout << endl;}
#endif
    }}

/*
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

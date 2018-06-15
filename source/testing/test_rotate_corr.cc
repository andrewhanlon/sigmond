#include <cstdio>
#include <ctime>
#include <algorithm>
#include "xml_handler.h"
#include "mcobs_handler.h"
#include "single_pivot.h"
#include "correlator_matrix_info.h"
#include "args_handler.h"

using namespace std;
using namespace LaphEnv;

double get_random_minusone_to_one()
{
 return double(int((rand() % 4096))-2048)/2048.0;
}


void read_improved_operators(ArgsHandler& xmlin, 
              const CorrelatorMatrixInfo& cormat_info,
              map<OperatorInfo,map<OperatorInfo,Scalar> >& trans)
{
 trans.clear();
 ArgsHandler xmli(xmlin,"ImprovedOperators");
 list<ArgsHandler> xmlo(xmli.getSubHandlers("ImprovedOperator"));
 for (list<ArgsHandler>::iterator it=xmlo.begin();it!=xmlo.end();it++){
    ArgsHandler xmlopname(*it,"OpName");
    OperatorInfo opimp;
    xmlopname.getItem("OpName",opimp);
    if (trans.find(opimp)!=trans.end())
       throw(std::invalid_argument("Repeated improved operator in SinglePivot"));
    it->insert(xmlopname);
    map<OperatorInfo,Scalar> opdef;
    list<ArgsHandler> xmlc(it->getSubHandlers("OpTerm"));
    for (list<ArgsHandler>::iterator ct=xmlc.begin();ct!=xmlc.end();ct++){
       OperatorInfo opterm;
       ct->getItem("OpTerm",opterm);
       if (opdef.find(opimp)!=opdef.end())
          throw(std::invalid_argument("Repeated operator term in ImprovedOperator in SinglePivot"));
       Scalar cf;
       ct->getScalar("Coefficient",cf);
       opdef.insert(make_pair(opterm,cf));
      it->insert(*ct);}
    trans.insert(make_pair(opimp,opdef));
    xmli.insert(*it);}

 XMLHandler xmlecho;
 xmli.echo(xmlecho);
 cout << "ReadImprovedOperators:"<<endl;
 //cout << xmlecho.output()<<endl;
 xmlin.insert(xmli);
}


void setup_improved_operators(const CorrelatorMatrixInfo& cormat_info,
                const std::map<OperatorInfo,
                std::map<OperatorInfo,Scalar> >& trans,
                CorrelatorMatrixInfo *& orig_cormat_info,
                bool subvev, TransMatrix& imp_transmat)
{
 delete orig_cormat_info;  // make sure set to zero first
 orig_cormat_info=0;
/* set<OperatorInfo> origset;
 for (map<OperatorInfo,map<OperatorInfo,Scalar> >::const_iterator 
     it=trans.begin();it!=trans.end();it++){
    for (map<OperatorInfo,Scalar>::const_iterator 
       ct=(it->second).begin();ct!=(it->second).end();ct++)
          origset.insert(ct->first);}
 orig_cormat_info=new CorrelatorMatrixInfo(origset,true,subvev);
 map<OperatorInfo,uint> origind;
 uint count=0;
 for (set<OperatorInfo>::iterator ot=origset.begin();ot!=origset.end();ot++,count++)
    origind.insert(make_pair(*ot,count));
 uint norig=origset.size();
 uint nops=trans.size();
 imp_transmat.resize(norig,nops);
 count=0;
 for (map<OperatorInfo,map<OperatorInfo,Scalar> >::const_iterator
      kt=trans.begin();kt!=trans.end();kt++,count++){
    for (map<OperatorInfo,Scalar>::const_iterator 
         ct=(kt->second).begin();ct!=(kt->second).end();ct++){
       Scalar cf=ct->second;
       const OperatorInfo& opt=ct->first;
       map<OperatorInfo,uint>::const_iterator ut=origind.find(opt);
       if (ut==origind.end())
          throw(std::invalid_argument("Something went wrong with indexing in setup_improved_operators"));
       uint oind=ut->second;
       imp_transmat(oind,count)=cf;}}
*/


 set<OperatorInfo> addset=cormat_info.getOperators();
    //  from improved operators, get original set of operators
 set<OperatorInfo> origset;
 for (map<OperatorInfo,map<OperatorInfo,Scalar> >::const_iterator 
     it=trans.begin();it!=trans.end();it++){
    if (addset.erase(it->first)==0)     // remove improved operators from "addset"
       throw(std::invalid_argument("Improved operator not in CorrelatorMatrixInfo"));
    for (map<OperatorInfo,Scalar>::const_iterator 
       ct=(it->second).begin();ct!=(it->second).end();ct++)
          origset.insert(ct->first);}
   // operators remaining in "addset" must be original operators (not improved)
 origset.insert(addset.begin(),addset.end());

 orig_cormat_info=new CorrelatorMatrixInfo(origset,true,cormat_info.subtractVEV());
 map<OperatorInfo,uint> origind;
 uint count=0;
 for (set<OperatorInfo>::iterator ot=origset.begin();ot!=origset.end();ot++,count++)
    origind.insert(make_pair(*ot,count));
 const set<OperatorInfo>& opset=cormat_info.getOperators();
 uint norig=origset.size();
 uint nops=opset.size();
 imp_transmat.resize(norig,nops);
 count=0;
 for (set<OperatorInfo>::const_iterator it=opset.begin();it!=opset.end();it++,count++){
    map<OperatorInfo,map<OperatorInfo,Scalar> >::const_iterator kt=trans.find(*it);
    if (kt==trans.end()){
       map<OperatorInfo,uint>::const_iterator ut=origind.find(*it);
       if (ut==origind.end())
          throw(std::invalid_argument("Something went wrong with indexing in setup_improved_operators"));
       uint oind=ut->second;
       imp_transmat(oind,count)=1.0;}
    else{
       for (map<OperatorInfo,Scalar>::const_iterator 
            ct=(kt->second).begin();ct!=(kt->second).end();ct++){
          Scalar cf=ct->second;
          const OperatorInfo& opt=ct->first;
          map<OperatorInfo,uint>::const_iterator ut=origind.find(opt);
          if (ut==origind.end())
             throw(std::invalid_argument("Something went wrong with indexing in setup_improved_operators"));
          uint oind=ut->second;
          imp_transmat(oind,count)=cf;}}}
}


#if defined COMPLEXNUMBERS

void make_fake_data_to_rotate(MCObsHandler& MC, CorrelatorMatrixInfo& cormat,
                              const string& maplefile, const TransMatrix& imp_transmat,
                              const string& binfile, uint tmin, uint tmax,
                              uint tauN, uint tau0, uint tauD, double mininvcondnum,
                              uint timeval, uint binindex)
{
 uint nops=cormat.getNumberOfOperators();
 bool vevs=cormat.subtractVEV();
 uint nbins=MC.getNumberOfBins();
 srand (time(NULL));

 set<MCObsInfo> obskeys;
 Vector<double> bins_re(nbins);
 Vector<double> bins_im(nbins);
 uint nlevels=nops;
 CMatrix Zcoef(nops,nlevels);
 vector<double> energies(nlevels);
 energies[0]=0.15;
 cout << "energies[0] := "<<energies[0]<<":"<<endl;
 for (uint level=1;level<nlevels;level++){
    energies[level]=energies[level-1]+0.08/double(level);
    cout << "energies["<<level<<"] := "<<energies[level]<<":"<<endl;}
 for (uint level=0;level<nlevels;level++){
    for (uint k=0;k<nops;k++){
       Zcoef(k,level)=complex<double>(10.0*get_random_minusone_to_one(),
                                      6.3*get_random_minusone_to_one());}}
 CVector Vevs(nops);
 for (uint k=0;k<nops;k++)
    Vevs[k]=complex<double>((k+3)*get_random_minusone_to_one(),(k+1)*get_random_minusone_to_one());

 ofstream fout(maplefile.c_str());
 fout.precision(16);
 if (vevs) fout << "VEVS:=true:"<<endl;
 else fout << "VEVS:=false:"<<endl;
 fout << "tmin:="<<tmin<<":"<<endl;
 fout << "tmax:="<<tmax<<":"<<endl<<endl;
 fout << "tauN:="<<tauN<<":"<<endl;
 fout << "tau0:="<<tau0<<":"<<endl;
 fout << "tauD:="<<tauD<<":"<<endl;
 fout << "mininvcondnum:="<<mininvcondnum<<":"<<endl;
 fout << "atimeval:="<<timeval<<":"<<endl;
 fout << "bin:="<<binindex<<":"<<endl;

 const set<OperatorInfo>& opset=cormat.getOperators();

 if (imp_transmat.size(0)==0){
    fout << "Nops:="<<nops<<":"<<endl;
    fout << "Norigops:="<<nops<<":"<<endl;}
 else{
    uint nkeep=imp_transmat.size(1);
    if (imp_transmat.size(0)!=nops)
       throw(std::invalid_argument("Bad transformation matrix"));
    fout << "Nops:="<<nkeep<<":"<<endl;
    fout << "Norigops:="<<nops<<":"<<endl;
    for (uint row=0;row<nops;row++)
    for (uint col=0;col<nkeep;col++)
       fout << "Trans["<<row+1<<","<<col+1<<"]:=Complex("
            <<imp_transmat(row,col).real()
            << ","<<imp_transmat(row,col).imag()<<"):"<<endl;}

 set<OperatorInfo>::const_iterator itrow,itcol;
 uint row=0;
 for (itrow=opset.begin();itrow!=opset.end();itrow++,row++){
    uint col=row;
    for (itcol=itrow;itcol!=opset.end();itcol++,col++){
       CorrelatorAtTimeInfo corrt(*itrow,*itcol,0,true,false);
       for (uint timeval=tmin;timeval<=tmax;timeval++){
          corrt.resetTimeSeparation(timeval);
          MCObsInfo obskey(corrt,RealPart);
          complex<double> ztemp(0.0,0.0);
          for (uint level=0;level<nlevels;level++)
             ztemp+=Zcoef(row,level)*conjugate(Zcoef(col,level))*exp(-energies[level]*timeval);
          if (vevs) ztemp+=Vevs[row]*conjugate(Vevs[col]);
          fout << "CorrData[["<<timeval<<","<<row<<","<<col<<"]]:=["<<endl;
          for (uint bin=0;bin<nbins;bin++){
             bins_re[bin]=ztemp.real()*(1.0+0.001*get_random_minusone_to_one());
             bins_im[bin]=ztemp.imag()*(1.0+0.001*get_random_minusone_to_one());
             fout <<"Complex("<<bins_re[bin]<<","<<bins_im[bin]<<"),"<<endl;}
          fout <<"NULL]:"<<endl;
          obskeys.insert(obskey);
          MC.putBins(obskey,bins_re);
          obskey.setToImaginaryPart();
          obskeys.insert(obskey);
          MC.putBins(obskey,bins_im);}}}
 if (vevs){
    uint row=0;
    for (itrow=opset.begin();itrow!=opset.end();itrow++,row++){
       MCObsInfo obskey(*itrow,RealPart);
       complex<double> ztemp=Vevs[row];
       fout << "VEVData["<<row<<"]:=["<<endl;
       for (uint bin=0;bin<nbins;bin++){
          bins_re[bin]=ztemp.real()*(1.0+0.001*get_random_minusone_to_one());
          bins_im[bin]=ztemp.imag()*(1.0+0.001*get_random_minusone_to_one());
          fout <<bins_re[bin]<<"+I*("<<bins_im[bin]<<"),"<<endl;}
       fout <<"NULL]:"<<endl;
       obskeys.insert(obskey);
       MC.putBins(obskey,bins_re);
       obskey.setToImaginaryPart();
       obskeys.insert(obskey);
       MC.putBins(obskey,bins_im);}}
 fout.close();
 XMLHandler xmlf;
 MC.writeBinsToFile(obskeys,binfile,xmlf,true);

 if (imp_transmat.size(0)==0){
    for (uint k=0;k<nops;k++){
       cout << endl;
       double rescale=0.0;
       for (uint level=0;level<nlevels;level++) 
          rescale+=std::pow(std::abs(Zcoef(k,level)),2);
       rescale=1.0/rescale;
       for (uint level=0;level<nlevels;level++){
          cout << "Zcoef["<<k<<","<<level<<"]:=Complex("<<Zcoef(k,level).real()
                                          <<","<<Zcoef(k,level).imag()<<"):  ";
          cout << "ZcoefMagSq["<<k<<","<<level<<"]:="<<rescale*std::pow(abs(Zcoef(k,level)),2)<<":"<<endl;}}}
 else{
    uint nkeep=imp_transmat.size(1);
    vector<complex<double> > zcf(nlevels);
    for (uint imp=0;imp<nkeep;imp++){
       double rescale=0.0;
       for (uint level=0;level<nlevels;level++){
          complex<double> res(0.0,0.0);
          for (uint k=0;k<nops;k++){
             res+=conj(imp_transmat(k,imp))*Zcoef(k,level);}
          zcf[level]=res;
          rescale+=std::pow(std::abs(res),2);}
       cout << endl;
       rescale=1.0/rescale;
       for (uint level=0;level<nlevels;level++){
          cout << "Zcoef["<<imp<<","<<level<<"]:=Complex("<<zcf[level].real()
                                          <<","<<zcf[level].imag()<<"):  ";
          cout << "ZcoefMagSq["<<imp<<","<<level<<"]:="<<rescale*std::pow(abs(zcf[level]),2)<<":"<<endl;}}}
}

#else

void make_fake_data_to_rotate(MCObsHandler& MC, CorrelatorMatrixInfo& cormat,
                              const string& maplefile, const TransMatrix& imp_transmat,
                              const string& binfile, uint tmin, uint tmax,
                              uint tauN, uint tau0, uint tauD, double mininvcondnum,
                              uint timeval, uint binindex)
{
 uint nops=cormat.getNumberOfOperators();
 bool vevs=cormat.subtractVEV();
 uint nbins=MC.getNumberOfBins();
 srand (time(NULL));

 set<MCObsInfo> obskeys;
 Vector<double> bins_re(nbins);
 uint nlevels=nops;
 RMatrix Zcoef(nops,nlevels);
 vector<double> energies(nlevels);
 energies[0]=0.15;
 for (uint level=1;level<nlevels;level++){
    energies[level]=energies[level-1]+0.08/double(level);}
 for (uint level=0;level<nlevels;level++){
    for (uint k=0;k<nops;k++)
       Zcoef(k,level)=10.0*get_random_minusone_to_one();}
 RVector Vevs(nops);
 for (uint k=0;k<nops;k++)
    Vevs[k]=(k+3)*get_random_minusone_to_one();

 ofstream fout(maplefile.c_str());
 fout.precision(16);
 if (vevs) fout << "VEVS:=true:"<<endl;
 else fout << "VEVS:=false:"<<endl;
 fout << "tmin:="<<tmin<<":"<<endl;
 fout << "tmax:="<<tmax<<":"<<endl<<endl;
 fout << "tauN:="<<tauN<<":"<<endl;
 fout << "tau0:="<<tau0<<":"<<endl;
 fout << "tauD:="<<tauD<<":"<<endl;
 fout << "mininvcondnum:="<<mininvcondnum<<":"<<endl;
 fout << "atimeval:="<<timeval<<":"<<endl;
 fout << "bin:="<<binindex<<":"<<endl;

 if (imp_transmat.size(0)==0){
    fout << "Nops:="<<nops<<":"<<endl;
    fout << "Norigops:="<<nops<<":"<<endl;}
 else{
    uint nkeep=imp_transmat.size(1);
    if (imp_transmat.size(0)!=nops)
       throw(std::invalid_argument("Bad transformation matrix"));
    fout << "Nops:="<<nkeep<<":"<<endl;
    fout << "Norigops:="<<nops<<":"<<endl;
    for (uint row=0;row<nops;row++)
    for (uint col=0;col<nkeep;col++)
       fout << "Trans["<<row+1<<","<<col+1<<"]:="<<imp_transmat(row,col)<<":"<<endl;}

 const set<OperatorInfo>& opset=cormat.getOperators();
 set<OperatorInfo>::const_iterator itrow,itcol;
 uint row=0;
 for (itrow=opset.begin();itrow!=opset.end();itrow++,row++){
    uint col=row;
    for (itcol=itrow;itcol!=opset.end();itcol++,col++){
       CorrelatorAtTimeInfo corrt(*itrow,*itcol,0,true,false);
       for (uint timeval=tmin;timeval<=tmax;timeval++){
          corrt.resetTimeSeparation(timeval);
          MCObsInfo obskey(corrt,RealPart);
          double ztemp=0.0;
          for (uint level=0;level<nlevels;level++)
             ztemp+=Zcoef(row,level)*conjugate(Zcoef(col,level))*exp(-energies[level]*timeval);
          if (vevs) ztemp+=Vevs[row]*conjugate(Vevs[col]);
          fout << "CorrData[["<<timeval<<","<<row<<","<<col<<"]]:=["<<endl;
          for (uint bin=0;bin<nbins;bin++){
             bins_re[bin]=ztemp*(1.0+0.001*get_random_minusone_to_one());
             fout <<bins_re[bin]<<","<<endl;}
          fout <<"NULL]:"<<endl;
          obskeys.insert(obskey);
          MC.putBins(obskey,bins_re);}}}
 if (vevs){
    uint row=0;
    for (itrow=opset.begin();itrow!=opset.end();itrow++,row++){
       MCObsInfo obskey(*itrow,RealPart);
       double ztemp=Vevs[row];
       fout << "VEVData["<<row<<"]:=["<<endl;
       for (uint bin=0;bin<nbins;bin++){
          bins_re[bin]=ztemp*(1.0+0.001*get_random_minusone_to_one());
          fout <<bins_re[bin]<<","<<endl;}
       fout <<"NULL]:"<<endl;
       obskeys.insert(obskey);
       MC.putBins(obskey,bins_re);}}
 fout.close();
 XMLHandler xmlf;
 MC.writeBinsToFile(obskeys,binfile,xmlf,true);
}


#endif


void make_fit_energy_files(MCObsHandler& MC, uint nlevels, const string& sampfile,
                           const string& energycommonname, const vector<MCObsInfo>& energy_names)
{
 vector<MCObsInfo> ekeys;
 if ((energycommonname.empty())&&(energy_names.size()==nlevels)){
    ekeys=energy_names;}
 else if ((!energycommonname.empty())&&(energy_names.empty())){
    MCObsInfo ecommonkey(energycommonname,0);
    for (uint level=0;level<nlevels;level++){
       ecommonkey.resetObsIndex(level);
       ekeys.push_back(ecommonkey);}}
 else{
    throw(std::invalid_argument("could not make fit energy files"));}

 double fitenergy;
 for (uint level=0;level<nlevels;level++){
    fitenergy=0.853*(1.2+get_random_minusone_to_one());
    cout << "fitenergy["<<level<<"] = "<<fitenergy<<endl;
    for (MC.begin();!MC.end();++MC){
       double value=fitenergy+0.008*get_random_minusone_to_one();
       MC.putCurrentSamplingValue(ekeys[level],value,true);}}
 XMLHandler xmllog;
 set<MCObsInfo> ekeyset(ekeys.begin(),ekeys.end());
 MC.writeSamplingValuesToFile(ekeyset,sampfile,xmllog,true);
}



   //  on first run, use <Setup> to make fake data bins
   //   This creates input XML file for second run and third run
   //  on second run, use sigmond_batch and do the correlator
   //   bin rotations, save to a file
   //  on third run, use <Finish> to read the data bins
   //   and do the checks

void testRotateCorrelator(XMLHandler& xml_in, int taskcount)
{
 if (xml_tag_count(xml_in,"TestRotateCorrelator")==0)
 return;

 cout << endl << "Starting TestRotateCorrelator"<<endl;
 XMLHandler xmlq(xml_in,"TestRotateCorrelator");
 int stage=0;
 if (xmlq.count_among_children("Setup")==1) stage=1;
 else if (xmlq.count_among_children("Finish")==1) stage=2;
 else throw(std::invalid_argument("Invalid XML in testRotateCorrelator"));

 if (stage==1){
 try{

 XMLHandler xmlr(xmlq,"Setup");
 MCBinsInfo bins_info(xmlr);
 MCSamplingInfo samp_info(xmlr);
 MCObsGetHandler MCOH(xmlr,bins_info,samp_info); 
 MCObsHandler MC(MCOH);
 string maplefile;
 xmlreadchild(xmlr,"MapleFile",maplefile);
 string binfile;
 xmlreadchild(xmlr,"BinFile",binfile);
 string xmlfile;
 xmlreadchild(xmlr,"InputXMLFile",xmlfile);
 string xmlfile2;
 xmlreadchild(xmlr,"InputXMLFile2",xmlfile2);
 string xmlfile3;
 xmlreadchild(xmlr,"InputXMLFile3",xmlfile3);
 string sampfile;
 xmlreadchild(xmlr,"SamplingFile",sampfile);
 string samptype;
 xmlreadchild(xmlr,"EnergyNamesTypeInSamplingFile",samptype);   // "common" or "specified"
 uint binvalue,timeindex;
 xmlreadchild(xmlr,"BinValue",binvalue);
 xmlreadchild(xmlr,"TimeIndex",timeindex);
 cout << endl<<endl;
 cout << "Number of measurements in ensemble = "<<MC.getNumberOfMeasurements()<<endl;
 cout << "Ensemble ID = "<<MCOH.getEnsembleId()<<endl;
 cout << "Number of bins = "<<MC.getNumberOfBins()<<endl;
 cout << "Maple file = "<<maplefile<<endl;
 cout << "Input XML file 1 = "<<xmlfile<<endl;
 cout << "Input XML file 2 = "<<xmlfile2<<endl;
 MC.setToBootstrapMode();

 XMLHandler xmlsetup(xmlq,"TaskSetup");

 XMLHandler xmlc(xmlsetup,"CorrelatorMatrixInfo");
 CorrelatorMatrixInfo cormat(xmlc);
 if (!cormat.isHermitian()){
    throw(std::invalid_argument("CorrelatorMatrix must be Hermitian for rotation"));}
 uint maxtimesep,mintimesep;
 xmlreadchild(xmlsetup,"MinTimeSep",mintimesep);
 xmlreadchild(xmlsetup,"MaxTimeSep",maxtimesep);
 cout << "min time sep = "<<mintimesep<<endl;
 cout << "max time sep = "<<maxtimesep<<endl;

 ArgsHandler xmlqq(xmlsetup,"SinglePivotInitiate");
 map<OperatorInfo,map<OperatorInfo,Scalar> > imp_trans;
 TransMatrix imp_transmat;
 CorrelatorMatrixInfo *orig_cormat_info=0;
 XMLHandler xmlimpspec;
 if (xmlqq.queryTag("ImprovedOperators")){
    cout << "Improved operators"<<endl;
    XMLHandler xmlnew(xmlq,"ImprovedOperators");
    xmlimpspec.set(xmlnew);
    ArgsHandler xmlimp(xmlqq,"ImprovedOperators");
    read_improved_operators(xmlimp,cormat,imp_trans);
    setup_improved_operators(cormat,imp_trans,orig_cormat_info,
                             cormat.subtractVEV(),imp_transmat);}
 uint tauN=xmlqq.getUInt("NormTime");
 uint tau0=xmlqq.getUInt("MetricTime");
 uint tauD=xmlqq.getUInt("DiagonalizeTime");
 double mininvcondnum=xmlqq.getReal("MinimumInverseConditionNumber");

             //   make the fake data and putBins into memory
 CorrelatorMatrixInfo *cmat=(orig_cormat_info!=0)?orig_cormat_info:&cormat;
 make_fake_data_to_rotate(MC,*cmat,maplefile,imp_transmat,binfile,
                          mintimesep,maxtimesep,tauN,tau0,tauD,
                          mininvcondnum,timeindex,binvalue);

             //   make the fit energy samplings file
 vector<MCObsInfo> ekeys;
 string energycommonname;
 if (samptype=="common"){
    energycommonname="FitEnergyCommon";}
 else{
    for (uint index=0;index<12;index++){
    ekeys.push_back(MCObsInfo("Gold",index)); ekeys.push_back(MCObsInfo("Silver",index)); 
    ekeys.push_back(MCObsInfo("Bronze",index)); ekeys.push_back(MCObsInfo("Iron",index)); 
    ekeys.push_back(MCObsInfo("Copper",index)); ekeys.push_back(MCObsInfo("Aluminum",index));}}
 make_fit_energy_files(MC,72,sampfile,energycommonname,ekeys);

             //  make the input XML file for the run

 XMLHandler xmlout("SigMonD");
 XMLHandler xmlinit("Initialize");
 xmlinit.put_child("ProjectName","TestRotate");
 xmlinit.put_child("LogFile","log_test_rotate_corr.xml");
 xmlinit.put_child("EchoXML");
 XMLHandler xmlt;  bins_info.output(xmlt); xmlinit.put_child(xmlt);
 samp_info.output(xmlt); xmlinit.put_child(xmlt);
 XMLHandler xmlb("MCObservables");
 XMLHandler xmlf("BinData"); xmlf.put_child("FileName",binfile);
 xmlb.put_child(xmlf); xmlinit.put_child(xmlb);
 xmlout.put_child(xmlinit);

 XMLHandler xmltask("TaskSequence");
 xmlsetup.seek_root();
 xmlsetup.rename_tag("Task");
 xmltask.put_child(xmlsetup);
 xmlout.put_child(xmltask);

 ofstream xout(xmlfile);
 xout << xmlout.output()<<endl;
 xout.close();
 delete orig_cormat_info;

             //  make the input XML file for the final stage

 xmlout.set_root("SigMonD");
 xmlout.put_child("EchoXML");
 xmltask.set_root("Task");
 xmltask.put_child("Name","SIGMOND_TEST");
 xmlinit.set_root("Finish");
 bins_info.output(xmlt); xmlinit.put_child(xmlt);
 samp_info.output(xmlt); xmlinit.put_child(xmlt);
 xmlinit.put_child("MCObservables");
 XMLHandler xmlrf(xmlsetup,"WriteRotatedCorrToFile");
 string xmlrotfile;
 xmlreadchild(xmlrf,"RotatedCorrFileName",xmlrotfile);
 xmlinit.put_child("RotatedCorrFile",xmlrotfile);  
 ArgsHandler xmlqqq(xmlqq,"RotatedCorrelator");
 OperatorInfo oprot;
 xmlqqq.getItem("Operator",oprot);
 XMLHandler xmlrot; oprot.output(xmlrot);
 XMLHandler xmlrott("RotatedCorrelator");
 xmlrott.put_child(xmlrot);
 xmlinit.put_child(xmlrott);  
 xmlinit.put_child("NumberLevels",make_string(cormat.getNumberOfOperators()));
 xmlinit.put_child("TimeIndex",make_string(timeindex));
 xmlinit.put_child("BinValue",make_string(binvalue));
 if (cormat.subtractVEV())
    xmlinit.put_child("SubtractVEV");
 XMLHandler xmltt("TestRotateCorrelator");
 xmltt.put_child(xmlinit);
 xmltask.put_child(xmltt);
 xmlout.put_child(xmltask);

 ofstream xout2(xmlfile2);
 xout2 << xmlout.output()<<endl;
 xout2.close();

             //  make the input XML file for the fit energy reordering

 xmlout.set_root("SigMonD");
 xmlinit.set_root("Initialize");
 xmlinit.put_child("ProjectName","TestRotate");
 xmlinit.put_child("LogFile","log_test_rotate_corr2.xml");
 xmlinit.put_child("EchoXML");
 bins_info.output(xmlt); xmlinit.put_child(xmlt);
 samp_info.output(xmlt); xmlinit.put_child(xmlt);
 xmlb.set_root("MCObservables");
 xmlf.set_root("SamplingData"); xmlf.put_child("FileName",sampfile);
 xmlb.put_child(xmlf); xmlinit.put_child(xmlb);
 xmlout.put_child(xmlinit);

 xmltask.set_root("TaskSequence");
 xmltask.put_child("Task");
 xmltask.seek_child("Task");
 xmltask.put_child("Action","DoRotCorrMatInsertFitInfos");
 xmltask.put_child("Type","SinglePivot");
 xmltask.put_child("SinglePivotInitiate");
 xmltask.seek_child("SinglePivotInitiate");
 XMLHandler xmlpivfile(xmlsetup,"WritePivotToFile");
 string pivfile;
 xmlread(xmlpivfile,"PivotFileName",pivfile,"Setup");
 xmltask.put_child("ReadPivotFromFile");
 xmltask.seek_child("ReadPivotFromFile");
 xmltask.put_child("PivotFileName",pivfile);            
 xmltask.seek_root();
 xmltask.seek_child("Task");
 xmltask.put_child("RotatedAmplitudeCommonName","CommonAmp");
 if (samptype=="common"){
    xmltask.put_child("EnergyFitCommonName",energycommonname);}
 else{
    uint nlev=cormat.getNumberOfOperators();
    for (uint k=0;k<nlev;k++){
       XMLHandler xmllev("EnergyFit");
       xmllev.put_child("Level",make_string(k));
       xmllev.put_child("Name",ekeys[k].getObsName());
       xmllev.put_child("IDIndex",make_string(ekeys[k].getObsIndex()));
       xmltask.put_child(xmllev);}}

 xmltask.put_child("ReorderByFitEnergy");
 xmlout.put_child(xmltask);

cout << xmltask.output()<<endl;

 ofstream xout3(xmlfile3);
 xout3 << xmlout.output()<<endl;
 xout3.close();

 }
  catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}}

 if (stage==2){
 try{

 ArgsHandler xmlr(xmlq,"Finish");
 XMLHandler xmlm(xmlq,"Finish");
 MCBinsInfo bins_info(xmlm);
 MCSamplingInfo samp_info(xmlm);
 MCObsGetHandler MCOH(xmlm,bins_info,samp_info); 
 MCObsHandler MC(MCOH);
 string rotcorrfile(xmlr.getString("RotatedCorrFile"));
 uint nlevels=xmlr.getUInt("NumberLevels");
 uint timeval=xmlr.getUInt("TimeIndex");
 uint binval=xmlr.getUInt("BinValue");
 ArgsHandler xmlg(xmlr,"RotatedCorrelator");
 OperatorInfo opinfo(xmlg.getItem<OperatorInfo>("RotatedCorrelator"));
 bool vevs=xmlr.getBool("SubtractVEV");

 cout.precision(14);
 XMLHandler xmlf;
 MC.readBinsFromFile(rotcorrfile,xmlf);
 cout <<xmlf.output()<<endl;
 for (uint col=0;col<nlevels;col++){
    opinfo.resetGenIrrepIDIndex(col);
    CorrelatorAtTimeInfo corr(opinfo,opinfo,timeval,true,false);
    MCObsInfo obskey(corr,RealPart);
    try{
    const RVector& bins=MC.getBins(obskey);
    cout << "corr[level="<<col<<"][t="<<timeval<<"][bin="<<binval<<"]="<<bins[binval]<<endl;}
    catch(const std::exception& xp){}}

#if defined COMPLEXNUMBERS
 if (vevs){
    cout << "VEVs for bin = "<<binval<<endl;
    for (uint col=0;col<nlevels;col++){
       opinfo.resetGenIrrepIDIndex(col);
       MCObsInfo obskey(opinfo,RealPart);
       try{
       const RVector& bins_re=MC.getBins(obskey);
     //  obskey.setToImaginaryPart();
     //  const RVector& bins_im=MC.getBins(obskey);
    //   cout << "("<<bins_re[binval]<<", "<<bins_im[binval]<<")"<<endl;}
       cout <<"vev[level="<<col<<"][bin="<<binval<<"]="<<bins_re[binval]<<endl;}
       catch(const std::exception& xp){}}
  //  cout << "VEVs mags for bin = "<<binval<<endl;
  //  for (uint col=0;col<nlevels;col++){
  //     opinfo.resetGenIrrepIDIndex(col);
  //     MCObsInfo obskey(opinfo,RealPart);
  //     const RVector& bins_re=MC.getBins(obskey);
  //     obskey.setToImaginaryPart();
  //     const RVector& bins_im=MC.getBins(obskey);
  //     cout << std::abs(complex<double>(bins_re[binval],bins_im[binval]))<<endl;}
       }
#else
 if (vevs){
    cout << "VEVs for bin = "<<binval<<endl;
    for (uint col=0;col<nlevels;col++){
       opinfo.resetGenIrrepIDIndex(col);
       MCObsInfo obskey(opinfo,RealPart);
       try{
       const RVector& bins=MC.getBins(obskey);
       cout <<"vev[level="<<col<<"][bin="<<binval<<"]="<<bins[binval]<<endl;}
       catch(const std::exception& xp){}}
   // cout << "VEVs mags for bin = "<<binval<<endl;
   // for (uint col=0;col<nlevels;col++){
   //    opinfo.resetGenIrrepIDIndex(col);
   //    MCObsInfo obskey(opinfo,RealPart);
   //    const RVector& bins=MC.getBins(obskey);
   //    cout << std::abs(bins[binval])<<endl;}
       }
#endif
 }
 catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}}
}


void testPivotCorrelator(XMLHandler& xml_in, int taskcount)
{
 if (xml_tag_count(xml_in,"TestPivotCorrelator")==0)
 return;

 cout << endl << "Starting TestPivotCorrelator"<<endl;
 

 try{

 XMLHandler xmlq(xml_in,"TestPivotCorrelator");
 XMLHandler xmlr(xmlq,"Setup");
 uint mintimesep,maxtimesep,tau0,nlevels,noperators;
 xmlreadchild(xmlr,"MinTimeSep",mintimesep);
 xmlreadchild(xmlr,"MaxTimeSep",maxtimesep);
 xmlreadchild(xmlr,"NumberLevels",nlevels);
 xmlreadchild(xmlr,"NumberOperators",noperators);
 xmlreadchild(xmlr,"MetricTime",tau0);
// double metric_min_inv_condnum;
 double min_inv_condnum;
 xmlreadchild(xmlr,"MinimumInverseConditionNumber",min_inv_condnum);
// xmlreadchild(xmlr,"MatrixMinimumInverseConditionNumber",matrix_min_inv_condnum);

 cout << endl<<endl<<"Starting complex Hermitian tests"<<endl<<endl;

 cout.precision(14);
 vector<ComplexHermitianMatrix> Cor(maxtimesep-mintimesep+1,nlevels);
 vector<double> energies(nlevels);
 energies[0]=0.2;
 for (uint n=1;n<nlevels;n++)
    energies[n]=energies[n-1]+0.08/sqrt(double(n));
 for (uint n=0;n<nlevels;n++)
    cout << "energy["<<n<<"] = "<<energies[n]<<endl;

 CMatrix Zcoef(noperators,nlevels);
 for (int j=0;j<int(noperators);j++)
 for (int n=0;n<int(nlevels);n++)
   Zcoef(j,n)=complex<double>( double((rand() % 4096)-2048)/237.0,
                               double((rand() % 4096)-2048)/841.0);

 for (uint t=mintimesep;t<=maxtimesep;t++){
    for (int i=0;i<int(noperators);i++)
    for (int j=0;j<int(noperators);j++){
      complex<double> tmp(0.0,0.0);
      for (int n=0;n<int(nlevels);n++)
         tmp+=Zcoef(i,n)*conjugate(Zcoef(j,n))*exp(-energies[n]*t);
      Cor[t-mintimesep].put(i,j,tmp);}}

 HermDiagonalizerWithMetric DM;
 DM.setExceptionsOff();
 DM.setMinInvCondNum(min_inv_condnum);
 cout << "Min Inv Cond Number = "<<DM.getMinInvCondNum()<<endl;
 cout << "Number of levels = "<<nlevels<<endl;
 cout << "Number of operators = "<<noperators<<endl;
 cout << "is metric set? "<<DM.isMetricSet()<<endl;
 cout << "metric rank = "<<DM.getMetricRank()<<endl;

 ComplexHermitianMatrix& B(Cor[tau0-mintimesep]);

 int info=DM.setMetric(B);
 cout << "info = "<<info<<endl;
 cout << "is metric set? "<<DM.isMetricSet()<<endl;

 RVector Beigvals;
 DM.getMetricEigenvalues(Beigvals);
 for (int i=0;i<int(Beigvals.size());i++)
    cout << "Beigval["<<i<<"] = "<<Beigvals[i]<<endl;
 cout << "metric rank = "<<DM.getMetricRank()<<endl;


 for (uint t=mintimesep;t<=maxtimesep;t++)
    if (t!=tau0){
       cout << endl<<endl<<"  t = "<<t<<endl<<endl;
       ComplexHermitianMatrix& A=Cor[t-mintimesep];
       int info=DM.setMatrix(A);
       cout << "info = "<<info<<endl;
       cout << " is matrix set? "<<DM.isMatrixSet()<<endl;
       cout << " rank of matrix = "<<DM.getMatrixRank()<<endl;
       cout << " null B in null A? "<<DM.isNullMetricInNullMatrix()<<endl;
       RVector eigvals;
       CMatrix eigvecs,orthovecs,Zmat;
       DM.getEigenvalues(eigvals);
       DM.getEigenvectors(eigvecs);
       DM.getOrthovectors(orthovecs);
       DM.getZMatrix(Zmat);

       for (int k=0;k<int(eigvals.size());k++)
          cout << "eigval["<<k<<"] = "<<eigvals[k]<<endl;


/*
       DM.getEigenvectors(Cor[t-mintimesep],Cor[tau0-mintimesep], 
                          Beigvals,eigvals,eigvecs,orthovecs,Zmat);

       cout << "Size of Beigvals = "<<Beigvals.size()<<endl;
       cout << "Size of eigvals= "<<eigvals.size()<<endl;
       cout << "eigvecs has "<<eigvecs.size(0)<<" rows and "<<eigvecs.size(1)<<" columns"<<endl;
       cout << "orthovecs has "<<orthovecs.size(0)<<" rows and "<<orthovecs.size(1)<<" columns"<<endl;
       cout << "Zmat has "<<Zmat.size(0)<<" rows and "<<Zmat.size(1)<<" columns"<<endl;
       Diagonalizer Dz;
       RVector LB; RMatrix Beigvecs;
       Dz.getEigenvectors(Cor[tau0-mintimesep],LB,Beigvecs);
       double lthreshold=LB[LB.size()-1]*DM.getMinInvCondNumOfMetric();
       cout << "lthreshold = "<<lthreshold<<endl;
       int Bdiscard=0;
       while ((Bdiscard<int(LB.size()))&&(LB[Bdiscard]<lthreshold)) Bdiscard++;
       cout << "Bdiscard = "<<Bdiscard<<endl;
       int n0=noperators-Bdiscard;
       RMatrix P0(int(noperators),n0);
       for (int k=0;k<int(noperators);k++)
       for (int l=0;l<n0;l++)
          P0(k,l)=Beigvecs(k,l+Bdiscard);
       RealSymmetricMatrix G(n0,n0);
       for (int k=0;k<n0;k++)
       for (int l=0;l<n0;l++){
          double tmp=0.0;
          for (int j=0;j<int(noperators);j++)
          for (int i=0;i<int(noperators);i++)
             tmp+=P0(i,k)*Cor[t-mintimesep](i,j)*P0(j,l);
          G(k,l)=tmp/sqrt(LB[k+Bdiscard]*LB[l+Bdiscard]);}

       RVector LA;
       Dz.getEigenvalues(G,LA);

       for (int k=0;k<int(LA.size());k++)
          cout << "LA["<<k<<"] = "<<LA[k]<<endl;


       lthreshold=LA[LA.size()-1]*DM.getMinInvCondNumOfMatrix();
       cout << "lthreshold = "<<lthreshold<<endl;
       int Gdiscard=0;
       while ((Gdiscard<int(LA.size()))&&(LA[Gdiscard]<lthreshold)) Gdiscard++;
       cout << "Gdiscard = "<<Gdiscard<<endl;

        //  test  G * orthovecs = orthovec * eigvals
       cout << "orthovecs has "<<orthovecs.size(0)<<" rows and "<<orthovecs.size(1)<<" columns"<<endl;
       int nn=orthovecs.size(1);
       if (int(eigvals.size())!=nn) cout << "PROBLEM"<<endl;
       if ((int(noperators)-Bdiscard-Gdiscard)!=nn) cout << "PROBLEM"<<endl;
       for (int k=0;k<n0;k++)
       for (int l=0;l<nn;l++){
          double tmp=0.0;
          for (int j=0;j<n0;j++)
             tmp+=G(k,j)*orthovecs(j,l);
          tmp-=orthovecs(k,l)*LA[l+Gdiscard];
          if (abs(tmp)>1e-10) cout << "WRONG!!"<<endl;}
      
       for (int k=0;k<int(Beigvals.size());k++)
          cout << "Beigval["<<k<<"] = "<<Beigvals[k]<<"  "<<LB[k]<<endl;
       cout <<endl;
       for (int k=0;k<int(eigvals.size());k++)
          cout << "eigval["<<k<<"] = "<<eigvals[k]<< "   "<<LA[k+Gdiscard]<<endl;

       cout << "checking rotated B"<<endl;
       for (int k=0;k<nn;k++)
       for (int l=0;l<nn;l++){
          double chk=0.0;
          for (int i=0;i<int(noperators);i++)
          for (int j=0;j<int(noperators);j++)
             chk+=eigvecs(i,k)*B(i,j)*eigvecs(j,l);
          if (k==l) chk-=1.0;
          if (abs(chk)>1e-10) cout << "WRONG!!"<<endl;}
       cout << "checking rotated A"<<endl;
       for (int k=0;k<nn;k++)
       for (int l=0;l<nn;l++){
          double chk=0.0;
          for (int i=0;i<int(noperators);i++)
          for (int j=0;j<int(noperators);j++)
             chk+=eigvecs(i,k)*A(i,j)*eigvecs(j,l);
          if (k==l) chk-=LA[k+Gdiscard];
          if (abs(chk)>1e-10) cout << "WRONG!!"<<endl;}

       if (nn!=int(nlevels))
          cout << "cannot check Zmat since retained levels less than initial"<<endl;
       else{
       cout << "checking Zmat"<<endl;
       for (int k=0;k<int(noperators);k++)
       for (int l=0;l<int(noperators);l++){
          double chkA=0.0;
          double chkB=0.0;
          for (int j=0;j<nn;j++){
             chkA+=Zmat(k,j)*LA[j+Gdiscard]*Zmat(l,j);
             chkB+=Zmat(k,j)*Zmat(l,j);}
          chkA-=A(k,l);
          chkB-=B(k,l);
          if (abs(chkA)>1e-10) cout << "WRONG Zmat A!!"<<endl;
          if (abs(chkB)>1e-10) cout << "WRONG Zmat B!!"<<endl;}}
*/
       }


/*

 cout << endl<<endl<<"Starting real symmetric tests"<<endl<<endl;

 cout.precision(14);
 vector<RealSymmetricMatrix> Cor(maxtimesep-mintimesep+1,nlevels);
 vector<double> energies(nlevels);
 energies[0]=0.2;
 for (uint n=1;n<nlevels;n++)
    energies[n]=energies[n-1]+0.08/sqrt(double(n));
 for (uint n=0;n<nlevels;n++)
    cout << "energy["<<n<<"] = "<<energies[n]<<endl;

 RMatrix Zcoef(noperators,nlevels);
 for (int j=0;j<int(noperators);j++)
 for (int n=0;n<int(nlevels);n++)
   Zcoef(j,n)=1.0/(1.0+0.6*double((j-n)*(j-n)));

 for (uint t=mintimesep;t<=maxtimesep;t++){
    for (int i=0;i<int(noperators);i++)
    for (int j=0;j<int(noperators);j++){
      double tmp=0.0;
      for (int n=0;n<int(nlevels);n++)
         tmp+=Zcoef(i,n)*Zcoef(j,n)*exp(-energies[n]*t);
      Cor[t-mintimesep](i,j)=tmp;}}

 for (uint t=mintimesep;t<=maxtimesep;t++)
    if (t!=tau0){
       cout << endl<<endl<<"  t = "<<t<<endl<<endl;
       RVector eigvals,Beigvals;
       RMatrix orthovecs,eigvecs,Zmat;
       RealSymmetricMatrix& A=Cor[t-mintimesep];
       RealSymmetricMatrix& B=Cor[tau0-mintimesep];
       DM.getEigenvectors(Cor[t-mintimesep],Cor[tau0-mintimesep], 
                          Beigvals,eigvals,eigvecs,orthovecs,Zmat);

       cout << "Size of Beigvals = "<<Beigvals.size()<<endl;
       cout << "Size of eigvals= "<<eigvals.size()<<endl;
       cout << "eigvecs has "<<eigvecs.size(0)<<" rows and "<<eigvecs.size(1)<<" columns"<<endl;
       cout << "orthovecs has "<<orthovecs.size(0)<<" rows and "<<orthovecs.size(1)<<" columns"<<endl;
       cout << "Zmat has "<<Zmat.size(0)<<" rows and "<<Zmat.size(1)<<" columns"<<endl;
       Diagonalizer Dz;
       RVector LB; RMatrix Beigvecs;
       Dz.getEigenvectors(Cor[tau0-mintimesep],LB,Beigvecs);
       double lthreshold=LB[LB.size()-1]*DM.getMinInvCondNumOfMetric();
       cout << "lthreshold = "<<lthreshold<<endl;
       int Bdiscard=0;
       while ((Bdiscard<int(LB.size()))&&(LB[Bdiscard]<lthreshold)) Bdiscard++;
       cout << "Bdiscard = "<<Bdiscard<<endl;
       int n0=noperators-Bdiscard;
       RMatrix P0(int(noperators),n0);
       for (int k=0;k<int(noperators);k++)
       for (int l=0;l<n0;l++)
          P0(k,l)=Beigvecs(k,l+Bdiscard);
       RealSymmetricMatrix G(n0,n0);
       for (int k=0;k<n0;k++)
       for (int l=0;l<n0;l++){
          double tmp=0.0;
          for (int j=0;j<int(noperators);j++)
          for (int i=0;i<int(noperators);i++)
             tmp+=P0(i,k)*Cor[t-mintimesep](i,j)*P0(j,l);
          G(k,l)=tmp/sqrt(LB[k+Bdiscard]*LB[l+Bdiscard]);}

       RVector LA;
       Dz.getEigenvalues(G,LA);

       for (int k=0;k<int(LA.size());k++)
          cout << "LA["<<k<<"] = "<<LA[k]<<endl;


       lthreshold=LA[LA.size()-1]*DM.getMinInvCondNumOfMatrix();
       cout << "lthreshold = "<<lthreshold<<endl;
       int Gdiscard=0;
       while ((Gdiscard<int(LA.size()))&&(LA[Gdiscard]<lthreshold)) Gdiscard++;
       cout << "Gdiscard = "<<Gdiscard<<endl;

        //  test  G * orthovecs = orthovec * eigvals
       cout << "orthovecs has "<<orthovecs.size(0)<<" rows and "<<orthovecs.size(1)<<" columns"<<endl;
       int nn=orthovecs.size(1);
       if (int(eigvals.size())!=nn) cout << "PROBLEM"<<endl;
       if ((int(noperators)-Bdiscard-Gdiscard)!=nn) cout << "PROBLEM"<<endl;
       for (int k=0;k<n0;k++)
       for (int l=0;l<nn;l++){
          double tmp=0.0;
          for (int j=0;j<n0;j++)
             tmp+=G(k,j)*orthovecs(j,l);
          tmp-=orthovecs(k,l)*LA[l+Gdiscard];
          if (abs(tmp)>1e-10) cout << "WRONG!!"<<endl;}
      
       for (int k=0;k<int(Beigvals.size());k++)
          cout << "Beigval["<<k<<"] = "<<Beigvals[k]<<"  "<<LB[k]<<endl;
       cout <<endl;
       for (int k=0;k<int(eigvals.size());k++)
          cout << "eigval["<<k<<"] = "<<eigvals[k]<< "   "<<LA[k+Gdiscard]<<endl;

       cout << "checking rotated B"<<endl;
       for (int k=0;k<nn;k++)
       for (int l=0;l<nn;l++){
          double chk=0.0;
          for (int i=0;i<int(noperators);i++)
          for (int j=0;j<int(noperators);j++)
             chk+=eigvecs(i,k)*B(i,j)*eigvecs(j,l);
          if (k==l) chk-=1.0;
          if (abs(chk)>1e-10) cout << "WRONG!!"<<endl;}
       cout << "checking rotated A"<<endl;
       for (int k=0;k<nn;k++)
       for (int l=0;l<nn;l++){
          double chk=0.0;
          for (int i=0;i<int(noperators);i++)
          for (int j=0;j<int(noperators);j++)
             chk+=eigvecs(i,k)*A(i,j)*eigvecs(j,l);
          if (k==l) chk-=LA[k+Gdiscard];
          if (abs(chk)>1e-10) cout << "WRONG!!"<<endl;}

       if (nn!=int(nlevels))
          cout << "cannot check Zmat since retained levels less than initial"<<endl;
       else{
       cout << "checking Zmat"<<endl;
       for (int k=0;k<int(noperators);k++)
       for (int l=0;l<int(noperators);l++){
          double chkA=0.0;
          double chkB=0.0;
          for (int j=0;j<nn;j++){
             chkA+=Zmat(k,j)*LA[j+Gdiscard]*Zmat(l,j);
             chkB+=Zmat(k,j)*Zmat(l,j);}
          chkA-=A(k,l);
          chkB-=B(k,l);
          if (abs(chkA)>1e-10) cout << "WRONG Zmat A!!"<<endl;
          if (abs(chkB)>1e-10) cout << "WRONG Zmat B!!"<<endl;}}

       }


    bool isMatrixPosDefRequired() const {return Aposdef;}


    void getEigenvalues(const RealSymmetricMatrix& A, 
                        const RealSymmetricMatrix& B,
                        RVector& eigvals);
    void getEigenvectors(const RealSymmetricMatrix& A, 
                         const RealSymmetricMatrix& B,
                         RVector& eigvals, RMatrix& eigvecs);
    void getOrthovectors(const RealSymmetricMatrix& A, 
                         const RealSymmetricMatrix& B,
                         RVector& eigvals, RMatrix& orthovecs);
    void getEigenvectors(const RealSymmetricMatrix& A, 
                         const RealSymmetricMatrix& B,
                         RVector& eigvals, RMatrix& eigvecs,
                         RMatrix& orthovecs);

    unsigned int getCurrentMetricRank() const {return rankB;}

    void getEigenvalues(const ComplexHermitianMatrix& A, 
                        const ComplexHermitianMatrix& B,
                        RVector& eigvals);
    void getEigenvectors(const ComplexHermitianMatrix& A, 
                         const ComplexHermitianMatrix& B,
                         RVector& eigvals, CMatrix& eigvecs);
    void getOrthovectors(const ComplexHermitianMatrix& A, 
                         const ComplexHermitianMatrix& B,
                         RVector& eigvals, CMatrix& orthovecs);
    void getEigenvectors(const ComplexHermitianMatrix& A, 
                         const ComplexHermitianMatrix& B,
                         RVector& eigvals, CMatrix& eigvecs,
                         CMatrix& orthovecs);

 private:

    void diagonalize(const RealSymmetricMatrix& A, 
                     const RealSymmetricMatrix& B,
                     RVector& eigvals, RMatrix& eigvecs, RMatrix& orthovecs,
                     bool calceigvecs, bool calcortho);

    void diagonalize(const ComplexHermitianMatrix& A, 
                     const ComplexHermitianMatrix& B,
                     RVector& eigvals, CMatrix& eigvecs, CMatrix& orthovecs,
                     bool calceigvecs, bool calcortho);
};









             //   make the fake data and putBins into memory

 make_fake_data_to_rotate(MC,cormat,maplefile,mintimesep,maxtimesep);

             //   now fire up the correlator rotator

 XMLHandler xmllog;
// RotatedCorrelatorMatrix RCM(MC,xmlq,xmllog);
// cout << xmllog.output()<<endl;

*/
 }
 catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}

}





void testPivotCorrelator0(XMLHandler& xml_in, int taskcount)
{
 if (xml_tag_count(xml_in,"TestPivotCorrelator0")==0)
 return;

 cout << endl << "Starting TestPivotCorrelator0"<<endl;
 
  /* initialize random seed: */
//  srand (time(NULL));

// try{

// XMLHandler xmlq(xml_in,"TestPivotCorrelator0");
// XMLHandler xmlr(xmlq,"Setup");
// int Noperators;
// xmlreadchild(xmlr,"NumberOperators",Noperators);
// string mode;
// xmlreadchild(xmlr,"Mode",mode);
// if ((mode=="real")||(mode=="Real")){

// cout << "Mode is Real"<<endl;
// cout << "Number of operators = "<<Noperators<<endl;

// RealSymmetricMatrix B(Noperators),A(Noperators);
// for (int i=0;i<int(Noperators);i++)
// for (int j=i;j<int(Noperators);j++){
//    B(i,j)=double((rand()%4096-2048))/2013.0;
//    A(i,j)=double((rand()%4096-2048))/2013.0;}

//    //  make A positive definite
// Diagonalizer dtemp;
// RVector Ltemp; RMatrix Utemp;
/* dtemp.getEigenvectors(A,Ltemp,Utemp);
// for (int k=0;k<Noperators;k++)
//    Ltemp[k]=fabs(Ltemp[k]);
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//       chk+=Utemp(i,k)*Ltemp[k]*Utemp(j,k);
//    A(i,j)=chk;} */
//    //  make B positive definite
/* dtemp.getEigenvectors(B,Ltemp,Utemp);
// for (int k=0;k<Noperators;k++)
//    Ltemp[k]=fabs(Ltemp[k]);
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//       chk+=Utemp(i,k)*Ltemp[k]*Utemp(j,k);
//    B(i,j)=chk;} */















// Diagonalizer dz;
// RVector L0; RMatrix U0;
// dz.getEigenvectors(B,L0,U0);

//    //  make Bsafe positive definite
// RealSymmetricMatrix Bsafe(Noperators);
// RVector L0safe(Noperators);
// for (int k=0;k<Noperators;k++)
//    L0safe[k]=(L0[k]>0.0?L0[k]:0.0);

// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//       chk+=U0(i,k)*L0safe[k]*U0(j,k);
//    Bsafe(i,j)=chk;}

// RealSymmetricMatrix Bsqrt(Noperators),Binvsqrt(Noperators);
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//       chk+=U0(i,k)*sqrt(L0safe[k])*U0(j,k);
//    Bsqrt(i,j)=chk;}
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//       if (L0safe[k]>0.0) chk+=U0(i,k)*U0(j,k)/sqrt(L0safe[k]);
//    Binvsqrt(i,j)=chk;}

//  // do some checks

// dz.getEigenvectors(B,L0,U0);
// for (int k=0;k<Noperators;k++)
//    cout << "B eigenvalue "<<k<<" = "<<L0[k]<<endl;
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//       chk+=U0(i,k)*L0[k]*U0(j,k);
//    chk-=B(i,j);
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}

// cout <<endl<< "Now diagonalizing with metric!!"<<endl<<endl;
// DiagonalizerWithMetric dzwm;
// RVector Beigvals, eigvals;
// RMatrix eigvecs,orthovecs;
// dzwm.removeMinInvCondNumOfMatrix();
// dzwm.getEigenvectors(A,B,Beigvals,eigvals,eigvecs,orthovecs);

// cout << "Noperators = "<<Noperators<<endl;
// int NP=eigvecs.size(1);
// cout << "Number of retained eigenvectors = "<<NP<<endl;
// cout << "Beigvals size = "<<Beigvals.size()<<endl;
// cout << "eigvals size = "<<eigvals.size()<<endl;
// cout << "eigvecs has "<<eigvecs.size(0)<<" rows and "<<eigvecs.size(1)<<" columns"<<endl;
// cout << "orthovecs has "<<orthovecs.size(0)<<" rows and "<<orthovecs.size(1)<<" columns"<<endl;
// for (int k=0;k<Noperators;k++)
//    cout << "B eigenvalue "<<k<<" = "<<Beigvals[k]<<endl;
// for (int k=0;k<int(eigvals.size());k++)
//    cout << "eigenvalue "<<k<<" = "<<eigvals[k]<<endl;
// 
// cout << "checking A * eigvecs = B * eigvecs * eigvals"<<endl;
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<NP;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//       chk+=A(i,k)*eigvecs(k,j);
//    double chk2=0.0;
//    for (int k=0;k<Noperators;k++)
//       chk2+=B(i,k)*eigvecs(k,j)*eigvals[j];
//    if (fabs(chk-chk2)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking A * eigvecs = Bsafe * eigvecs * eigvals"<<endl;
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<NP;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//       chk+=A(i,k)*eigvecs(k,j);
//    double chk2=0.0;
//    for (int k=0;k<Noperators;k++)
//       chk2+=Bsafe(i,k)*eigvecs(k,j)*eigvals[j];
//    if (fabs(chk-chk2)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking eigvecs^dag B eigvec = Id"<<endl;
// for (int i=0;i<NP;i++)
// for (int j=0;j<NP;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//    for (int l=0;l<Noperators;l++)
//       chk+=eigvecs(k,i)*B(k,l)*eigvecs(l,j);
//    if (i==j) chk-=1.0;
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking eigvecs^dag Bsafe eigvec = Id"<<endl;
// for (int i=0;i<NP;i++)
// for (int j=0;j<NP;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//    for (int l=0;l<Noperators;l++)
//       chk+=eigvecs(k,i)*Bsafe(k,l)*eigvecs(l,j);
//    if (i==j) chk-=1.0;
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking eigvecs^dag A eigvec = LambdaA_[NP,NP]"<<endl;
// for (int i=0;i<NP;i++)
// for (int j=0;j<NP;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//    for (int l=0;l<Noperators;l++)
//       chk+=eigvecs(k,i)*A(k,l)*eigvecs(l,j);
//    if (i==j) chk-=eigvals[i];
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking sqrt(Bsafe) eigvec = orthovecs"<<endl;
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<NP;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//       chk+=Bsqrt(i,k)*eigvecs(k,j);
//    chk-=orthovecs(i,j);
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking orthovec^dag orthovec = Id"<<endl;
// for (int i=0;i<NP;i++)
// for (int j=0;j<NP;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//       chk+=orthovecs(k,i)*orthovecs(k,j);
//    if (i==j) chk-=1.0;
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}


// cout << "checking orthovecs^dag * Bsafe^(-1/2)*A*Bsafe^(-1/2) * orthovecs = eigvals"<<endl;
// RealSymmetricMatrix G(Noperators);
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//    for (int l=0;l<Noperators;l++)
//       chk+=Binvsqrt(i,k)*A(k,l)*Binvsqrt(l,j);
//    G(i,j)=chk;}
// for (int i=0;i<NP;i++)
// for (int j=0;j<NP;j++){
//    double chk=0.0;
//    for (int k=0;k<Noperators;k++)
//    for (int l=0;l<Noperators;l++)
//       chk+=orthovecs(k,i)*G(k,l)*orthovecs(l,j);
//    if (i==j) chk-=eigvals[i];
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}

// }



// else{     //  complex version

// cout << "Mode is Complex"<<endl;
// cout << "Number of operators = "<<Noperators<<endl;

// ComplexHermitianMatrix B(Noperators),A(Noperators);
// for (int i=0;i<int(Noperators);i++){
//    B.put(i,i,complex<double>(double((rand()%4096-2048))/2013.0,0.0));
//    A.put(i,i,complex<double>(double((rand()%4096-2048))/2013.0,0.0));
//    for (int j=i+1;j<int(Noperators);j++){
//       B.put(i,j,complex<double>(double((rand()%4096-2048))/2013.0,
//                              double((rand()%4096-2048))/2013.0));
//       A.put(i,j,complex<double>(double((rand()%4096-2048))/2013.0,
//                               double((rand()%4096-2048))/2013.0));}}

// Diagonalizer dz;
// RVector L0; CMatrix U0;
// dz.getEigenvectors(B,L0,U0);

//    //  make Bsafe positive definite
// ComplexHermitianMatrix Bsafe(Noperators);
// RVector L0safe(Noperators);
// for (int k=0;k<Noperators;k++)
//    L0safe[k]=(L0[k]>0.0?L0[k]:0.0);
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//       chk+=U0(i,k)*L0[k]*conjugate(U0(j,k));
//    if (i==j){
//       if (fabs(chk.imag())>1e-10) cout << "B BAD"<<endl;
//       chk=complex<double>(chk.real(),0.0);}
//    Bsafe.put(i,j,chk);}

// ComplexHermitianMatrix Bsqrt(Noperators),Binvsqrt(Noperators);
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//       chk+=U0(i,k)*sqrt(L0safe[k])*conjugate(U0(j,k));
//    if (i==j){
//       if (fabs(chk.imag())>1e-10) cout << "Bsqrt BAD"<<endl;
//       chk=complex<double>(chk.real(),0.0);}
//    Bsqrt.put(i,j,chk);}
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//       if (L0safe[k]>0.0) chk+=U0(i,k)*conjugate(U0(j,k))/sqrt(L0safe[k]);
//    if (i==j){
//       if (fabs(chk.imag())>1e-10) cout << "Binvsqrt BAD"<<endl;
//       chk=complex<double>(chk.real(),0.0);}
//    Binvsqrt.put(i,j,chk);}

/* RVector LA; CMatrix UA;
// dz.getEigenvectors(AA,LA,UA);

//    //  make A positive definite
// for (int k=0;k<Noperators;k++)
//    LA[k]=fabs(LA[k]);
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//       chk+=UA(i,k)*LA[k]*conjugate(UA(j,k));
//    if (i==j){
//       if (fabs(chk.imag())>1e-10) cout << "AA BAD"<<endl;
//       chk=complex<double>(chk.real(),0.0);}
//    AA.put(i,j,chk);}
//   //  form  A = Bsqrt * AA * Bsqrt
// ComplexHermitianMatrix A(Noperators);
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//    for (int l=0;l<Noperators;l++)
//       chk+=Bsqrt(i,k)*AA(k,l)*Bsqrt(l,j);
//    if (i==j){
//       if (fabs(chk.imag())>1e-10) cout << "BAD"<<endl;
//       chk=complex<double>(chk.real(),0.0);}
//    A.put(i,j,chk);}
*/
//  // do some checks

// dz.getEigenvectors(B,L0,U0);
// for (int k=0;k<Noperators;k++)
//    cout << "B eigenvalue "<<k<<" = "<<L0[k]<<endl;
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//       chk+=U0(i,k)*L0[k]*conjugate(U0(j,k));
//    chk-=B(i,j);
//    if (abs(chk)>1e-10) cout << "PROBLEM"<<endl;}

// cout <<endl<< "Now diagonalizing with metric!!"<<endl<<endl;
// DiagonalizerWithMetric dzwm;
// RVector Beigvals, eigvals;
// CMatrix eigvecs,orthovecs;
// dzwm.removeMinInvCondNumOfMatrix();
// dzwm.getEigenvectors(A,B,Beigvals,eigvals,eigvecs,orthovecs);

// cout << "Noperators = "<<Noperators<<endl;
// int NP=eigvecs.size(1);
// cout << "Number of retained eigenvectors = "<<NP<<endl;
// cout << "Beigvals size = "<<Beigvals.size()<<endl;
// cout << "eigvals size = "<<eigvals.size()<<endl;
// cout << "eigvecs has "<<eigvecs.size(0)<<" rows and "<<eigvecs.size(1)<<" columns"<<endl;
// cout << "orthovecs has "<<orthovecs.size(0)<<" rows and "<<orthovecs.size(1)<<" columns"<<endl;
// for (int k=0;k<Noperators;k++)
//    cout << "B eigenvalue "<<k<<" = "<<Beigvals[k]<<endl;
// for (int k=0;k<int(eigvals.size());k++)
//    cout << "eigenvalue "<<k<<" = "<<eigvals[k]<<endl;

// cout << "checking A * eigvecs = B * eigvecs * eigvals"<<endl;
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<NP;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//       chk+=A(i,k)*eigvecs(k,j);
//    complex<double> chk2(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//       chk2+=B(i,k)*eigvecs(k,j)*eigvals[j];
//    if (abs(chk-chk2)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking A * eigvecs = Bsafe * eigvecs * eigvals"<<endl;
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<NP;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//       chk+=A(i,k)*eigvecs(k,j);
//    complex<double> chk2(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//       chk2+=Bsafe(i,k)*eigvecs(k,j)*eigvals[j];
//    if (abs(chk-chk2)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking eigvecs^dag B eigvec = Id"<<endl;
// for (int i=0;i<NP;i++)
// for (int j=0;j<NP;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//    for (int l=0;l<Noperators;l++)
//       chk+=conjugate(eigvecs(k,i))*B(k,l)*eigvecs(l,j);
//    if (i==j) chk-=1.0;
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking eigvecs^dag Bsafe eigvec = Id"<<endl;
// for (int i=0;i<NP;i++)
// for (int j=0;j<NP;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//    for (int l=0;l<Noperators;l++)
//       chk+=conjugate(eigvecs(k,i))*Bsafe(k,l)*eigvecs(l,j);
//    if (i==j) chk-=1.0;
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking eigvecs^dag A eigvec = LambdaA_[NP,NP]"<<endl;
// for (int i=0;i<NP;i++)
// for (int j=0;j<NP;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//    for (int l=0;l<Noperators;l++)
//       chk+=conjugate(eigvecs(k,i))*A(k,l)*eigvecs(l,j);
//    if (i==j) chk-=eigvals[i];
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking sqrt(Bsafe) eigvec = orthovecs"<<endl;
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<NP;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//       chk+=Bsqrt(i,k)*eigvecs(k,j);
//    chk-=orthovecs(i,j);
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}

// cout << "checking orthovec^dag orthovec = Id"<<endl;
// for (int i=0;i<NP;i++)
// for (int j=0;j<NP;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//       chk+=conjugate(orthovecs(k,i))*orthovecs(k,j);
//    if (i==j) chk-=1.0;
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}


// cout << "checking orthovecs^dag * Bsafe^(-1/2)*A*Bsafe^(-1/2) * orthovecs = eigvals"<<endl;
// ComplexHermitianMatrix G(Noperators);
// for (int i=0;i<Noperators;i++)
// for (int j=0;j<Noperators;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//    for (int l=0;l<Noperators;l++)
//       chk+=Binvsqrt(i,k)*A(k,l)*Binvsqrt(l,j);
//    if (i==j){
//       if (fabs(chk.imag())>1e-10) cout << "BAD"<<endl;
//       chk=complex<double>(chk.real(),0.0);}
//    G.put(i,j,chk);}
// for (int i=0;i<NP;i++)
// for (int j=0;j<NP;j++){
//    complex<double> chk(0.0,0.0);
//    for (int k=0;k<Noperators;k++)
//    for (int l=0;l<Noperators;l++)
//       chk+=conjugate(orthovecs(k,i))*G(k,l)*orthovecs(l,j);
//    if (i==j) chk-=eigvals[i];
//    if (fabs(chk)>1e-10) cout << "PROBLEM"<<endl;}







// }


// }
// catch(const std::exception& err){
//    cerr << "  Error: "<<err.what()<<endl;
//    cerr << "  ... exiting..."<<endl;
//    exit(1);}

}
// ***********************************************************************

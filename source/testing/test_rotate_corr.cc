#include <cstdio>
#include <ctime>
#include <algorithm>
#include "xml_handler.h"
#include "mcobs_handler.h"
#include "single_pivot.h"
#include "correlator_matrix_info.h"

using namespace std;
using namespace LaphEnv;

double get_random_minusone_to_one()
{
 return double(int((rand() % 4096))-2048)/2048.0;
}


#if defined COMPLEXNUMBERS

void make_fake_data_to_rotate(MCObsHandler& MC, CorrelatorMatrixInfo& cormat,
                              const string& maplefile, 
                              const string& binfile, uint tmin, uint tmax)
{
 uint nops=cormat.getNumberOfOperators();
 bool vevs=cormat.isVEVSubtracted();
 uint nbins=MC.getNumberOfBins();
 srand (time(NULL));

 set<MCObsInfo> obskeys;
 Vector<double> bins_re(nbins);
 Vector<double> bins_im(nbins);
 uint nlevels=nops;
 CMatrix Zcoef(nops,nlevels);
 vector<double> energies(nlevels);
 energies[0]=0.15;
 for (uint level=1;level<nlevels;level++){
    energies[level]=energies[level-1]+0.08/double(level);}
 for (uint level=0;level<nlevels;level++){
    for (uint k=0;k<nops;k++)
       Zcoef(k,level)=complex<double>(10.0*get_random_minusone_to_one(),
                                      6.3*get_random_minusone_to_one());}
 CVector Vevs(nops);
 for (uint k=0;k<nops;k++)
    Vevs[k]=complex<double>((k+3)*get_random_minusone_to_one(),(k+1)*get_random_minusone_to_one());

 ofstream fout(maplefile.c_str());
 fout.precision(16);
 fout << "Nops:="<<nops<<":"<<endl;
 if (vevs) fout << "VEVS:=true:"<<endl;
 else fout << "VEVS:=false:"<<endl;
 fout << "tmin:="<<tmin<<":"<<endl;
 fout << "tmax:="<<tmax<<":"<<endl<<endl;
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
          complex<double> ztemp(0.0,0.0);
          for (uint level=0;level<nlevels;level++)
             ztemp+=Zcoef(row,level)*conjugate(Zcoef(col,level))*exp(-energies[level]*timeval);
          if (vevs) ztemp+=Vevs[row]*conjugate(Vevs[col]);
          fout << "CorrData["<<timeval<<"]["<<row<<","<<col<<"]:=["<<endl;
          for (uint bin=0;bin<nbins;bin++){
             bins_re[bin]=ztemp.real()*(1.0+0.001*get_random_minusone_to_one());
             bins_im[bin]=ztemp.imag()*(1.0+0.001*get_random_minusone_to_one());
             fout <<bins_re[bin]<<"+I*("<<bins_im[bin]<<"),"<<endl;}
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
}

#else

void make_fake_data_to_rotate(MCObsHandler& MC, CorrelatorMatrixInfo& cormat,
                              const string& maplefile, 
                              const string& binfile, uint tmin, uint tmax)
{
 uint nops=cormat.getNumberOfOperators();
 bool vevs=cormat.isVEVSubtracted();
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
 fout << "Nops:="<<nops<<":"<<endl;
 if (vevs) fout << "VEVS:=true:"<<endl;
 else fout << "VEVS:=false:"<<endl;
 fout << "tmin:="<<tmin<<":"<<endl;
 fout << "tmax:="<<tmax<<":"<<endl<<endl;
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
          fout << "CorrData["<<timeval<<"]["<<row<<","<<col<<"]:=["<<endl;
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

   //  on first run, use <Setup> to make fake data bins
   //  on second run, use sigmond_batch and do the correlator
   //   bin rotations, save to a file
   //  on third run, use <Finish> to read the data bins
   //   and do the checks

void testRotateCorrelator(XMLHandler& xml_in, int taskcount)
{
 if (xml_tag_count(xml_in,"TestRotateCorrelator")==0)
 return;
/*
 cout << endl << "Starting TestRotateCorrelator"<<endl;
 XMLHandler xmlq(xml_in,"TestRotateCorrelator");
 int stage=0;
 if (xmlq.count_among_children("Setup")==1) stage=1;
 else if (xmlq.count_among_children("Finish")==1) stage=2;
 else throw(std::invalid_argument("Invalid XML in testRotateCorrelator"));

 if (stage==1){
 try{

 XMLHandler xmlr(xmlq,"Setup");
 MCObsGetHandler MCOH(xmlr); 
 MCObsHandler MC(MCOH);
 if (xmlr.count_among_children("TweakEnsemble")==1){
    XMLHandler xmlk(xmlr,"TweakEnsemble");
    int rebin;
    if (xmlreadifchild(xmlk,"Rebin",rebin)){
       if (rebin>1) MC.setRebin(rebin);}
    vector<int> ovec;
    if (xmlreadifchild(xmlk,"Omissions",ovec)){
       set<int> omissions(ovec.begin(),ovec.end());
       if (!omissions.empty()) MC.addOmissions(omissions);}}
 if (xmlr.count_among_children("Bootstrapper")==1){
    XMLHandler xmlb(xmlr,"Bootstrapper");
    int num_resamplings=1024;
    unsigned long bootseed=0, bootskip=64;
    bool precompute=false;
    xmlreadifchild(xmlb,"NumberResamplings",num_resamplings);
    xmlreadifchild(xmlb,"Seed",bootseed);
    xmlreadifchild(xmlb,"BootSkip",bootskip);
    if (xmlb.count_among_children("Precompute")==1) precompute=true;
    MC.setBootstrapper(num_resamplings,bootseed,bootskip,precompute);}
 string maplefile;
 xmlreadchild(xmlr,"MapleFile",maplefile);
 string binfile;
 xmlreadchild(xmlr,"BinFile",binfile);
 uint mintimesep,maxtimesep;
 xmlreadchild(xmlr,"MinTimeSep",mintimesep);
 xmlreadchild(xmlr,"MaxTimeSep",maxtimesep);

 cout << endl<<endl;
 cout << "Number of measurements in ensemble = "<<MC.getNumberOfMeasurements()<<endl;
 cout << "Ensemble ID = "<<MCOH.getEnsembleId()<<endl;
 cout << "Number of bins = "<<MC.getNumberOfBins()<<endl;
 cout << "Maple file = "<<maplefile<<endl;
 MC.setToBootstrapMode();

 XMLHandler xmld(xmlq,"RotatedCorrelatorMatrix");
 XMLHandler xmlc(xmld,"CorrelatorMatrixInfo");
 CorrelatorMatrixInfo cormat(xmlc);
 if (!cormat.isHermitian()){
    throw(std::invalid_argument("CorrelatorMatrix must be Hermitian for rotation"));}

             //   make the fake data and putBins into memory
 make_fake_data_to_rotate(MC,cormat,maplefile,binfile,mintimesep,maxtimesep);
 }
  catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}}

 if (stage==2){
 try{

 ArgsHandler xmlr(xmlq,"Finish");
 XMLHandler xmlm(xmlq,"Finish");
 MCObsGetHandler MCOH(xmlm); 
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
    const RVector& bins=MC.getBins(obskey);
    cout << bins[binval]<<endl;}

#if defined COMPLEXNUMBERS
 if (vevs){
    cout << "VEVs for bin = "<<binval<<endl;
    for (uint col=0;col<nlevels;col++){
       opinfo.resetGenIrrepIDIndex(col);
       MCObsInfo obskey(opinfo,RealPart);
       const RVector& bins_re=MC.getBins(obskey);
     //  obskey.setToImaginaryPart();
     //  const RVector& bins_im=MC.getBins(obskey);
    //   cout << "("<<bins_re[binval]<<", "<<bins_im[binval]<<")"<<endl;}
       cout <<bins_re[binval]<<endl;}
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
       const RVector& bins=MC.getBins(obskey);
       cout << bins[binval]<<endl;}
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
*/
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

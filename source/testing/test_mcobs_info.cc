#include "xml_handler.h"
#include "mcobs_info.h"
#include <ctime>
#include <map>

using namespace std;


void run2_an_obs(const MCObsInfo& mcobs)
{
 cout << endl<<endl<<"  ***************************** "<<endl<<endl;
 cout << "Long output:"<<endl;
 cout << mcobs.output(true)<<endl;
 cout << "Short output:"<<endl;
 string mcstr=mcobs.output(false);
 cout << mcstr<<endl;

 cout << " isVacuum(): "<<mcobs.isVacuum()<<endl;
 cout << " isVEV(): "<<mcobs.isVEV()<<endl;
 cout << " isCorrelatorAtTime: "<<mcobs.isCorrelatorAtTime()<<endl;
 cout << " isHermitianCorrelatorAtTime: "<<mcobs.isHermitianCorrelatorAtTime()<<endl;
 cout << " isRealPart: "<<mcobs.isRealPart()<<endl;
 cout << " isImaginaryPart: "<<mcobs.isImaginaryPart()<<endl;
 cout << " isSimple(): "<<mcobs.isSimple()<<endl;
 cout << " isNonSimple(): "<<mcobs.isNonSimple()<<endl;
 cout << " isPrimary(): "<<mcobs.isPrimary()<<endl;
 cout << " isSecondary(): "<<mcobs.isSecondary()<<endl;
 cout << " isBasicLapH(): "<<mcobs.isBasicLapH()<<endl;
 cout << " isGenIrrep(): "<<mcobs.isGenIrrep()<<endl;

 try{
    OperatorInfo vop(mcobs.getVEVInfo());
    cout << "getVEVInfo: "<<endl<<vop.output()<<endl;}
 catch(const std::exception& xx){
    cout << "getVEVInfo caught exception"<<endl;}

 try{
    OperatorInfo vop;
    mcobs.getVEVInfo(vop);
    cout << "getVEVInfo: "<<endl<<vop.output()<<endl;}
 catch(const std::exception& xx){
    cout << "getVEVInfo caught exception"<<endl;}

 try{
    CorrelatorAtTimeInfo ct(mcobs.getCorrelatorAtTimeInfo());
    cout << "getCorrelatorAtTimeInfo:"<<endl<<ct.output()<<endl;}
 catch(const std::exception& xx){
    cout << "getCorrelatorAtTimeInfo caught exception"<<endl;}

 try{
    OperatorInfo srcop(mcobs.getCorrelatorSourceInfo());
    cout << "getCorrelatorSourceInfo:" <<endl<<srcop.output()<<endl;}
 catch(const std::exception& xx){
    cout << "getCorrelatorSourceInfo caught exception"<<endl;}

 try{
    OperatorInfo snkop(mcobs.getCorrelatorSinkInfo());
    cout << "getCorrelatorSinkInfo:" <<endl<<snkop.output()<<endl;}
 catch(const std::exception& xx){
    cout << "getCorrelatorSinkInfo caught exception"<<endl;}

 try{
    unsigned int timeindex=mcobs.getCorrelatorTimeIndex();
    cout << "getCorrelatorTimeIndex:" <<timeindex<<endl;}
 catch(const std::exception& xx){
    cout << "getCorrelatorTimeIndex caught exception"<<endl;}

 try{
    CorrelatorInfo cinfo(mcobs.getCorrelatorInfo());
    cout << "getCorrelatorInfo:" <<endl<<cinfo.output()<<endl;}
 catch(const std::exception& xx){
    cout << "getCorrelatorInfo caught exception"<<endl;}

 try{
    string obsname(mcobs.getObsName());
    cout << "getObsName:" <<obsname<<endl;}
 catch(const std::exception& xx){
    cout << "getObsName caught exception"<<endl;}

 try{
    uint index=mcobs.getObsIndex();
    cout << "getObsIndex:" <<index<<endl;}
 catch(const std::exception& xx){
    cout << "getObsIndex caught exception"<<endl;}

}


void run_an_obs(const MCObsInfo& mcobs, const MCObsInfo& comparemcobs)
{
 run2_an_obs(mcobs);
 cout << "mcobs1==mcobs2? :"<<  (mcobs==comparemcobs) <<endl;
 cout << "mcobs1!=mcobs2? :"<<  (mcobs!=comparemcobs) <<endl; 
 cout << "mcobs1<mcobs2? :"<<  (mcobs<comparemcobs) <<endl; 
 cout << endl<<endl;
}







bool check_an_obs(XMLHandler& xmlobs)
{
 cout << endl<<endl<<"  ***************************** "<<endl<<endl;
 cout << "Input XML: "<<xmlobs.output()<<endl<<endl;
 try{
 MCObsInfo mcobs(xmlobs);
 cout << "Long output:"<<endl;
 cout << mcobs.output(true)<<endl;
 cout << "Short output:"<<endl;
 string mcstr=mcobs.output(false);
 cout << mcstr<<endl;
 cout << "isImagDiagOfHermCorr = "<<mcobs.isImagDiagOfHermCorr()<<endl;
 cout << "hasNoRelatedFlip = "<<mcobs.hasNoRelatedFlip()<<endl;
 cout << "time flipped: "<<mcobs.getTimeFlipped().output()<<endl;
 
 bool vac=(xmlobs.count("Vacuum")==1);
 bool vev=(xmlobs.count("VEV")==1);
 bool corr=(xmlobs.count("Correlator")==1);
 bool secondary=(xmlobs.count("ObsName")==1);
 bool primary=!secondary;
 bool hermcorr=(corr && xmlobs.count("HermitianMatrix")==1);
 bool subvev=(xmlobs.count("SubtractVEV")==1);
 bool simple=(vac)||(vev)||(corr && (!subvev))||(secondary && xmlobs.count("Simple")==1);
 bool nonsimple=!simple;
 bool realpart;
 if (xmlobs.count("Arg")==0)
    realpart=true;
 else{
    XMLHandler xmlarg(xmlobs,"Arg");
    string nv=xmlarg.get_text_content();
    realpart=(nv=="RealPart")||(nv=="Re");}
 bool imagpart=!realpart;
 uint blcount=xmlobs.count("Operator")+xmlobs.count("BLOperator")
             +xmlobs.count("OperatorString")+xmlobs.count("BLOperatorString");
 uint gicount=xmlobs.count("GIOperator")+xmlobs.count("GIOperatorString");
 bool basiclaph=(gicount==0)&&(blcount>0);
 bool genirrep=(blcount==0)&&(gicount>0);
 bool samesinksource=corr ? CorrelatorInfo(xmlobs).isSinkSourceSame() : false;

 bool result=true;
 bool checker;
 checker=mcobs.isVacuum(); 
 if (checker!=vac) {cout <<"vac mismatch"<<endl; result=false;}
 checker=mcobs.isVEV();
 if (checker!=vev) {cout <<"vev mismatch"<<endl; result=false;}
 checker=mcobs.isCorrelatorAtTime();
 if (checker!=corr) {cout <<"corr mismatch"<<endl; result=false;}
 checker=mcobs.isHermitianCorrelatorAtTime();
 if (checker!=hermcorr) {cout <<"hermcorr mismatch"<<endl; result=false;}
 checker=mcobs.isRealPart();
 if (checker!=realpart) {cout <<"realpart mismatch"<<endl; result=false;}
 checker=mcobs.isImaginaryPart();
 if (checker!=imagpart) {cout <<"imagpart mismatch"<<endl; result=false;}
 checker=mcobs.isSimple();
 if (checker!=simple) {cout <<"simple mismatch"<<endl; result=false;}
 checker=mcobs.isNonSimple(); 
 if (checker!=nonsimple) {cout <<"nonsimple mismatch"<<endl; result=false;}
 checker=mcobs.isPrimary();
 if (checker!=primary) {cout <<"primary mismatch"<<endl; result=false;}
 checker=mcobs.isSecondary();
 if (checker!=secondary) {cout <<"secondary mismatch"<<endl; result=false;}
 checker=mcobs.isBasicLapH();
 if (checker!=basiclaph) {cout <<"basiclaph mismatch"<<endl; result=false;}
 checker=mcobs.isGenIrrep();
 if (checker!=genirrep) {cout <<"genirrep mismatch"<<endl; result=false;}
 ComplexArg arg=(realpart)?RealPart:ImaginaryPart;

 try{
    OperatorInfo vop(mcobs.getVEVInfo());
    if (!vev) {cout << "problem getting vev"<<endl; result=false;}
    MCObsInfo mcheck(vop,arg);
    if (mcheck!=mcobs) {cout << "problem getting vev"<<endl; result=false;}}
 catch(const std::exception& xx){
    if (vev) {cout << "problem getting vev"<<endl; result=false;}}

 try{
    OperatorInfo vop;
    mcobs.getVEVInfo(vop);
    if (!vev) {cout << "problem getting vev"<<endl; result=false;}
    MCObsInfo mcheck(vop,arg);
    if (mcheck!=mcobs) {cout << "problem getting vev"<<endl; result=false;}}
 catch(const std::exception& xx){
    if (vev) {cout << "problem getting vev"<<endl; result=false;}}

 try{
    CorrelatorAtTimeInfo ct(mcobs.getCorrelatorAtTimeInfo());
    if (!corr) {cout << "problem getting corr"<<endl; result=false;}
    MCObsInfo mcheck(ct,arg);
    if (mcheck!=mcobs) {cout << "problem getting corr"<<endl; result=false;}}
 catch(const std::exception& xx){
    if (corr) {cout << "problem getting corr"<<endl; result=false;}}

 try{
    OperatorInfo srcop(mcobs.getCorrelatorSourceInfo());
    if (!corr) {cout << "problem getting corr source"<<endl; result=false;}}
 catch(const std::exception& xx){
    if (corr) {cout << "problem getting corr source"<<endl; result=false;}}

 try{
    OperatorInfo snkop(mcobs.getCorrelatorSinkInfo());
    if (!corr) {cout << "problem getting corr sink"<<endl; result=false;}}
 catch(const std::exception& xx){
    if (corr) {cout << "problem getting corr sink"<<endl; result=false;}}

 try{
    unsigned int timeindex=mcobs.getCorrelatorTimeIndex();
    if (!corr) {cout << "problem getting corr time"<<endl; result=false;}
    uint tcheck;
    xmlread(xmlobs,"TimeIndex",tcheck,"tester"); 
    if (timeindex!=tcheck) {cout << "problem getting time index"<<endl; result=false;}}
 catch(const std::exception& xx){
    if (corr) {cout << "problem getting corr time"<<endl; result=false;}}

 try{
    OperatorInfo srcop(mcobs.getCorrelatorSourceInfo());
    OperatorInfo snkop(mcobs.getCorrelatorSinkInfo());
    unsigned int timeindex=mcobs.getCorrelatorTimeIndex();
    if (!corr) {cout << "problem getting corr source"<<endl; result=false;}
    MCObsInfo mcheck(snkop,srcop,timeindex,hermcorr,arg,subvev);
    if (mcheck!=mcobs) {cout << "problem getting corr"<<endl; result=false;}}
 catch(const std::exception& xx){
    if (corr) {cout << "problem getting corr source/sink/time"<<endl; result=false;}}

 try{
    CorrelatorInfo cinfo(mcobs.getCorrelatorInfo());
    if (!corr) {cout << "problem getting corr info"<<endl; result=false;}
    unsigned int timeindex=mcobs.getCorrelatorTimeIndex();
    MCObsInfo mcheck(cinfo,timeindex,hermcorr,arg,subvev); 
    if (mcheck!=mcobs) {cout << "problem getting corr"<<endl; result=false;}}
 catch(const std::exception& xx){
    if (corr) {cout << "problem getting corr info"<<endl; result=false;}}

 try{
    string obsname(mcobs.getObsName());
    if (!secondary) {cout << "problem getting obsname"<<endl; result=false;}}
 catch(const std::exception& xx){
    if (secondary) {cout << "problem getting obsname"<<endl; result=false;}}

 try{
    uint index=mcobs.getObsIndex();
    if (!secondary) {cout << "problem getting obs index"<<endl; result=false;}
    uint tcheck;
    xmlread(xmlobs,"Index",tcheck,"tester");
    if (index!=tcheck) {cout << "problem getting obs index"<<endl; result=false;}}
 catch(const std::exception& xx){
    if (secondary) {cout << "problem getting obs index"<<endl; result=false;}}

 try{
    string obsname(mcobs.getObsName());
    uint index=mcobs.getObsIndex();
    if (!secondary) {cout << "problem getting obsname"<<endl; result=false;}
    MCObsInfo mcheck(obsname,index,simple,arg);
    if (mcheck!=mcobs) {cout << "problem getting obsname"<<endl; result=false;}}
 catch(const std::exception& xx){
    if (secondary) {cout << "problem getting obsname"<<endl; result=false;}}


 MCObsInfo mcheck2(mcobs);
 if (mcheck2!=mcobs) {cout << "constructor failed"<<endl; result=false;}
 if (!vac){
    if (realpart) mcheck2.setToImaginaryPart();
    else mcheck2.setToRealPart();
    if (mcheck2==mcobs) {cout << "deliberate mismatch"<<endl; result=false;}}
 mcheck2=mcobs;
 if (mcheck2!=mcobs) {cout << "equality failed"<<endl; result=false;}
 if (secondary){
    mcheck2.resetObsIndex(mcobs.getObsIndex()+8);
    if (mcheck2==mcobs) {cout << "deliberate mismatch"<<endl; result=false;}}
 mcheck2=mcobs;
 if (mcheck2!=mcobs) {cout << "equality failed"<<endl; result=false;}

 bool bcheck=samesinksource && imagpart && hermcorr;
 if (mcobs.isImagDiagOfHermCorr() != bcheck){ cout << "isImagDiagOfHermCorr failed"<<endl; result=false;}
 
 bcheck=hermcorr && (!samesinksource);
 if (mcobs.hasNoRelatedFlip() == bcheck){ cout << "hasNoRelatedFlip failed"<<endl; result=false;}

 MCObsInfo other(mcobs);
 if (corr){     
    uint tcheck;
    xmlread(xmlobs,"TimeIndex",tcheck,"tester"); 
    other=MCObsInfo(CorrelatorInfo(xmlobs),tcheck,hermcorr,realpart ? RealPart : ImaginaryPart, subvev);}
 if (mcobs!=other){ cout << "getTimeFlipped failed"<<endl; result=false;}

// cout << "numints = "<<mcobs.numints()<<endl;
// cout << "number bytes = "<<mcobs.numbytes()<<endl;

// unsigned int buf[24];
// mcobs.copyTo(buf);

// MCObsInfo mcobsBB(buf);
// cout << "Check copyTo and create: "<<(mcobs==mcobsBB)<<endl;


 if (result)
    cout << "Test OK"<<endl<<endl; 
 else
    cout << "Test FAILURE"<<endl<<endl;
 return result;}
 catch(const exception& xx){
    cout << "Exception in creation: test FAILURE"<<endl;
    cout << xx.what()<<endl;
    return false;}
}




void testMCObsInfo(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestMCObsInfo")==0)
 return;

 cout << endl<<endl<<"***************************************************"<<endl<<endl;
 cout << "Testing MCObsInfo"<<endl;

 MCObsInfo mc1;
 cout << mc1.output()<<endl;

 list<XMLHandler> mcobsxml=xml_in.find("DoATest");
 cout << "Found "<<mcobsxml.size()<<" MCObsInfo XML tags"<<endl<<endl;

 for (list<XMLHandler>::iterator it=mcobsxml.begin();it!=mcobsxml.end();it++){
    try{
       cout << endl<<endl<<" ********************"<<endl<<endl;
       cout << "Input XML: "<<it->output()<<endl;
       MCObsInfo mc2(*it);
       cout << mc2.output(true)<<endl;
       cout << mc2.output(false)<<endl;
       cout << "Set to RealPart"<<endl;
       mc2.setToRealPart();
       cout << mc2.output(false)<<endl;
       cout << "Set to ImaginaryPart"<<endl;
       mc2.setToImaginaryPart();
       cout << mc2.output(false)<<endl;
       }
    catch(const std::invalid_argument& errmsg){
       cout << "Whoops! Invalid XML: "<<errmsg.what()<<endl;}}

 cout << endl<<endl<<" *************run an obs(mc1,mcs1)*****"<<endl<<endl;
 run_an_obs(mc1,mc1);

 vector<OperatorInfo> qcdops(7);
 qcdops[0]=OperatorInfo("pion P=(0,1,0) A2m_1 DDL_8");  
 qcdops[1]=OperatorInfo("isosinglet_eta_eta A1gp_1 [P=(0,0,0) A1gp SD_2] [P=(0,0,0) A1gp SD_2]");   
 qcdops[2]=OperatorInfo("glueball P=(0,0,0) A1gp_1 TrEig");
 qcdops[3]=OperatorInfo("kaon P=(0,0,0) T1u_1 SS_0");
 qcdops[4]=OperatorInfo("nucleon P=(0,0,0) G1g_1 TDT_29");
 qcdops[5]=OperatorInfo("isodoublet_pion_nucleon G1g_1 [P=(0,0,1) A2m SS_1] [P=(0,0,-1) G1 SS_0]"); 
 qcdops[6]=OperatorInfo("isodoublet S=1 PSQ=4 G1g_1 Frappo",OperatorInfo::GenIrrep); 
 
 {CorrelatorAtTimeInfo corref(qcdops[2],qcdops[4],9);
 MCObsInfo mcref(corref);
 CorrelatorAtTimeInfo cortmp(qcdops[2],qcdops[5],12);
 MCObsInfo mc2(cortmp);
 mcref=mc2;

 for (int i=0;i<7;i++)
 for (int j=0;j<7;j++){
    CorrelatorAtTimeInfo corr(qcdops[i],qcdops[j],12);
    MCObsInfo mcobs(corr);
    run_an_obs(mcobs,mcref);
    MCObsInfo mcobsIm(corr);
    mcobsIm.setToImaginaryPart();
    run_an_obs(mcobsIm,mcref);

    MCObsInfo mcobs2(qcdops[i],qcdops[j],12);
    if (mcobs!=mcobs2) cout << "MISMATCH 1"<<endl;
    CorrelatorInfo cor(qcdops[i],qcdops[j]);
    mcobs2=MCObsInfo(cor,12);
    if (mcobs!=mcobs2) cout << "MISMATCH 2"<<endl;
    if (mcobs2.getCorrelatorAtTimeInfo()!=corr) cout << "MISMATCH 3"<<endl;
    if (mcobs2.getCorrelatorInfo()!=cor) cout << "MISMATCH 4"<<endl;
    if (mcobs2.getCorrelatorSourceInfo()!=qcdops[j]) cout << "MISMATCH 5"<<endl;
    if (mcobs2.getCorrelatorSinkInfo()!=qcdops[i]) cout << "MISMATCH 6"<<endl;
    if (mcobs2.getCorrelatorTimeIndex()!=12) cout << "MISMATCH 7"<<endl;

    CorrelatorAtTimeInfo ttt(mcobs);
    if (ttt!=corr) cout << "MISMATCH 12"<<endl;

    CorrelatorAtTimeInfo tmp1(corref);
    CorrelatorInfo tmp2(tmp1.getCorrelator());
    tmp1.resetTimeSeparation(45);
    mcobs.getCorrelatorAtTimeInfo(tmp1);
    if (tmp1!=corr) cout << "MISMATCH 8"<<endl;
    mcobs.getCorrelatorInfo(tmp2);
    if (tmp2!=cor) cout << "MISMATCH 9"<<endl;}
 }


 {MCObsInfo mcref(qcdops[2]);
  for (int i=0;i<6;i++){
    MCObsInfo mcobs(qcdops[i]);
    run_an_obs(mcobs,mcref);
    MCObsInfo mcobs2;
    mcobs2=mcobs;
    if (mcobs2.getVEVInfo()!=qcdops[i]) cout << "MISMATCH 10"<<endl;
    OperatorInfo vv;
    mcobs2.getVEVInfo(vv);
    if (vv!=qcdops[i]) cout << "MISMATCH 11"<<endl;}
 }


 list<XMLHandler> xmlops;
 XMLHandler xmlop;
 xmlop.set_root("BLOperatorString","pion P=(0,1,0) A2m_1 DDL_8");
 xmlops.push_back(xmlop); 
 xmlop.set_root("BLOperatorString","isosinglet_eta_eta A1gp_1 [P=(0,0,0) A1gp SD_2] [P=(0,0,0) A1gp SD_2]");   
 xmlops.push_back(xmlop); 
 xmlop.set_root("BLOperatorString","glueball P=(0,0,0) A1gp_1 TrEig");
 xmlops.push_back(xmlop); 
 xmlop.set_root("BLOperatorString","kaon P=(0,0,0) T1u_1 SS_0");
 xmlops.push_back(xmlop); 
 xmlop.set_root("BLOperatorString","nucleon P=(0,0,0) G1g_1 TDT_29");
 xmlops.push_back(xmlop); 
 xmlop.set_root("BLOperatorString","isodoublet_pion_nucleon G1g_1 [P=(0,0,1) A2m SS_1] [P=(0,0,-1) G1 SS_0]"); 
 xmlops.push_back(xmlop); 
 xmlop.set_root("GIOperatorString","isotriplet P=(0,0,0) A1um_1 Willow 24");
 xmlops.push_back(xmlop); 
 xmlop.set_root("GIOperatorString","isosinglet P=(0,1,0) B1p_1 Kincardine 3");
 xmlops.push_back(xmlop); 
 xmlop.set_root("GIOperatorString","isodoublet P=(0,0,1) E_1 ShellConch");
 xmlops.push_back(xmlop); 



 uint ncheck=0;
 uint nfail=0;
 cout <<endl<<endl<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl<<endl;
 XMLHandler xmlcheck;
 xmlcheck.set_root("MCObservable");
 xmlcheck.put_child("Vacuum");
 if (!check_an_obs(xmlcheck)) nfail++;
 ncheck++;

 for (list<XMLHandler>::iterator it=xmlops.begin();it!=xmlops.end();it++){
    xmlcheck.set_root("MCObservable");
    XMLHandler xmlv("VEV");
    xmlv.put_child(*it);
    xmlcheck.put_child(xmlv);
    xmlcheck.put_child("Arg","RealPart");
    if (!check_an_obs(xmlcheck)) nfail++;
    ncheck++;}

 for (list<XMLHandler>::iterator it=xmlops.begin();it!=xmlops.end();it++){
    xmlcheck.set_root("MCObservable");
    XMLHandler xmlv("VEV");
    xmlv.put_child(*it);
    xmlcheck.put_child(xmlv);
    xmlcheck.put_child("Arg","Im");
    if (!check_an_obs(xmlcheck)) nfail++;
    ncheck++;}

 for (uint hermcnt=0;hermcnt<2;hermcnt++)
 for (uint subvevcnt=0;subvevcnt<2;subvevcnt++)
 for (uint realcnt=0;realcnt<2;realcnt++)
 for (uint timesep=3;timesep<=25;timesep+=3){
 for (list<XMLHandler>::iterator itsrc=xmlops.begin();itsrc!=xmlops.end();itsrc++){
 for (list<XMLHandler>::iterator itsnk=xmlops.begin();itsnk!=xmlops.end();itsnk++){
    xmlcheck.set_root("MCObservable");
    XMLHandler xmlsrc("Source");
    xmlsrc.put_child(*itsrc);
    XMLHandler xmlsnk("Sink");
    xmlsnk.put_child(*itsnk);
    XMLHandler xmlc("Correlator");
    xmlc.put_child(xmlsrc);
    xmlc.put_child(xmlsnk);
    xmlc.put_child("TimeIndex",make_string(timesep));
    if (hermcnt>0) xmlc.put_child("HermitianMatrix");
    if (subvevcnt>0) xmlc.put_child("SubtractVEV");
    xmlcheck.put_child(xmlc);
    if (realcnt==0) xmlcheck.put_child("Arg","RealPart");
    else xmlcheck.put_child("Arg","Im");
    if (!check_an_obs(xmlcheck)) nfail++;
    ncheck++;}}}

 list<string> obsnames;
 obsnames.push_back("Kincardine");
 obsnames.push_back("Goderich");
 obsnames.push_back("Toronto834%!#");
 obsnames.push_back("NiagaraFalls");
 list<uint> indices;
 indices.push_back(4);
 indices.push_back(22);
 indices.push_back(11);
 indices.push_back(39);
 indices.push_back(2);
 indices.push_back(8);
 indices.push_back(1);

 for (uint simpcnt=0;simpcnt<2;simpcnt++)
 for (uint realcnt=0;realcnt<2;realcnt++)
 for (list<string>::iterator itn=obsnames.begin();itn!=obsnames.end();itn++){
 for (list<uint>::iterator it=indices.begin();it!=indices.end();it++){
    xmlcheck.set_root("MCObservable");
    xmlcheck.put_child("ObsName",*itn);
    xmlcheck.put_child("Index",make_string(*it));
    if (simpcnt>0) xmlcheck.put_child("Simple");
    if (realcnt==0) xmlcheck.put_child("Arg","RealPart");
    else xmlcheck.put_child("Arg","Im");
    if (!check_an_obs(xmlcheck)) nfail++;
    ncheck++;}}


 cout << endl<<endl<<"Total number of checks = "<<ncheck<<endl
      <<"Total number of check failures = "<<nfail<<endl<<endl;


}

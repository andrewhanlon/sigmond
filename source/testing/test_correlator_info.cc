#include "xml_handler.h"
#include "correlator_info.h"
#include "correlator_matrix_info.h"

using namespace std;

void run_a_corr(const CorrelatorInfo& corr, const CorrelatorInfo& comparecorr)
{
 cout << endl<<endl<<"  ***************************** "<<endl<<endl;
 cout << "Long output:"<<endl;
 cout << corr.output(true)<<endl;
 cout << "Short output:"<<endl;
 cout << corr.output(false)<<endl<<endl;
 cout << "SOURCE"<<endl;
 OperatorInfo temp;
 temp=corr.getSource();
 cout << temp.output()<<endl;
 cout << "SINK"<<endl;
 temp=corr.getSink();
 cout << temp.output()<<endl;

 CorrelatorInfo corrflip(corr.getTimeFlipped());
 cout << "CORR FLIPPED:"<<endl;
 cout << corrflip.output()<<endl;

 cout << "corr1==corr2? :"<<  (corr==comparecorr) <<endl;
 cout << "corr1!=corr2? :"<<  (corr!=comparecorr) <<endl; 
 cout << "corr1<corr2? :"<<  (corr<comparecorr) <<endl; 
 cout << endl<<endl;
}


void run_a_corrtime(const CorrelatorAtTimeInfo& corr, const CorrelatorAtTimeInfo& comparecorr)
{
 cout << endl<<endl<<"  ***************************** "<<endl<<endl;
 cout << "Long output:"<<endl;
 cout << corr.output(true)<<endl;
 cout << "Short output:"<<endl;
 cout << corr.output(false)<<endl<<endl;
 CorrelatorInfo cor(corr.getCorrelator());
 cout << "CORRELATOR:"<<endl;
 cout << cor.output()<<endl;
 cout << "SOURCE"<<endl;
 OperatorInfo temp;
 temp=corr.getSource();
 cout << temp.output()<<endl;
 cout << "SINK"<<endl;
 temp=corr.getSink();
 cout << temp.output()<<endl;

 cout << "Time = "<<corr.getTimeSeparation()<<endl;
 CorrelatorAtTimeInfo acorr(corr);
 acorr.resetTimeSeparation(3);
 cout << "Time reset to 3 = "<<acorr.getTimeSeparation()<<endl;
 acorr.resetTimeSeparation(17);
 cout << "Time reset to 17 = "<<acorr.getTimeSeparation()<<endl;
 acorr.resetTimeSeparation(41);
 cout << "Time reset to 41 = "<<acorr.getTimeSeparation()<<endl;
 cout << acorr.output()<<endl;
 acorr.resetVEVSubtracted(true);
 cout << "VEV reset to true = "<<acorr.isVEVsubtracted()<<endl;
 cout << acorr.output()<<endl;
 acorr.resetVEVSubtracted(false);
 cout << "VEV reset to false = "<<acorr.isVEVsubtracted()<<endl;
 cout << acorr.output()<<endl;

 cout << "corr1==corr2? :"<<  (corr==comparecorr) <<endl;
 cout << "corr1!=corr2? :"<<  (corr!=comparecorr) <<endl; 
 cout << "corr1<corr2? :"<<  (corr<comparecorr) <<endl; 
 cout << endl<<endl;
}


void testCorrelatorInfo(XMLHandler& xml_in)
{

 if (xml_tag_count(xml_in,"TestCorrelatorInfo")==0)
 return;

 vector<OperatorInfo> qcdops(14);
 qcdops[0]=OperatorInfo("pion P=(0,1,0) A2m_1 DDL_8");  
 qcdops[1]=OperatorInfo("kaon P=(0,0,0) T1u_1 DDL_18");
 qcdops[2]=OperatorInfo("isosinglet_eta_eta A1gp_1 [P=(0,0,0) A1gp SD_2] [P=(0,0,0) A1gp SD_2]");   
 qcdops[3]=OperatorInfo("glueball P=(0,0,0) A1gp_1 TrEig");
 qcdops[4]=OperatorInfo("isosinglet_eta_eta A1gp_1 [P=(0,0,0) A1up SS_0] [P=(0,0,0) A1up SS_0]");    
 qcdops[5]=OperatorInfo("isosinglet_kaon_kbar A1gp_1 [P=(0,0,0) A1g SD_0] [P=(0,0,0) A1g TDO_3]");
 qcdops[6]=OperatorInfo("pion P=(0,0,0) T1up_1 DDL_8");
 qcdops[7]=OperatorInfo("isoquintet_pion_pion Ep_1 CG_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B1p SS_1]");
 qcdops[8]=OperatorInfo("isoquintet_pion_pion Ep_1 CG_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B2p SS_2]");
 qcdops[9]=OperatorInfo("kaon P=(0,0,0) T1u_1 SS_0");
 qcdops[10]=OperatorInfo("isodoublet_kaon_pion T1u_1 [P=(2,0,0) A2 SS_1] [P=(-2,0,0) A2m SS_1]");
 qcdops[11]=OperatorInfo("pion P=(0,0,1) A2m_1 DDL_6");  
 qcdops[12]=OperatorInfo("nucleon P=(0,0,0) G1g_1 TDT_29");
 qcdops[13]=OperatorInfo("isodoublet_pion_nucleon G1g_1 [P=(0,0,1) A2m SS_1] [P=(0,0,-1) G1 SS_0]"); 
 
OperatorInfo tempop("pion P=(0,1,0) A2m_1 DDL_8");
cout << tempop.output(true)<<endl;

for (uint k=0;k<qcdops.size();k++){ cout << "op "<<k<<endl;cout<<qcdops[k].output()<<endl;}

 {CorrelatorInfo corref(qcdops[6],qcdops[3]);
 CorrelatorInfo cortmp(qcdops[3],qcdops[7]);
 corref=cortmp;

 cout << "Doing run_a_corr"<<endl;
 for (int i=0;i<14;i++)
 for (int j=0;j<14;j++){
    CorrelatorInfo corr(qcdops[i],qcdops[j]);
    run_a_corr(corr,corref);}}

 {CorrelatorAtTimeInfo corref(qcdops[6],qcdops[3],9);
 CorrelatorAtTimeInfo cortmp(qcdops[3],qcdops[7],12);
 corref=cortmp;

 cout << "Doing run_a_corrtime"<<endl;
 for (int i=0;i<14;i++)
 for (int j=0;j<14;j++){
    CorrelatorAtTimeInfo corr(qcdops[i],qcdops[j],12);
    run_a_corrtime(corr,corref);}}


 XMLHandler xmlop(xml_in,"TestInput1");
 CorrelatorInfo corr2(xmlop);
 cout << corr2.output()<<endl;
 cout << corr2.output(false)<<endl;

 XMLHandler xmlop2(xml_in,"TestInput2");
 CorrelatorAtTimeInfo corrtime2(xmlop2);
 cout << corrtime2.output()<<endl;
 cout << corrtime2.output(false)<<endl;



 vector<string> opstrings;
 opstrings.push_back("glueball P=(0,0,0) A1gp_1 TrEig");   
 opstrings.push_back("eta P=(0,0,0) A1gp_1 SD_2");         
 opstrings.push_back("phi P=(0,0,0) A1gp_1 TDO_5");        
 opstrings.push_back("isosinglet_pion_pion A1gp_1 [P=(0,0,0) A1um SS_0] [P=(0,0,0) A1um SS_0]");  
 opstrings.push_back("isosinglet_pion_pion A1gp_1 [P=(0,0,1) A2m SS_0] [P=(0,0,-1) A2m SS_0]");   
 opstrings.push_back("isosinglet_eta_eta A1gp_1 [P=(0,0,0) A1up SS_0] [P=(0,0,0) A1up SS_0]");    
 opstrings.push_back("isosinglet_eta_eta A1gp_1 [P=(0,0,0) A1gp SD_2] [P=(0,0,0) A1gp SD_2]");   
 opstrings.push_back("isosinglet_eta_eta A1gp_1 [P=(0,1,1) A2p TSD_2] [P=(0,-1,-1) A2p LSD_0]"); 
 opstrings.push_back("isosinglet_eta_phi A1gp_1 [P=(0,0,0) A1um SS_0] [P=(0,0,0) A1um SS_0]");
 opstrings.push_back("isosinglet_kaon_kbar A1gp_1 [P=(0,0,0) A1g SD_0] [P=(0,0,0) A1g TDO_3]");
 opstrings.push_back("eta P=(0,0,0) A1up_1 SS_0");  
 opstrings.push_back("eta P=(0,0,0) A1up_1 TDO_1"); 
 opstrings.push_back("phi P=(0,0,0) A1up_1 SS_0");  
 opstrings.push_back("pion P=(0,0,0) T1up_1 DDL_2");
 opstrings.push_back("pion P=(0,0,0) T1up_1 SS_0");
 opstrings.push_back("pion P=(0,0,0) T1up_1 SS_1");
 opstrings.push_back("pion P=(0,0,0) T1up_1 DDL_8");
 opstrings.push_back("pion P=(0,0,0) T1up_1 TDO_3");
 opstrings.push_back("pion P=(0,0,0) T1up_1 DDL_0");
 opstrings.push_back("pion P=(0,0,0) T1up_1 DDL_3");
 opstrings.push_back("pion P=(0,0,0) T1up_1 DDL_4");
 opstrings.push_back("pion P=(0,0,0) T1up_1 DDL_12");
 opstrings.push_back("pion P=(0,0,0) T1up_1 DDL_13");
 opstrings.push_back("pion P=(0,0,0) T1up_1 TDO_15");
 opstrings.push_back("pion P=(0,0,0) T1up_1 TDO_22");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(1,0,0) A2m SS_1] [P=(-1,0,0) A2m SS_1]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(1,0,1) A2m SS_0] [P=(-1,0,-1) A2m SS_0]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(0,0,0) T1gm SS_0] [P=(0,0,0) A1um SS_0]");
 opstrings.push_back("isotriplet_kaon_kbar T1up_1 [P=(1,0,0) A2 SS_1] [P=(-1,0,0) A2 SS_1]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(0,0,1) Em SS_1] [P=(0,0,-1) A2m LSD_1]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(0,0,0) A1um SS_0] [P=(0,0,0) T1gm SS_0]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(1,1,1) A2m SS_0] [P=(-1,-1,-1) A2m SS_0]");
 opstrings.push_back("isotriplet_phi_pion T1up_1 [P=(0,0,1) Em SS_1] [P=(0,0,-1) A2m SS_1]");
 opstrings.push_back("isotriplet_kaon_kbar T1up_1 [P=(1,0,1) A2 SS_0] [P=(-1,0,-1) A2 SS_0]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(0,0,1) A2p SS_1] [P=(0,0,-1) Ep SS_1]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(1,0,1) B1m SS_1] [P=(-1,0,-1) A2m SS_0]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(0,1,1) B2m SS_2] [P=(0,-1,-1) A2m SS_0]");
 opstrings.push_back("isotriplet_kaon_kbar T1up_1 [P=(0,0,1) E SS_2] [P=(0,0,-1) A2 SS_1]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(1,0,0) A2m LSD_3] [P=(-1,0,0) A2m SS_1]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(0,0,1) Em SS_2] [P=(0,0,-1) A2m SS_1]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(2,0,0) A2m SS_1] [P=(-2,0,0) A2m SS_1]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(0,0,1) A2m SS_1] [P=(0,0,-1) Em SS_0]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(1,0,0) A2m SS_1] [P=(-1,0,0) A2m SS_0]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(0,0,0) T1um SS_0] [P=(0,0,0) A1gm SS_0]");
 opstrings.push_back("isotriplet_kaon_kbar T1up_1 [P=(0,0,0) T1g SS_0] [P=(0,0,0) A1u SS_0]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(1,0,0) A1p SS_1] [P=(-1,0,0) A1p SS_1]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(0,0,1) A1p SS_1] [P=(0,0,-1) Ep SS_1]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(1,0,0) Ep SS_1] [P=(-1,0,0) Ep SS_1]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(0,0,0) A1up SS_0] [P=(0,0,0) T1gp SS_0]");
 opstrings.push_back("isotriplet_kaon_kbar T1up_1 [P=(1,1,1) A2 SS_0] [P=(-1,-1,-1) A2 SS_0]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(1,0,1) A2p SS_0] [P=(-1,0,-1) B1p SS_1]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(0,1,1) A2p SS_0] [P=(0,-1,-1) B2p SS_2]");
 opstrings.push_back("isotriplet_phi_pion T1up_1 [P=(1,0,1) B1m SS_1] [P=(-1,0,-1) A2m SS_0]");
 opstrings.push_back("isotriplet_phi_pion T1up_1 [P=(0,1,1) B2m SS_2] [P=(0,-1,-1) A2m SS_0]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(1,0,0) A2m SS_1] [P=(-1,0,0) A2m TSD_0]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(0,0,1) A2m SS_1] [P=(0,0,-1) Em TSD_1]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(1,1,1) Em SS_1] [P=(-1,-1,-1) A2m SS_0]");
 opstrings.push_back("isotriplet_kaon_kbar T1up_1 [P=(1,0,1) B1 SS_1] [P=(-1,0,-1) A2 SS_0]");
 opstrings.push_back("isotriplet_kaon_kbar T1up_1 [P=(0,1,1) B2 SS_3] [P=(0,-1,-1) A2 SS_0]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(0,0,1) A2m SS_1] [P=(0,0,-1) Em LSD_1]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(1,0,0) A2m TSD_2] [P=(-1,0,0) A2m TSD_1]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(1,0,1) A2m SS_1] [P=(-1,0,-1) A2m SS_0]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(0,0,0) T1gm SD_1] [P=(0,0,0) A1um SS_0]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(0,0,0) A1um TDO_1] [P=(0,0,0) T1gm SS_0]");
 opstrings.push_back("isotriplet_kaon_kbar T1up_1 [P=(1,0,1) A2 SS_0] [P=(-1,0,-1) A2 SS_1]");
 opstrings.push_back("isotriplet_kaon_kbar T1up_1 [P=(0,0,1) E SS_3] [P=(0,0,-1) A2 SS_0]");
 opstrings.push_back("isotriplet_pion_pion T1up_1 [P=(2,0,0) A2m SS_1] [P=(-2,0,0) A2m TSD_0]");
 opstrings.push_back("isotriplet_kaon_kbar T1up_1 [P=(1,1,1) A2 SS_0] [P=(-1,-1,-1) A2 SS_1]");
 opstrings.push_back("isotriplet_eta_pion T1up_1 [P=(1,1,1) Em SD_6] [P=(-1,-1,-1) A2m SS_0]");
 opstrings.push_back("pion P=(0,0,0) A1um_1 SS_0");   
 opstrings.push_back("pion P=(0,0,0) A1um_1 TDU_2");  
 opstrings.push_back("isoquintet_pion_pion A1gp_1 [P=(0,0,0) A1um SS_0] [P=(0,0,0) A1um SS_0]");    
 opstrings.push_back("isoquintet_pion_pion A1gp_1 [P=(0,0,1) A2m SS_1] [P=(0,0,-1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1gp_1 [P=(0,1,1) A2m SS_0] [P=(0,-1,-1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1gp_1 [P=(0,0,0) T1up SS_0] [P=(0,0,0) T1up SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1gp_1 [P=(1,1,1) A2m SS_0] [P=(-1,-1,-1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1gp_1 [P=(0,0,0) A1um SS_0] [P=(0,0,0) A1um SD_0]"); 
 opstrings.push_back("isoquintet_pion_pion A1gp_1 [P=(0,0,2) A2m SS_1] [P=(0,0,-2) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1gp_1 [P=(0,0,0) A1um SD_0] [P=(0,0,0) A1um SD_0]");
 opstrings.push_back("isoquintet_pion_pion A1gp_1 [P=(0,0,1) A2m SS_0] [P=(0,0,-1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1gp_1 [P=(0,1,1) A2m SS_1] [P=(0,-1,-1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1gp_1 [P=(0,0,0) T1up SS_1] [P=(0,0,0) T1up SS_1]");
 opstrings.push_back("isoquintet_pion_pion Egp_1 [P=(0,1,0) A2m SS_1] [P=(0,-1,0) A2m SS_1]");    
 opstrings.push_back("isoquintet_pion_pion Egp_1 [P=(0,1,1) A2m SS_0] [P=(0,-1,-1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Egp_1 [P=(0,0,0) T1up SS_0] [P=(0,0,0) T1up SS_0]");
 opstrings.push_back("isoquintet_pion_pion Egp_1 [P=(0,1,0) A1p SS_1] [P=(0,-1,0) A1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion Egp_1 [P=(0,0,1) Ep SS_1] [P=(0,0,-1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion Egp_1 CG_1 [P=(0,0,1) Ep SS_1] [P=(0,0,-1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion Egp_1 [P=(0,2,0) A2m SS_1] [P=(0,-2,0) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion Egp_1 [P=(0,1,0) A2m SS_0] [P=(0,-1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Egp_1 [P=(0,1,1) A2m SS_1] [P=(0,-1,-1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion Egp_1 [P=(0,0,0) T1up SS_1] [P=(0,0,0) T1up SS_1]");
 opstrings.push_back("isoquintet_pion_pion Egp_1 [P=(0,1,0) A1p SS_2] [P=(0,-1,0) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion T1up_1 [P=(1,0,0) A2m SS_1] [P=(-1,0,0) A2m SS_1]");  
 opstrings.push_back("isoquintet_pion_pion T1up_1 [P=(1,0,0) A2m SS_0] [P=(-1,0,0) A2m SS_1]");  
 opstrings.push_back("isoquintet_pion_pion T2gp_1 [P=(0,1,1) A2m SS_0] [P=(0,-1,-1) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion T2gp_1 [P=(0,0,0) T1up SS_0] [P=(0,0,0) T1up SS_0]");
 opstrings.push_back("isoquintet_pion_pion T2gp_1 [P=(1,1,1) A2m SS_0] [P=(-1,-1,-1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion T2gp_1 [P=(0,0,1) A1p SS_1] [P=(0,0,-1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion T2gp_1 [P=(1,0,0) Ep SS_1] [P=(-1,0,0) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion T2gp_1 [P=(1,0,1) A2m SS_0] [P=(-1,0,-1) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion T2gp_1 [P=(0,1,1) A2m SS_1] [P=(0,-1,-1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion T2gp_1 [P=(0,0,0) T1up SS_1] [P=(0,0,0) T1up SS_1]"); 
 opstrings.push_back("isoquintet_pion_pion T2gp_1 [P=(1,1,1) A2m SS_1] [P=(-1,-1,-1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion T2gp_1 [P=(0,0,1) A1p SS_2] [P=(0,0,-1) Ep SS_2]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SS_0] [P=(0,0,0) A2m SS_1]");  
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,0) A2m SS_1] [P=(1,0,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,-1) A2m SS_1] [P=(0,0,2) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,-1,0) A2m SS_0] [P=(1,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,-1) A2m SS_0] [P=(1,0,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,-1,-1) A2m SS_0] [P=(1,1,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-2,0,0) A2m SS_1] [P=(2,0,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SS_0] [P=(0,0,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SS_0] [P=(0,0,1) A2m TSD_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-2,-1,0) A2m SS_0] [P=(2,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) T1up SS_0] [P=(0,0,1) A1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) T1up SS_0] [P=(0,0,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) T1gm SS_0] [P=(0,0,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,-2) A2m SS_1] [P=(0,0,3) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um TDO_1] [P=(0,0,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-2,0,-1) A2m SS_0] [P=(2,0,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,0) A2m SS_1] [P=(1,0,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,0) A2m SS_1] [P=(1,0,1) B1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,0) A2m SS_1] [P=(1,0,1) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,0) T1up SS_0] [P=(0,0,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,-1) A2m SS_0] [P=(1,0,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,0) A1um SS_0] [P=(0,0,1) B2m TSD_1]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,0) A1p SS_1] [P=(1,0,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,0) A1p SS_1] [P=(1,0,1) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,0) Ep SS_1] [P=(1,0,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,-1,0) Ep SS_1] [P=(0,1,1) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,0) Ep SS_1] [P=(1,0,1) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-2,0,0) A2m SS_1] [P=(2,0,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,0) A2m SS_0] [P=(1,0,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,0) T1up DDL_8] [P=(0,0,1) Ep SS_2]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,-1) A2m SS_1] [P=(1,0,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,0) A1um TDO_1] [P=(0,0,1) B2m TSD_2]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,-1,0) A2m SS_0] [P=(1,1,1) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,0) T1up SS_0] [P=(0,0,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,0,0) A2m SS_1] [P=(1,0,1) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,0) A1um SS_0] [P=(0,0,1) B1m TSD_2]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,0,0) A1m SS_0] [P=(1,0,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,0,0) A1p SS_1] [P=(1,0,1) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,-1,0) Ep SS_1] [P=(0,1,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,0,0) Ep SS_1] [P=(1,0,1) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,-1,0) Ep SS_1] [P=(0,1,1) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,-1,-1) A2m SS_0] [P=(1,1,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,-1,0) A2m SS_1] [P=(1,1,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,0) T1up DDL_8] [P=(0,0,1) Ep SS_2]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,0,0) A2m SS_0] [P=(1,0,1) A1m TSD_1]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,0) A1um TDO_1] [P=(0,0,1) B1m TSD_3]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,0) A2m SS_1] [P=(1,0,1) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,-1,0) A2m SS_0] [P=(1,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,0) T1up SS_1] [P=(0,0,1) A1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,0) T1up SS_1] [P=(0,0,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,0) A1um SS_0] [P=(0,0,1) Em SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,-1,0) A2m SS_1] [P=(0,1,1) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,-1) A2m SS_0] [P=(1,0,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,0) A1um SS_0] [P=(0,0,1) Em TSD_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,-1,0) A1m SS_0] [P=(0,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,0) A1um SS_0] [P=(0,0,1) Em LSD_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,0) T1gm SS_0] [P=(0,0,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,0) A1p SS_1] [P=(1,0,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,-1,0) A1p SS_1] [P=(0,1,1) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,0) A1p SS_1] [P=(1,0,1) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,0) Ep SS_1] [P=(1,0,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,0) Ep SS_1] [P=(1,0,1) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,0) Ep SS_1] [P=(1,0,1) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 CG_1 [P=(0,-1,0) Ep SS_1] [P=(0,1,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 CG_1 [P=(0,-1,0) Ep SS_1] [P=(0,1,1) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 CG_1 [P=(0,-1,0) Ep SS_1] [P=(0,1,1) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-2,0,0) A2m SS_1] [P=(2,0,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,-1,-1) A2m SS_0] [P=(1,1,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,0) A2m SS_0] [P=(1,0,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,-1,0) A2m SS_1] [P=(1,1,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,0) T1up SS_1] [P=(0,0,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,0) A1um TDO_1] [P=(0,0,1) Em LSD_3]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SS_0] [P=(0,1,1) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) A2m SS_1] [P=(0,1,0) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,0) A2m SS_1] [P=(1,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,1) A2m SS_0] [P=(1,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,-1,0) A2m SS_1] [P=(0,2,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,2) A2m SS_1] [P=(0,1,-1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) T1up SS_0] [P=(0,1,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) T1up SS_0] [P=(0,1,1) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) T1up SS_0] [P=(0,1,1) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SS_0] [P=(0,1,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) A1p SS_1] [P=(0,1,0) A1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) A1p SS_1] [P=(0,1,0) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) Ep SS_1] [P=(0,1,0) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 CG_1 [P=(0,0,1) Ep SS_1] [P=(0,1,0) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SS_0] [P=(0,1,1) A2m LSD_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SS_0] [P=(0,1,1) A2m LSD_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,-1,0) A2m SS_0] [P=(1,2,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,1) A2m SS_0] [P=(1,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,-1,1) A2m SS_0] [P=(1,2,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SD_0] [P=(0,1,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) A2m SS_0] [P=(0,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,0) A2m SS_0] [P=(0,1,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,1) A2m SS_1] [P=(1,1,0) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,-1,0) A2m SS_0] [P=(0,2,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(-1,0,1) A2m SS_0] [P=(1,1,0) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) A1um SS_0] [P=(0,1,1) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) A2m SS_1] [P=(0,1,0) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) A1gm SS_0] [P=(0,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) T1up SS_0] [P=(0,1,1) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) T1up SS_0] [P=(0,1,1) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) A1p SS_1] [P=(0,1,0) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) Ep SS_1] [P=(0,1,0) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(-1,0,0) A2m SS_1] [P=(1,1,1) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) A1um SS_0] [P=(0,1,1) A1m LSD_7]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(-1,-1,0) A2m SS_0] [P=(1,2,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(-1,-1,1) A2m SS_0] [P=(1,2,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(-1,0,1) A2m SS_1] [P=(1,1,0) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) A1um SD_0] [P=(0,1,1) A1m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) A2m SS_0] [P=(0,1,0) A1m LSD_3]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,0) A2m SS_1] [P=(1,1,1) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,1) A2m SS_1] [P=(0,1,0) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,0) T1up SS_0] [P=(0,1,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,0) T1up SS_0] [P=(0,1,1) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,0) A1um SS_0] [P=(0,1,1) B2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,1) A1p SS_1] [P=(0,1,0) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,1) Ep SS_1] [P=(0,1,0) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,0) A1um SS_0] [P=(0,1,1) B2m LSD_7]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,-1,0) A2m SS_0] [P=(1,2,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,-1,1) A2m SS_0] [P=(1,2,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,0) A2m SS_0] [P=(1,1,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,1) A2m SS_0] [P=(0,1,0) A1m LSD_3]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,0) T1up DDL_8] [P=(0,1,1) A1p SS_0]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,0) T1up DDL_8] [P=(0,1,1) B1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,-1,0) A2m SS_1] [P=(0,2,1) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,2) A2m SS_1] [P=(0,1,-1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,0) T1up SS_0] [P=(0,1,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,0) T1up SS_0] [P=(0,1,1) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,0) A1um SS_0] [P=(0,1,1) B1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,1) A1p SS_1] [P=(0,1,0) A1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,1) A1p SS_1] [P=(0,1,0) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,1) Ep SS_1] [P=(0,1,0) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,0,0) A2m SS_1] [P=(1,1,1) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,0) A1um SS_0] [P=(0,1,1) B1m LSD_4]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,-1,0) A2m SS_0] [P=(1,2,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,-1,1) A2m SS_0] [P=(1,2,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,-1,0) A2m SS_0] [P=(0,2,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,2) A2m SS_0] [P=(0,1,-1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,0) T1up SS_1] [P=(0,1,1) A1p SS_0]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,0) A1um TDO_1] [P=(0,1,1) B1m TSD_2]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SS_0] [P=(1,1,1) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) A2m SS_1] [P=(1,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,0) A2m SS_1] [P=(2,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,1) A2m SS_0] [P=(2,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,2) A2m SS_1] [P=(1,1,-1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SS_0] [P=(1,1,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) T1up SS_0] [P=(1,1,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) T1up SS_0] [P=(1,1,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) A1p SS_1] [P=(1,1,0) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) A1p SS_1] [P=(1,1,0) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SD_0] [P=(1,1,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) A2m SS_0] [P=(1,1,0) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,0) A2m SS_0] [P=(2,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,1) A2m SS_1] [P=(2,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) A1um SS_0] [P=(1,1,1) A1m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) A2m SS_1] [P=(1,1,0) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(-1,0,1) A2m SS_0] [P=(2,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) A2m SS_0] [P=(1,1,0) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) A1gm SS_0] [P=(1,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) T1up SS_0] [P=(1,1,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) A1p SS_1] [P=(1,1,0) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B2p SS_1]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) A1um TDO_1] [P=(1,1,1) A1m SD_7]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) A2m SS_0] [P=(1,1,0) A1m TSD_1]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(-1,0,1) A2m SS_1] [P=(2,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) A2m SS_1] [P=(1,1,0) A1m LSD_6]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,1,0) A2m SS_1] [P=(1,0,1) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,0) A2m SS_1] [P=(2,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,1) A2m SS_1] [P=(1,1,0) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,1) A2m SS_0] [P=(2,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 CG_1 [P=(-1,0,1) A2m SS_0] [P=(2,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,2,0) A2m SS_1] [P=(1,-1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,1) A1m SS_0] [P=(1,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,0) A1um SS_0] [P=(1,1,1) Em SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,0) T1up SS_0] [P=(1,1,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,0) T1up SS_0] [P=(1,1,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 CG_1 [P=(0,0,0) T1up SS_0] [P=(1,1,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,1,0) A1p SS_1] [P=(1,0,1) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,1,0) A1p SS_1] [P=(1,0,1) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,1) A1p SS_1] [P=(1,1,0) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 CG_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) A1p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 CG_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 CG_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B2p SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,1,0) A2m SS_0] [P=(1,0,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,0) A2m SS_0] [P=(2,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,1) A2m SS_0] [P=(1,1,0) A1m TSD_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,1) A2m SS_1] [P=(2,1,0) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) A2m SS_1] [P=(0,0,1) A2m SS_1]");  
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SS_0] [P=(0,0,2) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,1) A2m SS_0] [P=(1,0,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,0) A2m SS_1] [P=(1,0,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,-1,1) A2m SS_0] [P=(1,1,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) A1p SS_1] [P=(0,0,1) A1p SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) Ep SS_1] [P=(0,0,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,-1,0) A2m SS_0] [P=(1,1,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,1) A2m SS_0] [P=(0,0,1) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(0,0,0) A1um SD_0] [P=(0,0,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,1) A2m SS_1] [P=(1,0,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion A1p_1 [P=(-1,0,0) A2m SS_0] [P=(1,0,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) A2m SS_1] [P=(0,0,1) A1m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) A1um SS_0] [P=(0,0,2) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) A2m SS_0] [P=(0,0,1) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,1) A2m SS_0] [P=(0,0,1) A1m LSD_3]");
 opstrings.push_back("isoquintet_pion_pion A2p_1 [P=(0,0,0) A1um TDO_1] [P=(0,0,2) A1m TSD_1]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,1) A2m SS_0] [P=(1,0,1) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,0) A2m SS_1] [P=(1,0,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,1) Ep SS_1] [P=(0,0,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,1) A2m SS_1] [P=(1,0,1) A2m SS_1]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(-1,0,0) A2m SS_0] [P=(1,0,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B1p_1 [P=(0,0,1) Ep SS_2] [P=(0,0,1) Ep SS_2]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,-1,1) A2m SS_0] [P=(1,1,1) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(0,0,1) Ep SS_1] [P=(0,0,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,-1,0) A2m SS_0] [P=(1,1,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion B2p_1 [P=(-1,0,1) A2m SS_0] [P=(1,0,1) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,0) A2m SS_1] [P=(1,0,2) A2m SS_0]");  
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,1) A1p SS_1] [P=(0,0,1) Ep SS_1]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,-1,0) A2m SS_0] [P=(1,1,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,-1,1) A2m SS_0] [P=(0,1,1) A1m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,0,0) A2m SS_0] [P=(1,0,2) A2m SS_0]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(0,0,1) A1p SS_2] [P=(0,0,1) Ep SS_2]");
 opstrings.push_back("isoquintet_pion_pion Ep_1 [P=(-1,-1,0) A2m SS_1] [P=(1,1,2) A2m SS_0]");
 opstrings.push_back("kaon P=(0,0,0) A1u_1 SS_0");
 opstrings.push_back("kaon P=(0,0,0) T1u_1 DDL_1");
 opstrings.push_back("kaon P=(0,0,0) T1u_1 SS_0");
 opstrings.push_back("kaon P=(0,0,0) T1u_1 SS_1");
 opstrings.push_back("kaon P=(0,0,0) T1u_1 DDL_5");
 opstrings.push_back("kaon P=(0,0,0) T1u_1 DDL_13");
 opstrings.push_back("kaon P=(0,0,0) T1u_1 TDO_42");
 opstrings.push_back("kaon P=(0,0,0) T1u_1 SD_3");
 opstrings.push_back("kaon P=(0,0,0) T1u_1 DDL_2");
 opstrings.push_back("kaon P=(0,0,0) T1u_1 DDL_6");
 opstrings.push_back("kaon P=(0,0,0) T1u_1 DDL_18");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(1,0,0) A2 SS_1] [P=(-1,0,0) A2m SS_1]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(1,0,0) A2 SS_1] [P=(-1,0,0) A2p SS_1]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(1,0,0) A2 SS_1] [P=(-1,0,0) A2p SS_1]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(1,0,1) A2 SS_0] [P=(-1,0,-1) A2m SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,0) T1g SS_0] [P=(0,0,0) A1um SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,1) E SS_2] [P=(0,0,-1) A2m SS_1]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,0) T1g DDL_6] [P=(0,0,0) A1um SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,1) A2 SS_1] [P=(0,0,-1) Ep SS_1]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(0,0,1) A2 SS_1] [P=(0,0,-1) Em SS_1]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(0,0,1) A2 SS_1] [P=(0,0,-1) Em SS_1]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(1,0,1) A2 SS_0] [P=(-1,0,-1) A2p SS_0]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(1,0,1) A2 SS_0] [P=(-1,0,-1) A2p SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(1,1,1) A2 SS_0] [P=(-1,-1,-1) A2m SS_0]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(0,0,0) A1u SS_0] [P=(0,0,0) T1gm SS_0]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(0,0,0) A1u SS_0] [P=(0,0,0) T1gm SS_0]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(0,0,1) E SS_2] [P=(0,0,-1) A2p SS_1]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(0,0,1) E SS_2] [P=(0,0,-1) A2p SS_1]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(1,0,1) B1 SS_1] [P=(-1,0,-1) A2m SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,1,1) B2 SS_3] [P=(0,-1,-1) A2m SS_0]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(1,0,0) A2 SS_1] [P=(-1,0,0) A2p SS_0]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(1,0,0) A2 SS_1] [P=(-1,0,0) A2p SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,0) A1u SS_0] [P=(0,0,0) T1gp SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,0) A1u SS_0] [P=(0,0,0) T1gm SS_0]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(0,0,1) A2 SS_1] [P=(0,0,-1) Em SS_2]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(0,0,1) A2 SS_1] [P=(0,0,-1) Em SS_2]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(1,0,1) A2 SS_0] [P=(-1,0,-1) B1p SS_1]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,1,1) A2 SS_0] [P=(0,-1,-1) B2p SS_2]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(1,0,1) A2 SS_0] [P=(-1,0,-1) B1m SS_1]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(0,1,1) A2 SS_0] [P=(0,-1,-1) B2m SS_2]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(1,0,1) A2 SS_0] [P=(-1,0,-1) B1m SS_1]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(0,1,1) A2 SS_0] [P=(0,-1,-1) B2m SS_2]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(0,0,0) A1u SS_0] [P=(0,0,0) T1gp SS_0]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(0,0,0) A1u SS_0] [P=(0,0,0) T1gp SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(1,0,0) A2 SS_0] [P=(-1,0,0) A2m SS_1]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,1) E SS_0] [P=(0,0,-1) A2m SS_1]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(0,0,0) T1g SS_0] [P=(0,0,0) A1up SS_0]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(0,0,0) T1g SS_0] [P=(0,0,0) A1up SS_0]");
 opstrings.push_back("isodoublet_kaon_eta T1u_1 [P=(1,1,1) A2 SS_0] [P=(-1,-1,-1) A2p SS_0]");
 opstrings.push_back("isodoublet_kaon_phi T1u_1 [P=(1,1,1) A2 SS_0] [P=(-1,-1,-1) A2p SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(2,0,0) A2 SS_1] [P=(-2,0,0) A2m SS_1]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,0) T1u SS_1] [P=(0,0,0) A1gm SS_0]"); 
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(1,0,0) A2 SS_1] [P=(-1,0,0) A2m SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(1,0,0) A2 SS_0] [P=(-1,0,0) A2m SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(1,0,1) A2 SS_1] [P=(-1,0,-1) A2m SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(1,0,1) A2 SS_0] [P=(-1,0,-1) A2m SS_1]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,1) E SS_3] [P=(0,0,-1) A2m SS_1]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,1) E SS_2] [P=(0,0,-1) A2m SS_0]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,1) A2 SS_0] [P=(0,0,-1) Ep SS_1]");
 opstrings.push_back("isodoublet_kaon_pion T1u_1 [P=(0,0,1) A2 SS_1] [P=(0,0,-1) Ep SS_2]");
 opstrings.push_back("pion P=(0,0,1) A2m_1 SS_0");   
 opstrings.push_back("pion P=(0,0,1) A2m_1 SS_1");   
 opstrings.push_back("pion P=(0,0,1) A2m_1 LSD_0");  
 opstrings.push_back("pion P=(0,0,1) A2m_1 LSD_1");  
 opstrings.push_back("pion P=(0,0,1) A2m_1 TSD_0");  
 opstrings.push_back("pion P=(0,0,1) A2m_1 TSD_1");  
 opstrings.push_back("pion P=(0,0,1) A2m_1 TSD_2");  
 opstrings.push_back("pion P=(0,0,1) A2m_1 DDL_0");  
 opstrings.push_back("pion P=(0,0,1) A2m_1 DDL_2");  
 opstrings.push_back("pion P=(0,0,1) A2m_1 DDL_8");  
 opstrings.push_back("pion P=(0,0,1) A2m_1 DDL_11"); 
 opstrings.push_back("pion P=(0,0,1) A2m_1 DDL_6");  
 opstrings.push_back("pion P=(0,1,0) A2m_1 SS_0");   
 opstrings.push_back("pion P=(0,1,0) A2m_1 SS_1");   
 opstrings.push_back("pion P=(0,1,0) A2m_1 LSD_0");  
 opstrings.push_back("pion P=(0,1,0) A2m_1 LSD_1");  
 opstrings.push_back("pion P=(0,1,0) A2m_1 TSD_0");  
 opstrings.push_back("pion P=(0,1,0) A2m_1 TSD_1");  
 opstrings.push_back("pion P=(0,1,0) A2m_1 TSD_2");  
 opstrings.push_back("pion P=(0,1,0) A2m_1 DDL_0");  
 opstrings.push_back("pion P=(0,1,0) A2m_1 DDL_2");  
 opstrings.push_back("pion P=(0,1,0) A2m_1 DDL_8");  
 opstrings.push_back("pion P=(0,1,0) A2m_1 DDL_11"); 
 opstrings.push_back("pion P=(0,1,0) A2m_1 DDL_6");  
 opstrings.push_back("pion P=(1,0,0) A2m_1 SS_0");   
 opstrings.push_back("pion P=(1,0,0) A2m_1 SS_1");   
 opstrings.push_back("pion P=(1,0,0) A2m_1 LSD_0");  
 opstrings.push_back("pion P=(1,0,0) A2m_1 LSD_1");  
 opstrings.push_back("pion P=(1,0,0) A2m_1 TSD_0");  
 opstrings.push_back("pion P=(1,0,0) A2m_1 TSD_1");  
 opstrings.push_back("pion P=(1,0,0) A2m_1 TSD_2");  
 opstrings.push_back("pion P=(1,0,0) A2m_1 DDL_0");  
 opstrings.push_back("pion P=(1,0,0) A2m_1 DDL_2");  
 opstrings.push_back("pion P=(1,0,0) A2m_1 DDL_8");  
 opstrings.push_back("pion P=(1,0,0) A2m_1 DDL_11"); 
 opstrings.push_back("pion P=(1,0,0) A2m_1 DDL_6");  
 opstrings.push_back("nucleon P=(0,0,0) G1g_1 TDT_29");
 opstrings.push_back("isodoublet_pion_nucleon G1g_1 [P=(0,0,1) A2m SS_1] [P=(0,0,-1) G1 SS_0]"); 

 int count=0;
 for (unsigned int k=0;k<opstrings.size();k++){
    OperatorInfo source(opstrings[k]);
    for (unsigned int j=0;j<opstrings.size();j++){
       OperatorInfo sink(opstrings[j]);
       cout <<endl<<endl<<" Test "<<count++<<endl<<endl;
       CorrelatorInfo corrtest(sink,source);
       string temp=corrtest.getSource().short_output();
       string temp2=corrtest.getSink().short_output();
       if (temp!=opstrings[k]) cout << "******MISMATCH ****"<<endl;
       if (temp2!=opstrings[j]) cout << "******MISMATCH2 ****"<<endl;}}

}





void testCorrelatorMatrixInfo(XMLHandler& xml_in)
{

 if (xml_tag_count(xml_in,"TestCorrelatorMatrixInfo")==0)
 return;

 set<OperatorInfo> qcdops;
 qcdops.insert(OperatorInfo("pion P=(0,1,0) A2m_1 DDL_8"));  
 qcdops.insert(OperatorInfo("kaon P=(0,0,0) T1u_1 DDL_18"));
 qcdops.insert(OperatorInfo("isoquintet_pion_pion Ep_1 CG_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B1p SS_1]"));

 
 {CorrelatorMatrixInfo cormat(qcdops,false,false);
 cout << cormat.output()<<endl;
 CorrelatorMatrixInfo cor2(cormat);
 cout << cor2.output()<<endl;
 qcdops.insert(OperatorInfo("isodoublet_kaon_pion T1u_1 [P=(1,0,0) A2 SS_0] [P=(-1,0,0) A2m SS_0]"));
 CorrelatorMatrixInfo cor3(qcdops,true,true);
 cormat=cor3;
 cout << cormat.output(false,0)<<endl;

 cout << "Iterating over the correlators: prefix"<<endl;
 for (cormat.begin();!cormat.end();++cormat)
    cout << cormat.getCurrentCorrelatorInfo().output()<<endl;

 cout << "Iterating over the correlators: postfix"<<endl;
 for (cormat.begin();!cormat.end();cormat++)
    cout << cormat.getCurrentCorrelatorInfo().output()<<endl;}

 cout <<endl<< "Testing input from XML:"<<endl<<endl;
 XMLHandler xmlop(xml_in,"TestInput1");
 CorrelatorMatrixInfo corrA(xmlop);
 cout << corrA.output()<<endl;
 cout << corrA.output(false)<<endl;

 XMLHandler xmlop2(xml_in,"TestInput2");
 CorrelatorMatrixInfo corrB(xmlop2);
 cout << corrB.output()<<endl;
 cout << corrB.output(false)<<endl;

 cout << "Equal?"<< (corrA==corrB)<<endl;
 cout << "Not Equal?"<< (corrA!=corrB) <<endl;
 cout << "A Less than B?"<< (corrA<corrB) <<endl;
 cout << "B Less than A?"<< (corrB<corrA) <<endl;

}

#include "xml_handler.h"
#include "basic_laph_operator_info.h"
#include "gen_irrep_operator_info.h"
#include "operator_info.h"

using namespace std;


void run_a_bl_op(const BasicLapHOperatorInfo& anop, const BasicLapHOperatorInfo& compareop)
{
 cout << endl<<endl<<"  ***************************** "<<endl<<endl;
 cout << "Long output:"<<endl;
 cout << anop.output(true)<<endl;
 cout << "Short output:"<<endl;
 cout << anop.short_output()<<endl;
 cout << "Short output again:"<<endl;
 cout << anop.output(false)<<endl<<endl;
 
 cout << "number of hadrons = "<<anop.getNumberOfHadrons()<<endl;
 cout << "is vacuum? "<< anop.isVacuum()<<endl;
 cout << "is glueball? "<< anop.isGlueball()<<endl;
 cout << "is meson? "<< anop.isMeson()<<endl;
 cout << "is baryon? "<< anop.isBaryon()<<endl;
 cout << "is meson-meson? "<< anop.isMesonMeson()<<endl;
 cout << "is meson-baryon? "<< anop.isMesonBaryon()<<endl;
// cout << "is rotated? "<< anop.isRotated()<<endl;
 cout << "Total X momentum = "<<anop.getXMomentum()<<endl;
 cout << "Total Y momentum = "<<anop.getYMomentum()<<endl;
 cout << "Total Z momentum = "<<anop.getZMomentum()<<endl;
 cout << "getLGIrrep: "<< anop.getLGIrrep()<<endl;
 cout << "getLGClebschGordonIdNum: "<<anop.getLGClebschGordonIdNum()<<endl;
 cout << "getLGIrrepRow: "<< anop.getLGIrrepRow()<<endl;
 cout << "getIsospin: "<< anop.getIsospin()<<endl; 
 cout << " getIsospinClebschGordonIdNum: "<< anop.getIsospinClebschGordonIdNum()<<endl; 
 cout << "getFlavor: "<< anop.getFlavor()<<endl; 
 cout << "getFlavorCode: "<< anop.getFlavorCode()<<endl;
   
 for (unsigned int h=0;h<4;h++){
    try {cout << "getFlavor("<<h<<"): "<<   anop.getFlavor(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "int getStrangeness("<<h<<"): "<< anop.getStrangeness(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "bool isGlueball("<<h<<"): "<< anop.isGlueball(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "bool isMeson("<<h<<"): "<< anop.isMeson(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "bool isBaryon("<<h<<"): "<< anop.isBaryon(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "bool isFermion("<<h<<"): "<< anop.isFermion(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "bool isBoson("<<h<<"): "<< anop.isBoson(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "std::string getLGIrrep("<<h<<"): "<< anop.getLGIrrep(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "std::string getSpatialType("<<h<<"): "<< anop.getSpatialType(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "unsigned int getSpatialIdNumber("<<h<<"): "<< anop.getSpatialIdNumber(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "unsigned int getDisplacementLength("<<h<<"): "<< anop.getDisplacementLength(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "Momentum getMomentum("<<h<<"): "<< anop.getMomentum(h).getMomentumString() <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "int getXMomentum("<<h<<"): "<< anop.getXMomentum(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "int getYMomentum("<<h<<"): "<< anop.getYMomentum(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}
    try {cout << "int getZMomentum("<<h<<"): "<< anop.getZMomentum(h) <<endl;} catch(const std::exception& xp){ cout << "exception "<<h<<" caught"<<endl;}}

 cout << "op1==op2? :"<<  (anop==compareop) <<endl;
 cout << "op1!=op2? :"<<  (anop!=compareop) <<endl; 
 cout << "op1<op2? :"<<  (anop<compareop) <<endl; 
 cout << endl<<endl;
}




void run_a_gi_op(const GenIrrepOperatorInfo& anop, const GenIrrepOperatorInfo& compareop)
{
 cout << endl<<endl<<"  ***************************** "<<endl<<endl;
 cout << "Long output:"<<endl;
 cout << anop.output(true)<<endl;
 cout << "Short output:"<<endl;
 cout << anop.short_output()<<endl;
 cout << "Short output again:"<<endl;
 cout << anop.output(false)<<endl<<endl;
 
 cout << "Total X momentum = "<<anop.getXMomentum()<<endl;
 cout << "Total Y momentum = "<<anop.getYMomentum()<<endl;
 cout << "Total Z momentum = "<<anop.getZMomentum()<<endl;
 cout << "getLGIrrep: "<< anop.getLGIrrep()<<endl;
 cout << "getLGIrrepRow: "<< anop.getLGIrrepRow()<<endl;
 cout << "getIsospin: "<< anop.getIsospin()<<endl; 
 cout << "getIDName: "<<anop.getIDName()<<endl;
 cout << "getIDIndex: "<<anop.getIDIndex()<<endl;

 cout << "op1==op2? :"<<  (anop==compareop) <<endl;
 cout << "op1!=op2? :"<<  (anop!=compareop) <<endl; 
 cout << "op1<op2? :"<<  (anop<compareop) <<endl; 
 cout << endl<<endl;
}




void testOperatorInfo(XMLHandler& xml_in)
{

 if (xml_tag_count(xml_in,"TestOperatorInfo")==0)
 return;

// cout << "codesize(0)="<<OperatorInfo::codesize(0)<<endl;
// cout << "codesize(1)="<<OperatorInfo::codesize(1)<<endl;
// cout << "codesize(2)="<<OperatorInfo::codesize(2)<<endl;
// cout << "codesize(3)="<<OperatorInfo::codesize(3)<<endl;

 BasicLapHOperatorInfo op1;
 run_a_bl_op(op1,op1);
 BasicLapHOperatorInfo op2("glueball P=(0,0,0) A1gp_1 TrEig");
 run_a_bl_op(op2,op2);
 run_a_bl_op(op2,op1);
 BasicLapHOperatorInfo op3("isosinglet_eta_eta A1gp_1 [P=(0,0,0) A1up SS_0] [P=(0,0,0) A1up SS_0]");    
 run_a_bl_op(op3,op3);
 BasicLapHOperatorInfo op4("isosinglet_kaon_kbar A1gp_1 [P=(0,0,0) A1g SD_0] [P=(0,0,0) A1g TDO_3]");
 run_a_bl_op(op4,op3);
 BasicLapHOperatorInfo op5("pion P=(0,0,0) T1up_1 DDL_8");
 run_a_bl_op(op5,op3);
 BasicLapHOperatorInfo op6("isoquintet_pion_pion Ep_1 CG_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B1p SS_1]");
 run_a_bl_op(op6,op3);
 BasicLapHOperatorInfo op7("isoquintet_pion_pion Ep_1 CG_1 [P=(0,0,1) Ep SS_1] [P=(1,1,0) B2p SS_2]");
 run_a_bl_op(op7,op7);
 BasicLapHOperatorInfo op8("kaon P=(0,0,0) T1u_1 SS_0");
 run_a_bl_op(op8,op6);
 BasicLapHOperatorInfo op9("isodoublet_kaon_pion T1u_1 [P=(2,0,0) A2 SS_1] [P=(-2,0,0) A2m SS_1]");
 run_a_bl_op(op9,op9);
 BasicLapHOperatorInfo op10("pion P=(0,0,1) A2m_1 DDL_6");  
 run_a_bl_op(op10,op7);

 BasicLapHOperatorInfo op11("nucleon P=(0,0,0) G1g_1 TDT_29");
 run_a_bl_op(op11,op11);
 BasicLapHOperatorInfo op12("isodoublet_pion_nucleon G1g_1 [P=(0,0,1) A2m SS_1] [P=(0,0,-1) G1 SS_0]"); 
 run_a_bl_op(op12,op12);
 run_a_bl_op(op12,op11);
 run_a_bl_op(op12,op10);


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

 opstrings.push_back("tquuuu1p P=(0,1,0) A2m_1 QDX_0");
 opstrings.push_back("tquudu3p P=(0,1,0) A1m_1 QDX_1");
 opstrings.push_back("tqdudu1p P=(0,1,0) A2m_1 QDX_2");
 opstrings.push_back("tqdudu3p P=(0,1,0) A1m_1 QDX_3");
 opstrings.push_back("tqdudu5p P=(0,1,0) A2m_1 QDX_4");
 opstrings.push_back("tqsuuu2p P=(0,1,0) A1m_1 QDX_5");
 opstrings.push_back("tqsudu2p P=(0,1,0) A2m_1 QDX_6");
 opstrings.push_back("tqsudu4p P=(0,1,0) A1m_1 QDX_7");
 opstrings.push_back("tqssdu3p P=(0,1,0) A2m_1 QDX_8");
 opstrings.push_back("tquuss1p P=(0,1,0) A1m_1 QDX_9");
 opstrings.push_back("tqsuss2p P=(0,1,0) A2m_1 QDX_10");
 opstrings.push_back("tqssss1p P=(0,1,0) A1m_1 QDX_11");

 opstrings.push_back("tquuuu1m P=(0,1,0) A2m_1 QDX_0");
 opstrings.push_back("tquudu3m P=(0,1,0) A1m_1 QDX_1");
 opstrings.push_back("tqdudu1m P=(0,1,0) A2m_1 QDX_2");
 opstrings.push_back("tqdudu3m P=(0,1,0) A1m_1 QDX_3");
 opstrings.push_back("tqdudu5m P=(0,1,0) A2m_1 QDX_4");
 opstrings.push_back("tqsuuu2m P=(0,1,0) A1m_1 QDX_5");
 opstrings.push_back("tqsudu2m P=(0,1,0) A2m_1 QDX_6");
 opstrings.push_back("tqsudu4m P=(0,1,0) A1m_1 QDX_7");
 opstrings.push_back("tqssdu3m P=(0,1,0) A2m_1 QDX_8");
 opstrings.push_back("tquuss1m P=(0,1,0) A1m_1 QDX_9");
 opstrings.push_back("tqsuss2m P=(0,1,0) A2m_1 QDX_10");
 opstrings.push_back("tqssss1m P=(0,1,0) A1m_1 QDX_11");


 for (unsigned int k=0;k<opstrings.size();k++){
    cout <<endl<<endl<<" Test "<<k<<endl<<endl;
    cout << "opstring = <"<<opstrings[k]<<">"<<endl<<endl;
    BasicLapHOperatorInfo opinfo(opstrings[k]);
    XMLHandler xmlt("OperatorString",opstrings[k]);
    BasicLapHOperatorInfo opinfo2(xmlt);
    cout << opinfo.output(true)<<endl;
    string temp=opinfo.short_output();
    if (temp!=opstrings[k]) cout << "******MISMATCH ****"<<endl;
    string temp2=opinfo2.short_output();
    if (temp2!=opstrings[k]) cout << "******MISMATCH2 ****"<<endl;
    cout << opinfo.output(false)<<endl;
    cout << temp<<endl;}

 cout << endl<<endl<<"Reading tests"<<endl<<endl;
 list<XMLHandler> tests=xml_in.find("Test");

 cout <<"Number of tests to do: "<<tests.size()<<endl;

 for (list<XMLHandler>::iterator it=tests.begin();it!=tests.end();it++){
    BasicLapHOperatorInfo opf(*it);
    cout << opf.output(true)<<endl;}


 XMLHandler xmlrot("GIOperator");
 xmlrot.put_child("Isospin","triplet");
 xmlrot.put_child("LGIrrep","T2gm");
 xmlrot.put_child("LGIrrepRow","2");
 xmlrot.put_child("Momentum","1 0 -1");
 xmlrot.put_child("IDName","RotatedName");
 xmlrot.put_child("IDIndex",make_string(8));
 xmlrot.put_child("Description","a very fragrant operator which smells nice");

// cout << xmlrot.output()<<endl;
 cout << "Test4"<<endl;
 GenIrrepOperatorInfo test4(xmlrot);
 cout << test4.output()<<endl;
 cout << test4.output(true)<<endl;


 xmlrot.set_root("GIOperator");
 xmlrot.put_child("Isospin","singlet");
 xmlrot.put_child("LGIrrep","Eu");
 xmlrot.put_child("LGIrrepRow","1");
 xmlrot.put_child("Momentum","0 0 -2");
 xmlrot.put_child("IDName","GhastlyFart");
 cout << "Test5"<<endl;
 GenIrrepOperatorInfo test5(xmlrot);
 cout << test5.output()<<endl;
 cout << test5.output(true)<<endl;

 cout << "equal? "<<(test4==test5)<<endl;
 cout << "not equal? "<<(test4!=test5)<<endl;

 xmlrot.set_root("GIOperator");
 xmlrot.put_child("IDName","BubbleGum");
 xmlrot.put_child("IDIndex",make_string(3));
 xmlrot.put_child("Description","a very nice operator which smells nice");
 xmlrot.put_child("Isospin","quintet");
 xmlrot.put_child("LGIrrep","G1g");
 xmlrot.put_child("LGIrrepRow","1");
 xmlrot.put_child("Momentum","1 1 1");
 cout << "Test6"<<endl;
 GenIrrepOperatorInfo test6(xmlrot);
 cout << test6.output()<<endl;
 cout << test6.output(true)<<endl;

 cout << "equal? "<<(test4==test6)<<endl;
 cout << "not equal? "<<(test4!=test6)<<endl;

 xmlrot.set_root("GIOperator");
 xmlrot.put_child("Isospin","singlet");
 xmlrot.put_child("LGIrrep","T1u");
 xmlrot.put_child("LGIrrepRow","1");
 xmlrot.put_child("Momentum","0 0 -2");
 xmlrot.put_child("IDName","GhastlyFart");
 cout << "Test5B"<<endl;
 GenIrrepOperatorInfo test5B(xmlrot);
 cout << test5B.output()<<endl;
 cout << test5B.output(true)<<endl;

 cout << "equal? "<<(test5==test5B)<<endl;
 cout << "not equal? "<<(test5!=test5B)<<endl;



 xmlrot.set_root("GIOperatorString","isosinglet P=(1,2,0) A1gm_2 Doggy ");
 cout << "Test7"<<endl;
 GenIrrepOperatorInfo test7(xmlrot);
 cout << test7.output()<<endl;
 cout << test7.short_output()<<endl;
 cout << op10.short_output()<<endl;
 cout << test7.output(true)<<endl;

 cout << "equal? "<<(test4==test7)<<endl;
 cout << "not equal? "<<(test4!=test7)<<endl;

 test7.resetIDIndex(3);
 cout << "reset index to 3 below"<<endl;
 cout << test7.output()<<endl;

 GenIrrepOperatorInfo gop11("isodoublet P=(0,0,0) G1g_1 phi=(0,0,0)-rho=(1,1,1) 7");
 run_a_gi_op(gop11,gop11);
 cout << gop11.output(true)<<endl;

 cout << endl<<endl<<" NOW to test OperatorInfo"<<endl<<endl;

 OperatorInfo tt1(gop11);
 cout << tt1.output()<<endl;
 cout << "is BasicLapH? "<<tt1.isBasicLapH()<<endl;
 cout << "is GenIrrep? "<<tt1.isGenIrrep()<<endl;
 tt1.resetGenIrrepIDIndex(4); cout << "reset index to 4 below"<<endl;
 cout << tt1.output()<<endl;


}

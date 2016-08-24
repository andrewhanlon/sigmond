#include "vev_data_handler.h"
#include "mcobs_handler.h"
#include <cstdio>
#include <ctime>
#include <map>

using namespace std;
using namespace LaphEnv;


void testVEVDataHandler(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestVEVDataHandler")==0)
 return;

 cout << endl << "Starting test_vev_handler"<<endl;
 /* 
 try{
 XMLHandler xmlr(xml_in,"TestVEVDataHandler");
 BLVEVDataHandler VH(xmlr); 
 cout << "Number of measurements in ensemble = "<<VH.getNumberOfMeasurements()<<endl;
 cout << "Ensemble ID = "<<VH.getEnsembleId()<<endl;

 set<OperatorInfo> opSet=VH.getOperatorSet();
 dcmplx result;
 for (set<OperatorInfo>::const_iterator 
      ct=opSet.begin();ct!=opSet.end();ct++){
    cout << "OperatorInfo: "<<endl;
    cout << ct->output()<<endl;
    for (unsigned int serial=0;serial<VH.getNumberOfMeasurements();serial++){
       try{VH.getData(*ct,serial,result);
       cout << " config="<<serial<<"  result: "<<result<<endl;}
       catch(const std::exception& xp){}}}

 MCDataHandler<BLVEVDataHandler,OperatorInfo,dcmplx> mch(VH);
 cout << endl<<endl<<" MC Data Handler NOW"<<endl<<endl;
 cout << "number of MC measurements = "<<mch.getNumberOfMeasurements()<<endl;

 for (set<OperatorInfo>::const_iterator 
      ct=opSet.begin();ct!=opSet.end();ct++){
    cout << "OperatorInfo: "<<endl;
    cout << ct->output()<<endl;
    result=mch.getMean(*ct);
    cout << " mean: "<<result<<endl;}

 }
 catch(const string& err){
    cerr << "  Error: "<<err<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}
*/
}




#include "corr_data_handler.h"
#include "mcobs_handler.h"
#include <cstdio>
#include <ctime>
#include <map>

using namespace std;
using namespace LaphEnv;


void testCorrDataHandler(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestCorrDataHandler")==0)
 return;

 cout << endl << "Starting test_corr_handler"<<endl;
/*  
 try{
 XMLHandler xmlr(xml_in,"TestCorrDataHandler");
 BLCorrelatorDataHandler CH(xmlr); 
 cout << "Number of measurements in ensemble = "<<CH.getNumberOfMeasurements()<<endl;
 cout << "Ensemble ID = "<<CH.getEnsembleId()<<endl;

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

 }
 catch(const std::invalid_argument& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}
*/
}




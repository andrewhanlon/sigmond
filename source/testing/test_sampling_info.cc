#include "sampling_info.h"
#include <cstdio>
#include <ctime>
#include <map>
#include <cstdlib>
#include <iostream>
#include <list>

using namespace std;


  // *************************************************************************


void testMCSamplingInfo(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestMCSamplingInfo")==0)
 return;

 cout << endl << "Starting TestMCSamplingInfo"<<endl;
 
 try{
 XMLHandler xmlb(xml_in,"TestMCSamplingInfo");  cout << xmlb.output()<<endl;
 MCBinsInfo binfo(xmlb);
 cout<< "Bins Info : "<<binfo.output()<<endl;

 list<XMLHandler> xmlr(xml_in.find("MCSamplingInfo"));
 for (list<XMLHandler>::iterator it=xmlr.begin();it!=xmlr.end();it++){
    MCSamplingInfo Ms(*it);
    cout << Ms.output()<<endl;
    cout << " is jackknife? "<<Ms.isJackknifeMode()<<endl;
    cout << " is bootstrap? "<<Ms.isBootstrapMode()<<endl;
    SamplingMode smode=Ms.getSamplingMode();
    cout << " sampling mode = "<<smode<<endl;
    cout << " number of samplings = "<<Ms.getNumberOfReSamplings(binfo)<<endl;
    cout << " seed = "<<Ms.getRNGSeed()<<endl;
    cout << " skip = "<<Ms.getSkipValue()<<endl;
    }

 for (list<XMLHandler>::iterator it1=xmlr.begin();it1!=xmlr.end();it1++){
    MCSamplingInfo Ms1(*it1);
    for (list<XMLHandler>::iterator it2=xmlr.begin();it2!=xmlr.end();it2++){
       MCSamplingInfo Ms2(*it2);
       cout << Ms1.str()<<" -- "<<Ms2.str() << endl;
       cout << "equal? "<< (Ms1==Ms2) <<endl;
       cout << "unequal? "<< (Ms1!=Ms2) <<endl<<" ________"<<endl;}}


 MCSamplingInfo s1;
 cout <<"s1: "<< s1.output()<<endl;

 MCSamplingInfo s2(s1);
 cout <<"s2: "<< s2.output()<<endl;

 Bootstrapper B(551,888,563,22,false);
 s1.setToBootstrapMode(B);
 cout << "s1 boot: "<<s1.output()<<endl;

 MCSamplingInfo s3(s1);
 cout <<"s3: "<< s3.output()<<endl;

 s2=s1;
 cout <<"s2: "<< s2.output()<<endl;

 s1.setToJackknifeMode();
 cout << "s1 "<<s1.output()<<endl;

 MCSamplingInfo s8(234,111,74);
 cout << "s8: "<<s8.output()<<endl;

 }
 catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}

}


// ***********************************************************************

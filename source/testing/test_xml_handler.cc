#include "xml_handler.h"
#include <cstdio>
#include <ctime>
#include <map>

using namespace std;



void testXMLHandler(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestXMLHandler")==0)
 return;

 cout << endl << "Starting test_xml_handler"<<endl;

 XMLHandler xmlt(xml_in,"TestXMLHandler");

 string sval;
 try{
    xmlread(xmlt,"TagA",sval,"testTaskHandler");
    cout << "TEST 1:  sval = "<<sval<<endl;}
 catch(const std::exception& xp){
    cout << "TEST 1:  Exception caught"<<endl;}
 try{
    xmlreadchild(xmlt,"TagA",sval,"testTaskHandler");
    cout << "TEST 2:  sval = "<<sval<<endl;}
 catch(const std::exception& xp){
    cout << "TEST 2:  Exception caught "<<xp.what()<<endl;}
 try{
    xmlreadchild(xmlt,"TagA",sval);
    cout << "TEST 3:  sval = "<<sval<<endl;}
 catch(const std::exception& xp){
    cout << "TEST 3:  Exception caught "<<xp.what()<<endl;}

 cout << "Current node = "<<xmlt.get_node_name()<<endl;
 cout << "TEST 4: count = "<<xml_tag_count(xmlt,"TagA")<<endl;
 cout << "TEST 5: count = "<<xml_child_tag_count(xmlt,"TagA")<<endl;
 cout << "TEST 5B: count = "<<xmlt.count_to_among_children("TagA")<<endl;

 xmlt.seek_first_child();
 cout << "Current node = "<<xmlt.get_node_name()<<endl;
 cout << "TEST 5C: count = "<<xmlt.count_to_among_children("TagA")<<endl;
 xmlt.seek_first_child();
 cout << "Current node = "<<xmlt.get_node_name()<<endl;
 cout << "TEST 5D: count = "<<xmlt.count_to_among_children("TagA")<<endl;
 xmlt.seek_first_child();
 cout << "Current node = "<<xmlt.get_node_name()<<endl;
 cout << "TEST 5E: count = "<<xmlt.count_to_among_children("TagA")<<endl;

 xmlt.seek_root();
 try{
    XMLHandler xmltt(xmlt);
    xml_tag_assert(xmltt,"Donkey");
    cout << "TEST 6:  Donkey assertion ok"<<endl;}
 catch(const std::exception& xp){
    cout << "TEST 6:  Exception caught "<<xp.what()<<endl;}

 try{
    XMLHandler xmltt(xmlt);
    xml_tag_assert(xmltt,"Spaceship","TestTaskHandler");
    cout << "TEST 7:  Spaceship assertion ok"<<endl;}
 catch(const std::exception& xp){
    cout << "TEST 7:  Exception caught "<<xp.what()<<endl;}

 try{
    XMLHandler xmltt(xmlt);
    xml_child_assert(xmltt,"Worm");
    cout << "TEST 8:  Worm child assertion ok"<<endl;}
 catch(const std::exception& xp){
    cout << "TEST 8:  Exception caught "<<xp.what()<<endl;}

 try{
    XMLHandler xmltt(xmlt);
    xml_child_assert(xmltt,"Bedbug","TestTaskHandler");
    cout << "TEST 9:  Bedbug child assertion ok"<<endl;}
 catch(const std::exception& xp){
    cout << "TEST 9:  Exception caught "<<xp.what()<<endl;}

 try{
    XMLHandler xmltt(xmlt,"Spaceship");
    cout << xmltt.output()<<endl;
    xml_root_assert(xmltt,"Spaceship","TestTaskHandler");
    cout << "TEST 10:  Spaceship root assertion ok"<<endl;}
 catch(const std::exception& xp){
    cout << "TEST 10:  Exception caught "<<xp.what()<<endl;}


}

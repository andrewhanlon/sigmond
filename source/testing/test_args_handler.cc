#include "args_handler.h"
#include "log_helper.h"
#include "operator_info.h"
#include "filelist_info.h"
using namespace std;



void testArgsHandler(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestArgsHandler")==0)
 return;

 XMLHandler xml1("Root1");
 xml1.put_child("cartan","     ");
 xml1.put_child("first_child",make_string(4));
 xml1.put_child("third_child","ernie and bert");
 xml1.put_child("second_child",make_string(7.321));
 xml1.put_child("kangaroo",make_string(9));
 xml1.put_child("dingo");
 xml1.put_child("flipper","true");
 xml1.put_child("flipper2","false");
 OperatorInfo opinfo("isotriplet P=(0,0,0) A1um_1 IDname 2", OperatorInfo::GenIrrep);
 XMLHandler xml2; opinfo.output(xml2);
 xml1.put_child(xml2);
 XMLHandler xml3("Root2");
 xml3.put_child("tagger1","frog");
 xml3.put_child("winer",make_string(-3.421));
 xml3.put_child("third_child","party");
 xml3.put_child("dingo");
 XMLHandler xml4("Root3");
 xml4.put_child("dingo");
 xml3.put_child(xml4);
 xml1.put_child(xml3);
 xml1.put_child("late_child","      ");
 xml1.put_child("later_child"," farm 788     ");
 xml1.put_child("latest_child"," %%%$##    ");
 xml1.put_child("fat_child"," Gretchen564_32at.12  ");
 xml1.put_child("intvec"," 2 -3 5");


 cout << "INPUT XML:"<<endl<<xml1.output()<<endl<<endl;
 
 ArgsHandler ah1(xml1);
 XMLHandler xmltemp;
 ah1.getInput(xmltemp);
 cout << "reshow input:"<<xmltemp.output()<<endl;
 cout << "reshow input again: "<<ah1.getInput()<<endl<<endl;

 try{
    ArgsHandler ah2(xml1,"Root3",true);
    cout << "input: "<<ah2.getInput()<<endl<<endl;}
 catch(const std::exception& errmsg){
    cout << "EXCEPTION CAUGHT CORRECTLY "<<errmsg.what()<<endl;}

 ArgsHandler ah3(xml1,"Root3",false);
 cout << "input: "<<ah3.getInput()<<endl<<endl;

 int k;
 double g;
 string str;
 int m=3;
 try{
   // xml1.seek_unique("third_child");
    ArgsHandler gin(xml1);
    gin.setOffEcho();

    cout <<endl<<endl<< "Testing reading an int"<<endl<<endl;
    gin.getInt("first_child",k);  cout << "k = "<<k<<" should be 4"<<endl;
    k=99;
    k=gin.getInt("first_child"); cout << "k = "<<k<<" should be 4"<<endl;
    k=88;
    gin.getOptionalInt("first_child",k); cout << "k = "<<k<<" should be 4"<<endl;
    k=88;
    gin.getOptionalInt("firstb_child",k); cout << "k = "<<k<<" should be 88"<<endl;
    try{
       gin.getInt("firstb_child",k);}
    catch(const std::exception& err){ cout << "CAUGHT CORRECTLY:  k = "<<k<<endl<<err.what()<<endl;}
    try{
       gin.getInt("second_child",k);}
    catch(const std::exception& err){ cout << "CAUGHT CORRECTLY:  k = "<<k<<endl<<err.what()<<endl;}
    try{
       gin.getInt("third_child",k);}
    catch(const std::exception& err){ cout << "CAUGHT CORRECTLY:  k = "<<k<<endl<<err.what()<<endl;}

    gin.getReal("second_child",g); cout << "g = "<<g<<" should be 7.321"<<endl;
    g=99.9;
    g=gin.getReal("second_child"); cout << "g = "<<g<<" should be 7.321"<<endl;
    g=88.8;
    gin.getOptionalReal("second_child",g); cout << "g = "<<g<<" should be 7.321"<<endl;
    g=88.8;
    gin.getOptionalReal("firstb_child",g); cout << "g = "<<g<<" should be 88.8"<<endl;
    try{
       gin.getReal("firstb_child",g);}
    catch(const std::exception& err){ cout << "CAUGHT CORRECTLY:  g = "<<g<<endl<<err.what()<<endl;}
    try{
       gin.getReal("second_child",g);}
    catch(const std::exception& err){ cout << "CAUGHT CORRECTLY:  g = "<<g<<endl<<err.what()<<endl;}
    try{
       gin.getReal("third_child",g);}
    catch(const std::exception& err){ cout << "CAUGHT CORRECTLY:  g = "<<g<<endl<<err.what()<<endl;}
    


    gin.setOnEcho();
    cout <<endl<<endl<<"Testing strings and names"<<endl<<endl;
    gin.getString("third_child",str);  cout << "str = <"<<str<<"> should be <ernie and bert>"<<endl;
    try{
    gin.getString("late_child",str);  cout << "str = <"<<str<<">"<<endl;}
    catch(const std::exception& errmsg){
      cout << "caught exception CORRECT since white space string: str = "<<str<<endl;}
    gin.getString("later_child",str);  cout << "str = <"<<str<<"> should be <farm 788>"<<endl;
    gin.getString("latest_child",str);  cout << "str = <"<<str<<"> should be <%%%$##>"<<endl;
    string str2(gin.getString("later_child")); cout << "str2 = <"<<str2<<"> should be <farm 788>"<<endl;
    try{
       string str3(gin.getString("late_child")); cout << "str3 = "<<str3<<endl;}
    catch(const std::exception& errmsg){
      cout << "caught exception CORRECT since white space string: for str3 "<<endl;}
 
    try{
    str=gin.getName("later_child");  cout << "name = <"<<str<<"> should be <farm 788>"<<endl;}
    catch(const std::exception& errmsg){
      cout << "caught exception CORRECT since invalid name: str = "<<str<<endl;}
    try{
    str=gin.getName("third_child");  cout << "name = <"<<str<<"> should be <farm 788>"<<endl;}
    catch(const std::exception& errmsg){
      cout << "caught exception CORRECT since invalid name: str = "<<str<<endl;}
    str=gin.getName("fat_child");  cout << "str = <"<<str<<"> should be <Gretchen564_32at.12>"<<endl;

    gin.getOptionalInt("kangaroo",m);
    double x=gin.getReal("second_child"); cout << "x = "<<x<<endl;
    int q=gin.getInt("first_child"); cout << "q = "<<q<<endl;
   // gin.getItem("dingo",dingo);  cout << "dingo = "<<dingo<<endl;

    OperatorInfo o1;
    gin.getItem("OperatorInfo",o1);
    cout << "read o1: "<<o1.output()<<endl;

//    OperatorInfo oop(gin.getItem<OperatorInfo>("OperatorInfo"));
//    cout << "oop: "<<oop.output()<<endl;

    try{
       string gg; gin.getString("cartan",gg);
       cout << "gg = <"<<gg<<">"<<endl;}
    catch(const std::exception& errmsg){
       cout << "CAUGHT EXCEPTION: "<<errmsg.what()<<endl;}

    ArgsHandler gin2(gin,"Root2");
    double w=gin2.getReal("winer"); cout <<" w = "<<w<<endl;
    gin.insert(gin2);


    cout << "test bools: "<<endl;

    bool vev=gin.getBool("VEV"); cout << "vev = "<<vev<<" should be 0"<<endl;
    bool dingo; gin.getBool("dingo",dingo); cout << "dingo = "<<dingo<<" should be 1"<<endl;
    bool h;
    h=gin.getBool("flipper"); cout << "flipper = "<<h<<" should be 1"<<endl;
    h=gin.getBool("flipper2"); cout << "flipper2 = "<<h<<" should be 0"<<endl;
    try{
    h=gin.getBool("first_child"); cout << " h = "<<h<<" should be 0"<<endl;  }
    catch(const std::exception& errmsg){
      cout << "caught exception CORRECT for h"<<endl;}

    bool bt=true;
    gin.getOptionalBool("Farter",bt); cout << "bt = "<<bt<<" should be 1"<<endl;
    gin.getOptionalBool("flipper",bt); cout << "bt = "<<bt<<" should be 1"<<endl;
    gin.getOptionalBool("flipper2",bt); cout << "bt = "<<bt<<" should be 0"<<endl;

    vector<int> p(gin.getIntVector("intvec"));
    for (uint ip=0;ip<p.size();ip++) cout << "p["<<ip<<"] = "<<p[ip]<<endl;
    vector<int> pp; gin.getIntVector("intvec",pp);
    for (uint ip=0;ip<pp.size();ip++) cout << "pp["<<ip<<"] = "<<pp[ip]<<endl;

    XMLHandler xmlout;
    gin.echo(xmlout);
    cout << xmlout.output()<<endl;
    cout << "Input root tag: "<<gin.getInputRootTag()<<endl;}
 catch(const std::exception& errmsg){
    cout << "Error: "<<errmsg.what()<<endl;}

 cout << endl<<endl<<"*********************************************************"<<endl<<endl;

 set<string> roottags;
 try{
    ArgsHandler gtest1(xml1,roottags);
    cout <<"ERROR"<<endl;}
 catch(const std::exception& errmsg){
    cout << "Caught Exception CORRECT "<<errmsg.what()<<endl;}

 roottags.insert("farmer");
 roottags.insert("gino");
 roottags.insert("first_pig");
 try{
    ArgsHandler gtest1(xml1,roottags);
    cout <<"ERROR"<<endl;}
 catch(const std::exception& errmsg){
    cout << "Caught Exception CORRECT "<<errmsg.what()<<endl;}

 roottags.insert("first_child");
 roottags.insert("sphaghetti");
 ArgsHandler gtest2(xml1,roottags);
 cout << "gtest2 input root = "<<gtest2.getInputRootTag()<<endl;

 roottags.insert("second_child");
 try{
    ArgsHandler gtest1(xml1,roottags);
    cout <<"ERROR"<<endl;}
 catch(const std::exception& errmsg){
    cout << "Caught Exception CORRECT "<<errmsg.what()<<endl;}

 cout << endl<<endl<<"*********************************************************"<<endl<<endl;

 LogHelper LL("helproot");
 LL.putInt("Int",-3);
 LL.putUInt("UInt",5);
 LL.putReal("Real",3.21534);
 LL.putString("String"," an unabashed string");
 LogHelper LL2("Root2");
 LL2.putInt("AnInt",-8);
 LL2.putString("AString","merry christmas");
 LL.put(LL2);
 LL.putBool("Bool",false);
 LL.putBool("Bool2",true);
 LL.putBoolAsEmpty("EmptyBool",true);
 LL.putBoolAsEmpty("EmptyBool2",false);

 LL.putItem(opinfo);
 LL.putItem("Source",opinfo);

 vector<int> p(3); p[0]=0; p[1]=-4; p[2]=9;
 LL.putIntVector("PVector",p);


 cout << LL.output()<<endl;

 cout << "additional tests:"<<endl;

 XMLHandler xmldum("Aroot");
 xmldum.put_child("Child1","fairy");
 xmldum.put_child("Child2","gnome");
 xmldum.seek_first_child();
 xmldum.seek_next_sibling();
 xmldum.put_child("Grandchild","orcs");

 LogHelper xmladd(xmldum);
 cout << xmladd.output()<<endl;

// cout << "a last test"<<endl;
// cout << xmldum.get_node_name()<<endl;
// xmldum.seek_next_node();
// cout << xmldum.output()<<endl;
// cout << xmldum.get_node_name()<<endl;
// XMLHandler xmlttt(xmldum);
// cout << xmlttt.output()<<endl;

 cout <<endl<<endl<<" *************Testing getMultiItems() for FileListInfo"<<endl<<endl;

 FileListInfo finfo("StubA",0,24);
 XMLHandler xmltestmulti("Root");
 XMLHandler xmlf; finfo.output(xmlf);
 xmltestmulti.put_child(xmlf);

 FileListInfo finfo2("StubB",6,16);
 finfo2.output(xmlf);
 xmltestmulti.put_child(xmlf);

 FileListInfo finfo3("StubC",2,22);
 finfo3.output(xmlf);
 xmltestmulti.put_child(xmlf);

 xmltestmulti.put_child("Stringer","ham");
 xmltestmulti.put_child("Stringer","beef");
 xmltestmulti.put_child("Stringer");

 ArgsHandler xmlaa(xmltestmulti);
 list<FileListInfo> alist;
 xmlaa.getMultiItems("FileListInfo",alist);
 for (list<FileListInfo>::iterator it=alist.begin();it!=alist.end();it++)
    cout << it->output()<<endl;

 list<string> slist;
 xmlaa.getMultiStrings("Stringer",slist);
 for (list<string>::iterator it=slist.begin();it!=slist.end();it++)
    cout << *it<<endl;

 cout << "ECHO"<< xmlaa.echo()<<endl;

 cout << endl<<endl<<"*********************************************************"<<endl<<endl;
}


// ******************************************************************************

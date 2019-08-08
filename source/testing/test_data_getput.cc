#include "xml_handler.h"
#include "args_handler.h"
#include "data_io_handler.h"
#include "bins_handler.h"
#include "samplings_handler.h"
#include <ctime>
#include <map>

using namespace std;


class HeadTemp
{
  int header_int;
 public:
   HeadTemp() : header_int(0) {}
   HeadTemp(int v) : header_int(v) {}
   bool checkHeader(XMLHandler& xmlr);
   void writeHeader(XMLHandler&);
};

void HeadTemp::writeHeader(XMLHandler& headerxml)
{
 headerxml.set_root("HeadTemp");
 headerxml.put_child("KeyType","Integer");
 headerxml.put_child("ValueType","Double");
 headerxml.put_child("Code",make_string(header_int));
}

bool HeadTemp::checkHeader(XMLHandler& xmlr)
{
 try{
    ArgsHandler xarg(xmlr,"HeadTemp",true);
    string check=xarg.getString("KeyType");
    if (check!="Integer") return false;
    check=xarg.getString("ValueType");
    if (check!="Double") return false;
    int icheck=xarg.getInt("Code");
    if (icheck!=header_int) return false;
    return true;}
 catch(const exception& xp){
    return false;}
}



void testDataGetPut(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestDataGetPut")==0)
 return;

 cout << endl<<endl<<"***************************************************"<<endl<<endl;
 cout << "Testing DataGetPut"<<endl;
 
 string filetype_id("DataGetPutTester");
 XMLHandler xmlr(xml_in,"TestDataGetPut");
 ArgsHandler xxarg(xmlr);
 bool overwrite=xxarg.getBool("Overwrite"); 
 WriteMode wmode=(overwrite) ? Overwrite : Update;
 bool use_checksums=xxarg.getBool("UseChecksums");
 int header_int=0;
 xxarg.getOptionalInt("HeaderInteger",header_int);
 cout << "overwrite = "<<overwrite<<endl;
 cout << "usechecksums = "<<use_checksums<<endl;
 cout << "header int = "<<header_int<<endl;
 HeadTemp theader(header_int);

 set<string> fputnames;
 list<XMLHandler> dhputs(xmlr.find("DoSomePuts"));
 for (list<XMLHandler>::iterator it=dhputs.begin();it!=dhputs.end();it++){
    cout <<endl<<endl<<"  ***** DO SOME PUTS *****"<<endl<<endl<< it->output()<<endl;
    ArgsHandler xarg(*it);
    string filename(xarg.getString("FileName"));
    fputnames.insert(filename);
    DataPutHandlerSF<HeadTemp,UIntKey,double> 
        DHput(theader,filename,filetype_id,wmode,use_checksums);
    list<XMLHandler> keyxmls(it->find("DoAPut"));
    for (list<XMLHandler>::iterator kt=keyxmls.begin();kt!=keyxmls.end();kt++){
       ArgsHandler karg(*kt);
       uint key; double value;
       key=karg.getInt("Key");
       value=karg.getReal("Value");
       cout << "Putting to filename "<<filename
            <<" key = "<<key<<" value = "<<value<<endl;
       DHput.putData(key,value);
   // DHput.flush();  // see if it compiles
   // DHput.queryData(key);
       }}

 cout << endl<<endl<<" *********************PUT SUMMARY*************"<<endl<<endl;
 for (set<string>::iterator it=fputnames.begin();it!=fputnames.end();it++){
    cout << "Filename: "<<*it<<endl;
    DataGetHandlerSF<HeadTemp,UIntKey,double> 
        DHget(theader,*it,filetype_id,use_checksums);
    std::set<UIntKey> keys(DHget.getKeys());
    for (set<UIntKey>::iterator kt=keys.begin();kt!=keys.end();kt++)
       cout << kt->getValue()<<endl;
    XMLHandler xmlkeysout;
    DHget.outputKeys(xmlkeysout);
    cout << "outputKeys:"<<xmlkeysout.output()<<endl;
    uint mapsize=DHget.size();
    cout << "size = "<<mapsize<<endl;}

 list<XMLHandler> dhgets(xmlr.find("DoSomeGets"));
 for (list<XMLHandler>::iterator it=dhgets.begin();it!=dhgets.end();it++){
    cout <<endl<<endl<<"  ***** DO SOME GETS *****"<<endl<<endl<< it->output()<<endl;
    ArgsHandler xarg(*it);
    string filename(xarg.getString("FileName"));
    DataGetHandlerSF<HeadTemp,UIntKey,double> 
        DHget(theader,filename,filetype_id,use_checksums);

/*
set<UIntKey> keep;
keep.insert(UIntKey(3));
keep.insert(UIntKey(5));
keep.insert(UIntKey(6));
keep.insert(UIntKey(7));
keep.insert(UIntKey(8));
keep.insert(UIntKey(9));
keep.insert(UIntKey(22));
cout <<endl<<"Keep = "<<DHget.keepKeys(keep)<<endl;
*/

    list<XMLHandler> keyxmls(it->find("Key"));
    for (list<XMLHandler>::iterator kt=keyxmls.begin();kt!=keyxmls.end();kt++){
       ArgsHandler karg(*kt);
       uint key; double value;
       key=karg.getInt("Key");
       cout <<endl<< "Getting from filename "<<filename<<" key = "<<key;
       if (DHget.queryData(key)){
          DHget.getData(key,value);
          cout<<" value = "<<value<<endl;}
       else{
          cout << "query was false"<<endl;}
       value=0;
       bool maybe=DHget.getDataMaybe(key,value);
       cout << "maybe = "<<maybe<<"   value = "<<value<<endl;}}


 list<XMLHandler> dhmultigets(xmlr.find("DoMultiGet"));
 for (list<XMLHandler>::iterator it=dhmultigets.begin();it!=dhmultigets.end();it++){
    cout <<endl<<endl<<"  ***** DO MULTI GET *****"<<endl<<endl<< it->output()<<endl;
    list<XMLHandler> fnamexml=it->find("FileName");
    set<string> fnames;
    for (list<XMLHandler>::iterator ft=fnamexml.begin();ft!=fnamexml.end();ft++){
       ArgsHandler xarg(*ft);
       string filename(xarg.getString("FileName"));
       fnames.insert(filename);}
    list<XMLHandler> keysxml=it->find("Key");
    list<UIntKey> gkeys;
    for (list<XMLHandler>::iterator kt=keysxml.begin();kt!=keysxml.end();kt++){
       ArgsHandler xarg(*kt);
       uint key=xarg.getUInt("Key");
       gkeys.push_back(key);}
    cout << "Getting from set of files: "<<endl;
    for (set<string>::iterator ft=fnames.begin();ft!=fnames.end();ft++)
       cout << "   "<<*ft<<endl;
    DataGetHandlerMF<HeadTemp,UIntKey,double> DHmultiget(
            theader,fnames,filetype_id,use_checksums);

/*
set<UIntKey> mkeep;
mkeep.insert(UIntKey(3));
mkeep.insert(UIntKey(5));
mkeep.insert(UIntKey(6));
mkeep.insert(UIntKey(7));
mkeep.insert(UIntKey(8));
mkeep.insert(UIntKey(9));
mkeep.insert(UIntKey(10));
mkeep.insert(UIntKey(53));
mkeep.insert(UIntKey(55));
mkeep.insert(UIntKey(56));
mkeep.insert(UIntKey(57));
mkeep.insert(UIntKey(58));
mkeep.insert(UIntKey(59));
//mkeep.insert(UIntKey(60));
cout <<endl<<"Keep = "<<DHmultiget.keepKeys(mkeep)<<endl;
*/
/*
set<UIntKey> mkeep;
mkeep.insert(UIntKey(50));
mkeep.insert(UIntKey(51));

    cout << "add file: "<<DHmultiget.addFile(string("test_put.2"),mkeep)<<endl;
    cout << "remove file: "<<endl; DHmultiget.removeFile(string("test_put.1"));
*/

    cout << "Total size = "<<DHmultiget.size()<<endl;
    std::set<UIntKey> keys(DHmultiget.getKeys());
    for (set<UIntKey>::iterator kt=keys.begin();kt!=keys.end();kt++)
       cout << kt->getValue()<<endl;
    XMLHandler xmlkeysout;
    DHmultiget.outputKeys(xmlkeysout);
    cout << "outputKeys:"<<xmlkeysout.output()<<endl;


    std::set<std::string> ff(DHmultiget.getFileNames());
    cout << "filename from getFileNames:"<<endl;
    for (set<string>::iterator xt=ff.begin();xt!=ff.end();xt++)
       cout << "   "<<*xt<<endl;
    for (list<UIntKey>::iterator kt=gkeys.begin();kt!=gkeys.end();kt++){
       cout <<endl<< "Do a get for key = "<<kt->getValue()<<endl;
       bool query=DHmultiget.queryData(*kt);
       cout << "query = "<<query<<endl;
       double value=0;
       if (query){
          DHmultiget.getData(*kt,value);
          cout << "value = "<<value<<endl;}
       value=0;
       if (DHmultiget.getDataMaybe(*kt,value)){
          cout << "from maybe value = "<<value<<endl;}
       else
          cout << "getDataMaybe returned false"<<endl;}
    }
/*    DataGetHandlerSF<HeadTemp,UIntKey,double> 
        DHget(theader,filename,filetype_id,use_checksums);
    if (DHget.queryData(key)){
       DHget.getData(key,value);
       cout<<" value = "<<value<<endl;}
    else{
       cout << "query was false"<<endl;}
    value=0;
    bool maybe=DHget.getDataMaybe(key,value);
    cout << "maybe = "<<maybe<<"   value = "<<value<<endl;}
*/

}



void testBinGetPut(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestBinsGetPut")==0)
 return;

 cout << endl<<endl<<"***************************************************"<<endl<<endl;
 cout << "Testing BinsGetPut"<<endl;
 
 XMLHandler xmlr(xml_in,"TestBinsGetPut");
 MCBinsInfo binfo(xmlr);
 cout << binfo.output()<<endl;
 int nbins=binfo.getNumberOfBins();
 cout << "nbins = "<<nbins<<endl;
 ArgsHandler xxarg(xmlr);
 bool overwrite=xxarg.getBool("Overwrite");
 WriteMode wmode=(overwrite) ? Overwrite : Update;
 bool use_checksums=xxarg.getBool("UseChecksums");

 set<string> fputnames;
 list<XMLHandler> dhputs(xmlr.find("DoPut"));
 for (list<XMLHandler>::iterator it=dhputs.begin();it!=dhputs.end();it++){
    cout <<endl<<endl<<"  ***** DO PUT *****"<<endl<<endl<< it->output()<<endl;
    ArgsHandler xarg(*it);
    string filename(xarg.getString("FileName"));
    fputnames.insert(filename);
    BinsPutHandler BH(binfo,filename,wmode,use_checksums);
    cout << "BH created"<<endl;

    list<XMLHandler> items(it->find("Item"));
    for (list<XMLHandler>::iterator kt=items.begin();kt!=items.end();kt++){
       MCObsInfo rkey(*kt);
       ArgsHandler karg(*kt);
       double value=karg.getReal("Value");
       cout << "Item:"<<endl;
       cout << rkey.output()<<endl;
       cout << "Value = "<<value<<endl;
       Vector<double> bins(nbins,value);
       BH.putData(rkey,bins);}
    }

 cout << endl<<endl<<" *********************PUT SUMMARY*************"<<endl<<endl;
 BinsGetHandler BG(binfo,fputnames,use_checksums);
 std::set<MCObsInfo> keys(BG.getKeys());
 for (set<MCObsInfo>::iterator kt=keys.begin();kt!=keys.end();kt++)
    cout << kt->output()<<endl;
 XMLHandler xmlkeysout;
 BG.outputKeys(xmlkeysout);
 cout << "outputKeys:"<<xmlkeysout.output()<<endl;
 uint mapsize=BG.size();
 cout << "size = "<<mapsize<<endl;


 list<XMLHandler> dhmultigets(xmlr.find("DoMultiGet"));
 for (list<XMLHandler>::iterator it=dhmultigets.begin();it!=dhmultigets.end();it++){
    cout <<endl<<endl<<"  ***** DO MULTI GET *****"<<endl<<endl<< it->output()<<endl;
    list<XMLHandler> fnamexml=it->find("FileName");
    set<string> fnames;
    for (list<XMLHandler>::iterator ft=fnamexml.begin();ft!=fnamexml.end();ft++){
       ArgsHandler xarg(*ft);
       string filename(xarg.getString("FileName"));
       fnames.insert(filename);}
    list<XMLHandler> keysxml=it->find("MCObservable");
    list<MCObsInfo> gkeys;
    for (list<XMLHandler>::iterator kt=keysxml.begin();kt!=keysxml.end();kt++){
       MCObsInfo akey(*kt);
       gkeys.push_back(akey);}
    cout << "Getting from set of files: "<<endl;
    for (set<string>::iterator ft=fnames.begin();ft!=fnames.end();ft++)
       cout << "   "<<*ft<<endl;
    BinsGetHandler BGG(binfo,fnames,use_checksums);

    for (list<MCObsInfo>::iterator kt=gkeys.begin();kt!=gkeys.end();kt++){
       cout << "Do a get for key = "<<kt->output()<<endl;
       bool query=BGG.queryData(*kt);
       cout << "query = "<<query<<endl;
       Vector<double> result;
       if (query){
          BGG.getData(*kt,result);
          cout << "result = "<<result[0]<<endl;}
       Vector<double> result2;
       if (BGG.getDataMaybe(*kt,result2)){
          cout << "from maybe value = "<<result2[0]<<endl;}
       else
          cout << "getDataMaybe returned false"<<endl;}
    }

}



void testSamplingGetPut(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestSamplingsGetPut")==0)
 return;

 cout << endl<<endl<<"***************************************************"<<endl<<endl;
 cout << "Testing SamplingsGetPut"<<endl;

 XMLHandler xmlr(xml_in,"TestSamplingsGetPut");
 MCBinsInfo binfo(xmlr);
 MCSamplingInfo bsamp(xmlr);
 cout << binfo.output()<<endl;
 cout << bsamp.output()<<endl;
 int nsamp=bsamp.getNumberOfReSamplings(binfo);
 cout << "nsamp = "<<nsamp<<endl;
 ArgsHandler xxarg(xmlr);
 bool overwrite=xxarg.getBool("Overwrite");
 WriteMode wmode=(overwrite) ? Overwrite : Update;
 bool use_checksums=xxarg.getBool("UseChecksums");

 set<string> fputnames;
 list<XMLHandler> dhputs(xmlr.find("DoPut"));
 for (list<XMLHandler>::iterator it=dhputs.begin();it!=dhputs.end();it++){
    cout <<endl<<endl<<"  ***** DO PUT *****"<<endl<<endl<< it->output()<<endl;
    ArgsHandler xarg(*it);
    string filename(xarg.getString("FileName"));
    fputnames.insert(filename);
    SamplingsPutHandler BH(binfo,bsamp,filename,wmode,use_checksums);
    cout << "BH created"<<endl;

    list<XMLHandler> items(it->find("Item"));
    for (list<XMLHandler>::iterator kt=items.begin();kt!=items.end();kt++){
       MCObsInfo rkey(*kt);
       ArgsHandler karg(*kt);
       double value=karg.getReal("Value");
       cout << "Item:"<<endl;
       cout << rkey.output()<<endl;
       cout << "Value = "<<value<<endl;
       Vector<double> bins(nsamp,value);
       BH.putData(rkey,bins);}
    }

 cout << endl<<endl<<" *********************PUT SUMMARY*************"<<endl<<endl;
 SamplingsGetHandler BG(binfo,bsamp,fputnames,use_checksums);
 std::set<MCObsInfo> keys(BG.getKeys());
 for (set<MCObsInfo>::iterator kt=keys.begin();kt!=keys.end();kt++)
    cout << kt->output()<<endl;
 XMLHandler xmlkeysout;
 BG.outputKeys(xmlkeysout);
 cout << "outputKeys:"<<xmlkeysout.output()<<endl;
 uint mapsize=BG.size();
 cout << "size = "<<mapsize<<endl;


 list<XMLHandler> dhmultigets(xmlr.find("DoMultiGet"));
 for (list<XMLHandler>::iterator it=dhmultigets.begin();it!=dhmultigets.end();it++){
    cout <<endl<<endl<<"  ***** DO MULTI GET *****"<<endl<<endl<< it->output()<<endl;
    list<XMLHandler> fnamexml=it->find("FileName");
    set<string> fnames;
    for (list<XMLHandler>::iterator ft=fnamexml.begin();ft!=fnamexml.end();ft++){
       ArgsHandler xarg(*ft);
       string filename(xarg.getString("FileName"));
       fnames.insert(filename);}
    list<XMLHandler> keysxml=it->find("MCObservable");
    list<MCObsInfo> gkeys;
    for (list<XMLHandler>::iterator kt=keysxml.begin();kt!=keysxml.end();kt++){
       MCObsInfo akey(*kt);
       gkeys.push_back(akey);}
    cout << "Getting from set of files: "<<endl;
    for (set<string>::iterator ft=fnames.begin();ft!=fnames.end();ft++)
       cout << "   "<<*ft<<endl;
    SamplingsGetHandler BGG(binfo,bsamp,fnames,use_checksums);

    for (list<MCObsInfo>::iterator kt=gkeys.begin();kt!=gkeys.end();kt++){
       cout << "Do a get for key = "<<kt->output()<<endl;
       bool query=BGG.queryData(*kt);
       cout << "query = "<<query<<endl;
       Vector<double> result;
       if (query){
          BGG.getData(*kt,result);
          cout << "result = "<<result[0]<<endl;}
       Vector<double> result2;
       if (BGG.getDataMaybe(*kt,result2)){
          cout << "from maybe value = "<<result2[0]<<endl;}
       else
          cout << "getDataMaybe returned false"<<endl;}
    }

}

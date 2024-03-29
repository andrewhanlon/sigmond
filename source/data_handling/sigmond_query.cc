#include <vector>
#include <list>
#include <set>
#include "xml_handler.h"
#include "io_map.h"
#include "mcobs_info.h"
#include "matrix.h"

using namespace std;

// ***************************************************************

void print_help()
{
 cout << endl;
 cout << " \"sigmond_query\" is used to display some of the information"<<endl;
 cout << "   contained in a Sigmond file, such as header info,"<<endl;
 cout << "   id string, a list of all record keys, etc."<<endl<<endl;
 cout << " Usage:  sigmond_query [options] file"<<endl;
 cout << " Options: -h, --help          display this help and exit"<<endl;
 cout << "          -i, --header        display the header XML info"<<endl;
 cout << "          -n, --numrec        display the number of records"<<endl;
 cout << "          -k, --keys          display all of the record key XMLs"<<endl;
 cout << "          -c, --checksums     display if checksums are stored in file"<<endl;
 cout << "          -e, --endian        display endian-ness of binary data in file"<<endl;
 cout << "          -v, --values        show values of records"<<endl;
 cout << "          -a, --all           display all above information"<<endl;
 cout << "  short form options can be put together, e.g.,  -nce"<<endl;
 cout << " "<<endl;
 cout << " Alternative: sigmond_query -x input.xml  file"<<endl;
 cout << "  where the input.xml file should have the format"<<endl;
 cout << "     <SigMonDQuery>"<<endl;
 cout << "          <Show>header</Show>"<<endl;
 cout << "          <Show>numrec</Show>"<<endl;
 cout << "          <Show>keys</Show>"<<endl;
 cout << "          <Show>checksums</Show>"<<endl;
 cout << "          <Show>endian</Show>"<<endl;
 cout << "          <Show>values</Show> (to show all values"<<endl;
 cout << "          or <ShowValues><MCObservable>...</MCObservable><ShowValues>...."<<endl;
 cout << "          <Show>all</Show>"<<endl;
 cout << "     </SigMonDQuery>"<<endl;
 cout << endl;
}


void printData(const vector<double>& data, char ftype)
{
 cout.precision(15);
 if (ftype=='B'){
    double mean=0.0;
    for (uint k=0;k<data.size();++k)
       mean+=data[k];
    cout << "Unweighted Mean Value of Bins = "<<mean/double(data.size())<<endl<<endl;}
 else if (ftype=='S')
    cout << "Full Sampling Mean Value = "<<data[0]<<endl<<endl;
 for (uint k=0;k<data.size();k++){
    cout << "["<<k<<"] = "<<data[k] << endl;}
}

void printData(const Array<double>& data, char ftype)
{
 cout.precision(15);
 if (ftype=='B'){
    double mean=0.0;
    for (uint k=0;k<data.size();++k)
       mean+=data[k];
    cout << "Unweighted Mean Value of Bins = "<<mean/double(data.size())<<endl<<endl;}
 else if (ftype=='S')
    cout << "Full Sampling Mean Value = "<<data[0]<<endl<<endl;
 for (uint k=0;k<data.size();k++){
    cout << "["<<k<<"] = "<<data[k] << endl;}
}

void printData(const Array<std::complex<double> >& data, char ftype)
{
 cout.precision(15);
 for (uint k=0;k<data.size();k++){
    cout << "["<<k<<"] = "<<data[k] << endl;}
}



template <typename K, typename D>
void squery_outputter(IOMap<K,D>& iom, char ftype, bool header, bool numrec, bool keys, bool csum,
                      bool endian, bool values, list<XMLHandler>& keyxmls)
{
 if (header){
    XMLHandler xmlo; xmlo.set_from_string(iom.getHeader());
    cout << xmlo.output()<<endl;}
 if (numrec)
    cout << "Number of records = "<<iom.size()<<endl;
 if (endian){
    if (iom.isFileLittleEndian())
       cout << "Binary data is in little endian format"<<endl;
    else if (iom.isFileBigEndian())
       cout << "Binary data is in big endian format"<<endl;}
 if (csum){
    if (iom.areChecksumsInFile())
       cout << "Check-sums are stored in file"<<endl;
    else
       cout << "Check-sums are not stored in file"<<endl;}
 if (keys){
    vector<K> thekeys;
    iom.getKeys(thekeys);
    if (keys){
       for (unsigned int k=0;k<thekeys.size();++k){
          XMLHandler xmlk;
          thekeys[k].output(xmlk);
          cout << "Record "<<k<<":"<<endl;
          cout << xmlk.output()<<endl;}}
    }
 if (values){
    vector<K> thekeys;
    iom.getKeys(thekeys);
    for (unsigned int k=0;k<thekeys.size();++k){
       D buffer; iom.get(thekeys[k],buffer);
       cout << "Record "<<k<<":"<<endl;
       printData(buffer,ftype);
       cout << endl;}
    }
 else{
    set<K> keys;
    for (list<XMLHandler>::iterator kt=keyxmls.begin();kt!=keyxmls.end();++kt){
       try{
          K key(*kt);
          keys.insert(key);}
       catch(const std::exception& xp){}}
    for (typename set<K>::const_iterator kt=keys.begin();kt!=keys.end();++kt){
       D buffer; iom.get(*kt,buffer);
       cout << "Record Key: "<<endl<<kt->output()<<endl;
       printData(buffer,ftype);
       cout << endl;}}
}

void calc_simple_jack_samples(const RVector& bins, RVector& samplings);

void calc_simple_boot_samples(const RVector& bins, RVector& samplings);



// *********************************************************


int main(int argc, const char* argv[])
{
 if (argc<2){
    print_help();
    return 0;}

     // convert arguments to C++ strings
 vector<string> tokens(argc-1);
 for (int k=1;k<argc;++k){
    tokens[k-1]=string(argv[k]);}

    // first, check to see if help was requested
 for (unsigned int k=0;k<tokens.size();++k){
    if ((tokens[k]==string("-h"))||(tokens[k]==string("--help"))){
       print_help();
       return 0;}}

    // next, check to see if -x was given
 string inputxmlfile;
 for (unsigned int k=0;k<tokens.size();++k){
    if (tokens[k]==string("-x")){
       if ((k!=0) || (tokens.size()!=3)){
          cout << "Invalid format"<<endl;
          return 0;}
       else{
          inputxmlfile=tokens[1];   
          if (inputxmlfile[0]=='-'){
             cout << "first argument after -x must be input xml file name"<<endl;
             return 0;}
          break;}}}

 bool valid=true,header=false,numrec=false,keys=false,
      csum=false,endian=false,values=false;
 list <XMLHandler> keyxmls;

 if (inputxmlfile.empty()){
    // now parse the command options and get the file name
    for (unsigned int k=0;k<tokens.size()-1;++k){
       if ((tokens[k][0]!='-')||(tokens[k].length()<2)){
          cout << "invalid argument "<<tokens[k]<<endl;
          valid=false;}
       if (tokens[k][1]=='-'){
          string longopt(tokens[k].substr(2));
          if (longopt==string("header")) header=true;
          else if (longopt==string("numrec")) numrec=true;
          else if (longopt==string("keys")) keys=true;
          else if (longopt==string("checksums")) csum=true;
          else if (longopt==string("endian")) endian=true;
          else if (longopt==string("values")) values=true;
          else if (longopt==string("all")){
             header=numrec=keys=csum=endian=values=true;}
          else if (longopt!=string("help")){
             cout << "invalid argument "<<tokens[k]<<endl;
             valid=false;}}
       else{
          for (unsigned int j=1;j<tokens[k].length();++j){
             char shortopt=tokens[k][j];
             if (shortopt=='i') header=true;
             else if (shortopt=='n') numrec=true;
             else if (shortopt=='k') keys=true;
             else if (shortopt=='c') csum=true;
             else if (shortopt=='e') endian=true;
             else if (shortopt=='v') values=true;
             else if (shortopt=='a'){
                header=numrec=keys=csum=endian=values=true;}
             else if (shortopt!='h'){
                cout << "invalid short-form option "<<shortopt<<endl;
                valid=false;}}}
       }}
 else{
    XMLHandler xmlin;
    xmlin.set_from_file(inputxmlfile);
    XMLHandler xmls(xmlin,"SigMonDQuery");
    list<XMLHandler> xmlss=xmls.find("Show");
    for (list<XMLHandler>::iterator xt=xmlss.begin();xt!=xmlss.end();++xt){
        string longopt;
        try{
           xmlread(*xt,"Show",longopt,"SigMonDQuery");
           if (longopt==string("header")) header=true;
           else if (longopt==string("numrec")) numrec=true;
           else if (longopt==string("keys")) keys=true;
           else if (longopt==string("checksums")) csum=true;
           else if (longopt==string("endian")) endian=true;
           else if (longopt==string("values")) values=true;
           else if (longopt==string("all")){
               header=numrec=keys=csum=endian=values=true;}}
        catch(const std::exception& xp){}}
    if (!values){
       keyxmls=xmls.find("ShowValues");}}

 string filename(tokens[tokens.size()-1]);
 if (filename[0]=='-'){
    cout << "last argument must be name of file: "<<filename<<endl;
    valid=false;}
 if (!valid) return 0;

 string ID;
 if (!IOMapPeekID(ID,filename)){
    cout << "Error with file "<<filename<<": could not extract ID string"<<endl;
    return 0;}

 IOMap<MCObsInfo,vector<double> > iom;
 string sID("Sigmond--SamplingsFile");
 string bID("Sigmond--BinsFile");
 IOMap<UIntKey,Array<std::complex<double> > > iopc;
 string spIDc("Sigmond--SinglePivotFile-CN");
// string rpIDc("Sigmond--RollingPivotFile-CN");
 IOMap<UIntKey,Array<double> > iopr;
 string spIDr("Sigmond--SinglePivotFile-RN");
// string rpIDr("Sigmond--RollingPivotFile-RN");

 try{
 if (ID==sID){
    iom.openReadOnly(filename,sID);
    cout <<endl<< "This is a Sigmond samplings file"<<iom.get_format()<<endl;
    squery_outputter(iom,'S',header,numrec,keys,csum,endian,values,keyxmls);}
 else if (ID==bID){
    iom.openReadOnly(filename,bID);
    cout <<endl<< "This is a Sigmond bins file"<<iom.get_format()<<endl;
    squery_outputter(iom,'B',header,numrec,keys,csum,endian,values,keyxmls);}
 else if (ID==spIDr){
    iopr.openReadOnly(filename,spIDr);
    cout <<endl<< "This is a Sigmond single pivot file with real numbers"<<iopr.get_format()<<endl;
    squery_outputter(iopr,'P',header,numrec,keys,csum,endian,values,keyxmls);}
// else if (ID==rpIDr){
//    cout <<endl<< "This is a Sigmond rolling pivot file with real numbers"<<endl;
//    iopr.openReadOnly(filename,rpIDr);
//    squery_outputter(iopr,'P',header,numrec,keys,csum,endian,values,keyxmls);}
 else if (ID==spIDc){
    iopc.openReadOnly(filename,spIDc);
    cout <<endl<< "This is a Sigmond single pivot file with complex numbers"<<iopc.get_format()<<endl;
    squery_outputter(iopc,'P',header,numrec,keys,csum,endian,values,keyxmls);}
// else if (ID==rpIDc){
//    cout <<endl<< "This is a Sigmond rolling pivot file with complex numbers"<<endl;
//    iopc.openReadOnly(filename,rpIDc);
//    squery_outputter(iopc,'P',header,numrec,keys,csum,endian,values,keyxmls);}
 else{
    cout <<endl<< "This file type is not known to Sigmond"<<endl;}}
 catch(const std::exception& msg){
    cout << "Error opening file "<<filename<<endl;
    return 0;}

 return 0;
}

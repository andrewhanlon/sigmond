#include <vector>
#include "xml_handler.h"
#include "io_map.h"
#include "mcobs_info.h"
#include "matrix.h"

using namespace std;

// **********************************************************

void print_help()
{
 cout << endl;
 cout << " \"sigmond_query\" is used to display some of the information"<<endl;
 cout << "   contained in a Sigmond file, such as header info,"<<endl;
 cout << "   id string, a list of all record keys, etc."<<endl<<endl;
 cout << " Usage:  laph_query [options] file"<<endl;
 cout << " Options: -h, --help          display this help and exit"<<endl;
 cout << "          -i, --header        display the header XML info"<<endl;
 cout << "          -n, --numrec        display the number of records"<<endl;
 cout << "          -k, --keys          display all of the record key XMLs"<<endl;
 cout << "          -c, --checksums     display if checksums are stored in file"<<endl;
 cout << "          -e, --endian        display endian-ness of binary data in file"<<endl;
 cout << "          -v, --values        show values of records"<<endl;
 cout << "          -a, --all           display all above information"<<endl;
 cout << "  short form options can be put together, e.g.,  -nce"<<endl;
 cout << endl;
}

void printData(vector<double> & data)
{
 for (std::vector<double>::iterator vec_it=data.begin(); vec_it!=data.end(); ++vec_it)
    cout << *vec_it << endl;
}

void printData(Array<std::complex<double> > & data)
{
 cout << "Not supported yet..." << endl;
}

void printData(Array<double> & data)
{
 cout << "Not supported yet..." << endl;
}


template <typename K, typename D>
void outputter(IOMap<K,D>& iom, bool header, bool numrec, bool keys, bool csum, bool endian, bool values)
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
 if (values){
    vector<K> thekeys;
    iom.getKeys(thekeys);
    for (unsigned int k=0;k<thekeys.size();++k){
       D buffer; iom.get(thekeys[k],buffer);
       cout << "Record "<<k<<":"<<endl;
       printData(buffer);
       cout << endl;}}
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
}

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

    // now parse the command options and get the file name
 bool valid=true,header=false,numrec=false,keys=false,
      csum=false,endian=false,values=false;
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
    }
 string filename(tokens[tokens.size()-1]);
 if (filename[0]=='-'){
    cout << "last argument must be name of file: "<<filename<<endl;
    valid=false;}
 if (!valid) return 0;

 ifstream fin(filename.c_str(),ios::binary);
 if (!fin){
    cout << "Error opening file "<<filename<<endl;
    return 0;}
 char idstring[33];
 if (!fin.read(idstring,33)){
    cout << "Error: could not extract ID string from file "<<filename<<endl;
    return 0;}
 string ID(&idstring[1],32); 
 ID=tidyString(ID); 

 IOMap<MCObsInfo,vector<double> > iom;
 string sID("Sigmond--SamplingsFile");
 string bID("Sigmond--BinsFile");
 IOMap<UIntKey,Array<std::complex<double> > > iopc;
 string spIDc("Sigmond--SinglePivotFile-CN");
 string rpIDc("Sigmond--RollingPivotFile-CN");
 IOMap<UIntKey,Array<double> > iopr;
 string spIDr("Sigmond--SinglePivotFile-RN");
 string rpIDr("Sigmond--RollingPivotFile-RN");

 try{
 if (ID==sID){
    cout <<endl<< "This is a Sigmond samplings file"<<endl;
    iom.openReadOnly(filename,sID);
    outputter(iom,header,numrec,keys,csum,endian,values);}
 else if (ID==bID){
    cout <<endl<< "This is a Sigmond bins file"<<endl;
    iom.openReadOnly(filename,bID);
    outputter(iom,header,numrec,keys,csum,endian,values);}
 else if (ID==spIDr){
    cout <<endl<< "This is a Sigmond single pivot file with real numbers"<<endl;
    iopr.openReadOnly(filename,spIDr);
    outputter(iopr,header,numrec,keys,csum,endian,values);}
 else if (ID==rpIDr){
    cout <<endl<< "This is a Sigmond rolling pivot file with real numbers"<<endl;
    iopr.openReadOnly(filename,rpIDr);
    outputter(iopr,header,numrec,keys,csum,endian,values);}
 else if (ID==spIDc){
    cout <<endl<< "This is a Sigmond single pivot file with complex numbers"<<endl;
    iopc.openReadOnly(filename,spIDc);
    outputter(iopc,header,numrec,keys,csum,endian,values);}
 else if (ID==rpIDc){
    cout <<endl<< "This is a Sigmond rolling pivot file with complex numbers"<<endl;
    iopc.openReadOnly(filename,rpIDc);
    outputter(iopc,header,numrec,keys,csum,endian,values);}
 else{
    cout <<endl<< "This file type is not known to Sigmond"<<endl;}}
 catch(const std::exception& msg){
    cout << "Error opening file "<<filename<<endl;
    return 0;}

 return 0;
}

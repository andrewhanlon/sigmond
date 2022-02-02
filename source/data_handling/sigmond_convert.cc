#include <cstdio>
#include <ctime>
#include <vector>
#include <map>
#include <iostream>
#include "xml_handler.h"
#include "io_map.h"
#include "mcobs_info.h"
#include "matrix.h"
#include "bins_info.h"
#include "bins_handler.h"
#include "samplings_handler.h"

using namespace std;

// *****************************************************************************
// *                                                                           *
// *   This program is used to convert Monte Carlo bins in textual (ascii)     *
// *   format into a format suitable for input to and subsequent analysis      *
// *   by SigMonD.  The output will be real double format.  No rebinning or    *
// *   other changes to the data are done.                                     *
// *                                                                           *
// *   Run this program with one argument, which is the name of an XML input   *
// *   file.  The XML input file must have content of the form below:          *
// *                                                                           *
// *    <SigmondConvert>                                                       *
// *     <LogFile>convert.log</LogFile>                                        *
// *     <EchoXML/>                                                            *
// *                                                                           *
// *      <Task>                                                               *
// *       <Action>ConvertBins</Action>                                        *
// *       <MCBinsInfo>                                                        *
// *         <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo>       *
// *         <TweakEnsemble>  (optional)                                       *
// *           <Rebin>2</Rebin>                                                *
// *           <Omissions>2 7 11</Omissions>                                   *
// *         </TweakEnsemble>                                                  *
// *       </MCBinsInfo>                                                       *
// *       <MCBinsData>                                                        *
// *         <Entry>                                                           *
// *           <MCObservable>...</MCObservable>                                *
// *           <Values>...</Values>  (real values, separated by                *
// *                                  spaces, in order of bin index,           *
// *                                  skipping omitted bins)                   *
// *         </Entry>                                                          *
// *            ... (other entries)                                            *
// *       </MCBinsData>                                                       *
// *       <OutputFileName>output.bins</OutputFileName>                        *
// *       <FileFormat>fstr</FileFormat> (or hdf5: default if absent)          *
// *      </Task>                                                              *
// *                                                                           *
// *         .... (other tasks)                                                *
// *                                                                           *
// *    </SigmondConvert>                                                      *
// *                                                                           *
// *                                                                           *
// *****************************************************************************


class ConvertHandler
{

   std::ofstream clog;
   std::list<XMLHandler> xmltasks;

   typedef void (ConvertHandler::*task_ptr)(XMLHandler&, XMLHandler&);
   std::map<std::string, task_ptr>    m_task_map;

       // Prevent copying ... handler might contain large
       // amounts of data

#ifndef NO_CXX11
   ConvertHandler() = delete;
   ConvertHandler(const ConvertHandler&) = delete;
   ConvertHandler& operator=(const ConvertHandler&) = delete;
#else
   ConvertHandler();
   ConvertHandler(const ConvertHandler&);
   ConvertHandler& operator=(const ConvertHandler&);
#endif

 public:
    
   ConvertHandler(XMLHandler& xmlin);
   ~ConvertHandler();

   void do_tasks();

 private:

   void do_task(XMLHandler& xml_in, XMLHandler& output);

   std::string get_date_time();

   void finish_log();


       // The important task subroutines

   void convertBins(XMLHandler &xml_in, XMLHandler& output);

};

// ***************************************************************

     // set up the known tasks, create logfile, open stream for logging

ConvertHandler::ConvertHandler(XMLHandler& xmlin)
{
 if (xmlin.get_node_name()!="SigmondConvert")
    throw(std::invalid_argument("Input file must have root tag <SigmondConvert>"));
 string nowstr=get_date_time();

 string logfile;
 if (xmlin.count_among_children("LogFile")==1)
    xmlread(xmlin,"LogFile",logfile,"ConvertHandler");
 else
    logfile=string("sigmondconvert_log_")+nowstr+".xml";

 clog.open(logfile.c_str());
 if (!clog.is_open()){
    cout << "Could not open log file "<<logfile<<" for output"<<endl;
    throw(std::invalid_argument("Could not write to log file"));}
 clog << "<LogSigmondConvert>"<<endl;

 if (xmlin.count_among_children("EchoXML")>=1){
    string input(xmlin.output());
    int pos=input.find("<SigmondConvert>");
    input.erase(0,pos+15);
    pos=input.find("</SigmondConvert>");
    input.erase(pos,string::npos);
    clog << " <InputXML>";
    clog << input<<"</InputXML>"<<endl;}

 xmltasks=xmlin.find_among_children("Task");
 
 m_task_map["ConvertBins"]=&ConvertHandler::convertBins;
// m_task_map["ConvertSamplings"]=&ConvertHandler::convertSamplings;

}


     // delete data, finish up logging, close log file

ConvertHandler::~ConvertHandler()
{
 finish_log();
 m_task_map.clear();
}


void ConvertHandler::finish_log()
{
 string nowstr=get_date_time();
 clog << " <FinishDateTime>"<<nowstr<<"</FinishDateTime>"<<endl;
 clog << "</SigmondConvert>"<<endl<<endl;
 clog.close();
}


void ConvertHandler::do_tasks()
{
 for (std::list<XMLHandler>::iterator it=xmltasks.begin();it!=xmltasks.end();++it){
    clog << endl<<"<BeginTask>****************************************</BeginTask>"<<endl;
    clog << endl<<"<Task>"<<endl;
    XMLHandler xmlout;
    do_task(*it,xmlout);
    clog << xmlout.output()<<endl;
    clog << "</Task>"<<endl;}
}


void ConvertHandler::do_task(XMLHandler& xml_task, XMLHandler& xml_out)
{
 try{
    XMLHandler xmlt(xml_task);
    if (xmlt.get_node_name()!="Task"){
       throw(std::invalid_argument("Input to do_task is not a Task tag"));}
    xmlt.seek_first_child();
    if (xmlt.get_node_name()!="Action"){
       throw(std::invalid_argument("Need Action tag as first child of Task tag"));}

    xmlt.seek_root();
    xml_child_assert(xmlt,"Action","do_task");
    if (!xmlt.is_simple_element()) 
       throw(std::invalid_argument("Action tag is not simple XML element"));
    string task_action=xmlt.get_text_content();

    map<string,task_ptr >::iterator taskit=m_task_map.find(task_action);
    if (taskit!=m_task_map.end()){
       (this->*(taskit->second))(xml_task,xml_out);}  // do the task!!
    else{
       throw(std::invalid_argument((string("Unknown task name: ")
               +task_action).c_str()));}}   // unknown task?
 catch(const std::exception& errmsg){
    if (xml_out.empty())
       xml_out.set_root("Error",string(errmsg.what()));
    else
       xml_out.put_child("Error",string(errmsg.what()));}
}


string ConvertHandler::get_date_time()
{
 time_t rawtime;
 struct tm *timeinfo;
 time(&rawtime);
 timeinfo=localtime(&rawtime);
 string dt=asctime(timeinfo);
 for (unsigned int k=0;k<dt.length();k++)
    if (dt[k]==' ') dt[k]='_';
 return tidyString(dt);
}


// ***************************************************************************************
 
void ConvertHandler::convertBins(XMLHandler& xml_in, XMLHandler& xmlout)
{
 MCBinsInfo binfo(xml_in);
 string binfile_name;
 xmlread(xml_in,"OutputFileName",binfile_name,"ConvertHandler");
 string fformat("default"); char ffmt='D';
 xmlreadifchild(xml_in,"FileFormat",fformat);
 if (fformat=="fstr") ffmt='F';
 else if (fformat=="hdf5") ffmt='H';
 else if (fformat=="default") ffmt='D';
 else throw(std::invalid_argument("<FileFormat> must be ftr or hdf5 or default in convertBins"));
 WriteMode wmode = Protect;  // protect mode
 BinsPutHandler BH(binfo,binfile_name,wmode=Protect,false,ffmt);
 XMLHandler xmld(xml_in,"MCBinsData");
 std::list<XMLHandler> entries=xmld.find("Entry");
 for (std::list<XMLHandler>::iterator it=entries.begin();it!=entries.end();++it){
    MCObsInfo rkey(*it);
    std::vector<double> data;
    xmlreadchild(*it,"Values",data,"ConvertHandler");
    BH.putData(rkey,data);}

 xmlout.set_root("ConvertedBins");
 xmlout.put_child("WroteToFile",binfile_name);
 xmlout.put_child("NumberOfEntries",make_string(entries.size()));
}

// ***************************************************************************************



int main(int argc, const char* argv[])
{
     // convert arguments to C++ strings
 vector<string> tokens(argc-1);
 for (int k=1;k<argc;++k){
    tokens[k-1]=string(argv[k]);}

 if (tokens.size()!=1){
    cout << "Error: requires a file name as the only argument"<<endl;
    return 1;}

 try{
    XMLHandler xmltask;  
    xmltask.set_exceptions_on();
    if (tokens.size()>0){
       string filename(tokens[0]);
       xmltask.set_from_file(filename);}

        // set up the task handler
    ConvertHandler tasker(xmltask);

        // do the tasks
    tasker.do_tasks();
    }
 catch(const std::exception& msg){
    cout << "Error: "<<msg.what()<<endl;
    return 1;}

 return 0;
}


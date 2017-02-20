#include "xml_handler.h"
#include "stopwatch.h"
#include <cstdio>
#include <ctime>
#include <map>
#include "testing.h"

using namespace std;


void output_datetime()
{
 time_t rawtime;
 struct tm *timeinfo;
 time(&rawtime);
 timeinfo=localtime(&rawtime);
 cout << "  Current date/time: "<<asctime(timeinfo);
}

void doSigmondTests(XMLHandler& xml_rdr,int taskcount) 
{
 testMatrix(xml_rdr);
 testBootstrapper(xml_rdr);
 testCorrDataHandler(xml_rdr);
 testCorrelatorInfo(xml_rdr);
 testCorrMatEstimates(xml_rdr);
 testVEVDataHandler(xml_rdr);
 testTaskHandler(xml_rdr);
 testXMLHandler(xml_rdr);
 testArgsHandler(xml_rdr);
 testOperatorInfo(xml_rdr);
 testMCObsInfo(xml_rdr);
 testMCObsGetHandler(xml_rdr);
 testMCObsGetHandlerFake(xml_rdr);
 testMCObsHandler(xml_rdr);
 testMCObsHandler2(xml_rdr);
 testMCObsHandlerIO(xml_rdr);
 testGracePlot(xml_rdr);
 testmulticompare(xml_rdr);
 testCorrelatorMatrixInfo(xml_rdr);
 testChiSquare(xml_rdr,taskcount);
 testGamma(xml_rdr);
 testEffEnergy(xml_rdr);
 testChisqTcorr(xml_rdr,taskcount);
 testRotateCorrelator(xml_rdr,taskcount);
 testPivotCorrelator(xml_rdr,taskcount);
 testPivotCorrelator0(xml_rdr,taskcount);
 testMCBinsInfo(xml_rdr);
 testMCSamplingInfo(xml_rdr);
 testDataGetPut(xml_rdr);
 testBinGetPut(xml_rdr);
 testSamplingGetPut(xml_rdr);
 testReorder(xml_rdr);
}



  //  Constructor sets up the known tasks.  Call
  //  member function "do_task" to perform the task.
  //  "TaskMap" is a map that associates a string
  //  containing a task name to a function pointer.

class Tasker
{
   typedef void (*task_ptr)(XMLHandler&,int);
   map<string, task_ptr> TaskMap;

 public:
    
   Tasker();
   ~Tasker(){}
   void do_task(XMLHandler& xml_in, bool echo, int taskcount);

};

     // set up the known tasks

Tasker::Tasker()
{
 TaskMap["SIGMOND_TEST"]=&doSigmondTests;
};


void Tasker::do_task(XMLHandler& xml_task, bool echo, int taskcount)
{
 if (echo){
    cout << "Input XML for this task:"<<endl
         <<xml_task.output()<<endl;}
 XMLHandler xmlt(xml_task);
 xml_child_assert(xmlt,"Name","do_task");
 if (!xmlt.is_simple_element()) 
    throw(std::invalid_argument("Name tag is not simple XML element"));
 string task_name=xmlt.get_text_content();
 cout << "  Task name = "<<task_name<<endl;

 map<string,task_ptr >::iterator taskit=TaskMap.find(task_name);
 if (taskit!=TaskMap.end()){
    (*(taskit->second))(xml_task,taskcount);}  // do the task!!
 else{
    throw(std::invalid_argument("Unknown task name"));}   // unknown task?
}


// ****************************************************************
// *                                                              *
// *          Main driver program to run all tasks                *
// *                                                              *
// *   Program takes a single argument that is the name of the    *
// *   input file, and standard output is used.  Use redirection  *
// *   for output to a log file.  Input file must contain a       *
// *   single XML document with root tag named "SigMonD".         *
// *   Inside the root tag should be one or more <Task> tags,     *
// *   and inside each <Task> element should be a <Name> tag      *
// *   whose content is the name of the task.  The name must      *
// *   be one of the allowed names specified in the "do_task"     *
// *   subroutine.  If a tag <EchoXML/> is present as a child of  *
// *   the root tag, then the XML input is echoed to standard     *
// *   output.                                                    *
// *                                                              *
// *   Sample input XML:                                          *
// *                                                              *
// *    <SigMonD>                                                 *
// *       <EchoXML/>                                             *
// *       <Task>                                                 *
// *          <Name>Task 1</Name>                                 *
// *       </Task>                                                *
// *       <Task>                                                 *
// *          <Name>Task 2</Name>                                 *
// *       </Task>                                                *
// *    </SigMonD>                                                *
// *                                                              *
// ****************************************************************


int main(int argc, const char* argv[])
{
 cout << endl << "Starting SigMonD tests"<<endl;
 output_datetime();
  
 if (argc!=2){
    cerr << "  SigMonD testing requires one argument: "
         << " the name of an input XML file"<<endl;
    cerr << "   ... exiting..."<<endl;
    return 1;}

 {StopWatch rolex,swatch;
 rolex.reset();
 rolex.start();
 string input_file = string( argv[1] );
 cout << "  Name of input XML file: "<<input_file<<endl<<endl;

 XMLHandler xml_in;
 xml_in.set_from_file(input_file);
 if (xml_in.fail()){
    cerr << "  Unable to read/parse XML content in input file"<<endl;
    cerr << "  ... exiting..."<<endl;
    return 1;}
 if (xml_in.get_node_name()!="SigMonD"){
    cerr << "  Root tag of input XML must be named \"SigMonD\""<<endl;
    cerr << "  ... exiting..."<<endl;
    return 1;}
 xml_in.set_exceptions_off();


 bool echo=false;
 int ntasks = 0;
 try{
    xml_in.seek_first_child();
    while (xml_in.good()){
       if (xml_in.get_node_name()=="Task") ntasks++;
       else if (xml_in.get_node_name()=="EchoXML") echo=true;
       else throw(std::invalid_argument("Invalid XML input"));
       xml_in.seek_next_sibling();}}
 catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    return 1;}

 cout << "  Number of tasks is "<<ntasks<<endl;
 xml_in.seek_root();
 xml_in.seek_first_child();
 Tasker T;

 for (int task=1;task<=ntasks;++task){
    swatch.reset();
    swatch.start();
    try{
       cout << endl<<endl<<"Starting Task "<<task<<endl<<endl;
       if (xml_in.fail()) throw(std::invalid_argument("XML input error"));
       if (xml_in.get_node_name()=="EchoXML")
          xml_in.seek_next_sibling();

       XMLHandler xml_task(xml_in);

           // the main task
       T.do_task(xml_task,echo,task);

       }
    catch(const std::exception& err){
       cerr << "Error on Task "<<task<<":  "<<err.what()<<endl;}
    swatch.stop();
    cout << "Task "<<task<<" done using time = " << swatch.getTimeInSeconds() 
         << " secs" << endl;
    xml_in.seek_next_sibling();
    }

 rolex.stop();
 cout <<endl<<endl;
 cout << "SigMonD testing: total time = "<< rolex.getTimeInSeconds() 
      << " secs" << endl;}
 cout << "SigMonD testing: completion" << endl;
 output_datetime();

 return 0;
}


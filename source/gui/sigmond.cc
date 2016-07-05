#include "xml_handler.h"
#include <cstdio>
#include <ctime>
#include <vector>
#include <map>
#include <iostream>
//#include "analysis_handler.h"

using namespace std;


void output_datetime()
{
 time_t rawtime;
 struct tm *timeinfo;
 time(&rawtime);
 timeinfo=localtime(&rawtime);
 cout << "  Current date/time: "<<asctime(timeinfo);
}

void print_help()
{
 cout << endl<<endl;
 cout << " \"sigmond\" is used for SIGnal extraction from MONte carlo Data."<<endl<<endl;
 cout << " Usage:"<<endl<<endl;
 cout << "    sigmond -h           display help and exit"<<endl;
 cout << "    sigmond -i [file]    prompt for info using command line interface (file optional)"<<endl;
 cout << "    sigmond -x [file]    prompt for info using graphical user interface (file optional)"<<endl;
 cout << "    sigmond -b file      batch mode using instructions in file"<<endl<<endl;
 cout << " If present, \"file\" must be the name of a file containing XML text of the form:"<<endl<<endl;
 cout << "     <SigMonD>"<<endl;
 cout << "        <ProjectType>...</ProjectType>"<<endl;
 cout << "        <Data>"<<endl;
 cout << "          <FileListInfo>...</FileListInfo>     <!-- data files --> "<<endl;
 cout << "              ....  "<<endl;                                          
 cout << "        </Data>  "<<endl;                                             
 cout << "        <Logfile>output.log</Logfile>"<<endl;
 cout << "        <Task>                <!-- batch mode only -->"<<endl;
 cout << "           <Name>Task 1</Name>"<<endl;
 cout << "        </Task> "<<endl;
 cout << "              .... "<<endl;
 cout << "     </SigMonD> "<<endl<<endl;
 cout << " <Task> tags are used only if the program runs in batch mode, and the tasks are"<<endl;
 cout << " done in the order specified in the XML file.  In interactive mode, all <Task> tags"<<endl;
 cout << " are ignored.  Read comments in the .h files of the source code for more information"<<endl;
 cout << " about the XML text needed to carry out the data analysis.  In interactive"<<endl;
 cout << " mode, limited help is sometimes available."<<endl<<endl;
}


// ****************************************************************
// *                                                              *
// *          Main driver program to run all tasks                *
// *                                                              *
// *   Program takes a single argument that is the name of the    *
// *   input file, and standard output is used.  Use redirection  *
// *   for output to a log file.  Input file must contain a       *
// *   single XML document with root tag named "LastLaph".        *
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
// *       <ProjectType>SingleRealCorrelator</ProjectType>        * 
// *       <Data>                                                 *
// *         <FileListInfo>...</FileListInfo>                     *
// *         <FileListInfo>...</FileListInfo>     <-- data files  *
// *             ....                                             *
// *       </Data>                                                *
// *       <Logfile>output.log</Logfile>                          *
// *       <Task>                                                 *
// *          <Name>Task 1</Name>                                 *
// *       </Task>                                                *
// *       <Task>                                                 *
// *          <Name>Task 2</Name>                                 *
// *       </Task>                                                *
// *    </SigMonD>                                                *
// *                                                              *
// *   <FileListInfo>                                             *
// *      <FileNameStub>  ...  </FileNameStub>                    *
// *      <MaxFileNumber> ...  </MaxFileNumber>                   *
// *      <MinFileNumber> ...  </MinFileNumber> (default=0)       *
// *   </FileListInfo>                                            *
// ****************************************************************


int main(int argc, const char* argv[])
{
 if (argc<2){
    print_help();
    return 0;}

     // convert arguments to C++ strings
 vector<string> tokens(argc-1);
 for (int k=1;k<argc;++k){
    tokens[k-1]=string(argv[k]);}

 if (tokens.size()>2){
    cout << "No more than two arguments allowed: use -h for help"<<endl<<endl;
    return 1;}

 else if (tokens[0]==string("-h")){

    print_help();
    return 0;}

 else if (tokens[0]==string("-i")){

    cout << "command line interface prompting"<<endl;
    XMLHandler xmlh;
    if (tokens.size()==1){
       cout << "Will prompt for info to make starting XML"<<endl;}
    else{
       string filename(tokens[1]);
       xmlh.set_from_file(filename);}
   // AnalysisHandler mcanal(xmlh);
    cout << "Will prompt for tasks"<<endl;}

 else if (tokens[0]==string("-x")){

    cout << "graphical user interface prompting"<<endl;
    XMLHandler xmlh;
    if (tokens.size()==1){
       cout << "Will prompt for info to make starting XML"<<endl;}
    else{
       string filename(tokens[1]);
       cout << "Filename is "<<filename<<endl;
       xmlh.set_from_file(filename);}
   // AnalysisHandler mcanal(xmlh);
    cout << "Will prompt for tasks"<<endl;}

 else if ((tokens[0]==string("-b"))&&(tokens.size()==1)){

    cout << "Error: batch mode requires a file name as last argument"<<endl;
    return 1;}

 else if (((tokens[0]==string("-b"))&&(tokens.size()==2))||
         ((tokens.size()==1)&&(tokens[0]!=string("-b")))){

    cout << "batch mode"<<endl;
    XMLHandler xmlh;
    string filename(tokens[tokens.size()-1]);
    cout << "Filename is "<<filename<<endl;
    xmlh.set_from_file(filename);
  //  AnalysisHandler mcanal(xmlh);
    }

 else{

    cout << "Invalid command line arguments"<<endl<<endl;
    return 1;}

 output_datetime();
  
 return 0;
}


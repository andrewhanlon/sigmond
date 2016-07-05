#include "task_handler.h"
#include <cstdio>
#include <ctime>
#include <vector>
#include <map>
#include <iostream>

#define RUNMODE batch
//#define RUNMODE "cli"
//#define RUNMODE "gui" 



using namespace std;

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
     // convert arguments to C++ strings
 vector<string> tokens(argc-1);
 for (int k=1;k<argc;++k){
    tokens[k-1]=string(argv[k]);}

#if (RUNMODE == batch)
 if (tokens.size()!=2){
    cout << "Error: batch mode requires a file name as the only argument"<<endl;
    return 1;}
#endif

 if (tokens.size()>1){
    cout << "No more than one argument allowed: input file name"<<endl<<endl;
    return 1;}

 XMLHandler xmltask;
 if (tokens.size()>0){
    string filename(tokens[1]);
    cout << "Filename is "<<filename<<endl;
    xmltask.set_from_file(filename);}

 TaskHandler tasker(xmltask);
  
 return 0;
}


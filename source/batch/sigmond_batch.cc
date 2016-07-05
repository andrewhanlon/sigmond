#include "task_handler.h"
#include <cstdio>
#include <ctime>
#include <vector>
#include <map>
#include <iostream>

using namespace std;

// ******************************************************************************
// *                                                                            *
// *          Main driver program to run "SigMonD" in batch mode                *
// *                                                                            *
// *   Program takes a single argument that is the name of the input file.      *
// *   Input file must contain a single XML document with root tag named        *
// *   <SigMonD>.  The input XML must have the form below:                      *
// *                                                                            *
// *    <SigMonD>                                                               *
// *       <ProjectName>NameOfProject</ProjectName>                             * 
// *       <Logfile>output.log</Logfile>                                        *
// *       <EchoXML/>                                                           *
// *       <MCObservables>  ...  </MCObservables>                               *
// *       <Bootstrapper>  ...  </Bootstrapper>                                 *
// *       <TweakEnsemble> ... </TweakEnsemble>                                 *
// *                                                                            *
// *       <Task><Action>...</Action> ...  </Task>                              *
// *       <Task><Action>...</Action> ...  </Task>                              *
// *           ....                                                             *
// *    </SigMonD>                                                              *
// *                                                                            *
// *                                                                            *
// *   (a) If <ProjectName> is missing, a default name will be created.         *
// *                                                                            *
// *   (b) If <Logfile> is missing, a default name for the log file is used.    *
// *                                                                            *
// *   (c) If <EchoXML> is missing, the input XML will not be written to the    *
// *       log file.                                                            *
// *                                                                            *
// *   (d) <MCObservables> describes the data to be input for analysis. See     *
// *       the class "MCObsGetHandler" in "source/laph_data/obs_get_handler.h"  *
// *       for a description of the XML needed in this tag.                     *
// *                                                                            *
// *   (e) The tag <Bootstrapper> is optional but controls how bootstrapping    *
// *       is done; default values are used if absent.  It has the form         *
// *       below, where each tag is optional.  See comments in the file         *
// *       "source/analysis/bootstrapper.h" for a desription of these tags.     *
// *                                                                            *
// *         <Bootstrapper>                                                     *
// *            <NumberResamplings>2048</NumberResamplings>                     *
// *            <Seed>6754</Seed>                                               *
// *            <BootSkip>127</BootSkip>                                        *
// *            <Precompute/>                                                   *
// *         </Bootstrapper>                                                    *
// *                                                                            *
// *   (f) The tag <TweakEnsemble> is optional: it controls rebinning the       *
// *       data, and possibly omitting certain configurations in the            *
// *       ensemble.  The XML must have the form below, where each tag          *
// *       is optional.                                                         *
// *                                                                            *
// *         <TweakEnsemble>                                                    *
// *           <Rebin>4</Rebin>                                                 *
// *           <Omissions> 6 9 88</Omissions>                                   *
// *         </TweakEnsemble>                                                   *
// *                                                                            *
// *   (g) The <Task> tags are needed in "batch" mode, but can be omitted in    *
// *   "cli" or "gui".  Each <Task> tag must begin with an <Action> tag.        *
// *                                                                            *
// *                                                                            *
// ******************************************************************************



int main(int argc, const char* argv[])
{
     // convert arguments to C++ strings
 vector<string> tokens(argc-1);
 for (int k=1;k<argc;++k){
    tokens[k-1]=string(argv[k]);}

 if (tokens.size()!=1){
    cout << "Error: batch mode requires a file name as the only argument"<<endl;
    return 1;}

 try{
    XMLHandler xmltask;  
    xmltask.set_exceptions_on();
    if (tokens.size()>0){
       string filename(tokens[0]);
       xmltask.set_from_file(filename);}

        // set up the task handler
    TaskHandler tasker(xmltask);

        // do the tasks in sequence
    tasker.do_batch_tasks(xmltask);
    }
 catch(const std::exception& msg){
    cout << "Error: "<<msg.what()<<endl;
    return 1;}

 return 0;
}


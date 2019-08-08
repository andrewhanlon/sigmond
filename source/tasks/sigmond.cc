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
// *                                                                            *
// *       <Initialize>                                                         *
// *         <ProjectName>NameOfProject</ProjectName>                           * 
// *         <Logfile>output.log</Logfile>                                      *
// *         <KnownEnsemblesFile>ensembles.xml</KnownEnsemblesFile> (optional)  *
// *         <EchoXML/>                                                         *
// *         <MCBinsInfo>  ...  </MCBinsInfo>                                   *
// *         <MCSamplingInfo> ... </MCSamplingInfo>                             *
// *         <MCObservables>  ...  </MCObservables>                             *
// *       </Initialize>                                                        *
// *                                                                            *
// *       <TaskSequence>                                                       *
// *         <Task><Action>...</Action> ...  </Task>                            *
// *         <Task><Action>...</Action> ...  </Task>                            *
// *           ....                                                             *
// *       </TaskSequence>                                                      *
// *                                                                            *
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
// *   (d) Various ensembles are made known to SigMonD in the ensembles XML     *
// *       file.  The absolute path to this file can be specified in            *
// *       the <KnownEnsemblesFile> tag.  If not given, a default location      *
// *       for this file has been stored during the compilation.                *
// *       This file must have information specified in the following XML       *
// *       format:                                                              *
// *                                                                            *
// *      <KnownEnsembles>                                                      *
// *        <Ensemble>...</Ensemble>                                            *
// *        <Ensemble>...</Ensemble>                                            *
// *           ....                                                             *
// *      </KnownEnsembles>                                                     *
// *                                                                            *
// *       with each ensemble specified by                                      *
// *                                                                            *
// *      <Ensemble>                                                            *
// *         <Id>clover_s24_t128_ud840_s743</Id>                                *
// *         <NStreams>4</NStreams>                                             *
// *         <NMeas>551</NMeas>                                                 *
// *         <NSpace>24</NSpace>                                                *
// *         <NTime>128</NTime>                                                 *
// *         <Weights> 0.999 0.998 ... </Weights> (optional)                    *
// *      </Ensemble>                                                           *
// *                                                                            *
// *   (e) The tag <MCBinsInfo> is mandatory: it specifies the ensemble,        *
// *       controls rebinning the data, and possibly omitting certain           *
// *       configurations in the ensemble.  The XML must have the form below:   *
// *                                                                            *
// *      <MCBinsInfo>                                                          *
// *        <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo>         *
// *        <TweakEnsemble>  (optional)                                         *
// *           <Rebin>2</Rebin>                                                 *
// *           <Omissions>2 7 11</Omissions>                                    *
// *        </TweakEnsemble>                                                    *
// *      </MCBinsInfo>                                                         *
// *                                                                            *
// *       Note that when reading from bin files (other than basic LapH files), *
// *       the omissions in the bin files MUST be the same as specified         *
// *       in <MCBinsInfo>.  The rebin value need NOT be the same.  The         *
// *       <Rebin> value must be an integer multiple of the rebin factors       *
// *       in the bin files.                                                    *
// *                                                                            *
// *   (f) The tag <MCSamplingInfo> is mandatory.  It controls the default      *
// *       resampling method:  jackknife or bootstrap.  This default method     *
// *       is assumed for all reading and writing sampling results to and       *
// *       from files.  Note that both jackknife and bootstrap resampling       *
// *       can be done in any program execution, but only one can be used       *
// *       for reading/writing to files.  This tag has the form below.  See     *
// *       comments for the MCSamplingInfo and Bootstrapper classes for more    *
// *       details about this tag.                                              *
// *                                                                            *
// *      <MCSamplingInfo>                                                      *
// *         <Jackknife/>                                                       *
// *      </MCSamplingInfo>                                                     *
// *                       OR                                                   *
// *      <MCSamplingInfo>                                                      *
// *         <Bootstrapper>                                                     *
// *            <NumberResamplings>2048</NumberResamplings>                     *
// *            <Seed>6754</Seed>                                               *
// *            <BootSkip>127</BootSkip>                                        *
// *            <Precompute/>  (optional)                                       *
// *         </Bootstrapper>                                                    *
// *      </MCSamplingInfo>                                                     *
// *                                                                            *
// *   (g) <MCObservables> describes the data to be input for analysis. See     *
// *       class "MCObsGetHandler" in "source/data_handling/obs_get_handler.h"  *
// *       for a description of the XML needed in this tag.  This handles       *
// *       input of only "standard" observables (see "mcobs_info.h").           *
// *       Only data for standard observables can be read through this tag.     *
// *       Data of "nonstandard" form, such as fit parameters, rotated          *
// *       correlators, and other user-defined observables, must be read        *
// *       from file in a <Task> tag.                                           *
// *                                                                            *
// *   (h) The <Task> tags are needed in "batch" mode, but can be omitted in    *
// *   "cli" or "gui".  Each <Task> tag must begin with an <Action> tag.        *
// *   The <Action> tag must be a string in the "m_task_map".  The remaining    *
// *   XML depends on the action being taken.                                   *
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


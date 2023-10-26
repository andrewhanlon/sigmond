#ifndef TASK_HANDLER_H
#define TASK_HANDLER_H

#include "xml_handler.h"
#include "scalar_defs.h"
#include "mcobs_handler.h"
#include "mcobs_info.h"
#include "obs_get_handler.h"
#include <map>
#include <iostream>
#include <fstream>
//#include "user_interface.h"
#include "minimizer.h"

class TaskHandlerData;  // base class for persistent data

// ******************************************************************************
// *                                                                            *
// *   "TaskHandler" is the workhorse class of SigMonD.  It manages the tasks   *
// *   to do, and maintains pointers to objects that get the data and analyze   *
// *   the data.  The constructor sets up the known tasks and creates the       *
// *   data handling objects.  An object of class "MCObsGetHandler" is created  *
// *   and is accessed through the pointer "m_getter"; this object is           *
// *   responsible for reading data from files.  An object of the key class     *
// *   "MCObsHandler" is created and is accessed through the pointer "m_obs";   *
// *   this object stores the Monte Carlo data and resampling estimates,        *
// *   and carries out the statistical analysis on these observables.           *
// *   Analysis tasks are performed using "do_batch_tasks" (in "batch" mode) or *
// *   "do_interactive_tasks" in interactive modes ("cli" or "gui").  Output    *
// *   is logged in XML format in a log file.                                   *
// *                                                                            *
// *   An object of this class keeps track of persistent data (data needed by   *
// *   several tasks) using a map "m_task_data_map" which associates an ID      *
// *   string with a pointer to the data.  The map uses pointers to a base      *
// *   class "TaskHandlerData", so any data needed that should be persistent    *
// *   must be defined using a class derived from "TaskHandlerData".  Any task  *
// *   creating such data is responsible for allocating the memory using "new", *
// *   but the destructor of this class will call "delete" to wipe it from      *
// *   memory at the end of program execution.  The base pointer returned by    *
// *   member function "get_task_data" can be dynamic_cast<Derived*> back to    *
// *   the derived class in subsequent tasks.                                   *
// *                                                                            *
// *   The "TaskHandler" constructor requires XML of the form                   *
// *                                                                            *
// *    <SigMonD>                                                               *
// *                                                                            *
// *       <Initialize>                                                         *
// *         <ProjectName>NameOfProject</ProjectName>                           * 
// *         <LogFile>output.log</LogFile>                                      *
// *         <KnownEnsemblesFile>/path/ensembles.xml</KnownEnsemblesFile> (optional)  *
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
// *   (b) If <LogFile> is missing, a default name for the log file is used.    *
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
// *        <Infos>                                                             *
// *          <EnsembleInfo>...</EnsembleInfo>                                  *
// *          <EnsembleInfo>...</EnsembleInfo>                                  *
// *           ....                                                             *
// *        </Infos>                                                            *
// *        <CLSEnsembleWeights>                                                *
// *          <Ensemble>...</Ensemble>                                          *
// *           ....                                                             *
// *        </CLSEnsembleWeights>                                               *
// *      </KnownEnsembles>                                                     *
// *                                                                            *
// *       with each ensemble in the <Infos> tags specified by                  *
// *                                                                            *
// *      <EnsembleInfo>                                                        *
// *         <Id>clover_s24_t128_ud840_s743</Id>                                *
// *         <NStreams>4</NStreams>                                             *
// *         <NMeas>551</NMeas>                                                 *
// *         <NSpace>24</NSpace>                                                *
// *         <NTime>128</NTime>                                                 *
// *         <Weighted/>  (if has CLS weights; omit otherwise)                  *
// *      </EnsembleInfo>                                                       *
// *                                                                            *
// *       The entries in the <CLSEnsembleWeights> tag must have the form:      *
// *                                                                            *
// *      <Ensemble>                                                            *
// *         <Id>cls21_D200_r000</Id>                                           *
// *         <Weights> 0.999 0.998 ... </Weights>                               *
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
// *       the class "MCObsGetHandler" in "source/laph_data/obs_get_handler.h"  *
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


class TaskHandler
{

   MCBinsInfo *m_bins_info;
   MCSamplingInfo *m_samp_info;
   MCObsGetHandler *m_getter;
   MCObsHandler *m_obs;
   //UserInterface *m_ui;
   std::ofstream clog;

   typedef void (TaskHandler::*task_ptr)(XMLHandler&, XMLHandler&, int);
   std::map<std::string, task_ptr>    m_task_map;
   std::map<std::string, TaskHandlerData*> m_task_data_map;

       // Prevent copying ... handler might contain large
       // amounts of data

#ifndef NO_CXX11
   TaskHandler() = delete;
   TaskHandler(const TaskHandler&) = delete;
   TaskHandler& operator=(const TaskHandler&) = delete;
#else
   TaskHandler();
   TaskHandler(const TaskHandler&);
   TaskHandler& operator=(const TaskHandler&);
#endif

 public:
    
   TaskHandler(XMLHandler& xmlin);
   ~TaskHandler();

   void do_batch_tasks(XMLHandler& xmlin);

   void clear_task_data();
   void erase_task_data(const std::string& tdname);
   void insert_task_data(const std::string& tdname, TaskHandlerData* tdata);
   TaskHandlerData* get_task_data(const std::string& tdname);
   MCObsHandler* getMCObsHandler() {return m_obs;}

 private:

   void do_task(XMLHandler& xml_in, XMLHandler& output, int taskcount);

   std::string get_date_time();

   void finish_log();


       // The important task subroutines

   void clearMemory(XMLHandler &xml_in, XMLHandler& output, int taskcount);
   void clearSamplings(XMLHandler &xml_in, XMLHandler& output, int taskcount);
   void eraseData(XMLHandler &xml_in, XMLHandler& output, int taskcount);
   void eraseSamplings(XMLHandler &xml_in, XMLHandler& output, int taskcount);
   void readFromFile(XMLHandler &xml_in, XMLHandler& output, int taskcount);
   void writeToFile(XMLHandler &xml_in, XMLHandler& output, int taskcount);
   void writeCorrMatToFile(XMLHandler &xmltask, XMLHandler& xmlout, int taskcount);

   void printXML(XMLHandler &xml_in, XMLHandler& output, int taskcount);
   void doPlot(XMLHandler &xml_in, XMLHandler& output, int taskcount);
   void doFit(XMLHandler &xml_in, XMLHandler& output, int taskcount);
   void doChecks(XMLHandler& xml_in, XMLHandler& output, int taskcount);
   void doObsFunction(XMLHandler& xml_in, XMLHandler& output, int taskcount);
   void doCorrMatrixRotation(XMLHandler& xml_in, XMLHandler& output, int taskcount);
   void doRotCorrMatrixInsertFitInfos(XMLHandler& xml_in, XMLHandler& output, int taskcount);
   void doRotCorrMatrixRelabelEnergyPlots(XMLHandler& xml_in, XMLHandler& xml_out, int taskcount);
   void doCorrMatrixZMagSquares(XMLHandler& xml_in, XMLHandler& output, int taskcount);

   void getFromPivot(XMLHandler& xml_in, XMLHandler& output, int taskcount);

       // Utility subroutines

   uint getLatticeTimeExtent() const;
   uint getLatticeXExtent() const;
   uint getLatticeYExtent() const;
   uint getLatticeZExtent() const;

};

// ***************************************************************

   // base class for persistent data

class TaskHandlerData 
{
 public:
   TaskHandlerData(){}
   virtual ~TaskHandlerData() {}
};

// ***************************************************************

//    Utility routine for getting user input about writing to file.
//      <WriteToFile> 
//         <FileName>name</FileName>
//         <FileType>bins</FileType> (or samplings)
//         <FileFormat>fstr</FileFormat> (or hdf5: default if absent)
//         <WriteMode>overwrite</WriteMode> (protect, update, overwrite)
//      </WriteToFile>
//    Returns false if no <WriteToFile> tag in xmlin.  If not file name,
//    throws an exception.  Default <WriteMode> is "protect".
//    "ftype" is output as 'N' (not specified), 'B' (bins), or 'S' (samplings).
//    "ffmt" is output as 'D' (default), 'F' (fstr), or 'H' (hdf5).
  
bool getWriteToFileInfo(XMLHandler& xml_in, std::string& fileName,
                        WriteMode& wmode, char& ftype, char& ffmt, 
                        XMLHandler& echo);
                        
// ****************************************************************
#endif  

#include "task_handler.h"
#include "stopwatch.h"
#include "correlator_matrix_info.h"
using namespace std;
using namespace LaphEnv;

// *************************************************************************

     // set up the known tasks, create logfile, open stream for logging

TaskHandler::TaskHandler(XMLHandler& xmlin)
                       : m_bins_info(0), m_samp_info(0), m_getter(0), m_obs(0)
{
 if (xmlin.get_node_name()!="SigMonD")
    throw(std::invalid_argument("Input file must have root tag <SigMonD>"));
 string nowstr=get_date_time();

 if (xmlin.count_among_children("Initialize")!=1)
    throw(std::invalid_argument("There must be one child <Initialize> tag"));
 XMLHandler xmli(xmlin,"Initialize");

 if (xmli.count_among_children("KnownEnsemblesFile")==1){
    string knownEnsFile;
    xmlread(xmli,"KnownEnsemblesFile",knownEnsFile,"TaskHandler");
    knownEnsFile=tidyString(knownEnsFile);
    if (!knownEnsFile.empty())
       MCEnsembleInfo::m_known_ensembles_filename=knownEnsFile;}

 if (xmli.count_among_children("MCBinsInfo")!=1)
    throw(std::invalid_argument("There must be one <MCBinsInfo> tag"));
 try{
    XMLHandler xmlb(xmli,"MCBinsInfo");
    m_bins_info=new MCBinsInfo(xmlb);}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument("Failure reading Bin information"));}

 if (xmli.count_among_children("MCSamplingInfo")!=1)
    throw(std::invalid_argument("There must be one <MCSamplingInfo> tag"));
 bool boot_precompute=false;
 try{
    XMLHandler xmls(xmli,"MCSamplingInfo");
    m_samp_info=new MCSamplingInfo(xmls);
    if (m_samp_info->isBootstrapMode()){
       if (xmls.count_among_children("Precompute")==1) boot_precompute=true;}}
 catch(const std::exception& errmsg){
    throw(std::invalid_argument("Failure reading sampling information"));}

 string logfile;
 if (xmli.count_among_children("LogFile")==1)
    xmlread(xmli,"LogFile",logfile,"TaskHandler");
 else
    logfile=string("sigmond_log_")+nowstr+".xml";

 clog.open(logfile.c_str());
 if (!clog.is_open()){
    cout << "Could not open log file "<<logfile<<" for output"<<endl;
    throw(std::invalid_argument("Could not write to log file"));}
 clog << "<LogSigMonD>"<<endl;

 string projname;
 if (xmli.count_among_children("ProjectName")==1){
    xmlread(xmli,"ProjectName",projname,"TaskHandler");}
 else
    projname=string("SigMonD Project ")+nowstr;
 clog << " <ProjectName>"<<projname<<"</ProjectName>"<<endl;
 clog << " <StartDateTime>"<<nowstr<<"</StartDateTime>"<<endl;

 if (xmli.count_among_children("EchoXML")>=1){
    string input(xmlin.output());
    int pos=input.find("<SigMonD>");
    input.erase(0,pos+9);
    pos=input.find("</SigMonD>");
    input.erase(pos,string::npos);
    clog << " <InputXML>";
    clog << input<<"</InputXML>"<<endl;}

 try{
    XMLHandler xmlr(xmli,"MCObservables");
    m_getter=new MCObsGetHandler(xmlr,*m_bins_info,*m_samp_info);}
 catch(const std::exception& errmsg){
    clog << endl<<"<ERROR>"<<errmsg.what()<<"</ERROR>"<<endl<<endl; 
    finish_log(); 
    throw(std::invalid_argument("Failure reading MCObservables and/or data"));}

 try{
    m_obs=new MCObsHandler(*m_getter,boot_precompute);}
 catch(const std::exception& errmsg){
    delete m_getter; m_getter=0;
    clog << endl<<"<ERROR>"<<errmsg.what()<<"</ERROR>"<<endl<<endl;
    finish_log(); 
    throw(std::invalid_argument("Bad MCObsHandler construction"));}

 m_task_map["ClearMemory"]=&TaskHandler::clearMemory;
 m_task_map["ClearSamplings"]=&TaskHandler::clearSamplings;
 m_task_map["EraseData"]=&TaskHandler::eraseData;
 m_task_map["EraseSamplings"]=&TaskHandler::eraseSamplings;
 m_task_map["ReadFromFile"]=&TaskHandler::readFromFile;
 m_task_map["WriteToFile"]=&TaskHandler::writeToFile;
 m_task_map["WriteCorrMatToFile"]=&TaskHandler::writeCorrMatToFile;

 m_task_map["PrintXML"]=&TaskHandler::printXML;
 m_task_map["DoPlot"]=&TaskHandler::doPlot;
 m_task_map["DoFit"]=&TaskHandler::doFit;
 m_task_map["DoChecks"]=&TaskHandler::doChecks;
 m_task_map["DoObsFunction"]=&TaskHandler::doObsFunction;
 m_task_map["DoCorrMatrixRotation"]=&TaskHandler::doCorrMatrixRotation;
 m_task_map["DoRotCorrMatInsertFitInfos"]=&TaskHandler::doRotCorrMatrixInsertFitInfos;
 m_task_map["DoCorrMatrixRelabelEnergyPlots"]=&TaskHandler::doRotCorrMatrixRelabelEnergyPlots;
 m_task_map["DoCorrMatrixZMagSquares"]=&TaskHandler::doCorrMatrixZMagSquares;

 m_task_map["GetFromPivot"]=&TaskHandler::getFromPivot;

// m_ui=new UserInterface;
}


     // delete data, finish up logging, close log file

TaskHandler::~TaskHandler()
{
 finish_log();
 delete m_bins_info;
 delete m_samp_info;
 delete m_obs;
 delete m_getter;
// delete m_ui;
 m_task_map.clear();
 clear_task_data();
}


void TaskHandler::finish_log()
{
 string nowstr=get_date_time();
 clog << " <FinishDateTime>"<<nowstr<<"</FinishDateTime>"<<endl;
 clog << "</LogSigMonD>"<<endl<<endl;
 clog.close();
}


void TaskHandler::do_batch_tasks(XMLHandler& xmlin)
{
 XMLHandler xmlt(xmlin,"TaskSequence");
 list<XMLHandler> taskxml=xmlt.find_among_children("Task");
 int count=0;
 clog << endl<<"<BeginTasks>****************************************</BeginTasks>"<<endl;
 for (list<XMLHandler>::iterator it=taskxml.begin();it!=taskxml.end();it++,count++){
    clog << endl<<"<Task>"<<endl;
    clog << " <Count>"<<count<<"</Count>"<<endl;
    XMLHandler xmlout;
    StopWatch rolex; rolex.start();
    do_task(*it,xmlout,count);
    rolex.stop();
    clog << xmlout.output()<<endl;
    clog << "<RunTimeInSeconds>"<<rolex.getTimeInSeconds()<<"</RunTimeInSeconds>"<<endl;
    clog << "</Task>"<<endl;}
}


void TaskHandler::do_task(XMLHandler& xml_task, XMLHandler& xml_out, int count)
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
       (this->*(taskit->second))(xml_task,xml_out,count);}  // do the task!!
    else{
       throw(std::invalid_argument(string("Unknown task name: ")
               +task_action));}}   // unknown task?
 catch(const std::exception& errmsg){
    if (xml_out.empty())
       xml_out.set_root("Error",string(errmsg.what()));
    else
       xml_out.put_child("Error",string(errmsg.what()));}
}


string TaskHandler::get_date_time()
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


void TaskHandler::clear_task_data()
{
 for (map<string,TaskHandlerData*>::iterator it
      =m_task_data_map.begin();it!=m_task_data_map.end();it++)
    delete it->second;
 m_task_data_map.clear();
}


void TaskHandler::erase_task_data(const string& tdname)
{
 string taskname=tidyName(tdname);
 map<string,TaskHandlerData*>::iterator it=m_task_data_map.find(taskname);
 if (it!=m_task_data_map.end()){
    delete it->second;
    m_task_data_map.erase(it);}
}


    //  TaskHandlerData should be created with "new" so is persistent

void TaskHandler::insert_task_data(const std::string& tdname, TaskHandlerData* tdata)
{
 string taskname=tidyName(tdname);
 if (taskname.empty()){
    throw(std::invalid_argument("Cannot insert task data since name is invalid"));}
 map<string,TaskHandlerData*>::iterator it=m_task_data_map.find(taskname);
 if (it!=m_task_data_map.end()){
    throw(std::invalid_argument("Cannot insert task data since name already in map"));}
 m_task_data_map.insert(make_pair(taskname,tdata));
}


   //  returns a base pointer, which can be dynamic_cast(...)
   //  Return null value if error.

TaskHandlerData* TaskHandler::get_task_data(const string& tdname)
{
 string taskname=tidyName(tdname);
 if (taskname.empty()) return 0;
 map<string,TaskHandlerData*>::iterator it=m_task_data_map.find(taskname);
 if (it==m_task_data_map.end()) return 0;
 return it->second;
}

  // *****************************************************************

   //   <Task>
   //     <Action>ClearMemory</Action>
   //   </Task>

void TaskHandler::clearMemory(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 m_obs->clearData();   // also clears all Samplings
 xmlout.set_root("ClearMemory","done");
}

   //   <Task>
   //     <Action>ClearSamplings</Action>
   //   </Task>

void TaskHandler::clearSamplings(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 m_obs->clearSamplings();   
 xmlout.set_root("ClearSamplings","done");
}

   //   <Task>
   //     <Action>EraseData</Action>
   //      <MCObservable>...</MCObservable> 
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::eraseData(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 xmlout.set_root("EraseData");
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       MCObsInfo obskey(*tt);
       m_obs->eraseData(obskey);   // also clears all Samplings
       XMLHandler xmlb; obskey.output(xmlb);
       xmlout.put_child(xmlb);}
}

   //   <Task>
   //     <Action>EraseSamplings</Action>
   //      <MCObservable>...</MCObservable> 
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::eraseSamplings(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 xmlout.set_root("EraseSamplings");
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       MCObsInfo obskey(*tt);
       m_obs->eraseSamplings(obskey);   // also clears all Samplings
       XMLHandler xmlb; obskey.output(xmlb);
       xmlout.put_child(xmlb);}
}


   //   <Task>
   //     <Action>ReadFromFile</Action>
   //      <FileType>bins</FileType> (or samplings)
   //      <FileName>name_of_file</FileName>
   //      <MCObservable>...</MCObservable>   (these are optional)
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::readFromFile(XMLHandler &xmltask, XMLHandler& xmlout, int taskcount)
{
 xmlout.set_root("ReadFromFile");
 string type;
 xmlreadchild(xmltask,"FileType",type,"ReadFromFile");
 if ((type!="bins")&&(type!="samplings"))
    throw(std::invalid_argument("<FileType> must be bins or samplings in ReadFromFile"));
 xmlout.put_child("FileType",type);
 string filename;
 xmlreadchild(xmltask,"FileName",filename,"TaskHandler");
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 set<MCObsInfo> obskeys;
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       obskeys.insert(MCObsInfo(*tt));}
 XMLHandler xmlf;
 if (type=="samplings"){
    if (obskeys.empty())
       m_obs->readSamplingValuesFromFile(filename,xmlf);
    else
       m_obs->readSamplingValuesFromFile(obskeys,filename,xmlf);}
 else{
    if (obskeys.empty())
       m_obs->readBinsFromFile(filename,xmlf);
    else
       m_obs->readBinsFromFile(obskeys,filename,xmlf);}
 xmlout.put_child(xmlf);
}

   //   <Task>
   //     <Action>WriteToFile</Action>   (uses samplings mode in <MCSamplingInfo> tag)
   //      <FileType>bins</FileType>  (or samplings)
   //      <FileName>name_of_file</FileName>
   //      <WriteMode>overwrite</WriteMode>   (optional: protect, or update, overwrite) 
   //      <MCObservable>...</MCObservable>   (these are needed)
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::writeToFile(XMLHandler &xmltask, XMLHandler& xmlout, int taskcount)
{
 xmlout.set_root("WriteToFile");
 string type;
 xmlreadchild(xmltask,"FileType",type,"WriteToFile");
 if ((type!="bins")&&(type!="samplings"))
    throw(std::invalid_argument("<FileType> must be bins or samplings in WriteToFile"));
 xmlout.put_child("FileType",type);
 string filename;
 xmlreadchild(xmltask,"FileName",filename,"WriteToFile");
 WriteMode wmode = Protect;  // protect mode
 if (xml_tag_count(xmltask,"WriteMode")==1){
    string fmode;
    xmlread(xmltask,"WriteMode",fmode,"WriteToFile");
    fmode=tidyString(fmode);
    if (fmode=="overwrite") wmode=Overwrite;
    else if (fmode=="update") wmode=Update;}
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 set<MCObsInfo> obskeys;
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       obskeys.insert(MCObsInfo(*tt));}
 XMLHandler xmlf;
 if (type=="bins")
    m_obs->writeBinsToFile(obskeys,filename,xmlf,wmode);
 else
    m_obs->writeSamplingValuesToFile(obskeys,filename,xmlf,wmode);
 xmlout.put_child(xmlf);
}

   //   <Task>
   //     <Action>WriteCorrMatToFile</Action>
   //       <FileType>bins</FileType> (or samplings)
   //      <FileName>name_of_file</FileName>
   //      <WriteMode>overwrite</WriteMode>   (optional: default protect, update, overwrite) 
   //   <CorrelatorMatrixInfo>
   //     <BLOperatorString>....</BLOperatorString>
   //      ....
   //     <HermitianMatrix/>
   //     <SubtractVEV/>
   //   </CorrelatorMatrixInfo>
   //   <MinTimeSep>3</MinTimeSep>
   //   <MaxTimeSep>25</MaxTimeSep>
   //   <SeparateVEVWrite/> (if type==samplings, the VEVs are NOT written to file by default;
   //   </Task>                include this tag if you still want the VEVs separately written)


void TaskHandler::writeCorrMatToFile(XMLHandler &xmltask, XMLHandler& xmlout, int taskcount)
{
 xmlout.set_root("WriteCorrMatToFile");
 string type;
 xmlreadchild(xmltask,"FileType",type,"WriteCorrMatToFile");
 if ((type!="bins")&&(type!="samplings"))
    throw(std::invalid_argument("<FileType> must be bins or samplings in WriteCorrMatToFile"));
 xmlout.put_child("FileType",type);
 string filename;
 xmlreadchild(xmltask,"FileName",filename,"WriteCorrMatToFile");
 xmlout.put_child("FileName",filename);
 WriteMode wmode = Protect;  // protect mode
 if (xml_tag_count(xmltask,"WriteMode")==1){
    string fmode;
    xmlread(xmltask,"WriteMode",fmode,"TaskHandler");
    fmode=tidyString(fmode);
    if (fmode=="overwrite") wmode=Overwrite;
    else if (fmode=="update") wmode=Update;}
 uint tmin,tmax;
 xmlreadchild(xmltask,"MinTimeSep",tmin,"WriteCorrMatToFile");
 xmlreadchild(xmltask,"MaxTimeSep",tmax,"WriteCorrMatToFile");
 CorrelatorMatrixInfo cormat(xmltask);
 bool herm=cormat.isHermitian();
 bool subvev=cormat.subtractVEV();
 bool vsep=subvev;
 if (subvev && type=="samplings"){
    if (xmltask.count("SeparateVEVWrite")>0) vsep=true;
    else vsep=false;}
 if (vsep) subvev=false;
 const set<OperatorInfo>& ops=cormat.getOperators();
 set<MCObsInfo> obskeys;
 for (set<OperatorInfo>::const_iterator it1=ops.begin();it1!=ops.end();++it1)
  for (set<OperatorInfo>::const_iterator it2=(herm)?it1:ops.begin();it2!=ops.end();++it2)
    for (uint tval=tmin;tval<=tmax;++tval){
       MCObsInfo key(*it1,*it2,tval,herm,RealPart,subvev);
       obskeys.insert(key);
#ifdef COMPLEXNUMBERS
       key.setToImaginaryPart();
       obskeys.insert(key);
#endif
       }
 if (vsep){
    for (set<OperatorInfo>::const_iterator it1=ops.begin();it1!=ops.end();++it1){
       MCObsInfo key(*it1,RealPart);
       obskeys.insert(key);
#ifdef COMPLEXNUMBERS
       key.setToImaginaryPart();
       obskeys.insert(key);
#endif
       }}
 XMLHandler xmlf;
 if (type=="bins")
    m_obs->writeBinsToFile(obskeys,filename,xmlf,wmode);
 else
    m_obs->writeSamplingValuesToFile(obskeys,filename,xmlf,wmode);
 xmlout.put_child(xmlf);
}


  // *****************************************************************

       // Utility subroutines


uint TaskHandler::getLatticeTimeExtent() const
{
 return m_obs->getLatticeTimeExtent();
}

uint TaskHandler::getLatticeXExtent() const
{
 return m_obs->getLatticeXExtent();
}

uint TaskHandler::getLatticeYExtent() const
{
 return m_obs->getLatticeYExtent();
}

uint TaskHandler::getLatticeZExtent() const
{
 return m_obs->getLatticeZExtent();
}


// ***************************************************************************************

//    Utility routine for getting user input about writing to file.
//      <WriteToFile> 
//         <FileName>name</FileName>
//         <FileType>bins</FileType> (or samplings)
//         <WriteMode>overwrite</WriteMode> (protect, update, overwrite)
//      </WriteToFile>
//    Returns false if no <WriteToFile> tag in xmlin.  If not file name,
//    throws an exception.  Default <WriteMode> is "protect".
//    "ftype" is output as 'N' (not specified), 'B' (bins), or 'S' (samplings).
  
bool getWriteToFileInfo(XMLHandler& xml_in, std::string& fileName,
                        WriteMode& wmode, char& ftype, XMLHandler& echo)
{
 wmode = Protect;  // protect mode
 ftype = 'N';
 fileName.clear();
 ArgsHandler xmlin(xml_in);
 if (xmlin.queryTag("WriteToFile")){
    ArgsHandler xmlf(xmlin,"WriteToFile");
    string filename(xmlf.getString("FileName"));
    if (filename.empty()) throw(std::invalid_argument("Empty file name so cannot WriteToFile"));
    string fmode("protect");
    xmlf.getOptionalString("WriteMode",fmode); 
    if (fmode=="overwrite") wmode=Overwrite;
    else if (fmode=="update") wmode=Update;
    string ftypestr;
    xmlf.getOptionalString("FileType",ftypestr);
    if (ftypestr=="bins") ftype='B';
    else if (ftypestr=="samplings") ftype='S';
    else throw(std::invalid_argument("Invalid file type for WriteToFile"));
    xmlf.echo(echo);
    return true;}
 return false;
}


// ***************************************************************************************

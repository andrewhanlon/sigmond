#include "task_handler.h"
#include "stopwatch.h"
using namespace std;
using namespace LaphEnv;

// *************************************************************************

     // set up the known tasks, create logfile, open stream for logging

TaskHandler::TaskHandler(XMLHandler& xmlin)
                       : m_getter(0), m_obs(0)
{
 if (xmlin.get_node_name()!="SigMonD")
    throw(std::invalid_argument("Input file must have root tag <SigMonD>"));
 string nowstr=get_date_time();

 if (xmlin.count_among_children("Initialize")!=1)
    throw(std::invalid_argument("There must be one child <Initialize> tag"));

 XMLHandler xmli(xmlin,"Initialize");
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
    m_getter=new MCObsGetHandler(xmlr);}
 catch(const std::exception& errmsg){
    clog << endl<<"<ERROR>"<<errmsg.what()<<"</ERROR>"<<endl<<endl; 
    finish_log(); 
    throw(std::invalid_argument("Failure reading MCObservables and/or data"));}
 try{
    m_obs=new MCObsHandler(*m_getter);}
 catch(const std::exception& errmsg){
    delete m_getter; m_getter=0;
    clog << endl<<"<ERROR>"<<errmsg.what()<<"</ERROR>"<<endl<<endl;
    finish_log(); 
    throw(std::invalid_argument("Bad MCObsHandler construction"));}

 if (xmli.count_among_children("TweakEnsemble")==1){
    XMLHandler xmlk(xmli,"TweakEnsemble");
    int rebin;
    if (xmlreadifchild(xmlk,"Rebin",rebin)){
       if (rebin>1) m_obs->setRebin(rebin);}
    vector<int> ovec;
    if (xmlreadifchild(xmlk,"Omissions",ovec)){
       set<int> omissions(ovec.begin(),ovec.end());
       if (!omissions.empty()) m_obs->addOmissions(omissions);}}

 if (xmli.count_among_children("Bootstrapper")==1){
    XMLHandler xmlb(xmli,"Bootstrapper");
    int num_resamplings=1024;
    unsigned long bootseed=0, bootskip=64;
    bool precompute=false;
    xmlreadifchild(xmlb,"NumberResamplings",num_resamplings);
    xmlreadifchild(xmlb,"Seed",bootseed);
    xmlreadifchild(xmlb,"BootSkip",bootskip);
    if (xmlb.count_among_children("Precompute")==1) precompute=true;
    m_obs->setBootstrapper(num_resamplings,bootseed,bootskip,precompute);}

 m_task_map["ClearMemory"]=&TaskHandler::clearMemory;
 m_task_map["ClearSamplings"]=&TaskHandler::clearSamplings;
 m_task_map["EraseData"]=&TaskHandler::eraseData;
 m_task_map["EraseSamplings"]=&TaskHandler::eraseSamplings;
 m_task_map["ReadSamplingsFromFile"]=&TaskHandler::readSamplingsFromFile;
 m_task_map["WriteSamplingsToFile"]=&TaskHandler::writeSamplingsToFile;
 m_task_map["ReadBinsFromFile"]=&TaskHandler::readBinsFromFile;
 m_task_map["WriteBinsToFile"]=&TaskHandler::writeBinsToFile;

 m_task_map["PrintXML"]=&TaskHandler::printXML;
 m_task_map["DoPlot"]=&TaskHandler::doPlot;
 m_task_map["DoFit"]=&TaskHandler::doFit;
 m_task_map["DoChecks"]=&TaskHandler::doChecks;
 m_task_map["DoObsFunction"]=&TaskHandler::doObsFunction;
 m_task_map["DoCorrMatrixRotation"]=&TaskHandler::doCorrMatrixRotation;
 m_task_map["DoCorrMatrixZMagSquares"]=&TaskHandler::doCorrMatrixZMagSquares;

 m_task_map["DoHamiltonian"]=&TaskHandler::doHamiltonian;

 m_ui=new UserInterface;
}


     // delete data, finish up logging, close log file

TaskHandler::~TaskHandler()
{
 finish_log();
 delete m_obs;
 delete m_getter;
 delete m_ui;
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
       throw(std::invalid_argument((string("Unknown task name: ")
               +task_action).c_str()));}}   // unknown task?
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
   //     <Action>ReadSamplingsFromFile</Action>
   //      <SamplingMode>Jackknife</SamplingMode>  (or Bootstrap or Current)
   //      <FileName>name_of_file</FileName>
   //      <MCObservable>...</MCObservable>   (these are optional)
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::readSamplingsFromFile(XMLHandler &xmltask, XMLHandler& xmlout, int taskcount)
{
 xmlout.set_root("ReadSamplingsFromFile");
 string smode;
 xmlreadchild(xmltask,"SamplingMode",smode,"TaskHandler");
 SamplingMode mode;
 if (smode=="Bootstrap") mode=Bootstrap;
 else if (smode=="Jackknife") mode=Jackknife;
 else if (smode=="Current") mode=m_obs->getCurrentSamplingMode();
 else throw(std::invalid_argument("Invalid Sampling Mode"));
 string filename;
 xmlreadchild(xmltask,"FileName",filename,"TaskHandler");
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 set<MCObsInfo> obskeys;
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       obskeys.insert(MCObsInfo(*tt));}
 XMLHandler xmlf;
 if (obskeys.empty())
    m_obs->readSamplingValuesFromFile(filename,mode,xmlf);
 else
    m_obs->readSamplingValuesFromFile(obskeys,filename,mode,xmlf);
 xmlout.put_child(xmlf);
}

   //   <Task>
   //     <Action>WriteSamplingsToFile</Action>
   //      <SamplingMode>Jackknife</SamplingMode>  (or Bootstrap or Current)
   //      <FileName>name_of_file</FileName>
   //      <FileMode>overwrite</FileMode>   (optional) 
   //      <MCObservable>...</MCObservable>   (these are needed)
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::writeSamplingsToFile(XMLHandler &xmltask, XMLHandler& xmlout, int taskcount)
{
 xmlout.set_root("WriteSamplingsToFile");
 string smode;
 xmlreadchild(xmltask,"SamplingMode",smode,"TaskHandler");
 SamplingMode mode;
 if (smode=="Bootstrap") mode=Bootstrap;
 else if (smode=="Jackknife") mode=Jackknife;
 else if (smode=="Current") mode=m_obs->getCurrentSamplingMode();
 else throw(std::invalid_argument("Invalid Sampling Mode"));
 string filename;
 xmlreadchild(xmltask,"FileName",filename,"TaskHandler");
 bool overwrite = false;  // protect mode
 if (xml_tag_count(xmltask,"FileMode")==1){
    string fmode;
    xmlread(xmltask,"FileMode",fmode,"FileListInfo");
    fmode=tidyString(fmode);
    if (fmode=="overwrite") overwrite=true;}
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 set<MCObsInfo> obskeys;
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       obskeys.insert(MCObsInfo(*tt));}
 XMLHandler xmlf;
 m_obs->writeSamplingValuesToFile(obskeys,filename,mode,xmlf,overwrite);
 xmlout.put_child(xmlf);
}



   //   <Task>
   //     <Action>ReadBinsFromFile</Action>
   //      <BinFileName>name_of_file</BinFileName>
   //      <MCObservable>...</MCObservable>   (these are optional)
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::readBinsFromFile(XMLHandler &xmltask, XMLHandler& xmlout, int taskcount)
{
 xmlout.set_root("ReadBinsFromFile");
 string filename;
 xmlreadchild(xmltask,"BinFileName",filename,"TaskHandler");
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 set<MCObsInfo> obskeys;
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       obskeys.insert(MCObsInfo(*tt));}
 XMLHandler xmlf;
 if (obskeys.empty())
    m_obs->readBinsFromFile(filename,xmlf);
 else
    m_obs->readBinsFromFile(obskeys,filename,xmlf);
 xmlout.put_child(xmlf);
}

   //   <Task>
   //     <Action>WriteBinsToFile</Action>
   //      <BinFileName>name_of_file</BinFileName>
   //      <FileMode>overwrite</FileMode>   (optional) 
   //      <MCObservable>...</MCObservable>   (these are needed)
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::writeBinsToFile(XMLHandler &xmltask, XMLHandler& xmlout, int taskcount)
{
 xmlout.set_root("WriteBinsToFile");
 string filename;
 xmlreadchild(xmltask,"BinFileName",filename,"TaskHandler");
 bool overwrite = false;  // protect mode
 if (xml_tag_count(xmltask,"FileMode")==1){
    string fmode;
    xmlread(xmltask,"FileMode",fmode,"FileListInfo");
    fmode=tidyString(fmode);
    if (fmode=="overwrite") overwrite=true;}
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 set<MCObsInfo> obskeys;
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       obskeys.insert(MCObsInfo(*tt));}
 XMLHandler xmlf;
 m_obs->writeBinsToFile(obskeys,filename,xmlf,overwrite);
 xmlout.put_child(xmlf);
}


  // *****************************************************************

       // Utility subroutines


uint TaskHandler::getLatticeTimeExtent() const
{
 return m_getter->getEnsembleInfo().getLatticeTimeExtent();
}

uint TaskHandler::getLatticeXExtent() const
{
 return m_getter->getEnsembleInfo().getLatticeXExtent();
}

uint TaskHandler::getLatticeYExtent() const
{
 return m_getter->getEnsembleInfo().getLatticeYExtent();
}

uint TaskHandler::getLatticeZExtent() const
{
 return m_getter->getEnsembleInfo().getLatticeZExtent();
}



// ***************************************************************************************
 

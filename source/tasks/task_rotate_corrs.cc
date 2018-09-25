#include "task_handler.h"
#include "task_utils.h"
#include "correlator_matrix_info.h"
#include "single_pivot.h"
#include "create_plots.h"
#include <tuple>

using namespace std;

// ***********************************************************************************
// *                                                                                 *
// *   Several tasks are defined in this file:                                       *
// *                                                                                 *
// *      (a)  <Action>DoCorrMatrixRotation</Action>                                 *
// *                                                                                 *
// *      (b)  <Action>DoRotCorrMatInsertFitInfos</Action>                           *
// *                                                                                 *
// *      (c)  <Action>DoCorrMatrixZMagSquares</Action>                              *
// *                                                                                 *
// *      (d)  <Action>DoCorrMatrixRelabelEnergyPlots</Action>                       *
// *                                                                                 *
// *   The task "DoCorrMatrixRotation" is done first, then fits to the diagonal      *
// *   element of the rotated correlation matrix are done.  From these fits,         *
// *   the level energies and amplitudes are obtained.  The task                     *
// *   "DoRotCorrMatInsertFitInfos" can be used to get insert the fit energies and   *
// *   fit amplitudes into memory, optionally reordering the level indices by        *
// *   ascending fit energy.  Then the amplitudes can be used finally to evaluate    *
// *   operator overlap factors in a "DoCorrMatrixZMagSquares" task.  The plots      *
// *   from the fits to the rotated diagonal correlators have an initial level       *
// *   ordering.  If you reorder the levels by ascending fit values, then the        *
// *   task "DoCorrMatrixRelabelEnergyPlots" can be used to modify the plots to      *
// *   agree with the new ordering.                                                  *
// *                                                                                 *
// *   To accomplish all of this, make sure to write the pivot to a file in the      *
// *   first task "DoCorrMatrixRotation" with a <WritePivotToFile> tag.  Then        *
// *   do all of the fits to the diagonal rotated correlators, and write these       *
// *   energy and amplitude fits into a samplings file(s).  Then do a                *
// *   "DoRotCorrMatInsertFitInfos" task, reading the pivot from file and inserting  *
// *   fit infos, specifying the sampling files containing the energies and          *
// *   amplitudes.  One can reordering the energies at this point.  Store the pivot  *
// *   in memory with an <AssignName> tag. In the final task                         *
// *   "DoCorrMatrixZMagSquares", initiate the pivot using a <GetFromMemory>.        *
// *                                                                                 *
// *                                                                                 *
// *   XML format for correlator matrix rotations:                                   *
// *                                                                                 *
// *     <Task>                                                                      *
// *        <Action>DoCorrMatrixRotation</Action>                                    *
// *        <MinTimeSep>3</MinTimeSep>                                               *
// *        <MaxTimeSep>30</MaxTimeSep>                                              *
// *        <RotateBy>bins</RotateBy>     (or samplings)                             *
// *        <Type>SinglePivot</Type>                                                 *
// *        <SinglePivotInitiate> ... </SinglePivotInitiate> (depends on type)       *
// *        <WriteRotatedCorrToFile>                                                 *
// *           <RotatedCorrFileName>rotated_corr_bins</RotatedCorrFileName>          *
// *           <Type>bins</Type>  (or samplings)                                     *
// *           <Overwrite/>                                                          *
// *        </WriteRotatedCorrToFile>                                                *
// *        <PlotRotatedEffectiveEnergies>   (optional)                              *
// *          <SamplingMode>Jackknife</SamplingMode> (or Bootstrap: optional)        *
// *          <EffEnergyType>TimeForward</EffEnergyType> (optional)                  *
// *              ( or TimeSymmetric, TimeSymmetricPlusConst, TimeForwardPlusConst)  *
// *          <TimeStep>3</TimeStep>  (default 1)                                    *
// *          <PlotFileStub>stubname</PlotFileStub> (produces several files with     *
// *                                            numerical suffices 0,1,...)          *
// *          <SymbolColor>blue</SymbolColor>                                        *
// *          <SymbolType>circle</SymbolType>                                        *
// *          <MaxErrorToPlot>1.0</MaxErrorToPlot>                                   *
// *        </PlotRotatedEffectiveEnergies>                                          *
// *        <PlotRotatedCorrelators>         (optional)                              *
// *          <SamplingMode>Jackknife</SamplingMode> (or Bootstrap: optional)        *
// *          <PlotFileStub>stubname</PlotFileStub> (produces several files with     *
// *                                            numerical suffices 0,1,...)          *
// *          <Arg>Re</Arg>  (or Im)                                                 *
// *          <SymbolColor>blue</SymbolColor>                                        *
// *          <SymbolType>circle</SymbolType>                                        *
// *          <Rescale>1.0</Rescale>          (optional)                             *
// *        </PlotRotatedEffectiveEnergies>                                          *
// *     </Task>                                                                     *
// *                                                                                 *
// *   The <Type> specifies which kind of pivot to create.  The tag                  *
// *   <SinglePivotInitiate> applies only for a type of SinglePivot.  For            *
// *   other types, there will be a different tag instead of <SinglePivotInitiate>.  *
// *                                                                                 *
// *                                                                                 *
// *   After fits to the diagonal elements of the rotated correlators are done       *
// *   to obtain the fit energies and the amplitudes, this information can be        *
// *   inserted into memory, optionally reordering the levels according to           *
// *   increasing fit energy.  The needed amplitudes can either be individually      *
// *   specified in the XML below, or if they are given the same name with the ID    *
// *   index equal to the level number, a short form is available.  The XML format   *
// *   for this task is                                                              *
// *                                                                                 *
// *     <Task>                                                                      *
// *       <Action>DoRotCorrMatInsertFitInfos</Action>                               *
// *        <Type>SinglePivot</Type>                                                 *
// *        <SinglePivotInitiate> ... </SinglePivotInitiate> (depends on type)       *
// *        <ReorderByFitEnergy/>  (optional)                                        *
// *        <EnergyFitCommonName>..</EnergyFitCommonName> (or individually specify)  *
// *        <EnergyFit>                                                              *
// *           <Level>0</Level>                                                      *
// *           <Name>A</Name><IDIndex>0</IDIndex> (name of fit energy observable)    *
// *        </EnergyFit>                                                             *  
// *            ... one for each level                                               *
// *                (specify amplitude using short form below...)                    *
// *        <RotatedAmplitudeCommonName>Amp</RotatedAmplitudeCommonName>             *
// *                (or specify each individually)                                   *
// *        <RotatedAmplitude>                                                       *
// *           <Level>0</Level>                                                      *
// *           <Name>A</Name><IDIndex>0</IDIndex>                                    *
// *        </RotatedAmplitude>                                                      *
// *            ...  (needed for all levels)                                         *
// *     </Task>                                                                     *
// *                                                                                 *
// *   Before computing the Zmag squares, the amplitudes from the rotated correlator *
// *   fits are needed.  These must be "inserted" into the pivot.  The needed        *
// *   amplitudes can either be individually specified in the XML below to compute   *
// *   the Zmag squares, or if they are given the same name with the ID index equal  *
// *   to the level number, a short form is available:                               *
// *                                                                                 *
// *   XML format for correlator matrix Zmag squares computation:                    *
// *                                                                                 *
// *     <Task>                                                                      *
// *        <Action>DoCorrMatrixZMagSquares</Action>                                 *
// *        <Type>SinglePivot</Type>                                                 *
// *        <SinglePivotInitiate> ... </SinglePivotInitiate> (depends on type)       *
// *           ...use <GetFromMemory> ...                                            *
// *        <DoPlots>  (optional)                                                    *
// *           <PlotFileStub> ... </PlotFileStub>                                    *
// *           <BarColor> ... </BarColor>           (optional: cyan default)         *
// *           <ZMagSqPlot>                                                          *
// *             <BLOperatorString>...</BLOperatorString>                            *
// *             <ObsName> ... </ObsName>                                            *
// *             <FileSuffix> ... </FileSuffix> (optional: default is index)         *
// *           </ZMagSqPlot>                                                         *
// *              ....                                                               *
// *        </DoPlots>                                                               *
// *     </Task>                                                                     *
// *                                                                                 *
// *         ... <ObsName>standard</ObsName> creates name for standard ops           *
// *   "DoPlots" produces an xmgrace file for each operator.                         *
// *                                                                                 *
// *                                                                                 *
// *   If you reorder the levels by ascending fit values, then the task              *
// *   "DoCorrMatrixRelabelEnergyPlots" can be used to modify the plots to           *
// *   agree with the new ordering.  The original plot files are specified either    *
// *   by a stub (stub_0.agr, stub_1.agr, ...) or the file names are specified       *
// *   in order.  The revised plot files are similarly specified.  If no revised     *
// *   plot files are given, the original files are overwritten.  The XML must       *
// *   have the form                                                                 *
// *                                                                                 *
// *     <Task>                                                                      *
// *        <Action>DoCorrMatrixRelabelEnergyPlots</Action>                          *
// *        <Type>SinglePivot</Type>                                                 *
// *        <SinglePivotInitiate> ... </SinglePivotInitiate> (depends on type)       *
// *           ...use <GetFromMemory> ...                                            *
// *        <OriginalPlotFiles>                                                      *
// *           <PlotFileStub> ... </PlotFileStub>                                    *
// *            or  <PlotFile>name0</PlotFile> ....                                  *
// *        </OriginalPlotFiles>                                                     *
// *        <RevisedPlotFiles>   (optional)                                          *
// *           <PlotFileStub> ... </PlotFileStub>                                    *
// *            or  <PlotFile>name0</PlotFile> ....                                  *
// *        </RevisedPlotFiles>                                                      *
// *     </Task>                                                                     *
// *                                                                                 *
// *                                                                                 *
// *   See the individual pivot header files for information about initiating        *
// *   the pivots.                                                                   *
// *                                                                                 *
// ***********************************************************************************




void TaskHandler::doCorrMatrixRotation(XMLHandler& xml_task, XMLHandler& xml_out, int taskcount)
{
 LogHelper xmlout;
 ArgsHandler xmltask(xml_task);
 xmlout.reset("DoCorrMatrixRotation");
 uint mintimesep,maxtimesep;
 xmltask.getUInt("MinTimeSep",mintimesep);
 xmltask.getUInt("MaxTimeSep",maxtimesep);

 string rotate_by(xmltask.getString("RotateBy"));
 if (rotate_by!="bins" && rotate_by!="samplings")
    throw(std::invalid_argument(string("Invalid entry in <RotateBy>")));
 string rotatetype(xmltask.getString("Type"));

 if (rotatetype=="SinglePivot"){
    ArgsHandler xmlpiv(xmltask,"SinglePivotInitiate");
    LogHelper xmllog;
    bool pkeep;
    SinglePivotOfCorrMat* pivoter=SinglePivotOfCorrMat::initiateSinglePivot(
                             *this,xmlpiv,xmllog,pkeep);
    xmlout.putItem(xmllog);
    if (pivoter==0){
       xmlout.output(xml_out);
       throw(std::runtime_error("Could not initiate Single Pivot"));}
    try{
    pivoter->doRotation(mintimesep,maxtimesep,rotate_by,xmllog);}
    catch(const std::exception& errmsg){
       xmlout.putItem(xmllog); xmlout.output(xml_out);
       throw(std::invalid_argument(string("Error in SinglePivotOfCorrMat::doRotation: ")
              +string(errmsg.what())));} 
    xmlout.putItem(xmllog);
       // save rotated correlators to file
    if (xmltask.queryTag("WriteRotatedCorrToFile")){
       try{
       ArgsHandler xmlf(xmltask,"WriteRotatedCorrToFile");
       string corrfile(xmlf.getString("RotatedCorrFileName"));
       bool overwrite=false;
       xmlf.getOptionalBool("Overwrite",overwrite); 
       if (corrfile.empty()) throw(std::invalid_argument("Empty file name"));
       string ftype("bins");
       xmlf.getOptionalString("Type",ftype);
       bool bins=(ftype=="samplings")?false:true;
       LogHelper xmlw;
       pivoter->writeRotated(mintimesep,maxtimesep,corrfile,overwrite,xmlw,bins);
       xmlout.putItem(xmlw);}
       catch(const std::exception& msg){
          xmlout.putString("Error",string(msg.what()));}}
    if (xmltask.queryTag("PlotRotatedEffectiveEnergies")){
       try{
       ArgsHandler xmlc(xmltask,"PlotRotatedEffectiveEnergies");
       LogHelper xmllog("PlotRotatedEffectiveEnergies");
       string instr("Jackknife");
       xmlc.getOptionalString("SamplingMode",instr);
       SamplingMode mode=Jackknife;
       if (instr=="Bootstrap") mode=Bootstrap;
       else if (instr=="Jackknife") mode=Jackknife;
       else throw(std::invalid_argument("Bad sampling mode"));
       instr="TimeForward";
       xmlc.getOptionalString("EffEnergyType",instr);
       uint efftype=0;
       if (instr=="TimeSymmetric") efftype=1;
       else if (instr=="TimeForward") efftype=0;
       else if (instr=="TimeSymmetricPlusConst") efftype=3;
       else if (instr=="TimeForwardPlusConst") efftype=2;
       else throw(std::invalid_argument("Bad effective energy type"));
       uint step=1;
       xmlc.getOptionalUInt("TimeStep",step);
       if ((step<1)||(step>getLatticeTimeExtent()/4))
          throw(std::invalid_argument("Bad effective energy time step"));
       string plotfilestub(xmlc.getString("PlotFileStub"));
       string color("blue"),symboltype("circle");
       xmlc.getOptionalString("SymbolColor",color);
       xmlc.getOptionalString("SymbolType",symboltype);
       double maxerror=0.0;
       xmlc.getOptionalReal("MaxErrorToPlot",maxerror);
       xmllog.putEcho(xmlc);
       uint nplots=pivoter->getNumberOfLevels();
       xmllog.putUInt("NumberOfPlots",nplots);
       bool herm=true;
       bool subvev=pivoter->subtractVEV();
       bool reweight=pivoter->reweight();
       GenIrrepOperatorInfo oprot(pivoter->getRotatedOperator());
       for (uint kp=0;kp<nplots;kp++){
          LogHelper xmlkp("EffEnergyPlot");
          xmlkp.putUInt("Index",kp);
          map<int,MCEstimate> results;
          oprot.resetIDIndex(kp);
          OperatorInfo opr(oprot);
          CorrelatorInfo corrinfo(opr,opr);
          getEffectiveEnergy(m_obs,corrinfo,herm,subvev,reweight,RealPart,mode,step,efftype,results);
          if (results.empty()){ 
             xmlkp.putString("Error","Could not make plot");
             xmllog.put(xmlkp);
             continue;}  // skip this plot
          if (maxerror>0.0){
             map<int,MCEstimate> raw(results);
             results.clear();
             for (map<int,MCEstimate>::const_iterator it=raw.begin();it!=raw.end();it++)
                if ((it->second).getSymmetricError()<std::abs(maxerror)) results.insert(*it);}
          vector<XYDYPoint> meffvals(results.size());
          uint k=0;
          for (map<int,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
             meffvals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(),
                                  (rt->second).getSymmetricError());}
          string plotfile(plotfilestub+"_"+make_string(kp)+".agr");
          string corrname("Corr");
          try{corrname=getCorrelatorStandardName(corrinfo);}
          catch(const std::exception& xp){}
          createEffEnergyPlot(meffvals,RealPart,corrname,plotfile,symboltype,color);
          xmlkp.putString("PlotStatus","Success");
          xmlkp.putString("PlotFile",plotfile);
          xmllog.put(xmlkp);}
       xmlout.putItem(xmllog); 
       xmlout.output(xml_out);}
       catch(const std::exception& msg){
          xmlout.putString("Error",string(msg.what()));}}
    
    if (xmltask.queryTag("PlotRotatedCorrelators")){
       try{
       ArgsHandler xmlc(xmltask,"PlotRotatedCorrelators");
       LogHelper xmllog("PlotRotatedCorrelators");
       string instr("Jackknife");
       xmlc.getOptionalString("SamplingMode",instr);
       SamplingMode mode=Jackknife;
       if (instr=="Bootstrap") mode=Bootstrap;
       else if (instr=="Jackknife") mode=Jackknife;
       else throw(std::invalid_argument("Bad sampling mode"));
       ComplexArg arg=RealPart;
       if (xmlc.queryTag("Arg")){
	 string arg_temp;
	 xmlc.getOptionalString("Arg",arg_temp);
	 if ((arg_temp=="Re")||(arg_temp=="RealPart")) arg=RealPart;
	 else if ((arg_temp=="Im")||(arg_temp=="ImaginaryPart")) arg=ImaginaryPart;
	 else throw(std::invalid_argument("Invalid Arg tag"));}
       string plotfilestub(xmlc.getString("PlotFileStub"));
       string color("blue"),symboltype("circle");
       xmlc.getOptionalString("SymbolColor",color);
       xmlc.getOptionalString("SymbolType",symboltype);
       xmllog.putEcho(xmlc);
       uint nplots=pivoter->getNumberOfLevels();
       xmllog.putUInt("NumberOfPlots",nplots);
       bool herm=true;
       bool subvev=pivoter->subtractVEV();
       bool reweight=pivoter->reweight();
       double rescale=1.0;
       xmlc.getOptionalReal("Rescale",rescale);
       GenIrrepOperatorInfo oprot(pivoter->getRotatedOperator());
       for (uint kp=0;kp<nplots;kp++){
          LogHelper xmlkp("CorrelatorPlot");
          xmlkp.putUInt("Index",kp);
          map<int,MCEstimate> results;
          oprot.resetIDIndex(kp);
          OperatorInfo opr(oprot);
          CorrelatorInfo corrinfo(opr,opr);
	  getCorrelatorEstimates(m_obs,corrinfo,herm,subvev,reweight,arg,mode,results);
          if (results.empty()){ 
             xmlkp.putString("Error","Could not make plot -- No correlator estimates could be obtained");
             xmllog.put(xmlkp);
             continue;}  // skip this plot
	  vector<XYDYPoint> corrvals(results.size());
	  uint k=0;
	  for (map<int,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
	    corrvals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(),
				  (rt->second).getSymmetricError());}
          string plotfile(plotfilestub+"_"+make_string(kp)+".agr");
          string corrname("Corr");
          try{corrname=getCorrelatorStandardName(corrinfo);}
          catch(const std::exception& xp){}
	  createCorrelatorPlot(corrvals,arg,corrname,plotfile,symboltype,color,rescale);
          xmlkp.putString("PlotStatus","Success");
          xmlkp.putString("PlotFile",plotfile);
	  if (arg==RealPart) xmlkp.putString("Arg","RealPart");
	  else xmlkp.putString("Arg","ImaginaryPart");
	  xmlkp.putBoolAsEmpty("HermitianMatrix", herm);
	  xmlkp.putBoolAsEmpty("SubtractVEV", subvev);
	  xmlkp.putBoolAsEmpty("Reweight", reweight);
	  if (mode==Jackknife) xmlkp.putString("SamplingMode","Jackknife");
	  else xmlkp.putString("SamplingMode","Bootstrap");
          xmllog.put(xmlkp);}
       xmlout.putItem(xmllog); 
       xmlout.output(xml_out);}
       catch(const std::exception& msg){
          xmlout.putString("Error",string(msg.what()));}}

       // delete pivoter if not put into persistent memory
    if (!pkeep) delete pivoter;}

 else{
    throw(std::invalid_argument("Unsupported rotation type"));}

 xmlout.output(xml_out);
}


void TaskHandler::doRotCorrMatrixInsertFitInfos(XMLHandler& xml_task, 
                        XMLHandler& xml_out, int taskcount)
{
 LogHelper xmlout;
 ArgsHandler xmltask(xml_task);
 xmlout.reset("DoRotCorrMatInsertFitInfos");
 map<uint,MCObsInfo> energyfits;
 string ecommon;
 xmltask.getOptionalString("EnergyFitCommonName",ecommon);
 if (!ecommon.empty()){
    xmlout.putString("EnergyFitCommonName",ecommon);}
 else{
    list<XMLHandler> xmlens=xml_task.find_among_children("EnergyFit"); 
    for (list<XMLHandler>::iterator it=xmlens.begin();it!=xmlens.end();it++){
       ArgsHandler xmle(*it);
       uint level=xmle.getUInt("Level");
       string name(xmle.getString("Name"));
       uint index=taskcount;
       xmle.getOptionalUInt("IDIndex",index);
       MCObsInfo energykey(name,index);
       energyfits.insert(make_pair(level,energykey));
       xmlout.putEcho(xmle);}}
 map<uint,MCObsInfo> ampfits;
 string common;
 xmltask.getOptionalString("RotatedAmplitudeCommonName",common);
 if (!common.empty()){
    xmlout.putString("RotatedAmplitudeCommonName",common);}
 else{
    list<XMLHandler> xmlamps=xml_task.find_among_children("RotatedAmplitude"); 
    for (list<XMLHandler>::iterator it=xmlamps.begin();it!=xmlamps.end();it++){
       ArgsHandler xmla(*it);
       uint level=xmla.getUInt("Level");
       string name(xmla.getString("Name"));
       uint index=taskcount;
       xmla.getOptionalUInt("IDIndex",index);
       MCObsInfo ampkey(name,index);
       ampfits.insert(make_pair(level,ampkey));
       xmlout.putEcho(xmla);}}
 bool reorder=false;
 xmltask.getOptionalBool("ReorderByFitEnergy",reorder);

 string rotatetype(xmltask.getString("Type"));

 if (rotatetype=="SinglePivot"){
    ArgsHandler xmlpiv(xmltask,"SinglePivotInitiate");
    LogHelper xmllog;
    bool pkeep;
    SinglePivotOfCorrMat* pivoter=SinglePivotOfCorrMat::initiateSinglePivot(
                             *this,xmlpiv,xmllog,pkeep);
    if (pivoter==0){
       xmlout.output(xml_out);
       throw(std::runtime_error("Could not initiate Single Pivot"));}

    if (!ecommon.empty()){
       MCObsInfo ecommonkey(ecommon,0);
       for (uint level=0;level<pivoter->getNumberOfLevels();level++){
          ecommonkey.resetObsIndex(level);
          pivoter->insertEnergyFitInfo(level,ecommonkey);
          LogHelper xmlinsert("Inserted"); xmlinsert.putUInt("Level",level);
          xmlinsert.putItem("EnergyFitInfo",ecommonkey);
          xmllog.put(xmlinsert);}}
    else{
       for (map<uint,MCObsInfo>::iterator it=energyfits.begin();it!=energyfits.end();it++){
          pivoter->insertEnergyFitInfo(it->first,it->second);
          LogHelper xmlinsert("Inserted"); xmlinsert.putUInt("Level",it->first);
          xmlinsert.putItem("EnergyFitInfo",it->second);
          xmllog.put(xmlinsert);}}

    if (!common.empty()){
       MCObsInfo commonkey(common,0);
       for (uint level=0;level<pivoter->getNumberOfLevels();level++){
          commonkey.resetObsIndex(level);
          pivoter->insertAmplitudeFitInfo(level,commonkey);
          LogHelper xmlinsert("Inserted"); xmlinsert.putUInt("Level",level);
          xmlinsert.putItem("AmplitudeFitInfo",commonkey);
          xmllog.put(xmlinsert);}}
    else{
       for (map<uint,MCObsInfo>::iterator it=ampfits.begin();it!=ampfits.end();it++){
          pivoter->insertAmplitudeFitInfo(it->first,it->second);
          LogHelper xmlinsert("Inserted"); xmlinsert.putUInt("Level",it->first);
          xmlinsert.putItem("AmplitudeFitInfo",it->second);
          xmllog.put(xmlinsert);}}

    try{
    if (reorder){ 
       LogHelper xmlreo;
       pivoter->reorderLevelsByFitEnergy(xmlreo);
       xmllog.put(xmlreo);}}
    catch(const std::exception& errmsg){
       xmlout.putItem(xmllog); xmlout.output(xml_out);
       throw(std::invalid_argument(string("Error in SinglePivotOfCorrMat::reorderLevelsByFitEnergy: ")
              +string(errmsg.what())));}
    xmlout.putItem(xmllog);

       // delete pivoter if not put into persistent memory
    if (!pkeep) delete pivoter;}

 xmlout.output(xml_out);
}


// *        <Action>DoCorrMatrixRelabelEnergyPlots</Action>                          *
// *        <Type>SinglePivot</Type>                                                 *
// *        <SinglePivotInitiate> ... </SinglePivotInitiate> (depends on type)       *
// *           ...use <GetFromMemory> ...                                            *
// *        <OriginalPlotFiles>                                                      *
// *           <PlotFileStub> ... </PlotFileStub>                                    *
// *            or  <PlotFile>name0</PlotFile> ....                                  *
// *        </OriginalPlotFiles>                                                     *
// *        <RevisedPlotFiles>   (optional)                                          *
// *           <PlotFileStub> ... </PlotFileStub>                                    *
// *            or  <PlotFile>name0</PlotFile> ....                                  *
// *        </RevisedPlotFiles>                                                      *
// *     </Task>                                                                     *

void TaskHandler::doRotCorrMatrixRelabelEnergyPlots(XMLHandler& xml_task, 
                        XMLHandler& xml_out, int taskcount)
{ /*
 LogHelper xmlout;
 ArgsHandler xmltask(xml_task);
 xmlout.reset("DoRotCorrMatrixRelabelEnergyPlots");

 ArgsHandler xmlof(xmltask,"OriginalPlotFiles");
 string stub;
 xmlof.getOptionalString("PlotFileStub",stub);
 if (!stub.empty()){
    xmlout.putString("PlotFileStub",stub);}
 else{
    list<XMLHandler> xmlf=xmlof.find_among_children("PlotFile"); 
    for (list<XMLHandler>::iterator it=xmlf.begin();it!=xmlf.end();it++){
       it->getString("PlotFile");
       uint level=xmla.getUInt("Level");
       string name(xmla.getString("Name"));
       uint index=taskcount;
       xmla.getOptionalUInt("IDIndex",index);
       MCObsInfo ampkey(name,index);
       ampfits.insert(make_pair(level,ampkey));
       xmlout.putEcho(xmla);}}  */
}

/*
void TaskHandler::doRotCorrMatrixReorderLevelsByEnergy(XMLHandler& xml_task, 
                        XMLHandler& xml_out, int taskcount)
{
 LogHelper xmlout;
 ArgsHandler xmltask(xml_task);
 xmlout.reset("DoRotCorrMatReorderLevelsByEnergy");
 map<uint,MCObsInfo> energyfits;
 list<XMLHandler> xmlens=xml_task.find_among_children("EnergyFit"); 
 for (list<XMLHandler>::iterator it=xmlens.begin();it!=xmlens.end();it++){
    ArgsHandler xmle(*it);
    uint level=xmle.getUInt("Level");
    string name(xmle.getString("Name"));
    uint index=taskcount;
    xmle.getOptionalUInt("IDIndex",index);
    MCObsInfo energykey(name,index);
    energyfits.insert(make_pair(level,energykey));
    xmlout.putEcho(xmle);}

 string rotatetype(xmltask.getString("Type"));

 if (rotatetype=="SinglePivot"){
    ArgsHandler xmlpiv(xmltask,"SinglePivotInitiate");
    LogHelper xmllog;
    bool pkeep;
    SinglePivotOfCorrMat* pivoter=SinglePivotOfCorrMat::initiateSinglePivot(
                             *this,xmlpiv,xmllog,pkeep);
    if (pivoter==0){
       xmlout.output(xml_out);
       throw(std::runtime_error("Could not initiate Single Pivot"));}

    for (map<uint,MCObsInfo>::iterator it=energyfits.begin();it!=energyfits.end();it++)
       pivoter->insertEnergyFitInfo(it->first,it->second);

    try{
    pivoter->reorderLevelsByFitEnergy();}
    catch(const std::exception& errmsg){
       xmlout.putItem(xmllog); xmlout.output(xml_out);
       throw(std::invalid_argument(string("Error in SinglePivotOfCorrMat::reorderLevelsByFitEnergy: ")
              +string(errmsg.what())));}
    xmlout.putItem(xmllog);

       // delete pivoter if not put into persistent memory
    if (!pkeep) delete pivoter;}

 xmlout.output(xml_out);
}*/



void TaskHandler::doCorrMatrixZMagSquares(XMLHandler& xml_task, 
                        XMLHandler& xml_out, int taskcount)
{
 LogHelper xmlout;
 ArgsHandler xmltask(xml_task);
 xmlout.reset("DoCorrMatrixZMagSquares");

 string rotatetype(xmltask.getString("Type"));

 if (rotatetype=="SinglePivot"){
    ArgsHandler xmlpiv(xmltask,"SinglePivotInitiate");
    LogHelper xmllog;
    bool pkeep;
    SinglePivotOfCorrMat* pivoter=SinglePivotOfCorrMat::initiateSinglePivot(
                             *this,xmlpiv,xmllog,pkeep);
    if (pivoter==0){
       xmlout.output(xml_out);
       throw(std::runtime_error("Could not initiate Single Pivot"));}

    Matrix<MCEstimate> ZMagSq;
    try{
    pivoter->computeZMagnitudesSquared(ZMagSq);}
    catch(const std::exception& errmsg){
       xmlout.putItem(xmllog); xmlout.output(xml_out);
       throw(std::invalid_argument(string("Error in SinglePivotOfCorrMat::computeZMagnitudesSquared: ")
              +string(errmsg.what())));}
    xmlout.putItem(xmllog);

    const set<OperatorInfo>& opset=pivoter->getOperators();
    uint nlevels=pivoter->getNumberOfLevels();
    uint opindex=0;
    for (set<OperatorInfo>::const_iterator opit=opset.begin();opit!=opset.end();opit++,opindex++){
       LogHelper xmlzop("OperatorZMagnitudeSquares");
       xmlzop.putInt("OperatorIndex",opindex);
       xmlzop.putItem(*opit);
       double rescale=0.0;
       for (uint level=0;level<nlevels;level++)
          rescale+=ZMagSq(opindex,level).getFullEstimate();
       rescale=1.0/rescale;
       for (uint level=0;level<nlevels;level++){
          ZMagSq(opindex,level).rescale(rescale);
          LogHelper xmlzcoef("ZMagSquare");
          xmlzcoef.putInt("OpIndex",opindex);
          xmlzcoef.putInt("Level",level);
          xmlzcoef.putItem("Value",ZMagSq(opindex,level));
          xmlzop.putItem(xmlzcoef);}
       xmlout.putItem(xmlzop);}

       // make plots of Z mag squares if requested

    if (xmltask.queryTag("DoPlots")){
       try{
       ArgsHandler xmlc(xmltask,"DoPlots");
       LogHelper xmllog("DoPlots");
       string plotfilestub(xmlc.getString("PlotFileStub"));
       string barcolor("cyan");
       xmlc.getOptionalString("BarColor",barcolor);
       list<ArgsHandler> zplots=xmlc.getSubHandlers("ZMagSqPlot");
       map<MCObsInfo,tuple<string,string,uint> > zplotinfos;
       for (list<ArgsHandler>::iterator zt=zplots.begin();zt!=zplots.end();zt++){
          OperatorInfo zop(zt->getItem<OperatorInfo>("ZMagSqPlot"));
          MCObsInfo obs(zop);
          string obsname,suffix;
          uint opindex=0;
          for (set<OperatorInfo>::const_iterator opit=opset.begin();opit!=opset.end();opit++,opindex++){
             if (obs==*opit){
                obsname=string("Operator index ")+make_string(opindex);
                suffix=make_string(opindex); break;}}
          if (opindex<opset.size()){
             zt->getOptionalString("ObsName",obsname);
             if (obsname=="standard") obsname=getOpStandardName(zop);
             zt->getOptionalString("FileSuffix",suffix);
             zplotinfos.insert(make_pair(obs,make_tuple(obsname,suffix,opindex)));}}
       for (map<MCObsInfo,tuple<string,string,uint> >::iterator it=zplotinfos.begin();it!=zplotinfos.end();it++){
          string obsname=get<0>(it->second);
          string suffix=get<1>(it->second);
          uint opindex=get<2>(it->second);
          string plotfile(plotfilestub+"_"+suffix+".agr");
          xmllog.putString("PlotFile",plotfile);
          vector<XYDYPoint> zmag_sqs(nlevels);
          for (uint level=0;level<nlevels;level++){
             zmag_sqs[level].xval=level;
             zmag_sqs[level].yval=ZMagSq(opindex,level).getFullEstimate();
             zmag_sqs[level].yerr=ZMagSq(opindex,level).getSymmetricError();}
          createCorrMatrixZMagSquaresPlot(zmag_sqs,obsname,plotfile,barcolor);}
       xmlout.putItem(xmllog);}
       catch(const std::exception& msg){
          xmlout.putString("Error",string(msg.what()));}} 

       // delete pivoter if not put into persistent memory
    if (!pkeep) delete pivoter;}

 xmlout.output(xml_out);
}


// ***************************************************************************************
 

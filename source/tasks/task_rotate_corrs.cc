#include "task_handler.h"
#include "task_utils.h"
#include "correlator_matrix_info.h"
#include "single_pivot.h"
#include "create_plots.h"

using namespace std;

// ***********************************************************************************
// *                                                                                 *
// *   XML format for correlator matrix rotations:                                   *
// *                                                                                 *
// *     <Task>                                                                      *
// *        <Action>DoCorrMatrixRotation</Action>                                    *
// *        <MinTimeSep>3</MinTimeSep>                                               *
// *        <MaxTimeSep>30</MaxTimeSep>                                              *
// *        <Type>SinglePivot</Type>                                                 *
// *        <SinglePivotInitiate> ... </SinglePivotInitiate> (depends on type)       *
// *        <WriteRotatedCorrToFile>                                                 *
// *           <RotatedCorrFileName>rotated_corr_bins</RotatedCorrFileName>          *
// *           <Overwrite/>                                                          *
// *        </WriteRotatedCorrToFile>                                                *
// *     </Task>                                                                     *
// *                                                                                 *
// *   The <Type> specifies which kind of pivot to create.  The tag                  *
// *   <SinglePivotInitiate> applies only for a type of SinglePivot.  For            *
// *   other types, there will be a different tag instead of <SinglePivotInitiate>.  *
// *                                                                                 *
// *                                                                                 *
// *   If the <Type> is SinglePivot, then the input XML to create the new pivot:     *
// *                                                                                 *
// *      <SinglePivotInitiate>                                                      *
// *         <RotatedCorrelator>                                                     *
// *           <GIOperator>...</GIOperator>                                          *
// *         </RotatedCorrelator>                                                    *
// *         <AssignName>PivTester</AssignName>  (optional)                          *
// *         <CorrelatorMatrixInfo> ... </CorrelatorMatrixInfo>                      *
// *         <NormTime>3</NormTime>                                                  *
// *         <MetricTime>6</MetricTime>                                              *
// *         <DiagonalizeTime>12</DiagonalizeTime>                                   *
// *         <MinimumInverseConditionNumber>0.01</MinimumInverseConditionNumber>     *
// *         <NegativeEigenvalueAlarm>-0.01</NegativeEigenvalueAlarm>  (optional)    *
// *         <CheckMetricErrors/>    (optional)                                      *
// *         <CheckCommonMetricMatrixNullSpace/>    (optional)                       *
// *         <WritePivotToFile>    (optional)                                        *
// *            <PivotFileName>pivot_test</PivotFileName>                            *
// *            <Overwrite/>                                                         *
// *         </WritePivotToFile>                                                     *
// *      </SinglePivotInitiate>                                                     *
// *                                                                                 *
// *   The <RotatedCorrelator> tag specifies the name to give the rotated            *
// *   operators.  Any integer index specified is ignored.  If the matrix to         *
// *   be rotated is "N" x "N", then the N rotated operators will all have the       *
// *   same isospin and irrep labels and ID name, but the ID index will vary         *
// *   from 0 to N-1.                                                                *
// *                                                                                 *
// *   The <CorrelatorMatrixInfo> tag specifies the original correlator matrix       *
// *   of operators to be rotated.                                                   *
// *                                                                                 *
// *                                                                                 *
// *   Input XML to set up a previously created pivot saved in a file:               *
// *                                                                                 *
// *      <SinglePivotInitiate>                                                      *
// *         <ReadPivotFromFile>                                                     *
// *            <PivotFileName>pivot_file</PivotFileName>                            *
// *         </ReadPivotFromFile>                                                    *
// *      </SinglePivotInitiate>                                                     *
// *                                                                                 *
// *   Input XML to set up a previously created pivot saved in memory:               *
// *                                                                                 *
// *      <SinglePivotInitiate>                                                      *
// *         <GetFromMemory>                                                         *
// *            <IDName>PivTester</IDName>                                           *
// *         </GetFromMemory>                                                        *
// *      </SinglePivotInitiate>                                                     *
// *                                                                                 *
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
    pivoter->doRotation(mintimesep,maxtimesep,xmllog);}
    catch(const std::exception& errmsg){
       xmlout.putItem(xmllog); xmlout.output(xml_out);
       throw(std::invalid_argument((string("Error in SinglePivotOfCorrMat::doRotation: ")
              +string(errmsg.what())).c_str()));} 
    xmlout.putItem(xmllog);
       // save rotated correlators to file
    if (xmltask.queryTag("WriteRotatedCorrToFile")){
       try{
       ArgsHandler xmlf(xmltask,"WriteRotatedCorrToFile");
       string corrfile(xmlf.getString("RotatedCorrFileName"));
       bool overwrite=false;
       xmlf.getOptionalBool("Overwrite",overwrite); 
       if (corrfile.empty()) throw(std::invalid_argument("Empty file name"));
       LogHelper xmlw;
       pivoter->writeRotated(mintimesep,maxtimesep,corrfile,overwrite,xmlw);
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
       bool subvev=pivoter->isVEVsubtracted();
       GenIrrepOperatorInfo oprot(pivoter->getRotatedOperator());
       for (uint kp=0;kp<nplots;kp++){
          LogHelper xmlkp("EffEnergyPlot");
          xmlkp.putUInt("Index",kp);
          map<int,MCEstimate> results;
          oprot.resetIDIndex(kp);
          OperatorInfo opr(oprot);
          CorrelatorInfo corrinfo(opr,opr);
          getEffectiveEnergy(m_obs,corrinfo,herm,subvev,RealPart,mode,step,efftype,results);
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

       // delete pivoter if not put into persistent memory
    if (!pkeep) delete pivoter;}

 xmlout.output(xml_out);
}




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
    xmlout.putItem(xmllog);
    if (pivoter==0){
       xmlout.output(xml_out);
       throw(std::runtime_error("Could not initiate Single Pivot"));}
    Matrix<MCEstimate> ZMagSq;
    try{
    pivoter->computeZMagnitudesSquared(ZMagSq);}
    catch(const std::exception& errmsg){
       xmlout.putItem(xmllog); xmlout.output(xml_out);
       throw(std::invalid_argument((string("Error in SinglePivotOfCorrMat::computeZMagnitudesSquared: ")
              +string(errmsg.what())).c_str()));}
    xmlout.putItem(xmllog);
       // make plots of Z mag squares
 /*   if (xmlpiv.queryTag("WriteRotatedCorrToFile")){
       try{
       ArgsHandler xmlf(xmlpiv,"WriteRotatedCorrToFile");
       string corrfile(xmlf.getString("FileName"));
       bool overwrite=false;
       xmlf.getOptionalBool("Overwrite",overwrite); 
       if (corrfile.empty()) throw(std::invalid_argument("Empty file name"));
       LogHelper xmlw;
       pivoter->writeRotated(mintimesep,maxtimesep,corrfile,overwrite,xmlw);
       xmlout.putItem(xmlw);}
       catch(const std::exception& msg){
          xmlout.putString("Error",string(msg.what()));}} */
       // delete pivoter if not put into persistent memory
    if (!pkeep) delete pivoter;}

 xmlout.output(xml_out);
}


// ***************************************************************************************
 

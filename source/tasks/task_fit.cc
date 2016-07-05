#include "task_handler.h"
#include "chisq_tcorr.h"
#include "chisq_fit.h"
#include "create_plots.h"
#include "task_utils.h"

using namespace std;

// *******************************************************************************
// *                                                                             *
// *    XML format for chi-square fitting:                                       *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoFit</Action>                                                  *
// *       <Type>TemporalCorrelator</Type>                                       *
// *       <MinimizerInfo>                 (optional)                            *
// *         <Method>Minuit2</Method>                                            *
// *         <ParameterRelTol>1e-6</ParameterRelTol>                             *
// *         <ChiSquareRelTol>1e-4</ChiSquareRelTol>                             *
// *         <MaximumIterations>1024</MaximumIterations>                         *
// *         <Verbosity>Low</Verbosity>                                          *
// *       </MinimizerInfo>                                                      *
// *       <SamplingMode>Bootstrap</SamplingMode>   (optional)                   *
// *       <TemporalCorrelatorFit>                                               *
// *         <Operator>.... </Operator>                                          *
// *         <MinimumTimeSeparation>3</MinimumTimeSeparation>                    *
// *         <MaximumTimeSeparation>12</MaximumTimeSeparation>                   *
// *         <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>                    *
// *         <Model>                                                             *
// *             <Type>TimeSymSingleExponential</Type>                           *
// *             <EnergyParameter>                                               *
// *                <Name>pion</Name><IDIndex>0</IDIndex> // default taskcount   *
// *             </EnergyParameter>                                              *
// *             <Amplitude>                                                     *
// *                <Name>A</Name><IDIndex>0</IDIndex>                           *
// *             </Amplitude>                                                    *
// *         </Model>                                                            *
// *       <DoEffectiveEnergyPlot>                                               *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <CorrName>standard</CorrName>   (optional)                         *
// *          <TimeStep>3</TimeStep>  (optional: 1 default)                      *
// *          <SymbolColor> ... </SymbolColor>                                   *
// *          <SymbolType> ... </SymbolType>                                     *
// *          <MaxErrorToPlot> ...</MaxErrorToPlot> (optional)                   *
// *          <Goodness>qual</Goodness>  "qual" or "chisq"                       *
// *          <ShowApproach/>   (optional)                                       *
// *       </DoEffectiveEnergyPlot>                                              *
// *       </TemporalCorrelatorFit>                                              *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *******************************************************************************


void TaskHandler::doFit(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{

 ChiSquareMinimizerInfo mz_info;  // default minimizer info
 if (xmltask.count_among_children("MinimizerInfo")>0){
    ChiSquareMinimizerInfo mz_user(xmltask);
    mz_info=mz_user;}
 string fittype;
 xmlreadchild(xmltask,"Type",fittype,"DoFit");
 xmlout.set_root("DoFit");
 XMLHandler xmlmz;
 mz_info.output(xmlmz);
 xmlout.put_child(xmlmz);
 xmlout.put_child("Type",fittype); 
 SamplingMode mode=Bootstrap;
 string instr;
 if (xmlreadifchild(xmltask,"SamplingMode",instr)){
    if (instr=="Bootstrap") mode=Bootstrap;
    else if (instr=="Jackknife") mode=Jackknife;
    else throw(std::invalid_argument("Bad sampling mode"));} 
 if (mode==Bootstrap){
    xmlout.put_child("SamplingMode","Bootstrap");
    m_obs->setToBootstrapMode();}
 else{
    xmlout.put_child("SamplingMode","Jackknife");
    m_obs->setToJackknifeMode();}
 vector<MCEstimate> bestfit_params;

 if (fittype=="TemporalCorrelator"){
    try{
    XMLHandler xmlf(xmltask,"TemporalCorrelatorFit");
    RealTemporalCorrelatorFit RTC(xmlf,*m_obs,taskcount);
    XMLHandler xmlof; RTC.output(xmlof);
    xmlout.put_child(xmlof);
    double chisq_dof,qual;
    doChiSquareFitting(RTC,mz_info,chisq_dof,qual,
                       bestfit_params,xmlout);

       // fit done, now do the plot if requested
    if (xmlf.count_among_children("DoEffectiveEnergyPlot")!=1) return;
    XMLHandler xmlp(xmlf,"DoEffectiveEnergyPlot");
    string plotfile;
    xmlreadifchild(xmlp,"PlotFile",plotfile);
    if (tidyString(plotfile).empty()){
       xmlout.put_child("Warning","No plot file but asked for plot!");
       return;}
    string symbolcolor("blue"),symboltype("circle");
    xmlreadifchild(xmlp,"SymbolColor",symbolcolor);
    xmlreadifchild(xmlp,"SymbolType",symboltype);
    string fitgood;
    xmlreadifchild(xmlp,"Goodness",fitgood);
    char goodtype='N';
    double goodness=qual;
    if (fitgood=="qual"){
       goodtype='Q'; }
    else if (fitgood=="chisq"){
       goodtype='X'; goodness=chisq_dof;}
    bool showapproach=(xml_child_tag_count(xmlp,"ShowApproach")>0);
    string corrname;
    xmlreadifchild(xmlp,"CorrName",corrname);
    uint step=1;
    if (xmlreadifchild(xmlp,"TimeStep",step)){
       if ((step<1)||(step>getLatticeTimeExtent()/4)){
          xmlout.put_child("PlotError","Bad effective energy time step");
          return;}}
    CorrelatorInfo corr(RTC.m_op,RTC.m_op);
    if (corrname=="standard") corrname=getCorrelatorStandardName(corr);
    bool hermitian=true;
    bool subvev=RTC.m_subt_vev;
    uint fit_tmin=RTC.m_tmin;
    uint fit_tmax=RTC.m_tmax;
    uint efftype=RTC.m_model_ptr->getEffMassType();
    double subt_const=0.0;
    if (efftype>1){    // subtract fit constant
       efftype-=2;     // efftypes 2 and 3 remove constant, but noisy
       subt_const=bestfit_params[bestfit_params.size()-1].getFullEstimate();}
    SamplingMode mode=m_obs->getCurrentSamplingMode();

    map<int,MCEstimate> results;
    getEffectiveEnergy(m_obs,corr,hermitian,subvev,RealPart,mode,step, 
                       efftype,results,subt_const);
    if (results.empty()){
       xmlout.put_child("PlotError","No effective energy estimates could be obtained");
       return;}
         // do some XML output
    xmlout.put_child("PlotFile",plotfile);
    XMLHandler xmlef;
    xmlef.set_root("EffectiveEnergy");
    xmlef.put_child("TimeStep",make_string(step));
    if (efftype==0) xmlef.put_child("EffEnergyType","TimeForward");
    else if (efftype==1) xmlef.put_child("EffEnergyType","TimeSymmetric");
    else if (efftype==2) xmlef.put_child("EffEnergyType","TimeForwardPlusConst");
    else if (efftype==3) xmlef.put_child("EffEnergyType","TimeSymmetricPlusConst");
    xmlef.seek_root();
    xmlef.seek_first_child();
    for (map<int,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++){
       XMLHandler xmlr("Estimate");
       xmlr.put_child("TimeSeparation",make_string(rt->first));
       xmlr.put_child("MeanValue",make_string((rt->second).getFullEstimate()));
       xmlr.put_child("SymmError",make_string((rt->second).getSymmetricError()));
       xmlef.put_sibling(xmlr);}
    xmlout.put_child(xmlef);
           // now prepare the plot
    double maxerror=0.0;
    if (xmlreadifchild(xmltask,"MaxErrorToPlot",maxerror)){
       map<int,MCEstimate> raw(results);
       results.clear();
       for (map<int,MCEstimate>::const_iterator it=raw.begin();it!=raw.end();it++)
          if ((it->second).getSymmetricError()<std::abs(maxerror)) results.insert(*it);}

    vector<XYDYPoint> meffvals(results.size());
    uint k=0;
    for (map<int,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
       meffvals[k]=XYDYPoint(rt->first, (rt->second).getFullEstimate(),
                            (rt->second).getSymmetricError());}

    TempCorrFitInfo fitinfo;
    RTC.m_model_ptr->setFitInfo(RTC.m_fitparam_info,bestfit_params,
                                fit_tmin,fit_tmax,step,
                                showapproach,step,chisq_dof,qual,fitinfo);

    createEffEnergyPlotWithFit(meffvals,RealPart,fitinfo,goodtype,goodness,corrname,
                               plotfile,symboltype,symbolcolor);
    }
    catch(const std::exception& errmsg){
       //xmlout.clear();
       xmlout.put_child("Error",string("DoFit with type TemporalCorrelator encountered an error: ")
               +string(errmsg.what()));
      // throw(std::invalid_argument("DoFit with type TemporalCorrelator encountered an error: ")+errmsg);}
    }}


}


// ***************************************************************************************
 

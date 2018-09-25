#include "task_handler.h"
#include "chisq_anisotropy.h"
#include "chisq_tcorr.h"
#include "chisq_fit.h"
#include "create_plots.h"
#include "task_utils.h"

using namespace std;

// *******************************************************************************
// *                                                                             *
// *    XML format for chi-square fitting:                                       *
// *                                                                             *
// *    Fit to a single real-valued temporal correlator.                         *
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
// *       <CovMatCalcSamplingMode>Bootstrap</CovMatCalcSamplingMode> (optional) *
// *       <TemporalCorrelatorFit>                                               *
// *         <Operator>.... </Operator>                                          *
// *         <SubtractVEV/>             (as appropriate)                         *
// *         <Reweight/>              (optional)                                 *
// *         <MinimumTimeSeparation>3</MinimumTimeSeparation>                    *
// *         <MaximumTimeSeparation>12</MaximumTimeSeparation>                   *
// *         <ExcludeTimes>4 8</ExcludeTimes>  (optional)                        *
// *         <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>                    *
// *         <Model>                                                             *
// *             <Type>TimeSymSingleExponential</Type>                           *
// *             <Energy>                                                        *
// *                <Name>pion</Name><IDIndex>0</IDIndex> // default taskcount   *
// *             </Energy>                                                       *
// *             <Amplitude>                                                     *
// *                <Name>A</Name><IDIndex>0</IDIndex>                           *
// *             </Amplitude>                                                    *
// *         </Model>                                                            *
// *       <DoEffectiveEnergyPlot> (optional)                                    *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <CorrName>standard</CorrName>   (optional)                         *
// *          <TimeStep>3</TimeStep>  (optional: 1 default)                      *
// *          <SymbolColor> ... </SymbolColor>                                   *
// *          <SymbolType> ... </SymbolType>                                     *
// *          <MaxErrorToPlot> ...</MaxErrorToPlot> (optional)                   *
// *          <Goodness>qual</Goodness>  "qual" or "chisq"                       *
// *          <ShowApproach/>   (optional)                                       *
// *          <ReferenceEnergy> (optional: includes energy ratio on plot)        *
// *            <Name>kaon</Name><IDIndex>0</IDIndex>                            *
// *          </ReferenceEnergy>                                                 *
// *       </DoEffectiveEnergyPlot>                                              *
// *       <InsertIntoPivot>  (optional)                                         *
// *         <Type>Single</Type>  (or Rolling or Principal)                      *
// *         <Name>the_pivot</Name>  (must already be in memory)                 *
// *         <Level>0</Level>                                                    *
// *       </InsertIntoPivot>                                                    *
// *       </TemporalCorrelatorFit>                                              *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *    Fit to two single real-valued temporal correlators. Plot is done         *
// *    for first correlator (second correlator is considered the reference      *
// *    correlator).  In plot, fit value is given as a RATIO of the first        *
// *    energy over the second energy.                                           *
// *                                                                             *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoFit</Action>                                                  *
// *       <Type>TwoTemporalCorrelator</Type>                                    *
// *       <MinimizerInfo>                 (optional)                            *
// *         <Method>Minuit2</Method>                                            *
// *         <ParameterRelTol>1e-6</ParameterRelTol>                             *
// *         <ChiSquareRelTol>1e-4</ChiSquareRelTol>                             *
// *         <MaximumIterations>1024</MaximumIterations>                         *
// *         <Verbosity>Low</Verbosity>                                          *
// *       </MinimizerInfo>                                                      *
// *       <SamplingMode>Bootstrap</SamplingMode>   (optional)                   *
// *       <CovMatCalcSamplingMode>Bootstrap</CovMatCalcSamplingMode> (optional) *
// *       <TwoTemporalCorrelatorFit>                                            *
// *         <CorrelatorOne>                                                     *
// *           <Operator>.... </Operator>                                        *
// *           <SubtractVEV/>             (as appropriate)                       *
// *           <Reweight/>              (optional)                               *
// *           <MinimumTimeSeparation>3</MinimumTimeSeparation>                  *
// *           <MaximumTimeSeparation>12</MaximumTimeSeparation>                 *
// *           <ExcludeTimes>4 8</ExcludeTimes>  (optional)                      *
// *           <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>                  *
// *           <Model>...</Model>                                                *
// *         </CorrelatorOne>                                                    *
// *         <CorrelatorTwo>                                                     *
// *           <Operator>.... </Operator>                                        *
// *           <SubtractVEV/>             (as appropriate)                       *
// *           <Reweight/>              (optional)                               *
// *           <MinimumTimeSeparation>3</MinimumTimeSeparation>                  *
// *           <MaximumTimeSeparation>12</MaximumTimeSeparation>                 *
// *           <ExcludeTimes>4 8</ExcludeTimes>  (optional)                      *
// *           <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>                  *
// *           <Model>...</Model>                                                *
// *         </CorrelatorTwo>                                                    *
// *         <EnergyRatio>                                                       *
// *            <Name>pion</Name><IDIndex>0</IDIndex>                            *
// *         </EnergyRatio>                                                      *
// *         <DoEffectiveEnergyPlot>     (plot correlator 1 only)                *
// *           <PlotFile> ... </PlotFile>                                        *
// *           <CorrName>standard</CorrName>   (optional)                        *
// *           <TimeStep>3</TimeStep>  (optional: 1 default)                     *
// *           <SymbolColor> ... </SymbolColor>                                  *
// *           <SymbolType> ... </SymbolType>                                    *
// *           <MaxErrorToPlot> ...</MaxErrorToPlot> (optional)                  *
// *           <Goodness>qual</Goodness>  "qual" or "chisq"                      *
// *           <ShowApproach/>   (optional)                                      *
// *         </DoEffectiveEnergyPlot>                                            *
// *         <InsertIntoPivot>  (optional, correlator 1 inserted)                *
// *           <Type>Single</Type>  (or Rolling or Principal)                    *
// *           <Name>the_pivot</Name>  (must already be in memory)               *
// *           <Level>0</Level>                                                  *
// *         </InsertIntoPivot>                                                  *
// *       </TwoTemporalCorrelatorFit>                                           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *    Fit the free-particle energies squared for various three-momenta squared *
// *    to estimate the lattice anisotropy  a_s/a_t.                             *
// *    The model used for the observables is                                    *
// *                                                                             *
// *      (a_t E)^2 = restmass_sq +  (2*Pi/Ns)^2 * nsq / xi^2                    *
// *                                                                             *
// *    where "restmass_sq" and "xi" are the two model parameters, and           *
// *    "Ns" is extent of the lattice in terms of number of sites                *
// *    in each of the three spatial directions, and "nsq" is the                *
// *    integer square of the three momentum.  Recall that the                   *
// *    (a_s P)^2 = (2*Pi/Ns)^2 * nsq.                                           *
// *                                                                             *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoFit</Action>                                                  *
// *       <Type>AnisotropyFromDispersion</Type>                                 *
// *       <MinimizerInfo>                 (optional)                            *
// *         <Method>Minuit2</Method>                                            *
// *         <ParameterRelTol>1e-6</ParameterRelTol>                             *
// *         <ChiSquareRelTol>1e-4</ChiSquareRelTol>                             *
// *         <MaximumIterations>1024</MaximumIterations>                         *
// *         <Verbosity>Low</Verbosity>                                          *
// *       </MinimizerInfo>                                                      *
// *       <SamplingMode>Bootstrap</SamplingMode>   (optional)                   *
// *       <CovMatCalcSamplingMode>Bootstrap</CovMatCalcSamplingMode> (optional) *
// *       <AnisotropyFromDispersionFit>                                         *
// *         <SpatialExtentNumSites>24</SpatialExtentNumSites>                   *
// *         <Energy>                                                            *
// *           <Name>pion</Name><IDIndex>0</IDIndex>                             *
// *           <IntMomSquared>0</IntMomSquared>                                  *
// *         </Energy>                                                           *
// *         <Energy>                                                            *
// *           <Name>pion</Name><IDIndex>1</IDIndex>                             *
// *           <IntMomSquared>1</IntMomSquared>                                  *
// *         </Energy>                                                           *
// *         <Energy>... </Energy>                                               *
// *         <Anisotropy>                                                        *
// *           <Name>PionXi</Name><IDIndex>0</IDIndex>                           *
// *         </Anisotropy>                                                       *
// *         <RestMassSquared>                                                   *
// *           <Name>PionRestMassSquared</Name><IDIndex>0</IDIndex>              *
// *         </RestMassSquared>                                                  *
// *         <DoPlot>                                                            *
// *           <PlotFile> ... </PlotFile>                                        *
// *           <ParticleName>pion</ParticleName>   (optional)                    *
// *           <SymbolColor> ... </SymbolColor>                                  *
// *           <SymbolType> ... </SymbolType>                                    *
// *           <Goodness>qual</Goodness>  "qual" or "chisq"                      *
// *         </DoPlot>                                                           *
// *       </AnisotropyFromDispersionFit>                                        *
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
 SamplingMode mode=Jackknife;
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
 SamplingMode covcalcmode=Jackknife;
 instr.clear();
 if (xmlreadifchild(xmltask,"CovMatCalcSamplingMode",instr)){
    if (instr=="Bootstrap") covcalcmode=Bootstrap;
    else if (instr=="Jackknife") covcalcmode=Jackknife;
    else throw(std::invalid_argument("Bad cov mat calc sampling mode"));} 
 if (covcalcmode==Bootstrap){
    xmlout.put_child("CovMatCalcSamplingMode","Bootstrap");
    m_obs->setCovMatToBootstrapMode();}
 else{
    xmlout.put_child("CovMatCalcSamplingMode","Jackknife");
    m_obs->setCovMatToJackknifeMode();}
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
    bool reweight=RTC.m_reweight;
    uint fit_tmin=RTC.getTmin();
    uint fit_tmax=RTC.getTmax();
    uint efftype=RTC.m_model_ptr->getEffMassType();
    double subt_const=0.0;
    if (efftype>1){    // subtract fit constant
       efftype-=2;     // efftypes 2 and 3 remove constant, but noisy
       subt_const=bestfit_params[bestfit_params.size()-1].getFullEstimate();}
    SamplingMode mode=m_obs->getCurrentSamplingMode();

    map<int,MCEstimate> results;
    getEffectiveEnergy(m_obs,corr,hermitian,subvev,reweight,RealPart,mode,step, 
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

    TCorrFitInfo fitinfo;
    RTC.m_model_ptr->setFitInfo(RTC.m_fitparam_info,bestfit_params,fit_tmin,fit_tmax,
                                showapproach,step,chisq_dof,qual,fitinfo);

    uint refcount=xmlp.count("ReferenceEnergy");
    if (refcount!=1){
       createEffEnergyPlotWithFit(meffvals,RealPart,fitinfo,goodtype,goodness,corrname,
                                  plotfile,symboltype,symbolcolor);}
    else if (refcount==1){
       XMLHandler xmlref(xmlp,"ReferenceEnergy");
       string refname; int refindex;
       xmlreadchild(xmlref,"Name",refname);
       if (refname.empty()) throw(std::invalid_argument("Must provide name for reference energy"));
       refindex=taskcount;
       xmlreadifchild(xmlref,"IDIndex",refindex);
       MCObsInfo refkey(refname,refindex);  // reference energy
       MCObsInfo enratio(string("TempEnergyRatioGwiqb"),taskcount);  // temporary name for ratio
       for (m_obs->setSamplingBegin();!m_obs->isSamplingEnd();m_obs->setSamplingNext()){
          double ratiovalue=m_obs->getCurrentSamplingValue(fitinfo.energy_key)
                           /m_obs->getCurrentSamplingValue(refkey);
          m_obs->putCurrentSamplingValue(enratio,ratiovalue);}
       MCEstimate ratioest=m_obs->getEstimate(enratio);
       XMLHandler xmlrat("EnergyRatioFitResult");
       XMLHandler xmlrr;
       ratioest.output(xmlrr); xmlrat.put_child(xmlrr);
       xmlout.put_child(xmlrat);
       createEffEnergyPlotWithFitAndEnergyRatio(meffvals,RealPart,fitinfo,
                           ratioest.getFullEstimate(),ratioest.getSymmetricError(),
                           goodtype,goodness,corrname,
                           plotfile,symboltype,symbolcolor);
       m_obs->eraseData(enratio);}
    }
    catch(const std::exception& errmsg){
       //xmlout.clear();
       xmlout.put_child("Error",string("DoFit with type TemporalCorrelator encountered an error: ")
               +string(errmsg.what()));
      // throw(std::invalid_argument("DoFit with type TemporalCorrelator encountered an error: ")+errmsg);}
    }}


 else if (fittype=="TwoTemporalCorrelator"){
    try{
    XMLHandler xmlf(xmltask,"TwoTemporalCorrelatorFit");
    TwoRealTemporalCorrelatorFit RTC(xmlf,*m_obs,taskcount);
    XMLHandler xmlof; RTC.output(xmlof);
    xmlout.put_child(xmlof);
    double chisq_dof,qual;
    doChiSquareFitting(RTC,mz_info,chisq_dof,qual,
                       bestfit_params,xmlout);

         // fit done, now evaluate energy ratio
    TCorrFitInfo fitinfo1, fitinfo2;
    int nparam1=RTC.m_model1_ptr->getNumberOfParams();
    vector<MCObsInfo> fitparam1_info(RTC.m_fitparam_info.begin(),RTC.m_fitparam_info.begin()+nparam1);
    vector<MCObsInfo> fitparam2_info(RTC.m_fitparam_info.begin()+nparam1,RTC.m_fitparam_info.end());
    vector<MCEstimate> bestfitparam1(bestfit_params.begin(),bestfit_params.begin()+nparam1);
    vector<MCEstimate> bestfitparam2(bestfit_params.begin()+nparam1,bestfit_params.end());
    RTC.m_model1_ptr->setFitInfo(fitparam1_info,bestfitparam1,RTC.getTmin1(),RTC.getTmax1(),
                                 false,1,chisq_dof,qual,fitinfo1);
    RTC.m_model2_ptr->setFitInfo(fitparam2_info,bestfitparam2,RTC.getTmin2(),RTC.getTmax2(),
                                 false,1,chisq_dof,qual,fitinfo2);
    for (m_obs->setSamplingBegin();!m_obs->isSamplingEnd();m_obs->setSamplingNext()){
       double ratiovalue=m_obs->getCurrentSamplingValue(fitinfo1.energy_key)
                        /m_obs->getCurrentSamplingValue(fitinfo2.energy_key);
       m_obs->putCurrentSamplingValue(RTC.m_energyratio,ratiovalue);}
    MCEstimate ratioest=m_obs->getEstimate(RTC.m_energyratio);
    XMLHandler xmlrat("EnergyRatioFitResult");
    XMLHandler xmlrr;
    ratioest.output(xmlrr); xmlrat.put_child(xmlrr);
    xmlout.put_child(xmlrat);
    
       // do the plot if requested for correlator 1
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
    CorrelatorInfo corr(RTC.m_op1,RTC.m_op1);
    if (corrname=="standard") corrname=getCorrelatorStandardName(corr);
    bool hermitian=true;
    bool subvev=RTC.m_subt_vev1;
    bool reweight=RTC.m_reweight1;
    uint efftype=RTC.m_model1_ptr->getEffMassType();
    double subt_const=0.0;
    if (efftype>1){    // subtract fit constant
       efftype-=2;     // efftypes 2 and 3 remove constant, but noisy
       subt_const=bestfit_params[bestfit_params.size()-1].getFullEstimate();}
    SamplingMode mode=m_obs->getCurrentSamplingMode();

    map<int,MCEstimate> results;
    getEffectiveEnergy(m_obs,corr,hermitian,subvev,reweight,RealPart,mode,step, 
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

    RTC.m_model1_ptr->setFitInfo(fitparam1_info,bestfitparam1,RTC.getTmin1(),RTC.getTmax1(),
                                 showapproach,step,chisq_dof,qual,fitinfo1);

    createEffEnergyPlotWithFitAndEnergyRatio(meffvals,RealPart,fitinfo1,
                               ratioest.getFullEstimate(),ratioest.getSymmetricError(),
                               goodtype,goodness,corrname,
                               plotfile,symboltype,symbolcolor);
    }
    catch(const std::exception& errmsg){
       //xmlout.clear();
       xmlout.put_child("Error",string("DoFit with type TwoTemporalCorrelator encountered an error: ")
               +string(errmsg.what()));
      // throw(std::invalid_argument("DoFit with type TemporalCorrelator encountered an error: ")+errmsg);}
    }}


 else if (fittype=="AnisotropyFromDispersion"){
    try{
    XMLHandler xmlf(xmltask,"AnisotropyFromDispersionFit");
    AnisotropyFromDispersionFit AFD(xmlf,*m_obs,taskcount);
    XMLHandler xmlof; AFD.output(xmlof);
    xmlout.put_child(xmlof);
    double chisq_dof,qual;
    doChiSquareFitting(AFD,mz_info,chisq_dof,qual,
                       bestfit_params,xmlout);

         // fit done, now make plot if requested
    if (xmlf.count_among_children("DoPlot")!=1) return;
    XMLHandler xmlp(xmlf,"DoPlot");
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
    MCEstimate xiest=m_obs->getEstimate(AFD.getAnisotropyKey());
         // do some XML output
    string particlename;
    xmlreadifchild(xmlp,"ParticleName",particlename);
    xmlout.put_child("PlotFile",plotfile);

    vector<XYDYPoint> Esq(AFD.m_nobs);
    vector<XYPoint> upperfit(2), lowerfit(2);
    uint kmin=0, kmax=0;
    for (uint k=0;k<AFD.m_nobs;k++){
       MCEstimate est=m_obs->getEstimate(AFD.m_obs_info[k]);
       Esq[k].xval=AFD.m_imomsq[k];
       Esq[k].yval=est.getFullEstimate();
       Esq[k].yerr=est.getSymmetricError();
       if (Esq[k].xval<Esq[kmin].xval) kmin=k;
       if (Esq[k].xval>Esq[kmax].xval) kmax=k;}
    MCObsInfo randtemp("RandomTemporary",0);
    doDispersionBySamplings(*m_obs,AFD.getAnisotropyKey(),AFD.getRestMassSquaredKey(), 
                            AFD.m_momsq_quantum*AFD.m_imomsq[kmin],randtemp);
    MCEstimate fit1=m_obs->getEstimate(randtemp);
    upperfit[0].xval=AFD.m_imomsq[kmin];
    upperfit[0].yval=fit1.getFullEstimate()+fit1.getSymmetricError();
    lowerfit[0].xval=AFD.m_imomsq[kmin];
    lowerfit[0].yval=fit1.getFullEstimate()-fit1.getSymmetricError();
    m_obs->eraseSamplings(randtemp);
    doDispersionBySamplings(*m_obs,AFD.getAnisotropyKey(),AFD.getRestMassSquaredKey(), 
                            AFD.m_momsq_quantum*AFD.m_imomsq[kmax],randtemp);
    MCEstimate fit2=m_obs->getEstimate(randtemp);
    upperfit[1].xval=AFD.m_imomsq[kmax];
    upperfit[1].yval=fit2.getFullEstimate()+fit2.getSymmetricError();
    lowerfit[1].xval=AFD.m_imomsq[kmax];
    lowerfit[1].yval=fit2.getFullEstimate()-fit2.getSymmetricError();
    m_obs->eraseSamplings(randtemp);

    createEnergyDispersionPlot(Esq,xiest.getFullEstimate(),xiest.getSymmetricError(),
                               goodtype,goodness,particlename,lowerfit,upperfit,
                               plotfile,symboltype,symbolcolor);
    }
    catch(const std::exception& errmsg){
       //xmlout.clear();
       xmlout.put_child("Error",string("DoFit with type AnisotropyFromDispersion encountered an error: ")
               +string(errmsg.what()));
      // throw(std::invalid_argument("DoFit with type TemporalCorrelator encountered an error: ")+errmsg);}
    }}
}


// ***************************************************************************************
 

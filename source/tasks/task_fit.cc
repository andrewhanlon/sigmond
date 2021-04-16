#include "task_handler.h"
#include "chisq_anisotropy.h"
#include "chisq_tcorr.h"
#include "chisq_fit.h"
#include "create_plots.h"
#include "plot_info.h"
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
// *       <Uncorrelated/>  (optional: performs an uncorrelated fit)             *
// *       <TemporalCorrelatorFit>                                               *
// *         <Operator>.... </Operator>                                          *
// *         <SubtractVEV/>             (as appropriate)                         *
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
// *                <Name>Amp</Name><IDIndex>0</IDIndex>                         *
// *             </Amplitude>                                                    *
// *         </Model>                                                            *
// *         <Priors>                                                            *
// *             <Energy>                                                        *
// *                <Name>pion_prior</Name><IDIndex>0</IDIndex>                  *
// *             </Energy>                                                       *
// *             <Amplitude>                                                     *
// *                <Name>Amp_prior</Name><IDIndex>0</IDIndex>                   *
// *             </Amplitude>                                                    *
// *         </Priors>                                                           *
// *       </TemporalCorrelatorFit>                                              *
// *       <DoEffectiveEnergyPlot> (optional)                                    *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <CorrName>standard</CorrName>   (optional)                         *
// *          <TimeStep>3</TimeStep>  (optional: 1 default)                      *
// *          <SymbolColor> ... </SymbolColor>                                   *
// *          <SymbolType> ... </SymbolType>                                     *
// *          <MaxRelativeErrorToPlot> ...</MaxRelativeErrorToPlot> (optional)   *
// *          <Goodness>qual</Goodness>  "qual" or "chisq"                       *
// *          <ShowApproach/>   (optional)                                       *
// *          <ReferenceEnergy> (optional: includes energy ratio on plot)        *
// *            <Name>kaon</Name><IDIndex>0</IDIndex>                            *
// *          </ReferenceEnergy>                                                 *
// *       </DoEffectiveEnergyPlot>                                              *
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
// *       <Uncorrelated/>  (optional: performs an uncorrelated fit)             *
// *       <TwoTemporalCorrelatorFit>                                            *
// *         <CorrelatorOne>                                                     *
// *           <Operator>.... </Operator>                                        *
// *           <SubtractVEV/>             (as appropriate)                       *
// *           <MinimumTimeSeparation>3</MinimumTimeSeparation>                  *
// *           <MaximumTimeSeparation>12</MaximumTimeSeparation>                 *
// *           <ExcludeTimes>4 8</ExcludeTimes>  (optional)                      *
// *           <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>                  *
// *           <Model>...</Model>                                                *
// *           <Priors>...</Priors>                                              *
// *         </CorrelatorOne>                                                    *
// *         <CorrelatorTwo>                                                     *
// *           <Operator>.... </Operator>                                        *
// *           <SubtractVEV/>             (as appropriate)                       *
// *           <MinimumTimeSeparation>3</MinimumTimeSeparation>                  *
// *           <MaximumTimeSeparation>12</MaximumTimeSeparation>                 *
// *           <ExcludeTimes>4 8</ExcludeTimes>  (optional)                      *
// *           <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>                  *
// *           <Model>...</Model>                                                *
// *           <Priors>...</Priors>                                              *
// *         </CorrelatorTwo>                                                    *
// *         <EnergyRatio>                                                       *
// *            <Name>pion</Name><IDIndex>0</IDIndex>                            *
// *         </EnergyRatio>                                                      *
// *       </TwoTemporalCorrelatorFit>                                           *
// *       <DoEffectiveEnergyPlot>     (plot correlator 1 only)                  *
// *         <PlotFile> ... </PlotFile>                                          *
// *         <CorrName>standard</CorrName>   (optional)                          *
// *         <TimeStep>3</TimeStep>  (optional: 1 default)                       *
// *         <SymbolColor> ... </SymbolColor>                                    *
// *         <SymbolType> ... </SymbolType>                                      *
// *         <MaxRelativeErrorToPlot> ...</MaxRelativeErrorToPlot> (optional)    *
// *         <Goodness>qual</Goodness>  "qual" or "chisq"                        *
// *         <ShowApproach/>   (optional)                                        *
// *       </DoEffectiveEnergyPlot>                                              *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoFit</Action>                                                  *
// *       <Type>TemporalCorrelatorTminVary</Type>                               *
// *       <MinimizerInfo>                 (optional)                            *
// *         <Method>Minuit2</Method>                                            *
// *         <ParameterRelTol>1e-6</ParameterRelTol>                             *
// *         <ChiSquareRelTol>1e-4</ChiSquareRelTol>                             *
// *         <MaximumIterations>1024</MaximumIterations>                         *
// *         <Verbosity>Low</Verbosity>                                          *
// *       </MinimizerInfo>                                                      *
// *       <SamplingMode>Bootstrap</SamplingMode>   (optional)                   *
// *       <CovMatCalcSamplingMode>Bootstrap</CovMatCalcSamplingMode> (optional) *
// *       <Uncorrelated/>  (optional) performs an uncorrelated fit              *
// *       <TemporalCorrelatorTminVaryFit>                                       *
// *         <Operator>.... </Operator>                                          *
// *         <SubtractVEV/>             (as appropriate)                         *
// *         <TminFirst>3</TminFirst>                                            *
// *         <TminLast>3</TminLast>                                              *
// *         <Tmax>12</Tmax>                                                     *
// *         <ExcludeTimes>4 8</ExcludeTimes>  (optional)                        *
// *         <Model>                                                             *
// *             <Type>TimeSymSingleExponential</Type>                           *
// *             <Energy>                                                        *
// *                <Name>pion</Name><IDIndex>0</IDIndex> // default taskcount   *
// *             </Energy>                                                       *
// *             <Amplitude>                                                     *
// *                <Name>Amp</Name><IDIndex>0</IDIndex>                         *
// *             </Amplitude>                                                    *
// *         </Model>                                                            *
// *         <Priors>...</Priors>                                                *
// *       </TemporalCorrelatorTminVaryFit>                                      *
// *       <Shift>                                                               *
// *         <Subtract/>    (optional: subtracts shift rather than adds)         *
// *         <ShiftInfo>                                                         *
// *           <Name>shift_obsname</Name>                                        *
// *           <IDIndex>0</IDIndex>                                              *
// *         </ShiftInfo>                                                        *
// *       </Shift>                                                              *
// *       <EnergyLevel>1</EnergyLevel> (0 default)                              *
// *       <DoPlot>                                                              *
// *         <PlotFile> ... </PlotFile>                                          *
// *         <CorrName>standard</CorrName>   (optional)                          *
// *         <SymbolType> ... </SymbolType>                                      *
// *         <GoodFitSymbolColor> ... </GoodFitSymbolColor>                      *
// *         <BadFitSymbolColor> ... </BadFitSymbolColor>                        *
// *         <CorrelatedFitSymbolHollow/>  (optional)                            *
// *         <UncorrelatedFitSymbolHollow/>  (optional)                          *
// *         <QualityThreshold>qual</QualityThreshold>  (0.1 default)            *
// *         <CorrelatedThreshold>1.2</CorrelatedThreshold>  (1.0 default)       *
// *         <ChosenFitInfo>           (optional)                                *
// *           <Name>fit_obsname</Name>                                          *
// *           <IDIndex>0</IDIndex>                                              *
// *         </ChosenFitInfo>                                                    *
// *       </DoPlot>                                                             *
// *       <DoShiftPlot>                                                         *
// *         <PlotFile> ... </PlotFile>                                          *
// *         <EnergyLevel>1</EnergyLevel>   (0 default)                          *
// *         <CorrName>standard</CorrName>   (optional)                          *
// *         <SymbolType> ... </SymbolType>                                      *
// *         <GoodFitSymbolColor> ... </GoodFitSymbolColor>                      *
// *         <BadFitSymbolColor> ... </BadFitSymbolColor>                        *
// *         <CorrelatedFitSymbolHollow/>  (optional)                            *
// *         <UncorrelatedFitSymbolHollow/>  (optional)                          *
// *         <QualityThreshold>qual</QualityThreshold>  (0.1 default)            *
// *         <CorrelatedThreshold>1.2</CorrelatedThreshold>  (1.0 default)       *
// *         <ChosenFitInfo>           (optional)                                *
// *           <Name>fit_obsname</Name>                                          *
// *           <IDIndex>0</IDIndex>                                              *
// *         </ChosenFitInfo>                                                    *
// *       </DoShiftPlot>                                                        *
// *    </Task>                                                                  *
// *                                                                             *
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
// *       <Uncorrelated/>  (optional) performs an uncorrelated fit              *
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
  string fittype;
  xmlreadchild(xmltask, "Type", fittype, "DoFit");
  xmlout.set_root("DoFit");
  xmlout.put_child("Type", fittype); 

  ChiSquareMinimizerInfo mz_info;  // default minimizer info
  if (xmltask.count_among_children("MinimizerInfo") > 0) {
    ChiSquareMinimizerInfo mz_user(xmltask);
    mz_info = mz_user;
  }
  XMLHandler xmlmz;
  mz_info.output(xmlmz);
  xmlout.put_child(xmlmz);
 
  m_obs->setToDefaultSamplingMode();
  m_obs->setCovMatToDefaultSamplingMode();
  SamplingMode mode = m_obs->getCurrentSamplingMode();
  string instr;
  if (xmlreadifchild(xmltask, "SamplingMode", instr)) {
    if (instr=="Bootstrap") mode = Bootstrap;
    else if (instr=="Jackknife") mode = Jackknife;
    else throw(std::invalid_argument("Bad sampling mode"));
    if (mode == Bootstrap) {
       xmlout.put_child("SamplingMode", "Bootstrap");
       m_obs->setToBootstrapMode();
       m_obs->setCovMatToBootstrapMode();
    }
    else {
       xmlout.put_child("SamplingMode", "Jackknife");
       m_obs->setToJackknifeMode();
       m_obs->setCovMatToJackknifeMode();
    }
  }
  SamplingMode covcalcmode = m_obs->getCovMatCurrentSamplingMode();
  instr.clear();
  if (xmlreadifchild(xmltask, "CovMatCalcSamplingMode", instr)) {
    if (instr == "Bootstrap") covcalcmode = Bootstrap;
    else if (instr == "Jackknife") covcalcmode = Jackknife;
    else throw(std::invalid_argument("Bad cov mat calc sampling mode"));
    if (covcalcmode == Bootstrap) {
      xmlout.put_child("CovMatCalcSamplingMode", "Bootstrap");
      m_obs->setCovMatToBootstrapMode();
    }
    else {
      xmlout.put_child("CovMatCalcSamplingMode", "Jackknife");
      m_obs->setCovMatToJackknifeMode();
    }
  }

  bool correlated = (xmltask.count_among_children("Uncorrelated") > 0) ? false: true;

  vector<MCEstimate> bestfit_params;

  if (fittype == "TemporalCorrelator") {
    try {
      XMLHandler xmlf(xmltask, "TemporalCorrelatorFit");
      RealTemporalCorrelatorFit RTC(xmlf, *m_obs, taskcount);
      XMLHandler xmlof;
      RTC.output(xmlof);
      xmlout.put_child(xmlof);
      FitResult fit_result = doChiSquareFitting(RTC, mz_info, correlated, xmlout);

      if (xmltask.count_among_children("DoEffectiveEnergyPlot") > 0) {
        XMLHandler xmlp(xmltask, "DoEffectiveEnergyPlot");
        FitEffEnergyPlotInfo plot_info(xmlp);
        makeFitPlot(plot_info, RTC, fit_result, m_obs, xmlout);
      }
    }
    catch(const std::exception& errmsg) {
      xmlout.put_child("Error", string("DoFit with type TemporalCorrelator encountered an error: ")
          +string(errmsg.what()));
    }
  }

  /*
  else if (fittype=="TwoTemporalCorrelator") {
    try{
      XMLHandler xmlf(xmltask,"TwoTemporalCorrelatorFit");
      TwoRealTemporalCorrelatorFit RTC(xmlf, *m_obs, taskcount);
      XMLHandler xmlof;
      RTC.output(xmlof);
      xmlout.put_child(xmlof);
      double chisq_dof,qual;
      FitResult fit_result = doChiSquareFitting(RTC, mz_info, correlated, xmlout);

      int nparam1 = RTC.m_model1_ptr->getNumberOfParams();
      doRatioBySamplings(*m_obs, RTC.m_fitparam_info[0], RTC.m_fitparam_info[nparam1], RTC.m_energyratio);
      MCEstimate ratioest=m_obs->getEstimate(RTC.m_energyratio);
      XMLHandler xmlrat("EnergyRatioFitResult");
      XMLHandler xmlrr;
      ratioest.output(xmlrr);
      xmlrat.put_child(xmlrr);
      xmlout.put_child(xmlrat);

      if (xmltask.count_among_children("DoEffectiveEnergyPlot") > 0) {
        XMLHandler xmlp(xmltask, "DoEffectiveEnergyPlot");
        FitEffEnergyPlotInfo plot_info(xmlp);
        makePlot(plot_info, RTC, fit_result, m_obs, xmlout);
      }
    }
    catch(const std::exception& errmsg) {
      xmlout.put_child("Error",string("DoFit with type TwoTemporalCorrelator encountered an error: ")
          +string(errmsg.what()));
    }
  }

  else if (fittype=="TemporalCorrelatorTminVary") {
    try {
      XMLHandler xmlf(xmltask,"TemporalCorrelatorTminVaryFit");
      uint tminfirst, tminlast, tmax;
      xmlread(xmlf,"TminFirst", tminfirst, "TemporalCorrelatorTminVary");
      xmlread(xmlf,"TminLast", tminlast, "TemporalCorrelatorTminVary");
      xmlread(xmlf,"Tmax", tmax, "TemporalCorrelatorTminVary");
      XMLHandler xmltf(xmlf, XMLHandler::subtree_copy);
      xmltf.rename_tag("TemporalCorrelatorFit");
      xmltf.put_child("MinimumTimeSeparation", make_string(tminfirst));
      xmltf.put_child("MaximumTimeSeparation", make_string(tmax));

      vector<FitResult> fit_results;
      for (uint tmin = tminfirst; tmin <= tminlast; ++tmin) {
        xmltf.seek_unique("MinimumTimeSeparation");
        xmltf.seek_next_node();       
        xmltf.set_text_content(make_string(tmin)); 
        try {
          RealTemporalCorrelatorFit RTC(xmltf, *m_obs, taskcount);
          const vector<uint>& tvalues = RTC.getTvalues();
          if (find(tvalues.begin(), tvalues.end(), tmin) == tvalues.end()) continue;
          const vector<MCObsInfo>& fitparam_infos = RTC.getFitParamInfos();
          for (uint k = 0; k < fitparam_infos.size(); ++k) {
            m_obs->eraseSamplings(fitparam_infos[k]);
          }
          XMLHandler xmlof;
          RTC.output(xmlof);
          xmlof.rename_tag("TemporalCorrelatorTminVaryFit");
          xmlout.put_child(xmlof);
          FitResult fit_result = doChiSquareFitting(RTC, mz_info, correlated, xmlout);
        }
        catch(const std::exception& errmsg) {
            xmlout.put_child("Error", string("DoFit within type TemporalCorrelator encountered an error: ")
                + string(errmsg.what()));
        }
      }

      if (xmltask.count_to_among_children("DoPlot") > 0) {
        XMLHandler xmlp(xmltask, "DoPlot");
        TminFitPlotInfo plot_info(xmlp);
        makePlot(plot_info, rtcs, fit_results, m_obs, xmlout);
      }
    }
    catch(const std::exception& errmsg) {
      xmlout.put_child("Error", string("DoFit with type TemporalCorrelatorTminVary encountered an error: ")
          + string(errmsg.what()));
    }
  }
  */

  /*
  else if (fittype=="AnisotropyFromDispersion"){
    try{
      XMLHandler xmlf(xmltask,"AnisotropyFromDispersionFit");
      AnisotropyFromDispersionFit AFD(xmlf,*m_obs,taskcount);
      XMLHandler xmlof; AFD.output(xmlof);
      xmlout.put_child(xmlof);
      FitResult fit_result = doChiSquareFitting(AFD, mz_info, correlated, xmlout);

      if (xmltask.count_among_children("DoPlot")>0) {
        XMLHandler xmlp(xmltask, "DoPlot");
        AnisotropyFitPlotInfo plot_info(xmlp, *m_obs, taskcount);
        make_plot(plot_info);
      }
    }
    catch(const std::exception& errmsg) {
      xmlout.put_child("Error",string("DoFit with type LogTemporalCorrelator encountered an error: ")
          +string(errmsg.what()));
    }
  }
  */
}
// ***************************************************************************************
 

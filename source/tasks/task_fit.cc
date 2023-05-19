#include "task_handler.h"
#include "chisq_anisotropy.h"
#include "chisq_disp.h"
#include "chisq_tcorr.h"
#include "chisq_fit.h"
#include "chisq_logtcorr.h"
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
// *       <Uncorrelated/>  (optional) performs an uncorrelated fit              *
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
// *       </TemporalCorrelatorFit>                                              *
// *    </Task>                                                                  *
// *                                                                             *
// *    Fit to the log of a single real-valued temporal correlator.              *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoFit</Action>                                                  *
// *       <Type>LogTemporalCorrelator</Type>                                    *
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
// *       <LogTemporalCorrelatorFit>                                            *
// *         <Operator>.... </Operator>                                          *
// *         <SubtractVEV/>             (as appropriate)                         *
// *         <MinimumTimeSeparation>3</MinimumTimeSeparation>                    *
// *         <MaximumTimeSeparation>12</MaximumTimeSeparation>                   *
// *         <ExcludeTimes>4 8</ExcludeTimes>  (optional)                        *
// *         <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>                    *
// *         <LogModel>                                                          *
// *             <Type>LogTimeSymSingleExponential</Type>                        *
// *             <Energy>                                                        *
// *                <Name>pion</Name><IDIndex>0</IDIndex> // default taskcount   *
// *             </Energy>                                                       *
// *             <LogAmplitude>                                                  *
// *                <Name>Amp</Name><IDIndex>0</IDIndex>                         *
// *             </LogAmplitude>                                                 *
// *         </LogModel>                                                         *
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
// *       </LogTemporalCorrelatorFit>                                           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *                                                                             *
// *      For extracting the `interaction energy' or energy difference           *
// *      between an interacting energy level and the nearest                    *
// *      non-interacting state, the ratio:                                      *
// *                               C_int(t)                                      *
// *               R(t) =  -------------------------                             *
// *                         prod_i C_non-int[i](t)                              *
// *      can be fit to the ansatz:  R(t) = A exp(-DeltaE*t), where:             *
// *        DeltaE          =  E_int - E_non-int,                                *
// *        C_int(t)        =  (diagonal) correlator for interacting level       *
// *        C_non-int[i](t) =  N correlators coresponding to closest             *
// *                           non-interacting level                             *
// *       (other fit ansatz are allowed as well, such as two exponentials)      *
// *                                                                             *
// *      Example:                                                               *
// *      if the closest non-interacting level to C_int is pi(0)K(1)K(1):        *
// *        C_non-int[0](t) = pion correlator for P^2=0                          *
// *        C_non-int[1](t) = kaon correlator for P^2=1                          *
// *        C_non-int[2](t) = kaon correlator for P^2=1                          *
// *                                                                             *
// *      Note: the resultant correlator will have (if specified) all            *
// *            VEVs subtracted in this task.                                    *
// *    <Task>                                                                   *
// *     <Action>DoFit</Action>                                                  *
// *       <Type>TemporalCorrelatorInteractionRatio</Type>                       *
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
// *       <TemporalCorrelatorInteractionRatioFit>                               *
// *         <Ratio>                                                             *
// *            <Operator>...</Operator>                                         *
// *         </Ratio>                                                            *
// *         <InteractingOperator>                                               *
// *            <Operator>...</Operator>                                         *
// *            <SubtractVEV />    (optional)                                    *
// *         </InteractingOperator>                                              *
// *         <NonInteractingOperator>                                            *
// *            <Operator>...</Operator>                                         *
// *            <SubtractVEV />    (optional)                                    *
// *         </NonInteractingOperator>                                           *
// *          ......                                                             *
// *         <NonInteractingOperator>                                            *
// *            <Operator>...</Operator>                                         *
// *            <SubtractVEV />    (optional)                                    *
// *         </NonInteractingOperator>                                           *
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
// *         <DoEffectiveEnergyPlot> (optional)                                  *
// *            <PlotFile> ... </PlotFile>                                       *
// *            <CorrName>standard</CorrName>   (optional)                       *
// *            <TimeStep>3</TimeStep>  (optional: 1 default)                    *
// *            <SymbolColor> ... </SymbolColor>                                 *
// *            <SymbolType> ... </SymbolType>                                   *
// *            <MaxRelativeErrorToPlot> ...</MaxRelativeErrorToPlot> (optional) *
// *            <Goodness>qual</Goodness>  "qual" or "chisq"                     *
// *            <ShowApproach/>   (optional)                                     *
// *            <ReferenceEnergy> (optional: includes energy ratio on plot)      *
// *              <Name>kaon</Name><IDIndex>0</IDIndex>                          *
// *            </ReferenceEnergy>                                               *
// *         </DoEffectiveEnergyPlot>                                            *
// *       </TemporalCorrelatorInteractionRatioFit>                              *
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
// *       <Uncorrelated/>  (optional) performs an uncorrelated fit              *
// *       <TwoTemporalCorrelatorFit>                                            *
// *         <CorrelatorOne>                                                     *
// *           <Operator>.... </Operator>                                        *
// *           <SubtractVEV/>             (as appropriate)                       *
// *           <MinimumTimeSeparation>3</MinimumTimeSeparation>                  *
// *           <MaximumTimeSeparation>12</MaximumTimeSeparation>                 *
// *           <ExcludeTimes>4 8</ExcludeTimes>  (optional)                      *
// *           <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>                  *
// *           <Model>...</Model>                                                *
// *         </CorrelatorOne>                                                    *
// *         <CorrelatorTwo>                                                     *
// *           <Operator>.... </Operator>                                        *
// *           <SubtractVEV/>             (as appropriate)                       *
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
// *           <MaxRelativeErrorToPlot> ...</MaxRelativeErrorToPlot> (optional)  *
// *           <Goodness>qual</Goodness>  "qual" or "chisq"                      *
// *           <ShowApproach/>   (optional)                                      *
// *         </DoEffectiveEnergyPlot>                                            *
// *       </TwoTemporalCorrelatorFit>                                           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoFit</Action>                                                  *
// *       <Type>NSimTemporalCorrelator</Type>                                   *
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
// *       <NSimTemporalCorrelatorFit>                                           *
// *         <Fits>                                                              *
// *         <TemporalCorrelatorFit>                                             *
// *           <Operator>.... </Operator>                                        *
// *           <SubtractVEV/>             (as appropriate)                       *
// *           <MinimumTimeSeparation>3</MinimumTimeSeparation>                  *
// *           <MaximumTimeSeparation>12</MaximumTimeSeparation>                 *
// *           <ExcludeTimes>4 8</ExcludeTimes>  (optional)                      *
// *           <LargeTimeNoiseCutoff>1.0</LargeTimeNoiseCutoff>                  *
// *           <Model> ... </Model>                                              *
// *           <DoEffectiveEnergyPlot> (optional)                                *
// *             <PlotFile> ... </PlotFile>                                      *
// *             <CorrName>standard</CorrName>   (optional)                      *
// *             <TimeStep>3</TimeStep>  (optional: 1 default)                   *
// *             <SymbolColor> ... </SymbolColor>                                *
// *             <SymbolType> ... </SymbolType>                                  *
// *             <MaxRelativeErrorToPlot> ...</MaxRelativeErrorToPlot> (optional)*
// *             <Goodness>qual</Goodness>  "qual" or "chisq"                    *
// *             <ShowApproach/>   (optional)                                    *
// *             <ReferenceEnergy> (optional: includes energy ratio on plot)     *
// *               <Name>kaon</Name><IDIndex>0</IDIndex>                         *
// *             </ReferenceEnergy>                                              *
// *           </DoEffectiveEnergyPlot>                                        *
// *         </TemporalCorrelatorFit>                                            *
// *         <TemporalCorrelatorFit>                                             *
// *              ...                                                            *
// *         </TemporalCorrelatorFit>                                            *
// *              ... include as many TemporalCorrelatorFits here as desired     *
// *         </Fits>                                                             *
// *       </NSimTemporalCorrelatorFit>                                          *
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
// *       </TemporalCorrelatorTminVaryFit>                                      *
// *       <DoEnergyDifference>                                                  *
// *         <SpatialExtentNumSites>32</SpatialExtentNumSites>                   *
// *         <Anisotropy>       (optional)                                       *
// *            <Name>aniso_fit_name</Name>                                      *
// *            <IDIndex>0</IDIndex>                                             *
// *         </Anisotropy>                                                       *
// *         <ScatteringParticleEnergyFit>                                       *
// *            <IntMomSquared>4</IntMomSquared>                                 *
// *            <Name>scatting_part_atrest_energy_fit_name</Name>                *
// *            <IDIndex>0</IDIndex>                                             *
// *         </ScatteringParticleEnergyFit>                                      *
// *         <ScatteringParticleEnergyFit>                                       *
// *            <IntMomSquared>2</IntMomSquared>                                 *
// *            <Name>scatting_part_atrest_energy_fit_name</Name>                *
// *            <IDIndex>0</IDIndex>                                             *
// *         </ScatteringParticleEnergyFit>                                      *
// *            ...                                                              *
// *       </DoEnergyDifference>                                                 *
// *       <ChosenFitInfo>                                                       *
// *         <Name>fit_obsname</Name>                                            *
// *         <IDIndex>0</IDIndex>                                                *
// *       </ChosenFitInfo>                                                      *
// *       <PlotInfo>                                                            *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <CorrName>standard</CorrName>   (optional)                         *
// *          <SymbolType> ... </SymbolType>                                     *
// *          <GoodFitSymbolColor> ... </GoodFitSymbolColor>                     *
// *          <BadFitSymbolColor> ... </BadFitSymbolColor>                       *
// *          <CorrelatedFitSymbolHollow/>  (optional)                           *
// *          <UncorrelatedFitSymbolHollow/>  (optional)                         *
// *          <QualityThreshold>qual</QualityThreshold>  (0.1 default)           *
// *          <CorrelatedThreshold>1.2</CorrelatedThreshold>  (1.0 default)      *
// *       </PlotInfo>                                                           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoFit</Action>                                                  *
// *       <Type>LogTemporalCorrelatorTminVary</Type>                            *
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
// *       <LogTemporalCorrelatorTminVaryFit>                                    *
// *         <Operator>.... </Operator>                                          *
// *         <SubtractVEV/>             (as appropriate)                         *
// *         <TminFirst>3</TminFirst>                                            *
// *         <TminLast>3</TminLast>                                              *
// *         <Tmax>12</Tmax>                                                     *
// *         <ExcludeTimes>4 8</ExcludeTimes>  (optional)                        *
// *         <LogModel>                                                          *
// *             <Type>TimeSymSingleExponential</Type>                           *
// *             <Energy>                                                        *
// *                <Name>pion</Name><IDIndex>0</IDIndex> // default taskcount   *
// *             </Energy>                                                       *
// *             <LogAmplitude>                                                  *
// *                <Name>Amp</Name><IDIndex>0</IDIndex>                         *
// *             </LogAmplitude>                                                 *
// *         </LogModel>                                                         *
// *       </LogTemporalCorrelatorTminVaryFit>                                   *
// *       <DoEnergyDifference>                                                  *
// *         <SpatialExtentNumSites>32</SpatialExtentNumSites>                   *
// *         <Anisotropy>       (optional)                                       *
// *            <Name>aniso_fit_name</Name>                                      *
// *            <IDIndex>0</IDIndex>                                             *
// *         </Anisotropy>                                                       *
// *         <ScatteringParticleEnergyFit>                                       *
// *            <IntMomSquared>4</IntMomSquared>                                 *
// *            <Name>scatting_part_atrest_energy_fit_name</Name>                *
// *            <IDIndex>0</IDIndex>                                             *
// *         </ScatteringParticleEnergyFit>                                      *
// *         <ScatteringParticleEnergyFit>                                       *
// *            <IntMomSquared>2</IntMomSquared>                                 *
// *            <Name>scatting_part_atrest_energy_fit_name</Name>                *
// *            <IDIndex>0</IDIndex>                                             *
// *         </ScatteringParticleEnergyFit>                                      *
// *            ...                                                              *
// *       </DoEnergyDifference>                                                 *
// *       <ChosenFitInfo>                                                       *
// *         <Name>fit_obsname</Name>                                            *
// *         <IDIndex>0</IDIndex>                                                *
// *       </ChosenFitInfo>                                                      *
// *       <PlotInfo>                                                            *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <CorrName>standard</CorrName>   (optional)                         *
// *          <SymbolType> ... </SymbolType>                                     *
// *          <GoodFitSymbolColor> ... </GoodFitSymbolColor>                     *
// *          <BadFitSymbolColor> ... </BadFitSymbolColor>                       *
// *          <CorrelatedFitSymbolHollow/>  (optional)                           *
// *          <UncorrelatedFitSymbolHollow/>  (optional)                         *
// *          <QualityThreshold>qual</QualityThreshold>  (0.1 default)           *
// *          <CorrelatedThreshold>1.2</CorrelatedThreshold>  (1.0 default)      *
// *       </PlotInfo>                                                           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoFit</Action>                                                  *
// *       <Type>TemporalCorrelatorInteractionRatioTminVary</Type>               *
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
// *       <TemporalCorrelatorInteractionRatioTminVaryFit>                       *
// *         <Ratio>                                                             *
// *            <Operator>...</Operator>                                         *
// *         </Ratio>                                                            *
// *         <InteractingOperator>                                               *
// *            <Operator>...</Operator>                                         *
// *            <SubtractVEV />    (optional)                                    *
// *         </InteractingOperator>                                              *
// *         <NonInteractingOperator>                                            *
// *            <Operator>...</Operator>                                         *
// *            <SubtractVEV />    (optional)                                    *
// *         </NonInteractingOperator>                                           *
// *          ......                                                             *
// *         <NonInteractingOperator>                                            *
// *            <Operator>...</Operator>                                         *
// *            <SubtractVEV />    (optional)                                    *
// *         </NonInteractingOperator>                                           *
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
// *       </TemporalCorrelatorTminVaryFit>                                      *
// *       <DoReconstructEnergy>                                                 *
// *         <SpatialExtentNumSites>32</SpatialExtentNumSites>                   *
// *         <Anisotropy>       (optional)                                       *
// *            <Name>aniso_fit_name</Name>                                      *
// *            <IDIndex>0</IDIndex>                                             *
// *         </Anisotropy>                                                       *
// *         <ScatteringParticleEnergyFit>                                       *
// *            <IntMomSquared>4</IntMomSquared>                                 *
// *            <Name>scatting_part_atrest_energy_fit_name</Name>                *
// *            <IDIndex>0</IDIndex>                                             *
// *         </ScatteringParticleEnergyFit>                                      *
// *         <ScatteringParticleEnergyFit>                                       *
// *            <IntMomSquared>2</IntMomSquared>                                 *
// *            <Name>scatting_part_atrest_energy_fit_name</Name>                *
// *            <IDIndex>0</IDIndex>                                             *
// *         </ScatteringParticleEnergyFit>                                      *
// *            ...                                                              *
// *       </DoReconstructEnergy>                                                *
// *       <ChosenFitInfo>                                                       *
// *         <Name>fit_obsname</Name>                                            *
// *         <IDIndex>0</IDIndex>                                                *
// *       </ChosenFitInfo>                                                      *
// *       <PlotInfo>                                                            *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <CorrName>standard</CorrName>   (optional)                         *
// *          <SymbolType> ... </SymbolType>                                     *
// *          <GoodFitSymbolColor> ... </GoodFitSymbolColor>                     *
// *          <BadFitSymbolColor> ... </BadFitSymbolColor>                       *
// *          <CorrelatedFitSymbolHollow/>  (optional)                           *
// *          <UncorrelatedFitSymbolHollow/>  (optional)                         *
// *          <QualityThreshold>qual</QualityThreshold>  (0.1 default)           *
// *          <CorrelatedThreshold>1.2</CorrelatedThreshold>  (1.0 default)      *
// *       </PlotInfo>                                                           *
// *    </Task>                                                                  *
// *                                                                             *
// *                                                                             *
// *    Fit the free-particle energies squared for various three-momenta squared *
// *    to estimate the lattice anisotropy  a_s/a_t.                             *
// *    The model used for the observables is                                    *
// *                                                                             *
// *      (a_t E)^2 = restmass_sq + c * (2*Pi/Ns)^2 * nsq / xi^2                    *
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
// *       <Type>Dispersion</Type>                                               *
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
// *       <DispersionFit>                                                       *
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
// *         <Coefficient>                                                       *
// *           <Name>PionC</Name><IDIndex>0</IDIndex>                            *
// *         </Coefficient>                                                      *
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
// *       </DispersionFit>                                                      *
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
// *    Fit the free-particle energies squared for various three-momenta squared *
// *    to the lattice dispersion relation                                       *
// *                                                                             *
// *     sinh^2(a_t E) = sinh^2(restmass) + disp * sum_j sin^2(2*Pi n_j/Ns)      *
// *                                                                             *
// *    where "restmass" and "disp" are the two model parameters, and            *
// *    "Ns" is extent of the lattice in terms of number of sites                *
// *    in each of the three spatial directions, and "n" is the                  *
// *    vector of integers describing the momentum direction.                    *
// *    Recall that P_i = 2*Pi*n_i/L.                                            *
// *                                                                             *
// *    <Task>                                                                   *
// *     <Action>DoFit</Action>                                                  *
// *       <Type>LatticeDispersionRelation</Type>                                *
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
// *       <LatticeDispersionRelationFit>                                        *
// *         <SpatialExtentNumSites>24</SpatialExtentNumSites>                   *
// *         <Energy>                                                            *
// *           <Name>pion</Name><IDIndex>0</IDIndex>                             *
// *           <IntMom>0 0 0</IntMom>                                            *
// *         </Energy>                                                           *
// *         <Energy>                                                            *
// *           <Name>pion</Name><IDIndex>1</IDIndex>                             *
// *           <IntMom>0 0 1</IntMom>                                            *
// *         </Energy>                                                           *
// *         <Energy>... </Energy>                                               *
// *         <RestMass>                                                          *
// *           <Name>PionRestMassSquared</Name><IDIndex>0</IDIndex>              *
// *         </RestMass>                                                         *
// *         <Disp>                                                              *
// *           <Name>DispTerm</Name><IDIndex>0</IDIndex>                         *
// *         </Disp>                                                             *
// *         <DoPlot>                                                            *
// *           <PlotFile> ... </PlotFile>                                        *
// *           <ParticleName>pion</ParticleName>   (optional)                    *
// *           <SymbolColor> ... </SymbolColor>                                  *
// *           <SymbolType> ... </SymbolType>                                    *
// *           <Goodness>qual</Goodness>  "qual" or "chisq"                      *
// *         </DoPlot>                                                           *
// *       </LatticeDispersionRelationFit>                                       *
// *    </Task>                                                                  *
// *                                                                             *
// *******************************************************************************


void TaskHandler::doFit(XMLHandler& xmltask, XMLHandler& xmlout, int taskcount)
{
 string fittype;
 xmlreadchild(xmltask,"Type",fittype,"DoFit");
 xmlout.set_root("DoFit");
 xmlout.put_child("Type",fittype); 

 ChiSquareMinimizerInfo mz_info;  // default minimizer info
 if (xmltask.count_among_children("MinimizerInfo")>0){
    ChiSquareMinimizerInfo mz_user(xmltask);
    mz_info=mz_user;}
 XMLHandler xmlmz;
 mz_info.output(xmlmz);
 xmlout.put_child(xmlmz);
 
 m_obs->setToDefaultSamplingMode();
 m_obs->setCovMatToDefaultSamplingMode();
 SamplingMode mode=m_obs->getCurrentSamplingMode();
 string instr;
 if (xmlreadifchild(xmltask,"SamplingMode",instr)){
    if (instr=="Bootstrap") mode=Bootstrap;
    else if (instr=="Jackknife") mode=Jackknife;
    else throw(std::invalid_argument("Bad sampling mode"));
    if (mode==Bootstrap){
       xmlout.put_child("SamplingMode","Bootstrap");
       m_obs->setToBootstrapMode();
       m_obs->setCovMatToBootstrapMode();}
    else{
       xmlout.put_child("SamplingMode","Jackknife");
       m_obs->setToJackknifeMode();
       m_obs->setCovMatToJackknifeMode();}}
 SamplingMode covcalcmode=m_obs->getCovMatCurrentSamplingMode();
 instr.clear();
 if (xmlreadifchild(xmltask,"CovMatCalcSamplingMode",instr)){
    if (instr=="Bootstrap") covcalcmode=Bootstrap;
    else if (instr=="Jackknife") covcalcmode=Jackknife;
    else throw(std::invalid_argument("Bad cov mat calc sampling mode"));
    if (covcalcmode==Bootstrap){
       xmlout.put_child("CovMatCalcSamplingMode","Bootstrap");
       m_obs->setCovMatToBootstrapMode();}
    else{
       xmlout.put_child("CovMatCalcSamplingMode","Jackknife");
       m_obs->setCovMatToJackknifeMode();}}

 bool uncorrelated=(xmltask.count_among_children("Uncorrelated")>0) ? true: false;
 if (uncorrelated){
   m_obs->setToUnCorrelated();
   xmlout.put_child("Uncorrelated");}
 else
   m_obs->setToCorrelated();

 vector<MCEstimate> bestfit_params;

 if (fittype=="TemporalCorrelator"){
    try{
    XMLHandler xmlf(xmltask,"TemporalCorrelatorFit");
    RealTemporalCorrelatorFit RTC(xmlf,*m_obs,taskcount);
    XMLHandler xmlof; RTC.output(xmlof);
    xmlout.put_child(xmlof);
    double chisq_dof,qual;
    doChiSquareFitting(RTC,mz_info,chisq_dof,qual,bestfit_params,xmlout);

       // fit done, now do the plot if requested
    if (xmlf.count_among_children("DoEffectiveEnergyPlot")!=1) return;
    uint lattice_time_extent = getLatticeTimeExtent();
    RTC.plot( xmlf, taskcount, qual, chisq_dof, lattice_time_extent, bestfit_params, xmlout);
        
    }
    catch(const std::exception& errmsg){
       xmlout.put_child("Error",string("DoFit with type TemporalCorrelator encountered an error: ")
               +string(errmsg.what()));
    }}


 else if (fittype=="LogTemporalCorrelator"){
    try{
    XMLHandler xmlf(xmltask,"LogTemporalCorrelatorFit");
    LogRealTemporalCorrelatorFit RTC(xmlf,*m_obs,taskcount);
    XMLHandler xmlof; RTC.output(xmlof);
    xmlout.put_child(xmlof);
    double chisq_dof,qual;
    doChiSquareFitting(RTC,mz_info,chisq_dof,qual,bestfit_params,xmlout);

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
    string corrname;
    bool showapproach=(xml_child_tag_count(xmlp,"ShowApproach")>0);
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
    uint fit_tmin=RTC.getTmin();
    uint fit_tmax=RTC.getTmax();
    uint efftype=RTC.m_model_ptr->getEffMassType();
    double subt_const=0.0;
   // if (efftype>1){    // subtract fit constant
   //    efftype-=2;     // efftypes 2 and 3 remove constant, but noisy
   //    subt_const=bestfit_params[bestfit_params.size()-1].getFullEstimate();}
    SamplingMode mode=m_obs->getCurrentSamplingMode();

    map<double,MCEstimate> results;
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
    for (map<double,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++){
       XMLHandler xmlr("Estimate");
       xmlr.put_child("TimeSeparation",make_string(rt->first));
       xmlr.put_child("MeanValue",make_string((rt->second).getFullEstimate()));
       xmlr.put_child("SymmError",make_string((rt->second).getSymmetricError()));
       xmlef.put_sibling(xmlr);}
    xmlout.put_child(xmlef);
           // now prepare the plot
    //double maxrelerror=0.0;
    /* TODO: probably should change this
    if (xmlreadifchild(xmlp,"MaxRelativeErrorToPlot",maxrelerror)){
       map<double,MCEstimate> raw(results);
       results.clear();
       for (map<double,MCEstimate>::const_iterator it=raw.begin();it!=raw.end();it++)
          if ((it->second).getRelativeError()<std::abs(maxrelerror)) results.insert(*it);}
    */

    vector<XYDYPoint> meffvals(results.size());
    uint k=0;
    for (map<double,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
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
       xmlout.put_child("Error",string("DoFit with type LogTemporalCorrelator encountered an error: ")
               +string(errmsg.what()));
    }}


 else if (fittype=="TemporalCorrelatorInteractionRatio"){
    try{
     XMLHandler xmlf(xmltask,"TemporalCorrelatorInteractionRatioFit");
       
     XMLHandler xmltf(xmlf,XMLHandler::subtree_copy);
     bool erase_orig=true;
     setUpRatioFit( *m_obs, xmlf, xmltf, true, xmlout, erase_orig );
        
     RealTemporalCorrelatorFit RTC(xmltf,*m_obs,taskcount);
     XMLHandler xmlof; RTC.output(xmlof);
     xmlof.rename_tag("TemporalCorrelatorInteractionRatioFit");
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
     uint fit_tmin=RTC.getTmin();
     uint fit_tmax=RTC.getTmax();
     uint efftype=RTC.m_model_ptr->getEffMassType();
     double subt_const=0.0;
   // if (efftype>1){    // subtract fit constant
   //    efftype-=2;     // efftypes 2 and 3 remove constant, but noisy
   //    subt_const=bestfit_params[bestfit_params.size()-1].getFullEstimate();}
     SamplingMode mode=m_obs->getCurrentSamplingMode();

     map<double,MCEstimate> results;
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
     for (map<double,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++){
        XMLHandler xmlr("Estimate");
        xmlr.put_child("TimeSeparation",make_string(rt->first));
        xmlr.put_child("MeanValue",make_string((rt->second).getFullEstimate()));
        xmlr.put_child("SymmError",make_string((rt->second).getSymmetricError()));
        xmlef.put_sibling(xmlr);}
     xmlout.put_child(xmlef);
           // now prepare the plot
     //double maxrelerror=0.0;
    /* TODO: probably should change this
     if (xmlreadifchild(xmlp,"MaxRelativeErrorToPlot",maxrelerror)){
        map<double,MCEstimate> raw(results);
        results.clear();
        for (map<double,MCEstimate>::const_iterator it=raw.begin();it!=raw.end();it++)
           if ((it->second).getRelativeError()<std::abs(maxrelerror)) results.insert(*it);}
    */

     vector<XYDYPoint> meffvals(results.size());
     uint k=0;
     for (map<double,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
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
       xmlout.put_child("Error",string("DoFit with type TemporalCorrelator encountered an error: ")
               +string(errmsg.what()));
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
    uint efftype=RTC.m_model1_ptr->getEffMassType();
    double subt_const=0.0;
   // if (efftype>1){    // subtract fit constant
   //    efftype-=2;     // efftypes 2 and 3 remove constant, but noisy
   //    subt_const=bestfit_params[bestfit_params.size()-1].getFullEstimate();}
    SamplingMode mode=m_obs->getCurrentSamplingMode();

    map<double,MCEstimate> results;
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
    for (map<double,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++){
       XMLHandler xmlr("Estimate");
       xmlr.put_child("TimeSeparation",make_string(rt->first));
       xmlr.put_child("MeanValue",make_string((rt->second).getFullEstimate()));
       xmlr.put_child("SymmError",make_string((rt->second).getSymmetricError()));
       xmlef.put_sibling(xmlr);}
    xmlout.put_child(xmlef);
           // now prepare the plot
    //double maxrelerror=0.0;
    /* TODO: probably should change this
    if (xmlreadifchild(xmlp,"MaxRelativeErrorToPlot",maxrelerror)){
       map<double,MCEstimate> raw(results);
       results.clear();
       for (map<double,MCEstimate>::const_iterator it=raw.begin();it!=raw.end();it++)
          if ((it->second).getRelativeError()<std::abs(maxrelerror)) results.insert(*it);}
    */

    vector<XYDYPoint> meffvals(results.size());
    uint k=0;
    for (map<double,MCEstimate>::const_iterator rt=results.begin();rt!=results.end();rt++,k++){
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
       xmlout.put_child("Error",string("DoFit with type TwoTemporalCorrelator encountered an error: ")
               +string(errmsg.what()));
    }}

 //N Simultaneous Temporal Correlator Fit
 else if (fittype=="NSimTemporalCorrelator"){
    try{
    XMLHandler xmlf(xmltask,"NSimTemporalCorrelatorFit");
    NSimRealTemporalCorrelatorFit NSimRTC(xmlf,*m_obs,taskcount);
    XMLHandler xmlof; NSimRTC.output(xmlof);
    xmlout.put_child(xmlof);
    double chisq_dof,qual;
    doChiSquareFitting(NSimRTC,mz_info,chisq_dof,qual,
                       bestfit_params,xmlout);
        
    uint lattice_time_extent = getLatticeTimeExtent();
    
    XMLHandler xmlc(xmlf,"Fits");
    list<XMLHandler> xmlccs = xmlc.find("TemporalCorrelatorFit");
    list<XMLHandler> xmlccs2 = xmlc.find("TemporalCorrelatorInteractionRatioFit");
    uint i = 0;
    for(list<XMLHandler>::iterator it = xmlccs.begin(); it != xmlccs.end(); ++it){
        if (it->count_among_children("DoEffectiveEnergyPlot")==1){
            NSimRTC.plot( i, *it, taskcount, qual, chisq_dof, lattice_time_extent, bestfit_params, xmlout);
        }
        i++;
    }
    for(list<XMLHandler>::iterator it = xmlccs2.begin(); it != xmlccs2.end(); ++it){
        if (it->count_among_children("DoEffectiveEnergyPlot")==1){
            NSimRTC.plot( i, *it, taskcount, qual, chisq_dof, lattice_time_extent, bestfit_params, xmlout);
        }
        i++;
    }
        
    }
    catch(const std::exception& errmsg){
       xmlout.put_child("Error",string("DoFit with type NSimTemporalCorrelator encountered an error: ")
               +string(errmsg.what()));
    }}


 else if (fittype=="TemporalCorrelatorTminVary"){
    try{
    XMLHandler xmlf(xmltask,"TemporalCorrelatorTminVaryFit");
    uint tminfirst,tminlast,tmax;
    xmlread(xmlf,"TminFirst",tminfirst,"TemporalCorrelatorTminVary");
    xmlread(xmlf,"TminLast",tminlast,"TemporalCorrelatorTminVary");
    xmlread(xmlf,"Tmax",tmax,"TemporalCorrelatorTminVary");
    XMLHandler xmltf(xmlf,XMLHandler::subtree_copy);
    xmltf.rename_tag("TemporalCorrelatorFit");
    xmltf.put_child("MaximumTimeSeparation",make_string(tmax));
    xmltf.put_child("MinimumTimeSeparation",make_string(tminfirst));

    list<pair<MCObsInfo,double> > scattering_particles;
    MCObsInfo aniso_obsinfo;
    if (xmltask.count_among_children("DoEnergyDifference")==1){
       XMLHandler xmled(xmltask,"DoEnergyDifference");
       uint num_spatial_sites=0;
       xmlread(xmled,"SpatialExtentNumSites",num_spatial_sites,"TemporalCorrelatorTminVary");
       double m_momsq_quantum=6.2831853071795864770/double(num_spatial_sites);
       m_momsq_quantum*=m_momsq_quantum;
       if (xmled.count_to_among_children("Anisotropy")==1){
         XMLHandler xmlani(xmled,"Anisotropy");
         string name;
         xmlread(xmlani,"Name",name,"Anisotropy");
         int index=0;
         xmlreadifchild(xmlani,"IDIndex",index);
         aniso_obsinfo = MCObsInfo(name,index);}
       list<XMLHandler> scattering_xml=xmled.find_among_children("ScatteringParticleEnergyFit");
       for (list<XMLHandler>::iterator st=scattering_xml.begin();st!=scattering_xml.end();++st){
          uint psq;
          xmlread(*st,"IntMomSquared",psq,"ScatteringParticleEnergyFit");
          double psqfactor=psq*m_momsq_quantum;
          string name;
          xmlread(*st,"Name",name,"ScatteringParticleEnergyFit");
          int index=0;
          xmlreadifchild(*st,"IDIndex",index);
          scattering_particles.push_back(make_pair(MCObsInfo(name,index),psqfactor));}}

    MCObsInfo chosen_fit_info;
    if (xmltask.count_among_children("ChosenFitInfo")==1){
       XMLHandler xmlchosen(xmltask,"ChosenFitInfo");
       string name;
       xmlread(xmlchosen,"Name",name,"ChosenFitInfo");
       int index=0;
       xmlreadifchild(xmlchosen,"IDIndex",index);
       chosen_fit_info = MCObsInfo(name,index);}

    XMLHandler xmlp(xmltask,"PlotInfo");
    string plotfile;
    xmlread(xmlp,"PlotFile",plotfile,"TemporalCorrelatorTminVary");
    if (plotfile.empty()) throw(std::invalid_argument("Must have plot file name"));
    string corrname("standard");
    xmlreadif(xmlp,"CorrName",corrname,"TemporalCorrelatorTminVary");
    string symbol("circle");
    xmlreadif(xmlp,"SymbolType",symbol,"TemporalCorrelatorTminVary");
    double qualthreshold=0.1;
    xmlreadif(xmlp,"QualityThreshold",qualthreshold,"TemporalCorrelatorTminVary");
    string goodfitcolor("blue");
    xmlreadif(xmlp,"GoodFitSymbolColor",goodfitcolor,"TemporalCorrelatorTminVary");
    string badfitcolor("red");
    xmlreadif(xmlp,"BadFitSymbolColor",badfitcolor,"TemporalCorrelatorTminVary");
    double correlatedthreshold=1.0;
    xmlreadif(xmlp,"CorrelatedThreshold",correlatedthreshold,"TemporalCorrelatorTminVary");
    bool correlatedfit_hollow=false;
    if (xml_child_tag_count(xmlp,"CorrelatedFitSymbolHollow")>0) correlatedfit_hollow=true;
    bool uncorrelatedfit_hollow=false;
    if (xml_child_tag_count(xmlp,"UncorrelatedFitSymbolHollow")>0) uncorrelatedfit_hollow=true;
    vector<XYDYDYPoint> goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits;
    for (uint tmin=tminfirst;tmin<=tminlast;++tmin){
       xmltf.seek_unique("MinimumTimeSeparation");
       xmltf.seek_next_node();       
       xmltf.set_text_content(make_string(tmin)); 
       try{
       RealTemporalCorrelatorFit RTC(xmltf,*m_obs,taskcount);
       CorrelatorInfo corr(RTC.m_op,RTC.m_op);
       if (corrname=="standard") corrname=getCorrelatorStandardName(corr);
       const vector<uint>& tvalues=RTC.getTvalues();
       if (find(tvalues.begin(),tvalues.end(),tmin)==tvalues.end()) continue;
       int dof = tvalues.size() - RTC.m_model_ptr->getNumberOfParams();
       if (dof < 1) continue;
       const vector<MCObsInfo>& fitparam_infos=RTC.getFitParamInfos();
       for (uint k=0;k<fitparam_infos.size();++k)
          m_obs->eraseSamplings(fitparam_infos[k]);
       XMLHandler xmlof; RTC.output(xmlof);
       xmlof.rename_tag("TemporalCorrelatorTminVaryFit");
       xmlout.put_child(xmlof);
       double chisq_dof,qual;
       doChiSquareFitting(RTC,mz_info,chisq_dof,qual,bestfit_params,xmlout);
       TCorrFitInfo fitinfo;
       uint meff_tstep=1; bool showapproach=false;
       RTC.m_model_ptr->setFitInfo(RTC.m_fitparam_info,bestfit_params,tmin,tmax,
                                   showapproach,meff_tstep,chisq_dof,qual,fitinfo);
       MCObsInfo energy_key=fitinfo.energy_key;
       if (scattering_particles.size()>0){
          if (aniso_obsinfo.isVacuum()) // no anisotropy
            doEnergyDifferenceBySamplings(*m_obs,energy_key,scattering_particles,energy_key);
          else
            doEnergyDifferenceBySamplings(*m_obs,energy_key,aniso_obsinfo,scattering_particles,energy_key);}

       bool correlated=false;
       if (!chosen_fit_info.isVacuum()){
          MCObsInfo diff_obs;
          doCorrelatedDifferenceBySamplings(*m_obs,chosen_fit_info,energy_key,diff_obs);
          MCEstimate diff_est=m_obs->getEstimate(diff_obs);
          double diff_val=diff_est.getFullEstimate();
          double diff_up,diff_down;
          if (diff_est.isJackknifeMode())
             diff_up=diff_down=correlatedthreshold*diff_est.getSymmetricError();
          else{
             diff_up=correlatedthreshold*(diff_est.getUpperConfLimit()-diff_val);
             diff_down=correlatedthreshold*(diff_val-diff_est.getLowerConfLimit());}

          double upper_limit=diff_up+diff_val;
          double lower_limit=diff_val-diff_down;
          correlated = upper_limit >= 0. && lower_limit <= 0.;}

       MCEstimate energy=m_obs->getEstimate(energy_key);
       double y=energy.getFullEstimate();
       double dyup,dydn;
       if (energy.isJackknifeMode()) 
          dyup=dydn=energy.getSymmetricError();
       else{
          dyup=energy.getUpperConfLimit()-y;
          dydn=y-energy.getLowerConfLimit();}
       if (qual>=0.1 && correlated) goodcorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else if (qual>=0.1 && !correlated) gooduncorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else if (qual<0.1 && correlated) badcorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else baduncorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));}
       catch(const std::exception& errmsg){
          xmlout.put_child("Error",string("DoFit within type TemporalCorrelatorTminVary encountered an error: ")
                 +string(errmsg.what()));
       }}
    XMLHandler xmlplog("TminPlot");
    xmlplog.put_child("PlotFile",plotfile);
    xmlplog.put_child("QualityThreshold",make_string(qualthreshold));
    xmlplog.put_child("CorrelatedThreshold",make_string(correlatedthreshold));
    xmlplog.put_child("NumberOfGoodCorrelatedFitPoints",make_string(goodcorrelatedfits.size()));
    xmlplog.put_child("NumberOfGoodUncorrelatedFitPoints",make_string(gooduncorrelatedfits.size()));
    xmlplog.put_child("NumberOfBadCorrelatedFitPoints",make_string(badcorrelatedfits.size()));
    xmlplog.put_child("NumberOfBadUncorrelatedFitPoints",make_string(baduncorrelatedfits.size()));
    xmlout.put_child(xmlplog);

    XYDYDYPoint chosen_fit(0,0,0,0);
    if (!chosen_fit_info.isVacuum()){
       MCEstimate chosen_fit_estimate=m_obs->getEstimate(chosen_fit_info);
       double y=chosen_fit_estimate.getFullEstimate();
       double dyup,dydn;
       if (chosen_fit_estimate.isJackknifeMode()) 
          dyup=dydn=chosen_fit_estimate.getSymmetricError();
       else{
          dyup=chosen_fit_estimate.getUpperConfLimit()-y;
          dydn=y-chosen_fit_estimate.getLowerConfLimit();}
       chosen_fit = XYDYDYPoint(1,y,dyup,dydn);}
    createTMinPlot(goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits,
                   corrname,plotfile,symbol,goodfitcolor,badfitcolor,correlatedfit_hollow,
                   uncorrelatedfit_hollow,chosen_fit);}
    catch(const std::exception& errmsg){
       xmlout.put_child("Error",string("DoFit with type TemporalCorrelatorTminVary encountered an error: ")
               +string(errmsg.what()));
    }}

 else if (fittype=="LogTemporalCorrelatorTminVary"){
    try{
    XMLHandler xmlf(xmltask,"LogTemporalCorrelatorTminVaryFit");
    uint tminfirst,tminlast,tmax;
    xmlread(xmlf,"TminFirst",tminfirst,"LogTemporalCorrelatorTminVary");
    xmlread(xmlf,"TminLast",tminlast,"LogTemporalCorrelatorTminVary");
    xmlread(xmlf,"Tmax",tmax,"LogTemporalCorrelatorTminVary");
    XMLHandler xmltf(xmlf,XMLHandler::subtree_copy);
    xmltf.rename_tag("LogTemporalCorrelatorFit");
    xmltf.put_child("MaximumTimeSeparation",make_string(tmax));
    xmltf.put_child("MinimumTimeSeparation",make_string(tminfirst));

    list<pair<MCObsInfo,double> > scattering_particles;
    MCObsInfo aniso_obsinfo;
    if (xmltask.count_among_children("DoEnergyDifference")==1){
       XMLHandler xmled(xmltask,"DoEnergyDifference");
       uint num_spatial_sites=0;
       xmlread(xmled,"SpatialExtentNumSites",num_spatial_sites,"TemporalCorrelatorTminVary");
       double m_momsq_quantum=6.2831853071795864770/double(num_spatial_sites);
       m_momsq_quantum*=m_momsq_quantum;
       if (xmled.count_to_among_children("Anisotropy")==1){
         XMLHandler xmlani(xmled,"Anisotropy");
         string name;
         xmlread(xmlani,"Name",name,"Anisotropy");
         int index=0;
         xmlreadifchild(xmlani,"IDIndex",index);
         aniso_obsinfo = MCObsInfo(name,index);}
       list<XMLHandler> scattering_xml=xmled.find_among_children("ScatteringParticleEnergyFit");
       for (list<XMLHandler>::iterator st=scattering_xml.begin();st!=scattering_xml.end();++st){
          uint psq;
          xmlread(*st,"IntMomSquared",psq,"ScatteringParticleEnergyFit");
          double psqfactor=psq*m_momsq_quantum;
          string name;
          xmlread(*st,"Name",name,"ScatteringParticleEnergyFit");
          int index=0;
          xmlreadifchild(*st,"IDIndex",index);
          scattering_particles.push_back(make_pair(MCObsInfo(name,index),psqfactor));}}

    MCObsInfo chosen_fit_info;
    if (xmltask.count_among_children("ChosenFitInfo")==1){
       XMLHandler xmlchosen(xmltask,"ChosenFitInfo");
       string name;
       xmlread(xmlchosen,"Name",name,"ChosenFitInfo");
       int index=0;
       xmlreadifchild(xmlchosen,"IDIndex",index);
       chosen_fit_info = MCObsInfo(name,index);}

    XMLHandler xmlp(xmltask,"PlotInfo");
    string plotfile;
    xmlread(xmlp,"PlotFile",plotfile,"LogTemporalCorrelatorTminVary");
    if (plotfile.empty()) throw(std::invalid_argument("Must have plot file name"));
    string corrname("standard");
    xmlreadif(xmlp,"CorrName",corrname,"LogTemporalCorrelatorTminVary");
    string symbol("circle");
    xmlreadif(xmlp,"SymbolType",symbol,"LogTemporalCorrelatorTminVary");
    double qualthreshold=0.1;
    xmlreadif(xmlp,"QualityThreshold",qualthreshold,"LogTemporalCorrelatorTminVary");
    string goodfitcolor("blue");
    xmlreadif(xmlp,"GoodFitSymbolColor",goodfitcolor,"LogTemporalCorrelatorTminVary");
    string badfitcolor("red");
    xmlreadif(xmlp,"BadFitSymbolColor",badfitcolor,"LogTemporalCorrelatorTminVary");
    double correlatedthreshold=1.0;
    xmlreadif(xmlp,"CorrelatedThreshold",correlatedthreshold,"TemporalCorrelatorTminVary");
    bool correlatedfit_hollow=false;
    if (xml_child_tag_count(xmlp,"CorrelatedFitSymbolHollow")>0) correlatedfit_hollow=true;
    bool uncorrelatedfit_hollow=false;
    if (xml_child_tag_count(xmlp,"UncorrelatedFitSymbolHollow")>0) uncorrelatedfit_hollow=true;
    vector<XYDYDYPoint> goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits;
    for (uint tmin=tminfirst;tmin<=tminlast;++tmin){
       xmltf.seek_unique("MinimumTimeSeparation");
       xmltf.seek_next_node();       
       xmltf.set_text_content(make_string(tmin)); 
       try{
       LogRealTemporalCorrelatorFit RTC(xmltf,*m_obs,taskcount);
       CorrelatorInfo corr(RTC.m_op,RTC.m_op);
       if (corrname=="standard") corrname=getCorrelatorStandardName(corr);
       const vector<uint>& tvalues=RTC.getTvalues();
       if (find(tvalues.begin(),tvalues.end(),tmin)==tvalues.end()) continue;
       int dof = tvalues.size() - RTC.m_model_ptr->getNumberOfParams();
       if (dof < 1) continue;
       const vector<MCObsInfo>& fitparam_infos=RTC.getFitParamInfos();
       for (uint k=0;k<fitparam_infos.size();++k)
          m_obs->eraseSamplings(fitparam_infos[k]);
       XMLHandler xmlof; RTC.output(xmlof);
       xmlof.rename_tag("LogTemporalCorrelatorTminVaryFit");
       xmlout.put_child(xmlof);
       double chisq_dof,qual;
       doChiSquareFitting(RTC,mz_info,chisq_dof,qual,bestfit_params,xmlout);
       TCorrFitInfo fitinfo;
       uint meff_tstep=1; bool showapproach=false;
       RTC.m_model_ptr->setFitInfo(RTC.m_fitparam_info,bestfit_params,tmin,tmax,
                                   showapproach,meff_tstep,chisq_dof,qual,fitinfo);
       MCObsInfo energy_key=fitinfo.energy_key;
       if (scattering_particles.size()>0){
          if (aniso_obsinfo.isVacuum()) // no anisotropy
            doEnergyDifferenceBySamplings(*m_obs,energy_key,scattering_particles,energy_key);
          else
            doEnergyDifferenceBySamplings(*m_obs,energy_key,aniso_obsinfo,scattering_particles,energy_key);}

       bool correlated=false;
       if (!chosen_fit_info.isVacuum()){
          MCObsInfo diff_obs;
          doCorrelatedDifferenceBySamplings(*m_obs,chosen_fit_info,energy_key,diff_obs);
          MCEstimate diff_est=m_obs->getEstimate(diff_obs);
          double diff_val=diff_est.getFullEstimate();
          double diff_up,diff_down;
          if (diff_est.isJackknifeMode())
             diff_up=diff_down=correlatedthreshold*diff_est.getSymmetricError();
          else{
             diff_up=correlatedthreshold*(diff_est.getUpperConfLimit()-diff_val);
             diff_down=correlatedthreshold*(diff_val-diff_est.getLowerConfLimit());}

          double upper_limit=diff_up+diff_val;
          double lower_limit=diff_val-diff_down;
          correlated = upper_limit >= 0. && lower_limit <= 0.;}

       MCEstimate energy=m_obs->getEstimate(energy_key);
       double y=energy.getFullEstimate();
       double dyup,dydn;
       if (energy.isJackknifeMode()) 
          dyup=dydn=energy.getSymmetricError();
       else{
          dyup=energy.getUpperConfLimit()-y;
          dydn=y-energy.getLowerConfLimit();}
       if (qual>=0.1 && correlated) goodcorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else if (qual>=0.1 && !correlated) gooduncorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else if (qual<0.1 && correlated) badcorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else baduncorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));}
       catch(const std::exception& errmsg){
          xmlout.put_child("Error",string("DoFit within type LogTemporalCorrelatorTminVary encountered an error: ")
                 +string(errmsg.what()));
       }}
    XMLHandler xmlplog("TminPlot");
    xmlplog.put_child("PlotFile",plotfile);
    xmlplog.put_child("QualityThreshold",make_string(qualthreshold));
    xmlplog.put_child("CorrelatedThreshold",make_string(correlatedthreshold));
    xmlplog.put_child("NumberOfGoodCorrelatedFitPoints",make_string(goodcorrelatedfits.size()));
    xmlplog.put_child("NumberOfGoodUncorrelatedFitPoints",make_string(gooduncorrelatedfits.size()));
    xmlplog.put_child("NumberOfBadCorrelatedFitPoints",make_string(badcorrelatedfits.size()));
    xmlplog.put_child("NumberOfBadUncorrelatedFitPoints",make_string(baduncorrelatedfits.size()));
    xmlout.put_child(xmlplog);

    XYDYDYPoint chosen_fit(0,0,0,0);
    if (!chosen_fit_info.isVacuum()){
       MCEstimate chosen_fit_estimate=m_obs->getEstimate(chosen_fit_info);
       double y=chosen_fit_estimate.getFullEstimate();
       double dyup,dydn;
       if (chosen_fit_estimate.isJackknifeMode()) 
          dyup=dydn=chosen_fit_estimate.getSymmetricError();
       else{
          dyup=chosen_fit_estimate.getUpperConfLimit()-y;
          dydn=y-chosen_fit_estimate.getLowerConfLimit();}
       chosen_fit = XYDYDYPoint(1,y,dyup,dydn);}
    createTMinPlot(goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits,
                   corrname,plotfile,symbol,goodfitcolor,badfitcolor,correlatedfit_hollow,
                   uncorrelatedfit_hollow,chosen_fit);}
    catch(const std::exception& errmsg){
       xmlout.put_child("Error",string("DoFit with type LogTemporalCorrelatorTminVary encountered an error: ")
               +string(errmsg.what()));
    }}



 else if (fittype=="TemporalCorrelatorInteractionRatioTminVary"){
    try{
    XMLHandler xmlf(xmltask,"TemporalCorrelatorInteractionRatioTminVaryFit");
    XMLHandler xmlres(xmlf,"Ratio");
    OperatorInfo ratio_op(xmlres);
    XMLHandler xmlint(xmlf,"InteractingOperator");
    bool numvev=(xmlint.count("SubtractVEV")>0) ? true: false;
    pair<OperatorInfo,bool> numerator=make_pair(OperatorInfo(xmlint),numvev);
    vector<pair<OperatorInfo,bool> > denominator;
    list<XMLHandler> denomxml=xmlf.find_among_children("NonInteractingOperator");
    for (list<XMLHandler>::iterator it=denomxml.begin();it!=denomxml.end();++it){
      OperatorInfo opinfo(*it);
      bool subvev=(it->count("SubtractVEV")>0) ? true: false;
      denominator.push_back(make_pair(opinfo,subvev));}
    uint nterms=denominator.size();
    if (nterms<2) throw(std::invalid_argument("Two or more NonInteractingOperators required"));
    XMLHandler xmlo, xmldp;
    xmlo.set_root("Ratio");
    ratio_op.output(xmldp);
    xmlo.put_child(xmldp);
    xmlout.put_child(xmlo);
    xmlo.set_root("InteractingOperator");
    numerator.first.output(xmldp);
    xmlo.put_child(xmldp);
    xmlout.put_child(xmlo);
    xmlo.set_root("NonInteractingOperators");
    for (vector<pair<OperatorInfo,bool> >::const_iterator
           it=denominator.begin();it!=denominator.end();it++){
       XMLHandler xmloo; it->first.output(xmloo); xmlo.put_child(xmloo);}
    xmlout.put_child(xmlo);
    set<MCObsInfo> obskeys;
    bool erase_orig=true;
    uint tminfirst,tminlast,tmax;
    xmlread(xmlf,"TminFirst",tminfirst,"TemporalCorrelatorInteractionRatioTminVary");
    xmlread(xmlf,"TminLast",tminlast,"TemporalCorrelatorInteractionRatioTminVary");
    xmlread(xmlf,"Tmax",tmax,"TemporalCorrelatorInteractionRatioTminVary");

    doCorrelatorInteractionRatioBySamplings(*m_obs,numerator,denominator,
                                            0,(tmax<64)?64:tmax,ratio_op,obskeys,erase_orig);


    XMLHandler xmltf(xmlf,XMLHandler::subtree_copy);
    xmltf.rename_tag("TemporalCorrelatorFit");
    xmltf.put_child("MaximumTimeSeparation",make_string(tmax));
    xmltf.put_child("MinimumTimeSeparation",make_string(tminfirst));
    XMLHandler xmlro; ratio_op.output(xmlro);
    xmltf.put_child(xmlro); 

    list<pair<MCObsInfo,double> > scattering_particles;
    MCObsInfo aniso_obsinfo;
    if (xmltask.count_among_children("DoReconstructEnergy")==1){
       XMLHandler xmled(xmltask,"DoReconstructEnergy");
       uint num_spatial_sites=0;
       xmlread(xmled,"SpatialExtentNumSites",num_spatial_sites,"TemporalCorrelatorTminVary");
       double m_momsq_quantum=6.2831853071795864770/double(num_spatial_sites);
       m_momsq_quantum*=m_momsq_quantum;
       if (xmled.count_to_among_children("Anisotropy")==1){
         XMLHandler xmlani(xmled,"Anisotropy");
         string name;
         xmlread(xmlani,"Name",name,"Anisotropy");
         int index=0;
         xmlreadifchild(xmlani,"IDIndex",index);
         aniso_obsinfo = MCObsInfo(name,index);}
       list<XMLHandler> scattering_xml=xmled.find_among_children("ScatteringParticleEnergyFit");
       for (list<XMLHandler>::iterator st=scattering_xml.begin();st!=scattering_xml.end();++st){
          uint psq;
          xmlread(*st,"IntMomSquared",psq,"ScatteringParticleEnergyFit");
          double psqfactor=psq*m_momsq_quantum;
          string name;
          xmlread(*st,"Name",name,"ScatteringParticleEnergyFit");
          int index=0;
          xmlreadifchild(*st,"IDIndex",index);
          scattering_particles.push_back(make_pair(MCObsInfo(name,index),psqfactor));}}

    MCObsInfo chosen_fit_info;
    if (xmltask.count_among_children("ChosenFitInfo")==1){
       XMLHandler xmlchosen(xmltask,"ChosenFitInfo");
       string name;
       xmlread(xmlchosen,"Name",name,"ChosenFitInfo");
       int index=0;
       xmlreadifchild(xmlchosen,"IDIndex",index);
       chosen_fit_info = MCObsInfo(name,index);}

    XMLHandler xmlp(xmltask,"PlotInfo");
    string plotfile;
    xmlread(xmlp,"PlotFile",plotfile,"TemporalCorrelatorInteractionRatioTminVary");
    if (plotfile.empty()) throw(std::invalid_argument("Must have plot file name"));
    string corrname("standard");
    xmlreadif(xmlp,"CorrName",corrname,"TemporalCorrelatorInteractionRatioTminVary");
    string symbol("circle");
    xmlreadif(xmlp,"SymbolType",symbol,"TemporalCorrelatorInteractionRatioTminVary");
    double qualthreshold=0.1;
    xmlreadif(xmlp,"QualityThreshold",qualthreshold,"TemporalCorrelatorInteractionRatioTminVary");
    string goodfitcolor("blue");
    xmlreadif(xmlp,"GoodFitSymbolColor",goodfitcolor,"TemporalCorrelatorInteractionRatioTminVary");
    string badfitcolor("red");
    xmlreadif(xmlp,"BadFitSymbolColor",badfitcolor,"TemporalCorrelatorInteractionRatioTminVary");
    double correlatedthreshold=1.0;
    xmlreadif(xmlp,"CorrelatedThreshold",correlatedthreshold,"TemporalCorrelatorTminVary");
    bool correlatedfit_hollow=false;
    if (xml_child_tag_count(xmlp,"CorrelatedFitSymbolHollow")>0) correlatedfit_hollow=true;
    bool uncorrelatedfit_hollow=false;
    if (xml_child_tag_count(xmlp,"UncorrelatedFitSymbolHollow")>0) uncorrelatedfit_hollow=true;
    vector<XYDYDYPoint> goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits;
    for (uint tmin=tminfirst;tmin<=tminlast;++tmin){
       xmltf.seek_unique("MinimumTimeSeparation");
       xmltf.seek_next_node();       
       xmltf.set_text_content(make_string(tmin)); 
       try{
       RealTemporalCorrelatorFit RTC(xmltf,*m_obs,taskcount);
       CorrelatorInfo corr(RTC.m_op,RTC.m_op);
       if (corrname=="standard") corrname=getCorrelatorStandardName(corr);
       const vector<uint>& tvalues=RTC.getTvalues();
       if (find(tvalues.begin(),tvalues.end(),tmin)==tvalues.end()) continue;
       int dof = tvalues.size() - RTC.m_model_ptr->getNumberOfParams();
       if (dof < 1) continue;
       const vector<MCObsInfo>& fitparam_infos=RTC.getFitParamInfos();
       for (uint k=0;k<fitparam_infos.size();++k)
          m_obs->eraseSamplings(fitparam_infos[k]);
       XMLHandler xmlof; RTC.output(xmlof);
       xmlof.rename_tag("TemporalCorrelatorInteractionRatioTminVaryFit");
       xmlout.put_child(xmlof);
       double chisq_dof,qual;
       doChiSquareFitting(RTC,mz_info,chisq_dof,qual,bestfit_params,xmlout);
       TCorrFitInfo fitinfo;
       uint meff_tstep=1; bool showapproach=false;
       RTC.m_model_ptr->setFitInfo(RTC.m_fitparam_info,bestfit_params,tmin,tmax,
                                   showapproach,meff_tstep,chisq_dof,qual,fitinfo);
       MCObsInfo energy_key=fitinfo.energy_key;
       if (scattering_particles.size()>0){
          if (aniso_obsinfo.isVacuum()) // no anisotropy
            doReconstructEnergyBySamplings(*m_obs,energy_key,scattering_particles,energy_key);
          else
            doReconstructEnergyBySamplings(*m_obs,energy_key,aniso_obsinfo,scattering_particles,energy_key);}

       bool correlated=false;
       if (!chosen_fit_info.isVacuum()){
          MCObsInfo diff_obs;
          doCorrelatedDifferenceBySamplings(*m_obs,chosen_fit_info,energy_key,diff_obs);
          MCEstimate diff_est=m_obs->getEstimate(diff_obs);
          double diff_val=diff_est.getFullEstimate();
          double diff_up,diff_down;
          if (diff_est.isJackknifeMode())
             diff_up=diff_down=correlatedthreshold*diff_est.getSymmetricError();
          else{
             diff_up=correlatedthreshold*(diff_est.getUpperConfLimit()-diff_val);
             diff_down=correlatedthreshold*(diff_val-diff_est.getLowerConfLimit());}

          double upper_limit=diff_up+diff_val;
          double lower_limit=diff_val-diff_down;
          correlated = upper_limit >= 0. && lower_limit <= 0.;}

       MCEstimate energy=m_obs->getEstimate(energy_key);
       double y=energy.getFullEstimate();
       double dyup,dydn;
       if (energy.isJackknifeMode()) 
          dyup=dydn=energy.getSymmetricError();
       else{
          dyup=energy.getUpperConfLimit()-y;
          dydn=y-energy.getLowerConfLimit();}
       if (qual>=0.1 && correlated) goodcorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else if (qual>=0.1 && !correlated) gooduncorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else if (qual<0.1 && correlated) badcorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else baduncorrelatedfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));}
       catch(const std::exception& errmsg){
          xmlout.put_child("Error",string("DoFit within type TemporalCorrelatorInteractionRatioTminVary encountered an error: ")
                 +string(errmsg.what()));
       }}
    XMLHandler xmlplog("TminPlot");
    xmlplog.put_child("PlotFile",plotfile);
    xmlplog.put_child("QualityThreshold",make_string(qualthreshold));
    xmlplog.put_child("CorrelatedThreshold",make_string(correlatedthreshold));
    xmlplog.put_child("NumberOfGoodCorrelatedFitPoints",make_string(goodcorrelatedfits.size()));
    xmlplog.put_child("NumberOfGoodUncorrelatedFitPoints",make_string(gooduncorrelatedfits.size()));
    xmlplog.put_child("NumberOfBadCorrelatedFitPoints",make_string(badcorrelatedfits.size()));
    xmlplog.put_child("NumberOfBadUncorrelatedFitPoints",make_string(baduncorrelatedfits.size()));
    xmlout.put_child(xmlplog);

    XYDYDYPoint chosen_fit(0,0,0,0);
    if (!chosen_fit_info.isVacuum()){
       MCEstimate chosen_fit_estimate=m_obs->getEstimate(chosen_fit_info);
       double y=chosen_fit_estimate.getFullEstimate();
       double dyup,dydn;
       if (chosen_fit_estimate.isJackknifeMode()) 
          dyup=dydn=chosen_fit_estimate.getSymmetricError();
       else{
          dyup=chosen_fit_estimate.getUpperConfLimit()-y;
          dydn=y-chosen_fit_estimate.getLowerConfLimit();}
       chosen_fit = XYDYDYPoint(1,y,dyup,dydn);}
    createTMinPlot(goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits,
                   corrname,plotfile,symbol,goodfitcolor,badfitcolor,correlatedfit_hollow,
                   uncorrelatedfit_hollow,chosen_fit);}
    catch(const std::exception& errmsg){
       xmlout.put_child("Error",string("DoFit with type TemporalCorrelatorInteractionRatioTminVary encountered an error: ")
               +string(errmsg.what()));
    }}

 else if (fittype=="NSimTemporalCorrelatorTminVary"){
    try{
    XMLHandler xmlf(xmltask,"NSimTemporalCorrelatorTminVaryFit");
    uint tminfirst,tminlast,tmax;
    xmlread(xmlf,"TminFirst",tminfirst,"NSimTemporalCorrelatorTminVary");
    xmlread(xmlf,"TminLast",tminlast,"NSimTemporalCorrelatorTminVary");
    xmlread(xmlf,"Tmax",tmax,"NSimTemporalCorrelatorTminVary");
        
    XMLHandler xmlfit(xmlf,"NSimTemporalCorrelatorFit");
    XMLHandler xmltf(xmlfit,XMLHandler::subtree_copy);
        
    XMLHandler xmlc(xmltf,"Fits");
    list<XMLHandler> xmlccs = xmlc.find("TemporalCorrelatorFit");
    list<XMLHandler> xmlccs2 = xmlc.find("TemporalCorrelatorInteractionRatioFit");
    uint nfits = xmlccs.size()+xmlccs2.size();
           
    vector<MCObsInfo> chosen_fit_info;
    chosen_fit_info.resize(nfits);
        
    uint ii = 0;
    for(list<XMLHandler>::iterator it = xmlccs.begin(); it != xmlccs.end(); ++it){
        if(!it->count_among_children("Fixed")){
            it->put_child("MaximumTimeSeparation",make_string(tmax));
            it->put_child("MinimumTimeSeparation",make_string(tminfirst));
        }
        if (it->count_among_children("ChosenFitInfo")==1){
           XMLHandler xmlchosen(*it,"ChosenFitInfo");
           string name;
           xmlread(xmlchosen,"Name",name,"ChosenFitInfo");
           int index=0;
           xmlreadifchild(xmlchosen,"IDIndex",index);
           chosen_fit_info[ii] = MCObsInfo(name,index);}
        ii++;
    }
    for(list<XMLHandler>::iterator it = xmlccs2.begin(); it != xmlccs2.end(); ++it){
        if(!it->count_among_children("Fixed")){
            it->put_child("MaximumTimeSeparation",make_string(tmax));
            it->put_child("MinimumTimeSeparation",make_string(tminfirst));
        }
        if (it->count_among_children("ChosenFitInfo")==1){
           XMLHandler xmlchosen(*it,"ChosenFitInfo");
           string name;
           xmlread(xmlchosen,"Name",name,"ChosenFitInfo");
           int index=0;
           xmlreadifchild(xmlchosen,"IDIndex",index);
           chosen_fit_info[ii] = MCObsInfo(name,index);}
        ii++;
    }
        
//     XMLHandler xmltf(xmlf,XMLHandler::subtree_copy);
//     xmltf.rename_tag("TemporalCorrelatorFit");
//     list<pair<MCObsInfo,double> > scattering_particles;
//     MCObsInfo aniso_obsinfo;
//     if (xmltask.count_among_children("DoEnergyDifference")==1){
//        XMLHandler xmled(xmltask,"DoEnergyDifference");
//        uint num_spatial_sites=0;
//        xmlread(xmled,"SpatialExtentNumSites",num_spatial_sites,"TemporalCorrelatorTminVary");
//        double m_momsq_quantum=6.2831853071795864770/double(num_spatial_sites);
//        m_momsq_quantum*=m_momsq_quantum;
//        if (xmled.count_to_among_children("Anisotropy")==1){
//          XMLHandler xmlani(xmled,"Anisotropy");
//          string name;
//          xmlread(xmlani,"Name",name,"Anisotropy");
//          int index=0;
//          xmlreadifchild(xmlani,"IDIndex",index);
//          aniso_obsinfo = MCObsInfo(name,index);}
//        list<XMLHandler> scattering_xml=xmled.find_among_children("ScatteringParticleEnergyFit");
//        for (list<XMLHandler>::iterator st=scattering_xml.begin();st!=scattering_xml.end();++st){
//           uint psq;
//           xmlread(*st,"IntMomSquared",psq,"ScatteringParticleEnergyFit");
//           double psqfactor=psq*m_momsq_quantum;
//           string name;
//           xmlread(*st,"Name",name,"ScatteringParticleEnergyFit");
//           int index=0;
//           xmlreadifchild(*st,"IDIndex",index);
//           scattering_particles.push_back(make_pair(MCObsInfo(name,index),psqfactor));}}

//     if (xmltask.count_among_children("ChosenFitInfo")==1){
//        XMLHandler xmlchosen(xmltask,"ChosenFitInfo");
//        string name;
//        xmlread(xmlchosen,"Name",name,"ChosenFitInfo");
//        int index=0;
//        xmlreadifchild(xmlchosen,"IDIndex",index);
//        chosen_fit_info = MCObsInfo(name,index);}

//     string plotfile;
//     xmlread(xmlp,"PlotFile",plotfile,"NSimTemporalCorrelatorTminVary");
//     if (plotfile.empty()) throw(std::invalid_argument("Must have plot file name"));
//     string corrname("standard");
//     xmlreadif(xmlp,"CorrName",corrname,"NSimTemporalCorrelatorTminVary");
//     string symbol("circle");
//     xmlreadif(xmlp,"SymbolType",symbol,"TemporalCorrelatorTminVary");
//     double qualthreshold=0.1;
//     xmlreadif(xmlp,"QualityThreshold",qualthreshold,"TemporalCorrelatorTminVary");
//     string goodfitcolor("blue");
//     xmlreadif(xmlp,"GoodFitSymbolColor",goodfitcolor,"TemporalCorrelatorTminVary");
//     string badfitcolor("red");
//     xmlreadif(xmlp,"BadFitSymbolColor",badfitcolor,"TemporalCorrelatorTminVary");
    

//     bool correlatedfit_hollow=false;
//     if (xml_child_tag_count(xmlp,"CorrelatedFitSymbolHollow")>0) correlatedfit_hollow=true;
//     bool uncorrelatedfit_hollow=false;
//     if (xml_child_tag_count(xmlp,"UncorrelatedFitSymbolHollow")>0) uncorrelatedfit_hollow=true;
    vector<vector<XYDYDYPoint>> goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits;
    goodcorrelatedfits.resize(nfits);
    gooduncorrelatedfits.resize(nfits);
    badcorrelatedfits.resize(nfits);
    baduncorrelatedfits.resize(nfits);
    for (uint tmin=tminfirst;tmin<=tminlast;++tmin){
        for(list<XMLHandler>::iterator it = xmlccs.begin(); it != xmlccs.end(); ++it){
            if(!it->count_among_children("Fixed")){
                it->seek_unique("MinimumTimeSeparation");
                it->seek_next_node();
                it->set_text_content(make_string(tmin));
            }
        }
        for(list<XMLHandler>::iterator it = xmlccs2.begin(); it != xmlccs2.end(); ++it){
            if(!it->count_among_children("Fixed")){
                it->seek_unique("MinimumTimeSeparation");
                it->seek_next_node();
                it->set_text_content(make_string(tmin));
            }
        }
       try{
        NSimRealTemporalCorrelatorFit NSimRTC(xmltf,*m_obs,taskcount);
//        CorrelatorInfo corr(RTC.m_op,RTC.m_op);
//        if (corrname=="standard") corrname=getCorrelatorStandardName(corr);
//        const vector<uint>& tvalues=RTC.getTvalues();
//        if (find(tvalues.begin(),tvalues.end(),tmin)==tvalues.end()) continue;
//        int dof = tvalues.size() - RTC.m_model_ptr->getNumberOfParams();
//        if (dof < 1) continue;
       const vector<MCObsInfo>& fitparam_infos=NSimRTC.getFitParamInfos();
       for (uint k=0;k<fitparam_infos.size();++k)
          m_obs->eraseSamplings(fitparam_infos[k]);
       XMLHandler xmlof; NSimRTC.output(xmlof);
//        xmlof.rename_tag("TemporalCorrelatorTminVaryFit");
       xmlout.put_child(xmlof);
       double chisq_dof,qual;
       doChiSquareFitting(NSimRTC,mz_info,chisq_dof,qual,bestfit_params,xmlout);
       for( uint i=0;i<nfits;i++){
           TCorrFitInfo fitinfo;
           double correlatedthreshold=1.0;
           uint meff_tstep=1; bool showapproach=false;
           //this is not the correct way to set fitinfo.energy_mean, fitinfo.energy_err, but we don't use those here
           NSimRTC.m_fits[i]->m_model_ptr->setFitInfo(NSimRTC.m_fits[i]->m_fitparam_info,bestfit_params,tmin,tmax,
                                       showapproach,meff_tstep,chisq_dof,qual,fitinfo); 
           MCObsInfo energy_key=fitinfo.energy_key;

           bool correlated=false;
           
           if (!chosen_fit_info[i].isVacuum()){
              MCObsInfo diff_obs;
              doCorrelatedDifferenceBySamplings(*m_obs,chosen_fit_info[i],energy_key,diff_obs);
              MCEstimate diff_est=m_obs->getEstimate(diff_obs);
              double diff_val=diff_est.getFullEstimate();
              double diff_up,diff_down;
              if (diff_est.isJackknifeMode())
                 diff_up=diff_down=correlatedthreshold*diff_est.getSymmetricError();
              else{
                 diff_up=correlatedthreshold*(diff_est.getUpperConfLimit()-diff_val);
                 diff_down=correlatedthreshold*(diff_val-diff_est.getLowerConfLimit());}

              double upper_limit=diff_up+diff_val;
              double lower_limit=diff_val-diff_down;
              correlated = upper_limit >= 0. && lower_limit <= 0.;}

           MCEstimate energy=m_obs->getEstimate(energy_key);
           double y=energy.getFullEstimate();
           double dyup,dydn;
           if (energy.isJackknifeMode()) 
              dyup=dydn=energy.getSymmetricError();
           else{
              dyup=energy.getUpperConfLimit()-y;
              dydn=y-energy.getLowerConfLimit();}
           if (qual>=0.1 && correlated) goodcorrelatedfits[i].push_back(XYDYDYPoint(tmin,y,dyup,dydn));
           else if (qual>=0.1 && !correlated) gooduncorrelatedfits[i].push_back(XYDYDYPoint(tmin,y,dyup,dydn));
           else if (qual<0.1 && correlated) badcorrelatedfits[i].push_back(XYDYDYPoint(tmin,y,dyup,dydn));
           else baduncorrelatedfits[i].push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       }
       }catch(const std::exception& errmsg){
          xmlout.put_child("Error",string("DoFit within type NSimTemporalCorrelatorTminVary encountered an error: ")
                 +string(errmsg.what()));
       }
    }
    XMLHandler xmlplog("TminPlots");
    
     uint i = 0;
    for(list<XMLHandler>::iterator it = xmlccs.begin(); it != xmlccs.end(); ++it){
//         if(!it->count_among_children("Fixed")){
            string plotfile;
            XMLHandler xmlp(*it,"PlotInfo");
            xmlread(xmlp,"PlotFile",plotfile,"TemporalCorrelatorFit");
            xmlplog.put_child("PlotFile",plotfile);
            if (plotfile.empty()) continue;
            
            string corrname("standard");
            xmlreadif(xmlp,"CorrName",corrname,"TemporalCorrelatorFit");
            string symbol("circle");
            xmlreadif(xmlp,"SymbolType",symbol,"TemporalCorrelatorFit");
            double qualthreshold=0.1;
            xmlreadif(xmlp,"QualityThreshold",qualthreshold,"TemporalCorrelatorFit");
            string goodfitcolor("blue");
            xmlreadif(xmlp,"GoodFitSymbolColor",goodfitcolor,"TemporalCorrelatorFit");
            string badfitcolor("red");
            xmlreadif(xmlp,"BadFitSymbolColor",badfitcolor,"TemporalCorrelatorFit");
            
            bool correlatedfit_hollow=false;
            if (xml_child_tag_count(xmlp,"CorrelatedFitSymbolHollow")>0) correlatedfit_hollow=true;
            bool uncorrelatedfit_hollow=false;
            if (xml_child_tag_count(xmlp,"UncorrelatedFitSymbolHollow")>0) uncorrelatedfit_hollow=true;
            
            XYDYDYPoint chosen_fit(0,0,0,0);
            if (!chosen_fit_info[i].isVacuum()){
               MCEstimate chosen_fit_estimate=m_obs->getEstimate(chosen_fit_info[i]);
               double y=chosen_fit_estimate.getFullEstimate();
               double dyup,dydn;
               if (chosen_fit_estimate.isJackknifeMode()) 
                  dyup=dydn=chosen_fit_estimate.getSymmetricError();
               else{
                  dyup=chosen_fit_estimate.getUpperConfLimit()-y;
                  dydn=y-chosen_fit_estimate.getLowerConfLimit();}
               chosen_fit = XYDYDYPoint(1,y,dyup,dydn);}
            createTMinPlot(goodcorrelatedfits[i],gooduncorrelatedfits[i],
                   badcorrelatedfits[i],baduncorrelatedfits[i],corrname,plotfile,symbol,
                   goodfitcolor,badfitcolor,correlatedfit_hollow,
                   uncorrelatedfit_hollow,chosen_fit);
//         }
        i++;
    }
    for(list<XMLHandler>::iterator it = xmlccs2.begin(); it != xmlccs2.end(); ++it){
//         if(!it->count_among_children("Fixed")){
            string plotfile;
            XMLHandler xmlp(*it,"PlotInfo");
            xmlread(xmlp,"PlotFile",plotfile,"TemporalCorrelatorFit");
//             xmlplog.put_child("PlotFile",plotfile);
            if (plotfile.empty()) continue;
            
            string corrname("standard");
            xmlreadif(xmlp,"CorrName",corrname,"TemporalCorrelatorFit");
            string symbol("circle");
            xmlreadif(xmlp,"SymbolType",symbol,"TemporalCorrelatorFit");
            double qualthreshold=0.1;
            xmlreadif(xmlp,"QualityThreshold",qualthreshold,"TemporalCorrelatorFit");
            string goodfitcolor("blue");
            xmlreadif(xmlp,"GoodFitSymbolColor",goodfitcolor,"TemporalCorrelatorFit");
            string badfitcolor("red");
            xmlreadif(xmlp,"BadFitSymbolColor",badfitcolor,"TemporalCorrelatorFit");
            
            bool correlatedfit_hollow=false;
            if (xml_child_tag_count(xmlp,"CorrelatedFitSymbolHollow")>0) correlatedfit_hollow=true;
            bool uncorrelatedfit_hollow=false;
            if (xml_child_tag_count(xmlp,"UncorrelatedFitSymbolHollow")>0) uncorrelatedfit_hollow=true;
            
        
            XYDYDYPoint chosen_fit(0,0,0,0);
            if (!chosen_fit_info[i].isVacuum()){
               MCEstimate chosen_fit_estimate=m_obs->getEstimate(chosen_fit_info[i]);
               double y=chosen_fit_estimate.getFullEstimate();
               double dyup,dydn;
               if (chosen_fit_estimate.isJackknifeMode()) 
                  dyup=dydn=chosen_fit_estimate.getSymmetricError();
               else{
                  dyup=chosen_fit_estimate.getUpperConfLimit()-y;
                  dydn=y-chosen_fit_estimate.getLowerConfLimit();}
               chosen_fit = XYDYDYPoint(1,y,dyup,dydn);}
        
            createTMinPlot(goodcorrelatedfits[i],gooduncorrelatedfits[i],
                   badcorrelatedfits[i],baduncorrelatedfits[i],corrname,plotfile,symbol,
                   goodfitcolor,badfitcolor,correlatedfit_hollow,
                   uncorrelatedfit_hollow,chosen_fit);
//         }
        i++;
    }
        
//     xmlplog.put_child("QualityThreshold",make_string(qualthreshold));
//     xmlplog.put_child("CorrelatedThreshold",make_string(correlatedthreshold));
//     xmlplog.put_child("NumberOfGoodCorrelatedFitPoints",make_string(goodcorrelatedfits.size()));
//     xmlplog.put_child("NumberOfGoodUncorrelatedFitPoints",make_string(gooduncorrelatedfits.size()));
//     xmlplog.put_child("NumberOfBadCorrelatedFitPoints",make_string(badcorrelatedfits.size()));
//     xmlplog.put_child("NumberOfBadUncorrelatedFitPoints",make_string(baduncorrelatedfits.size()));
    xmlout.put_child(xmlplog);

//     XYDYDYPoint chosen_fit(0,0,0,0);
//     if (!chosen_fit_info.isVacuum()){
//        MCEstimate chosen_fit_estimate=m_obs->getEstimate(chosen_fit_info);
//        double y=chosen_fit_estimate.getFullEstimate();
//        double dyup,dydn;
//        if (chosen_fit_estimate.isJackknifeMode()) 
//           dyup=dydn=chosen_fit_estimate.getSymmetricError();
//        else{
//           dyup=chosen_fit_estimate.getUpperConfLimit()-y;
//           dydn=y-chosen_fit_estimate.getLowerConfLimit();}
//        chosen_fit = XYDYDYPoint(1,y,dyup,dydn);}
//     createTMinPlot(goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits,
//                    corrname,plotfile,symbol,goodfitcolor,badfitcolor,correlatedfit_hollow,
//                    uncorrelatedfit_hollow,chosen_fit);}
    }catch(const std::exception& errmsg){
       xmlout.put_child("Error",string("DoFit with type NSimTemporalCorrelatorTminVary encountered an error: ")
               +string(errmsg.what()));
    }
 }

 else if (fittype=="TemporalCorrelatorTmaxVary"){
    try{
    XMLHandler xmlf(xmltask,"TemporalCorrelatorTmaxVaryFit");
    uint tmaxfirst,tmaxlast,tmin;
    xmlread(xmlf,"TmaxFirst",tmaxfirst,"TemporalCorrelatorTmaxVary");
    xmlread(xmlf,"TmaxLast",tmaxlast,"TemporalCorrelatorTmaxVary");
    xmlread(xmlf,"Tmin",tmin,"TemporalCorrelatorTmaxVary");
    XMLHandler xmltf(xmlf,XMLHandler::subtree_copy);
    xmltf.rename_tag("TemporalCorrelatorFit");
    xmltf.put_child("MaximumTimeSeparation",make_string(tmaxfirst));
    xmltf.put_child("MinimumTimeSeparation",make_string(tmin));

    list<pair<MCObsInfo,double> > scattering_particles;
    MCObsInfo aniso_obsinfo;
    if (xmltask.count_among_children("DoEnergyDifference")==1){
       XMLHandler xmled(xmltask,"DoEnergyDifference");
       uint num_spatial_sites=0;
       xmlread(xmled,"SpatialExtentNumSites",num_spatial_sites,"TemporalCorrelatorTmaxVary");
       double m_momsq_quantum=6.2831853071795864770/double(num_spatial_sites);
       m_momsq_quantum*=m_momsq_quantum;
       if (xmled.count_to_among_children("Anisotropy")==1){
         XMLHandler xmlani(xmled,"Anisotropy");
         string name;
         xmlread(xmlani,"Name",name,"Anisotropy");
         int index=0;
         xmlreadifchild(xmlani,"IDIndex",index);
         aniso_obsinfo = MCObsInfo(name,index);}
       list<XMLHandler> scattering_xml=xmled.find_among_children("ScatteringParticleEnergyFit");
       for (list<XMLHandler>::iterator st=scattering_xml.begin();st!=scattering_xml.end();++st){
          uint psq;
          xmlread(*st,"IntMomSquared",psq,"ScatteringParticleEnergyFit");
          double psqfactor=psq*m_momsq_quantum;
          string name;
          xmlread(*st,"Name",name,"ScatteringParticleEnergyFit");
          int index=0;
          xmlreadifchild(*st,"IDIndex",index);
          scattering_particles.push_back(make_pair(MCObsInfo(name,index),psqfactor));}}

    MCObsInfo chosen_fit_info;
    if (xmltask.count_among_children("ChosenFitInfo")==1){
       XMLHandler xmlchosen(xmltask,"ChosenFitInfo");
       string name;
       xmlread(xmlchosen,"Name",name,"ChosenFitInfo");
       int index=0;
       xmlreadifchild(xmlchosen,"IDIndex",index);
       chosen_fit_info = MCObsInfo(name,index);}

    XMLHandler xmlp(xmltask,"PlotInfo");
    string plotfile;
    xmlread(xmlp,"PlotFile",plotfile,"TemporalCorrelatorTmaxVary");
    if (plotfile.empty()) throw(std::invalid_argument("Must have plot file name"));
    string corrname("standard");
    xmlreadif(xmlp,"CorrName",corrname,"TemporalCorrelatorTmaxVary");
    string symbol("circle");
    xmlreadif(xmlp,"SymbolType",symbol,"TemporalCorrelatorTmaxVary");
    double qualthreshold=0.1;
    xmlreadif(xmlp,"QualityThreshold",qualthreshold,"TemporalCorrelatorTmaxVary");
    string goodfitcolor("blue");
    xmlreadif(xmlp,"GoodFitSymbolColor",goodfitcolor,"TemporalCorrelatorTmaxVary");
    string badfitcolor("red");
    xmlreadif(xmlp,"BadFitSymbolColor",badfitcolor,"TemporalCorrelatorTmaxVary");
    double correlatedthreshold=1.0;
    xmlreadif(xmlp,"CorrelatedThreshold",correlatedthreshold,"TemporalCorrelatorTmaxVary");
    bool correlatedfit_hollow=false;
    if (xml_child_tag_count(xmlp,"CorrelatedFitSymbolHollow")>0) correlatedfit_hollow=true;
    bool uncorrelatedfit_hollow=false;
    if (xml_child_tag_count(xmlp,"UncorrelatedFitSymbolHollow")>0) uncorrelatedfit_hollow=true;
    vector<XYDYDYPoint> goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits;
    for (uint tmax=tmaxfirst;tmax<=tmaxlast;++tmax){
       xmltf.seek_unique("MaximumTimeSeparation");
       xmltf.seek_next_node();       
       xmltf.set_text_content(make_string(tmax)); 
       try{
       RealTemporalCorrelatorFit RTC(xmltf,*m_obs,taskcount);
       CorrelatorInfo corr(RTC.m_op,RTC.m_op);
       if (corrname=="standard") corrname=getCorrelatorStandardName(corr);
       const vector<uint>& tvalues=RTC.getTvalues();
       if (find(tvalues.begin(),tvalues.end(),tmin)==tvalues.end()) continue;
       int dof = tvalues.size() - RTC.m_model_ptr->getNumberOfParams();
       if (dof < 1) continue;
       const vector<MCObsInfo>& fitparam_infos=RTC.getFitParamInfos();
       for (uint k=0;k<fitparam_infos.size();++k)
          m_obs->eraseSamplings(fitparam_infos[k]);
       XMLHandler xmlof; RTC.output(xmlof);
       xmlof.rename_tag("TemporalCorrelatorTmaxVaryFit");
       xmlout.put_child(xmlof);
       double chisq_dof,qual;
       doChiSquareFitting(RTC,mz_info,chisq_dof,qual,bestfit_params,xmlout);
       TCorrFitInfo fitinfo;
       uint meff_tstep=1; bool showapproach=false;
       RTC.m_model_ptr->setFitInfo(RTC.m_fitparam_info,bestfit_params,tmin,tmax,
                                   showapproach,meff_tstep,chisq_dof,qual,fitinfo);
       MCObsInfo energy_key=fitinfo.energy_key;
       if (scattering_particles.size()>0){
          if (aniso_obsinfo.isVacuum()) // no anisotropy
            doEnergyDifferenceBySamplings(*m_obs,energy_key,scattering_particles,energy_key);
          else
            doEnergyDifferenceBySamplings(*m_obs,energy_key,aniso_obsinfo,scattering_particles,energy_key);}

       bool correlated=false;
       if (!chosen_fit_info.isVacuum()){
          MCObsInfo diff_obs;
          doCorrelatedDifferenceBySamplings(*m_obs,chosen_fit_info,energy_key,diff_obs);
          MCEstimate diff_est=m_obs->getEstimate(diff_obs);
          double diff_val=diff_est.getFullEstimate();
          double diff_up,diff_down;
          if (diff_est.isJackknifeMode())
             diff_up=diff_down=correlatedthreshold*diff_est.getSymmetricError();
          else{
             diff_up=correlatedthreshold*(diff_est.getUpperConfLimit()-diff_val);
             diff_down=correlatedthreshold*(diff_val-diff_est.getLowerConfLimit());}

          double upper_limit=diff_up+diff_val;
          double lower_limit=diff_val-diff_down;
          correlated = upper_limit >= 0. && lower_limit <= 0.;}

       MCEstimate energy=m_obs->getEstimate(energy_key);
       double y=energy.getFullEstimate();
       double dyup,dydn;
       if (energy.isJackknifeMode()) 
          dyup=dydn=energy.getSymmetricError();
       else{
          dyup=energy.getUpperConfLimit()-y;
          dydn=y-energy.getLowerConfLimit();}
       if (qual>=0.1 && correlated) goodcorrelatedfits.push_back(XYDYDYPoint(tmax,y,dyup,dydn));
       else if (qual>=0.1 && !correlated) gooduncorrelatedfits.push_back(XYDYDYPoint(tmax,y,dyup,dydn));
       else if (qual<0.1 && correlated) badcorrelatedfits.push_back(XYDYDYPoint(tmax,y,dyup,dydn));
       else baduncorrelatedfits.push_back(XYDYDYPoint(tmax,y,dyup,dydn));}
       catch(const std::exception& errmsg){
          xmlout.put_child("Error",string("DoFit within type TemporalCorrelatorTmaxVary encountered an error: ")
                 +string(errmsg.what()));
       }}
    XMLHandler xmlplog("TmaxPlot");
    xmlplog.put_child("PlotFile",plotfile);
    xmlplog.put_child("QualityThreshold",make_string(qualthreshold));
    xmlplog.put_child("CorrelatedThreshold",make_string(correlatedthreshold));
    xmlplog.put_child("NumberOfGoodCorrelatedFitPoints",make_string(goodcorrelatedfits.size()));
    xmlplog.put_child("NumberOfGoodUncorrelatedFitPoints",make_string(gooduncorrelatedfits.size()));
    xmlplog.put_child("NumberOfBadCorrelatedFitPoints",make_string(badcorrelatedfits.size()));
    xmlplog.put_child("NumberOfBadUncorrelatedFitPoints",make_string(baduncorrelatedfits.size()));
    xmlout.put_child(xmlplog);

    XYDYDYPoint chosen_fit(0,0,0,0);
    if (!chosen_fit_info.isVacuum()){
       MCEstimate chosen_fit_estimate=m_obs->getEstimate(chosen_fit_info);
       double y=chosen_fit_estimate.getFullEstimate();
       double dyup,dydn;
       if (chosen_fit_estimate.isJackknifeMode()) 
          dyup=dydn=chosen_fit_estimate.getSymmetricError();
       else{
          dyup=chosen_fit_estimate.getUpperConfLimit()-y;
          dydn=y-chosen_fit_estimate.getLowerConfLimit();}
       chosen_fit = XYDYDYPoint(1,y,dyup,dydn);}
    createTMinPlot(goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits,
                   corrname,plotfile,symbol,goodfitcolor,badfitcolor,correlatedfit_hollow,
                   uncorrelatedfit_hollow,chosen_fit,false,true);}
    catch(const std::exception& errmsg){
       xmlout.put_child("Error",string("DoFit with type TemporalCorrelatorTmaxVary encountered an error: ")
               +string(errmsg.what()));
    }}
    
 else if (fittype=="TemporalCorrelatorInteractionRatioTmaxVary"){
    try{
    XMLHandler xmlf(xmltask,"TemporalCorrelatorInteractionRatioTmaxVaryFit");
    XMLHandler xmlres(xmlf,"Ratio");
    OperatorInfo ratio_op(xmlres);
    XMLHandler xmlint(xmlf,"InteractingOperator");
    bool numvev=(xmlint.count("SubtractVEV")>0) ? true: false;
    pair<OperatorInfo,bool> numerator=make_pair(OperatorInfo(xmlint),numvev);
    vector<pair<OperatorInfo,bool> > denominator;
    list<XMLHandler> denomxml=xmlf.find_among_children("NonInteractingOperator");
    for (list<XMLHandler>::iterator it=denomxml.begin();it!=denomxml.end();++it){
      OperatorInfo opinfo(*it);
      bool subvev=(it->count("SubtractVEV")>0) ? true: false;
      denominator.push_back(make_pair(opinfo,subvev));}
    uint nterms=denominator.size();
    if (nterms<2) throw(std::invalid_argument("Two or more NonInteractingOperators required"));
    XMLHandler xmlo, xmldp;
    xmlo.set_root("Ratio");
    ratio_op.output(xmldp);
    xmlo.put_child(xmldp);
    xmlout.put_child(xmlo);
    xmlo.set_root("InteractingOperator");
    numerator.first.output(xmldp);
    xmlo.put_child(xmldp);
    xmlout.put_child(xmlo);
    xmlo.set_root("NonInteractingOperators");
    for (vector<pair<OperatorInfo,bool> >::const_iterator
           it=denominator.begin();it!=denominator.end();it++){
       XMLHandler xmloo; it->first.output(xmloo); xmlo.put_child(xmloo);}
    xmlout.put_child(xmlo);
    set<MCObsInfo> obskeys;
    bool erase_orig=true;
    uint tmaxfirst,tmaxlast,tmin;
    xmlread(xmlf,"TmaxFirst",tmaxfirst,"TemporalCorrelatorInteractionRatioTmaxVary");
    xmlread(xmlf,"TmaxLast",tmaxlast,"TemporalCorrelatorInteractionRatioTmaxVary");
    xmlread(xmlf,"Tmin",tmin,"TemporalCorrelatorInteractionRatioTmaxVary");

    doCorrelatorInteractionRatioBySamplings(*m_obs,numerator,denominator,
                                            0,(tmaxlast<64)?64:tmaxlast,ratio_op,obskeys,erase_orig);


    XMLHandler xmltf(xmlf,XMLHandler::subtree_copy);
    xmltf.rename_tag("TemporalCorrelatorFit");
    xmltf.put_child("MaximumTimeSeparation",make_string(tmaxfirst));
    xmltf.put_child("MinimumTimeSeparation",make_string(tmin));
    XMLHandler xmlro; ratio_op.output(xmlro);
    xmltf.put_child(xmlro); 

    list<pair<MCObsInfo,double> > scattering_particles;
    MCObsInfo aniso_obsinfo;
    if (xmltask.count_among_children("DoReconstructEnergy")==1){
       XMLHandler xmled(xmltask,"DoReconstructEnergy");
       uint num_spatial_sites=0;
       xmlread(xmled,"SpatialExtentNumSites",num_spatial_sites,"TemporalCorrelatorTmaxVary");
       double m_momsq_quantum=6.2831853071795864770/double(num_spatial_sites);
       m_momsq_quantum*=m_momsq_quantum;
       if (xmled.count_to_among_children("Anisotropy")==1){
         XMLHandler xmlani(xmled,"Anisotropy");
         string name;
         xmlread(xmlani,"Name",name,"Anisotropy");
         int index=0;
         xmlreadifchild(xmlani,"IDIndex",index);
         aniso_obsinfo = MCObsInfo(name,index);}
       list<XMLHandler> scattering_xml=xmled.find_among_children("ScatteringParticleEnergyFit");
       for (list<XMLHandler>::iterator st=scattering_xml.begin();st!=scattering_xml.end();++st){
          uint psq;
          xmlread(*st,"IntMomSquared",psq,"ScatteringParticleEnergyFit");
          double psqfactor=psq*m_momsq_quantum;
          string name;
          xmlread(*st,"Name",name,"ScatteringParticleEnergyFit");
          int index=0;
          xmlreadifchild(*st,"IDIndex",index);
          scattering_particles.push_back(make_pair(MCObsInfo(name,index),psqfactor));}}

    MCObsInfo chosen_fit_info;
    if (xmltask.count_among_children("ChosenFitInfo")==1){
       XMLHandler xmlchosen(xmltask,"ChosenFitInfo");
       string name;
       xmlread(xmlchosen,"Name",name,"ChosenFitInfo");
       int index=0;
       xmlreadifchild(xmlchosen,"IDIndex",index);
       chosen_fit_info = MCObsInfo(name,index);}

    XMLHandler xmlp(xmltask,"PlotInfo");
    string plotfile;
    xmlread(xmlp,"PlotFile",plotfile,"TemporalCorrelatorInteractionRatioTmaxVary");
    if (plotfile.empty()) throw(std::invalid_argument("Must have plot file name"));
    string corrname("standard");
    xmlreadif(xmlp,"CorrName",corrname,"TemporalCorrelatorInteractionRatioTmaxVary");
    string symbol("circle");
    xmlreadif(xmlp,"SymbolType",symbol,"TemporalCorrelatorInteractionRatioTmaxVary");
    double qualthreshold=0.1;
    xmlreadif(xmlp,"QualityThreshold",qualthreshold,"TemporalCorrelatorInteractionRatioTmaxVary");
    string goodfitcolor("blue");
    xmlreadif(xmlp,"GoodFitSymbolColor",goodfitcolor,"TemporalCorrelatorInteractionRatioTmaxVary");
    string badfitcolor("red");
    xmlreadif(xmlp,"BadFitSymbolColor",badfitcolor,"TemporalCorrelatorInteractionRatioTmaxVary");
    double correlatedthreshold=1.0;
    xmlreadif(xmlp,"CorrelatedThreshold",correlatedthreshold,"TemporalCorrelatorTmaxVary");
    bool correlatedfit_hollow=false;
    if (xml_child_tag_count(xmlp,"CorrelatedFitSymbolHollow")>0) correlatedfit_hollow=true;
    bool uncorrelatedfit_hollow=false;
    if (xml_child_tag_count(xmlp,"UncorrelatedFitSymbolHollow")>0) uncorrelatedfit_hollow=true;
    vector<XYDYDYPoint> goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits;
    for (uint tmax=tmaxfirst;tmax<=tmaxlast;++tmax){
       xmltf.seek_unique("MaximumTimeSeparation");
       xmltf.seek_next_node();       
       xmltf.set_text_content(make_string(tmax)); 
       try{
       RealTemporalCorrelatorFit RTC(xmltf,*m_obs,taskcount);
       CorrelatorInfo corr(RTC.m_op,RTC.m_op);
       if (corrname=="standard") corrname=getCorrelatorStandardName(corr);
       const vector<uint>& tvalues=RTC.getTvalues();
       if (find(tvalues.begin(),tvalues.end(),tmin)==tvalues.end()) continue;
       int dof = tvalues.size() - RTC.m_model_ptr->getNumberOfParams();
       if (dof < 1) continue;
       const vector<MCObsInfo>& fitparam_infos=RTC.getFitParamInfos();
       for (uint k=0;k<fitparam_infos.size();++k)
          m_obs->eraseSamplings(fitparam_infos[k]);
       XMLHandler xmlof; RTC.output(xmlof);
       xmlof.rename_tag("TemporalCorrelatorInteractionRatioTmaxVaryFit");
       xmlout.put_child(xmlof);
       double chisq_dof,qual;
       doChiSquareFitting(RTC,mz_info,chisq_dof,qual,bestfit_params,xmlout);
       TCorrFitInfo fitinfo;
       uint meff_tstep=1; bool showapproach=false;
       RTC.m_model_ptr->setFitInfo(RTC.m_fitparam_info,bestfit_params,tmin,tmax,
                                   showapproach,meff_tstep,chisq_dof,qual,fitinfo);
       MCObsInfo energy_key=fitinfo.energy_key;
       if (scattering_particles.size()>0){
          if (aniso_obsinfo.isVacuum()) // no anisotropy
            doReconstructEnergyBySamplings(*m_obs,energy_key,scattering_particles,energy_key);
          else
            doReconstructEnergyBySamplings(*m_obs,energy_key,aniso_obsinfo,scattering_particles,energy_key);}

       bool correlated=false;
       if (!chosen_fit_info.isVacuum()){
          MCObsInfo diff_obs;
          doCorrelatedDifferenceBySamplings(*m_obs,chosen_fit_info,energy_key,diff_obs);
          MCEstimate diff_est=m_obs->getEstimate(diff_obs);
          double diff_val=diff_est.getFullEstimate();
          double diff_up,diff_down;
          if (diff_est.isJackknifeMode())
             diff_up=diff_down=correlatedthreshold*diff_est.getSymmetricError();
          else{
             diff_up=correlatedthreshold*(diff_est.getUpperConfLimit()-diff_val);
             diff_down=correlatedthreshold*(diff_val-diff_est.getLowerConfLimit());}

          double upper_limit=diff_up+diff_val;
          double lower_limit=diff_val-diff_down;
          correlated = upper_limit >= 0. && lower_limit <= 0.;}

       MCEstimate energy=m_obs->getEstimate(energy_key);
       double y=energy.getFullEstimate();
       double dyup,dydn;
       if (energy.isJackknifeMode()) 
          dyup=dydn=energy.getSymmetricError();
       else{
          dyup=energy.getUpperConfLimit()-y;
          dydn=y-energy.getLowerConfLimit();}
       if (qual>=0.1 && correlated) goodcorrelatedfits.push_back(XYDYDYPoint(tmax,y,dyup,dydn));
       else if (qual>=0.1 && !correlated) gooduncorrelatedfits.push_back(XYDYDYPoint(tmax,y,dyup,dydn));
       else if (qual<0.1 && correlated) badcorrelatedfits.push_back(XYDYDYPoint(tmax,y,dyup,dydn));
       else baduncorrelatedfits.push_back(XYDYDYPoint(tmax,y,dyup,dydn));}
       catch(const std::exception& errmsg){
          xmlout.put_child("Error",string("DoFit within type TemporalCorrelatorInteractionRatioTmaxVary encountered an error: ")
                 +string(errmsg.what()));
       }}
    XMLHandler xmlplog("TmaxPlot");
    xmlplog.put_child("PlotFile",plotfile);
    xmlplog.put_child("QualityThreshold",make_string(qualthreshold));
    xmlplog.put_child("CorrelatedThreshold",make_string(correlatedthreshold));
    xmlplog.put_child("NumberOfGoodCorrelatedFitPoints",make_string(goodcorrelatedfits.size()));
    xmlplog.put_child("NumberOfGoodUncorrelatedFitPoints",make_string(gooduncorrelatedfits.size()));
    xmlplog.put_child("NumberOfBadCorrelatedFitPoints",make_string(badcorrelatedfits.size()));
    xmlplog.put_child("NumberOfBadUncorrelatedFitPoints",make_string(baduncorrelatedfits.size()));
    xmlout.put_child(xmlplog);

    XYDYDYPoint chosen_fit(0,0,0,0);
    if (!chosen_fit_info.isVacuum()){
       MCEstimate chosen_fit_estimate=m_obs->getEstimate(chosen_fit_info);
       double y=chosen_fit_estimate.getFullEstimate();
       double dyup,dydn;
       if (chosen_fit_estimate.isJackknifeMode()) 
          dyup=dydn=chosen_fit_estimate.getSymmetricError();
       else{
          dyup=chosen_fit_estimate.getUpperConfLimit()-y;
          dydn=y-chosen_fit_estimate.getLowerConfLimit();}
       chosen_fit = XYDYDYPoint(1,y,dyup,dydn);}
    createTMinPlot(goodcorrelatedfits,gooduncorrelatedfits,badcorrelatedfits,baduncorrelatedfits,
                   corrname,plotfile,symbol,goodfitcolor,badfitcolor,correlatedfit_hollow,
                   uncorrelatedfit_hollow,chosen_fit,false,true);}
    catch(const std::exception& errmsg){
       xmlout.put_child("Error",string("DoFit with type TemporalCorrelatorTmaxVary encountered an error: ")
               +string(errmsg.what()));
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
    doAnisoDispersionBySamplings(*m_obs,AFD.getAnisotropyKey(),AFD.getRestMassSquaredKey(), 
                            AFD.m_momsq_quantum*AFD.m_imomsq[kmin],randtemp);
    MCEstimate fit1=m_obs->getEstimate(randtemp);
    upperfit[0].xval=AFD.m_imomsq[kmin];
    upperfit[0].yval=fit1.getFullEstimate()+fit1.getSymmetricError();
    lowerfit[0].xval=AFD.m_imomsq[kmin];
    lowerfit[0].yval=fit1.getFullEstimate()-fit1.getSymmetricError();
    m_obs->eraseSamplings(randtemp);
    doAnisoDispersionBySamplings(*m_obs,AFD.getAnisotropyKey(),AFD.getRestMassSquaredKey(), 
                            AFD.m_momsq_quantum*AFD.m_imomsq[kmax],randtemp);
    MCEstimate fit2=m_obs->getEstimate(randtemp);
    upperfit[1].xval=AFD.m_imomsq[kmax];
    upperfit[1].yval=fit2.getFullEstimate()+fit2.getSymmetricError();
    lowerfit[1].xval=AFD.m_imomsq[kmax];
    lowerfit[1].yval=fit2.getFullEstimate()-fit2.getSymmetricError();
    m_obs->eraseSamplings(randtemp);

    createAnisoEnergyDispersionPlot(Esq,xiest.getFullEstimate(),xiest.getSymmetricError(),
                               goodtype,goodness,particlename,lowerfit,upperfit,
                               plotfile,symboltype,symbolcolor);
    }
    catch(const std::exception& errmsg){
       xmlout.put_child("Error",string("DoFit with type AnisotropyFromDispersion encountered an error: ")
               +string(errmsg.what()));
    }}

 else if (fittype=="Dispersion"){
    try{
    XMLHandler xmlf(xmltask,"DispersionFit");
    DispersionFit DF(xmlf,*m_obs,taskcount);
    XMLHandler xmlof; DF.output(xmlof);
    xmlout.put_child(xmlof);
    double chisq_dof,qual;
    doChiSquareFitting(DF,mz_info,chisq_dof,qual,
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
    MCEstimate cest=m_obs->getEstimate(DF.getCoefficientKey());
         // do some XML output
    string particlename;
    xmlreadifchild(xmlp,"ParticleName",particlename);
    xmlout.put_child("PlotFile",plotfile);

    vector<XYDYPoint> Esq(DF.m_nobs);
    vector<XYPoint> upperfit(2), lowerfit(2);
    uint kmin=0, kmax=0;
    for (uint k=0;k<DF.m_nobs;k++){
       MCEstimate est=m_obs->getEstimate(DF.m_obs_info[k]);
       Esq[k].xval=DF.m_imomsq[k];
       Esq[k].yval=est.getFullEstimate();
       Esq[k].yerr=est.getSymmetricError();
       if (Esq[k].xval<Esq[kmin].xval) kmin=k;
       if (Esq[k].xval>Esq[kmax].xval) kmax=k;}
    MCObsInfo randtemp("RandomTemporary",0);
    doCoeffDispersionBySamplings(*m_obs,DF.getCoefficientKey(),DF.getRestMassSquaredKey(), 
                            DF.m_momsq_quantum*DF.m_imomsq[kmin],randtemp);
    MCEstimate fit1=m_obs->getEstimate(randtemp);
    upperfit[0].xval=DF.m_imomsq[kmin];
    upperfit[0].yval=fit1.getFullEstimate()+fit1.getSymmetricError();
    lowerfit[0].xval=DF.m_imomsq[kmin];
    lowerfit[0].yval=fit1.getFullEstimate()-fit1.getSymmetricError();
    m_obs->eraseSamplings(randtemp);
    doCoeffDispersionBySamplings(*m_obs,DF.getCoefficientKey(),DF.getRestMassSquaredKey(), 
                            DF.m_momsq_quantum*DF.m_imomsq[kmax],randtemp);
    MCEstimate fit2=m_obs->getEstimate(randtemp);
    upperfit[1].xval=DF.m_imomsq[kmax];
    upperfit[1].yval=fit2.getFullEstimate()+fit2.getSymmetricError();
    lowerfit[1].xval=DF.m_imomsq[kmax];
    lowerfit[1].yval=fit2.getFullEstimate()-fit2.getSymmetricError();
    m_obs->eraseSamplings(randtemp);

    createCoeffEnergyDispersionPlot(Esq,cest.getFullEstimate(),cest.getSymmetricError(),
                               goodtype,goodness,particlename,lowerfit,upperfit,
                               plotfile,symboltype,symbolcolor);
    }
    catch(const std::exception& errmsg){
       xmlout.put_child("Error",string("DoFit with type Dispersion encountered an error: ")
               +string(errmsg.what()));
    }}

  /*
 else if (fittype=="LatticeDispersionRelation"){
    try{
    XMLHandler xmlf(xmltask,"LatticeDispersionRelationFit");
    LatticeDispersionFit LDF(xmlf,*m_obs,taskcount);
    XMLHandler xmlof; LDF.output(xmlof);
    xmlout.put_child(xmlof);
    double chisq_dof,qual;
    doChiSquareFitting(LDF,mz_info,chisq_dof,qual,
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
    doAnisoDispersionBySamplings(*m_obs,AFD.getAnisotropyKey(),AFD.getRestMassSquaredKey(), 
                            AFD.m_momsq_quantum*AFD.m_imomsq[kmin],randtemp);
    MCEstimate fit1=m_obs->getEstimate(randtemp);
    upperfit[0].xval=AFD.m_imomsq[kmin];
    upperfit[0].yval=fit1.getFullEstimate()+fit1.getSymmetricError();
    lowerfit[0].xval=AFD.m_imomsq[kmin];
    lowerfit[0].yval=fit1.getFullEstimate()-fit1.getSymmetricError();
    m_obs->eraseSamplings(randtemp);
    doAnisoDispersionBySamplings(*m_obs,AFD.getAnisotropyKey(),AFD.getRestMassSquaredKey(), 
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
       xmlout.put_child("Error",string("DoFit with type AnisotropyFromDispersion encountered an error: ")
               +string(errmsg.what()));
    }}
    */


}
// ***************************************************************************************
 

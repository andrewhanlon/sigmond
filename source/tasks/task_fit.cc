#include "task_handler.h"
#include "chisq_anisotropy.h"
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
// *       <PlotInfo>                                                            *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <ObsName>standard</ObsName>   (optional)                           *
// *          <SymbolType> ... </SymbolType>                                     *
// *          <GoodFitSymbolColor> ... </GoodFitSymbolColor>                     *
// *          <BadFitSymbolColor> ... </BadFitSymbolColor>                       *
// *          <GoodFitSymbolHollow/>  (optional)                                 *
// *          <BadFitSymbolHollow/>  (optional)                                  *
// *          <QualityThreshold>qual</QualityThreshold>  (0.1 default)           *
// *          <ChosenFitTmin>22</ChosenFitTmin> (optional)                       *
// *          <ChosenFitSymbolColor>black</ChosenFitSymbolColor> (optional)      *
// *          <ChosenFitDrawLines/>  (optional)                                  *
// *          <PrintChosenValue/> (optional)                                     *
// *       </PlotInfo>                                                           *
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
// *       <PlotInfo>                                                            *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <ObsName>standard</ObsName>   (optional)                           *
// *          <SymbolType> ... </SymbolType>                                     *
// *          <GoodFitSymbolColor> ... </GoodFitSymbolColor>                     *
// *          <BadFitSymbolColor> ... </BadFitSymbolColor>                       *
// *          <GoodFitSymbolHollow/>  (optional)                                 *
// *          <BadFitSymbolHollow/>  (optional)                                  *
// *          <QualityThreshold>qual</QualityThreshold>  (0.1 default)           *
// *          <ChosenFitTmin>22</ChosenFitTmin> (optional)                       *
// *          <ChosenFitSymbolColor>black</ChosenFitSymbolColor> (optional)      *
// *          <ChosenFitDrawLines/>  (optional)                                  *
// *          <PrintChosenValue/> (optional)                                     *
// *       </PlotInfo>                                                           *
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
// *    <Task>                                                                   *
// *     <Action>DoFit</Action>                                                  *
// *       <Type>TemporalCorrelatorInteractiongRatioTminVary</Type>              *
// *       <MinimizerInfo>                 (optional)                            *
// *         <Method>Minuit2</Method>                                            *
// *         <ParameterRelTol>1e-6</ParameterRelTol>                             *
// *         <ChiSquareRelTol>1e-4</ChiSquareRelTol>                             *
// *         <MaximumIterations>1024</MaximumIterations>                         *
// *         <Verbosity>Low</Verbosity>                                          *
// *       </MinimizerInfo>                                                      *
// *       <SamplingMode>Bootstrap</SamplingMode>   (optional)                   *
// *       <CovMatCalcSamplingMode>Bootstrap</CovMatCalcSamplingMode> (optional) *
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
// *       <PlotInfo>                                                            *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <ObsName>standard</ObsName>   (optional)                           *
// *          <SymbolType> ... </SymbolType>                                     *
// *          <GoodFitSymbolColor> ... </GoodFitSymbolColor>                     *
// *          <BadFitSymbolColor> ... </BadFitSymbolColor>                       *
// *          <GoodFitSymbolHollow/>  (optional)                                 *
// *          <BadFitSymbolHollow/>  (optional)                                  *
// *          <QualityThreshold>qual</QualityThreshold>  (0.1 default)           *
// *          <ChosenFitTmin>22</ChosenFitTmin> (optional)                       *
// *          <ChosenFitSymbolColor>black</ChosenFitSymbolColor> (optional)      *
// *          <ChosenFitDrawLines/>  (optional)                                  *
// *          <PrintChosenValue/> (optional)                                     *
// *       </PlotInfo>                                                           *
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
    double maxrelerror=0.0;
    if (xmlreadifchild(xmltask,"MaxRelativeErrorToPlot",maxrelerror)){
       map<double,MCEstimate> raw(results);
       results.clear();
       for (map<double,MCEstimate>::const_iterator it=raw.begin();it!=raw.end();it++)
          if ((it->second).getRelativeError()<std::abs(maxrelerror)) results.insert(*it);}

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
    uint efftype=RTC.m_model1_ptr->getEffMassType();
    double subt_const=0.0;
    if (efftype>1){    // subtract fit constant
       efftype-=2;     // efftypes 2 and 3 remove constant, but noisy
       subt_const=bestfit_params[bestfit_params.size()-1].getFullEstimate();}
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
    double maxrelerror=0.0;
    if (xmlreadifchild(xmltask,"MaxRelativeErrorToPlot",maxrelerror)){
       map<double,MCEstimate> raw(results);
       results.clear();
       for (map<double,MCEstimate>::const_iterator it=raw.begin();it!=raw.end();it++)
          if ((it->second).getRelativeError()<std::abs(maxrelerror)) results.insert(*it);}

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
       //xmlout.clear();
       xmlout.put_child("Error",string("DoFit with type TwoTemporalCorrelator encountered an error: ")
               +string(errmsg.what()));
      // throw(std::invalid_argument("DoFit with type TemporalCorrelator encountered an error: ")+errmsg);}
    }}


 else if (fittype=="LogTemporalCorrelator"){
    try{
    XMLHandler xmlf(xmltask,"LogTemporalCorrelatorFit");
    LogRealTemporalCorrelatorFit RTC(xmlf,*m_obs,taskcount);
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
    double subt_const=0.0;
    SamplingMode mode=m_obs->getCurrentSamplingMode();

    map<double,MCEstimate> results;
    getEffectiveEnergy(m_obs,corr,hermitian,subvev,RealPart,mode,step, 
                       0,results,subt_const);
    if (results.empty()){
       xmlout.put_child("PlotError","No effective energy estimates could be obtained");
       return;}
         // do some XML output
    xmlout.put_child("PlotFile",plotfile);
    XMLHandler xmlef;
    xmlef.set_root("EffectiveEnergy");
    xmlef.put_child("TimeStep",make_string(step));
    xmlef.put_child("EffEnergyType","TimeForward");
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
    double maxrelerror=0.0;
    if (xmlreadifchild(xmltask,"MaxRelativeErrorToPlot",maxrelerror)){
       map<double,MCEstimate> raw(results);
       results.clear();
       for (map<double,MCEstimate>::const_iterator it=raw.begin();it!=raw.end();it++)
          if ((it->second).getRelativeError()<std::abs(maxrelerror)) results.insert(*it);}

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
    XMLHandler xmlp(xmltask,"PlotInfo");
    string plotfile;
    xmlread(xmlp,"PlotFile",plotfile,"TemporalCorrelatorTminVary");
    if (plotfile.empty()) throw(std::invalid_argument("Must have plot file name"));
    string obsname("standard");
    xmlreadif(xmlp,"ObsName",obsname,"TemporalCorrelatorTminVary");
    string symbol("circle");
    xmlreadif(xmlp,"SymbolType",symbol,"TemporalCorrelatorTminVary");
    double qualthreshold=0.1;
    xmlreadif(xmlp,"QualityThreshold",qualthreshold,"TemporalCorrelatorTminVary");
    string goodfitcolor("blue");
    xmlreadif(xmlp,"GoodFitSymbolColor",goodfitcolor,"TemporalCorrelatorTminVary");
    string badfitcolor("red");
    xmlreadif(xmlp,"BadFitSymbolColor",badfitcolor,"TemporalCorrelatorTminVary");
    bool badfit_hollow=false;
    if (xml_child_tag_count(xmlp,"BadFitSymbolHollow")>0) badfit_hollow=true;
    bool goodfit_hollow=false;
    if (xml_child_tag_count(xmlp,"GoodFitSymbolHollow")>0) goodfit_hollow=true;
    int tmin_chosen_fit=-1;
    xmlreadif(xmlp,"ChosenFitTmin",tmin_chosen_fit,"TemporalCorrelatorTminVary");
    string chosenfitcolor="black";
    xmlreadif(xmlp,"ChosenFitSymbolColor",chosenfitcolor,"TemporalCorrelatorTminVary");
    bool chosen_fit_lines=false;
    if (xml_child_tag_count(xmlp,"ChosenFitDrawLines")>0) chosen_fit_lines=true;
    bool print_chosen_value=false;
    if (xml_child_tag_count(xmlp,"PrintChosenValue")>0) print_chosen_value=true;
    vector<XYDYDYPoint> goodfits,badfits;
    for (uint tmin=tminfirst;tmin<=tminlast;++tmin){
       xmltf.seek_unique("MinimumTimeSeparation");
       xmltf.seek_next_node();       
       xmltf.set_text_content(make_string(tmin)); 
       RealTemporalCorrelatorFit RTC(xmltf,*m_obs,taskcount);
       if (obsname=="standard"){
         CorrelatorInfo corr(RTC.m_op,RTC.m_op);
         obsname=getCorrelatorStandardName(corr);}
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
       try{
       doChiSquareFitting(RTC,mz_info,chisq_dof,qual,
                          bestfit_params,xmlout);
       TCorrFitInfo fitinfo;
       uint meff_tstep=1; bool showapproach=false;
       RTC.m_model_ptr->setFitInfo(RTC.m_fitparam_info,bestfit_params,tmin,tmax,
                                   showapproach,meff_tstep,chisq_dof,qual,fitinfo);
       MCEstimate energy=m_obs->getEstimate(fitinfo.energy_key);
       double y=energy.getFullEstimate();
       double dyup,dydn;
       if (energy.isJackknifeMode()) 
          dyup=dydn=energy.getSymmetricError();
       else{
          dyup=energy.getUpperConfLimit()-y;
          dydn=y-energy.getLowerConfLimit();}
       if (qual>0.1) goodfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else badfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));}
       catch(const std::exception& xp){}}
    XMLHandler xmlplog("TminPlot");
    xmlplog.put_child("PlotFile",plotfile);
    xmlplog.put_child("QualityThreshold",make_string(qualthreshold));
    xmlplog.put_child("NumberOfBadFitPoints",make_string(badfits.size()));
    xmlplog.put_child("NumberOfGoodFitPoints",make_string(goodfits.size()));
    if (tmin_chosen_fit>0) xmlplog.put_child("ChosenFitTmin",make_string(tmin_chosen_fit));
    xmlout.put_child(xmlplog);
    createTMinPlot(goodfits,badfits,obsname,plotfile,symbol,goodfitcolor,
                   badfitcolor,goodfit_hollow,badfit_hollow,tmin_chosen_fit,
                   chosenfitcolor,chosen_fit_lines,print_chosen_value);}
    catch(const std::exception& errmsg){
       //xmlout.clear();
       xmlout.put_child("Error",string("DoFit with type TemporalCorrelatorTminVary encountered an error: ")
               +string(errmsg.what()));
      // throw(std::invalid_argument("DoFit with type TemporalCorrelator encountered an error: ")+errmsg);}
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
    XMLHandler xmlp(xmltask,"PlotInfo");
    string plotfile;
    xmlread(xmlp,"PlotFile",plotfile,"LogTemporalCorrelatorTminVary");
    if (plotfile.empty()) throw(std::invalid_argument("Must have plot file name"));
    string obsname("standard");
    xmlreadif(xmlp,"ObsName",obsname,"LogTemporalCorrelatorTminVary");
    string symbol("circle");
    xmlreadif(xmlp,"SymbolType",symbol,"LogTemporalCorrelatorTminVary");
    double qualthreshold=0.1;
    xmlreadif(xmlp,"QualityThreshold",qualthreshold,"LogTemporalCorrelatorTminVary");
    string goodfitcolor("blue");
    xmlreadif(xmlp,"GoodFitSymbolColor",goodfitcolor,"LogTemporalCorrelatorTminVary");
    string badfitcolor("red");
    xmlreadif(xmlp,"BadFitSymbolColor",badfitcolor,"LogTemporalCorrelatorTminVary");
    bool badfit_hollow=false;
    if (xml_child_tag_count(xmlp,"BadFitSymbolHollow")>0) badfit_hollow=true;
    bool goodfit_hollow=false;
    if (xml_child_tag_count(xmlp,"GoodFitSymbolHollow")>0) goodfit_hollow=true;
    int tmin_chosen_fit=-1;
    xmlreadif(xmlp,"ChosenFitTmin",tmin_chosen_fit,"LogTemporalCorrelatorTminVary");
    string chosenfitcolor="black";
    xmlreadif(xmlp,"ChosenFitSymbolColor",chosenfitcolor,"LogTemporalCorrelatorTminVary");
    bool chosen_fit_lines=false;
    if (xml_child_tag_count(xmlp,"ChosenFitDrawLines")>0) chosen_fit_lines=true;
    bool print_chosen_value=false;
    if (xml_child_tag_count(xmlp,"PrintChosenValue")>0) print_chosen_value=true;
    vector<XYDYDYPoint> goodfits,badfits;
    for (uint tmin=tminfirst;tmin<=tminlast;++tmin){
       xmltf.seek_unique("MinimumTimeSeparation");
       xmltf.seek_next_node();       
       xmltf.set_text_content(make_string(tmin)); 
       LogRealTemporalCorrelatorFit RTC(xmltf,*m_obs,taskcount);
       if (obsname=="standard"){
         CorrelatorInfo corr(RTC.m_op,RTC.m_op);
         obsname=getCorrelatorStandardName(corr);}
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
       try{
       doChiSquareFitting(RTC,mz_info,chisq_dof,qual,
                          bestfit_params,xmlout);
       TCorrFitInfo fitinfo;
       uint meff_tstep=1; bool showapproach=false;
       RTC.m_model_ptr->setFitInfo(RTC.m_fitparam_info,bestfit_params,tmin,tmax,
                                   showapproach,meff_tstep,chisq_dof,qual,fitinfo);
       MCEstimate energy=m_obs->getEstimate(fitinfo.energy_key);
       double y=energy.getFullEstimate();
       double dyup,dydn;
       if (energy.isJackknifeMode()) 
          dyup=dydn=energy.getSymmetricError();
       else{
          dyup=energy.getUpperConfLimit()-y;
          dydn=y-energy.getLowerConfLimit();}
       if (qual>0.1) goodfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else badfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));}
       catch(const std::exception& xp){}}
    XMLHandler xmlplog("TminPlot");
    xmlplog.put_child("PlotFile",plotfile);
    xmlplog.put_child("QualityThreshold",make_string(qualthreshold));
    xmlplog.put_child("NumberOfBadFitPoints",make_string(badfits.size()));
    xmlplog.put_child("NumberOfGoodFitPoints",make_string(goodfits.size()));
    if (tmin_chosen_fit>0) xmlplog.put_child("ChosenFitTmin",make_string(tmin_chosen_fit));
    xmlout.put_child(xmlplog);
    createTMinPlot(goodfits,badfits,obsname,plotfile,symbol,goodfitcolor,
                   badfitcolor,goodfit_hollow,badfit_hollow,tmin_chosen_fit,
                   chosenfitcolor,chosen_fit_lines,print_chosen_value);}
    catch(const std::exception& errmsg){
       //xmlout.clear();
       xmlout.put_child("Error",string("DoFit with type LogTemporalCorrelatorTminVary encountered an error: ")
               +string(errmsg.what()));
      // throw(std::invalid_argument("DoFit with type TemporalCorrelator encountered an error: ")+errmsg);}
    }}


 else if (fittype=="TemporalCorrelatorInteractionRatio"){
    try{
     XMLHandler xmlf(xmltask,"TemporalCorrelatorInteractionRatioFit");
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
     uint tmin,tmax;
     xmlreadchild(xmlf,"MinimumTimeSeparation",tmin);
     xmlreadchild(xmlf,"MaximumTimeSeparation",tmax);
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
     doCorrelatorInteractionRatioBySamplings(*m_obs,numerator,denominator,
                                             0,(tmax<64)?64:tmax,ratio_op,obskeys,erase_orig);
     XMLHandler xmltf(xmlf,XMLHandler::subtree_copy);
     xmltf.rename_tag("TemporalCorrelatorFit");
     XMLHandler xmlro; ratio_op.output(xmlro);
     xmltf.put_child(xmlro); 
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
     double maxrelerror=0.0;
     if (xmlreadifchild(xmltask,"MaxRelativeErrorToPlot",maxrelerror)){
        map<double,MCEstimate> raw(results);
        results.clear();
        for (map<double,MCEstimate>::const_iterator it=raw.begin();it!=raw.end();it++)
           if ((it->second).getRelativeError()<std::abs(maxrelerror)) results.insert(*it);}

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
       //xmlout.clear();
       xmlout.put_child("Error",string("DoFit with type TemporalCorrelator encountered an error: ")
               +string(errmsg.what()));
      // throw(std::invalid_argument("DoFit with type TemporalCorrelator encountered an error: ")+errmsg);}
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

    XMLHandler xmlp(xmltask,"PlotInfo");
    string plotfile;
    xmlread(xmlp,"PlotFile",plotfile,"TemporalCorrelatorInteractionRatioTminVary");
    if (plotfile.empty()) throw(std::invalid_argument("Must have plot file name"));
    string obsname("standard");
    xmlreadif(xmlp,"ObsName",obsname,"TemporalCorrelatorInteractionRatioTminVary");
    string symbol("circle");
    xmlreadif(xmlp,"SymbolType",symbol,"TemporalCorrelatorInteractionRatioTminVary");
    double qualthreshold=0.1;
    xmlreadif(xmlp,"QualityThreshold",qualthreshold,"TemporalCorrelatorInteractionRatioTminVary");
    string goodfitcolor("blue");
    xmlreadif(xmlp,"GoodFitSymbolColor",goodfitcolor,"TemporalCorrelatorInteractionRatioTminVary");
    string badfitcolor("red");
    xmlreadif(xmlp,"BadFitSymbolColor",badfitcolor,"TemporalCorrelatorInteractionRatioTminVary");
    bool badfit_hollow=false;
    if (xml_child_tag_count(xmlp,"BadFitSymbolHollow")>0) badfit_hollow=true;
    bool goodfit_hollow=false;
    if (xml_child_tag_count(xmlp,"GoodFitSymbolHollow")>0) goodfit_hollow=true;
    int tmin_chosen_fit=-1;
    xmlreadif(xmlp,"ChosenFitTmin",tmin_chosen_fit,"TemporalCorrelatorInteractionRatioTminVary");
    string chosenfitcolor="black";
    xmlreadif(xmlp,"ChosenFitSymbolColor",chosenfitcolor,"TemporalCorrelatorInteractionRatioTminVary");
    bool chosen_fit_lines=false;
    if (xml_child_tag_count(xmlp,"ChosenFitDrawLines")>0) chosen_fit_lines=true;
    bool print_chosen_value=false;
    if (xml_child_tag_count(xmlp,"PrintChosenValue")>0) print_chosen_value=true;
    vector<XYDYDYPoint> goodfits,badfits;
    for (uint tmin=tminfirst;tmin<=tminlast;++tmin){
       xmltf.seek_unique("MinimumTimeSeparation");
       xmltf.seek_next_node();       
       xmltf.set_text_content(make_string(tmin)); 
       RealTemporalCorrelatorFit RTC(xmltf,*m_obs,taskcount);
       if (obsname=="standard"){
         CorrelatorInfo corr(RTC.m_op,RTC.m_op);
         obsname=getCorrelatorStandardName(corr);}
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
       try{
       doChiSquareFitting(RTC,mz_info,chisq_dof,qual,
                          bestfit_params,xmlout);
       TCorrFitInfo fitinfo;
       uint meff_tstep=1; bool showapproach=false;
       RTC.m_model_ptr->setFitInfo(RTC.m_fitparam_info,bestfit_params,tmin,tmax,
                                   showapproach,meff_tstep,chisq_dof,qual,fitinfo);
       MCEstimate energy=m_obs->getEstimate(fitinfo.energy_key);
       double y=energy.getFullEstimate();
       double dyup,dydn;
       if (energy.isJackknifeMode()) 
          dyup=dydn=energy.getSymmetricError();
       else{
          dyup=energy.getUpperConfLimit()-y;
          dydn=y-energy.getLowerConfLimit();}
       if (qual>0.1) goodfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));
       else badfits.push_back(XYDYDYPoint(tmin,y,dyup,dydn));}
       catch(const std::exception& xp){}}
    XMLHandler xmlplog("TminPlot");
    xmlplog.put_child("PlotFile",plotfile);
    xmlplog.put_child("QualityThreshold",make_string(qualthreshold));
    xmlplog.put_child("NumberOfBadFitPoints",make_string(badfits.size()));
    xmlplog.put_child("NumberOfGoodFitPoints",make_string(goodfits.size()));
    if (tmin_chosen_fit>0) xmlplog.put_child("ChosenFitTmin",make_string(tmin_chosen_fit));
    xmlout.put_child(xmlplog);
    createTMinPlot(goodfits,badfits,obsname,plotfile,symbol,goodfitcolor,
                   badfitcolor,goodfit_hollow,badfit_hollow,tmin_chosen_fit,
                   chosenfitcolor,chosen_fit_lines,print_chosen_value);}
    catch(const std::exception& errmsg){
       //xmlout.clear();
       xmlout.put_child("Error",string("DoFit with type TemporalCorrelatorTminVary encountered an error: ")
               +string(errmsg.what()));
      // throw(std::invalid_argument("DoFit with type TemporalCorrelator encountered an error: ")+errmsg);}
    }}


}


// ***************************************************************************************
 

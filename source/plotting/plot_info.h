#ifndef PLOT_INFO_H
#define PLOT_INFO_H

#include <string>
#include "xml_handler.h"
#include "mcobs_info.h"

// *******************************************************************************
// *                                                                             *
// *        <FitEffEnergyPlotInfo>                                               *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <CorrName>standard</CorrName>   (optional)                         *
// *          <TimeStep>3</TimeStep>  (optional: 1 default)                      *
// *          <SymbolColor> ... </SymbolColor>                                   *
// *          <SymbolType> ... </SymbolType>                                     *
// *          <MaxErrorToPlot> ...</MaxErrorToPlot> (optional)                   *
// *          <Goodness>qual</Goodness>  "qual" or "chisq"                       *
// *          <ReferenceEnergy> (optional: includes energy ratio on plot)        *
// *            <Name>kaon</Name><IDIndex>0</IDIndex>                            *
// *          </ReferenceEnergy>                                                 *
// *        </FitEffEnergyPlotInfo>                                              *
// *                                                                             *
// *******************************************************************************

struct FitEffEnergyPlotInfo
{
  std::string plotfile;
  std::string corrname = "standard";
  uint timestep = 1;
  std::string symbolcolor = "blue";
  std::string symboltype = "circle";
  double maxerror = 0.;
  std::string goodness = "chisq";
  MCObsInfo ref_energy;

  FitEffEnergyPlotInfo(XMLHandler& xmlin);

  FitEffEnergyPlotInfo(const std::string& in_plotfile) : plotfile(in_plotfile) {}
};


// *******************************************************************************
// *                                                                             *
// *        <DataFitRatioPlotInfo>                                               *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <CorrName>standard</CorrName>   (optional)                         *
// *          <SymbolColor> ... </SymbolColor>                                   *
// *          <SymbolType> ... </SymbolType>                                     *
// *          <MaxErrorToPlot> ...</MaxErrorToPlot> (optional)                   *
// *          <Goodness>qual</Goodness>  "qual" or "chisq"                       *
// *        </DataFitRatioPlotInfo>                                              *
// *                                                                             *
// *******************************************************************************

struct DataFitRatioPlotInfo
{
  std::string plotfile;
  std::string corrname = "standard";
  std::string symbolcolor = "blue";
  std::string symboltype = "circle";
  double maxerror = 0.;
  std::string goodness = "chisq";

  DataFitRatioPlotInfo(XMLHandler& xmlin);

  DataFitRatioPlotInfo(const std::string& in_plotfile) : plotfile(in_plotfile) {}
};


// *******************************************************************************
// *                                                                             *
// *        <TminFitPlotInfo>                                                    *
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
// *         <SubtractShift/>     (optional: subtracts shift rather than adds)   *
// *         <ShiftInfo>                                                         *
// *           <Name>shift_obsname</Name>                                        *
// *           <IDIndex>0</IDIndex>                                              *
// *         </ShiftInfo>                                                        *
// *         <ChosenFitInfo>           (optional)                                *
// *           <Name>fit_obsname</Name>                                          *
// *           <IDIndex>0</IDIndex>                                              *
// *         </ChosenFitInfo>                                                    *
// *        </TminFitPlotInfo>                                                   *
// *                                                                             *
// *******************************************************************************

struct TminFitPlotInfo
{
  std::string plotfile;
  uint energy_level = 0;
  std::string corrname = "standard";
  std::string symboltype = "circle";
  std::string goodfit_symbolcolor = "blue";
  std::string badfit_symbolcolor = "red";
  bool correlatedfit_hollow = true;
  bool uncorrelatedfit_hollow = true;
  double quality_threshold = 0.1;
  double correlated_threshold = 1.;
  bool subtract_shift = false;
  MCObsInfo shift_info;
  MCObsInfo chosen_fit;

  TminFitPlotInfo(XMLHandler& xmlin);
};


// *******************************************************************************
// *                                                                             *
// *        <AnisotropyFitPlotInfo>                                              *
// *           <PlotFile> ... </PlotFile>                                        *
// *           <ParticleName>pion</ParticleName>   (optional)                    *
// *           <SymbolColor> ... </SymbolColor>                                  *
// *           <SymbolType> ... </SymbolType>                                    *
// *           <Goodness>qual</Goodness>  "qual" or "chisq"                      *
// *        </AnisotropyFitPlotInfo>                                             *
// *                                                                             *
// *******************************************************************************

struct AnisotropyFitPlotInfo
{
  std::string plotfile;
  std::string particle_name;
  std::string symbolcolor = "blue";
  std::string symboltype = "circle";
  std::string goodness = "chisq";

  AnisotropyFitPlotInfo(XMLHandler& xmlin);
};

#endif

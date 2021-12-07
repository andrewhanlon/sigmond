#ifndef PLOT_INFO_H
#define PLOT_INFO_H

#include <string>
#include "xml_handler.h"
#include "mcobs_info.h"

// *******************************************************************************
// *                                                                             *
// *        <EffEnergyWithFitPlotInfo>                                           *
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
// *        </EffEnergyWithFitPlotInfo>                                          *
// *                                                                             *
// *******************************************************************************

struct EffEnergyWithFitPlotInfo
{
  std::string plotfile;
  std::string corrname = "standard";
  uint timestep = 1;
  std::string symbolcolor = "blue";
  std::string symboltype = "circle";
  double maxerror = 0.;
  std::string goodness = "chisq";
  MCObsInfo ref_energy;

  EffEnergyWithFitPlotInfo(XMLHandler& xmlin);

  EffEnergyWithFitPlotInfo(const std::string& in_plotfile) : plotfile(in_plotfile) {}
};

// *******************************************************************************
// *                                                                             *
// *        <ThreePointCorrelatorPlotInfo>                                       *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <PlotLable>Plot Label!</PlotLabel>   (optional)                    *
// *          <ComplexArg>RealPart</ComplexArg>   (default: RealPart)            *
// *          <Labels>label 1 | label 2</Labels>    (optional)                   *
// *          <TimeSeparations>6 8 10</TimeSeparations>    (optional)            *
// *          <SymbolColors> ... </SymbolColors>                                 *
// *          <SymbolTypes> ... </SymbolTypes>                                   *
// *        </ThreePointCorrelatorPlotInfo>                                      *
// *                                                                             *
// *******************************************************************************

struct ThreePointCorrelatorPlotInfo
{
  std::string plotfile;
  std::string plotlabel = "";
  std::vector<uint> time_seps = {};
  ComplexArg complex_arg = RealPart;
  std::vector<std::string> labels = {""};
  std::vector<std::string> symbolcolors = {"blue"};
  std::vector<std::string> symboltypes = {"circle"};

  ThreePointCorrelatorPlotInfo(XMLHandler& xmlin);

  ThreePointCorrelatorPlotInfo(const std::string& in_plotfile,
                               const std::vector<uint>& in_tseps) 
      : plotfile(in_plotfile), time_seps(in_tseps) {}
};

// *******************************************************************************
// *                                                                             *
// *        <ThreePointCorrelatorWithFitPlotInfo>                                *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <PlotLable>Plot Label!</PlotLabel>   (optional)                    *
// *          <Labels>label 1 | label 2</Labels>    (optional)                   *
// *          <TimeSeparations>6 8 10</TimeSeparations>    (optional)            *
// *          <Goodness>qual</Goodness>  "qual" or "chisq"                       *
// *          <SymbolColors> ... </SymbolColors>                                 *
// *          <SymbolTypes> ... </SymbolTypes>                                   *
// *        </ThreePointCorrelatorWithFitPlotInfo>                               *
// *                                                                             *
// *******************************************************************************

struct ThreePointCorrelatorWithFitPlotInfo
{
  std::string plotfile;
  std::string plotlabel = "";
  std::vector<uint> time_seps = {};
  std::vector<std::string> labels = {""};
  std::string goodness = "chisq";
  std::vector<std::string> symbolcolors = {"blue"};
  std::vector<std::string> symboltypes = {"circle"};

  ThreePointCorrelatorWithFitPlotInfo(XMLHandler& xmlin);

  ThreePointCorrelatorWithFitPlotInfo(const std::string& in_plotfile,
                                      const std::vector<uint>& in_tseps) 
      : plotfile(in_plotfile), time_seps(in_tseps) {}
};


// *******************************************************************************
// *                                                                             *
// *        <DataFitRatioPlotInfo>                                               *
// *          <PlotFile> ... </PlotFile>                                         *
// *          <PlotLable>Plot Label!</PlotLabel>   (optional)                    *
// *          <Labels>label 1 | label 2</Labels>    (optional)                   *
// *          <SymbolColors> black blue </SymbolColor>  (default: blue)          *
// *          <SymbolTypes> circle square </SymbolTypes>  (default: circle)      *
// *          <MaxErrorToPlot> ...</MaxErrorToPlot> (optional)                   *
// *        </DataFitRatioPlotInfo>                                              *
// *                                                                             *
// *******************************************************************************

struct DataFitRatioPlotInfo
{
  std::string plotfile;
  std::string plotlabel = "";
  std::vector<std::string> labels = {""};
  std::vector<std::string> symbolcolors = {"blue"};
  std::vector<std::string> symboltypes = {"circle"};
  double maxerror = 0.;

  DataFitRatioPlotInfo(XMLHandler& xmlin);

  DataFitRatioPlotInfo(const std::string& in_plotfile) : plotfile(in_plotfile) {}
};


// *******************************************************************************
// *                                                                             *
// *        <TminFitPlotInfo>                                                    *
// *         <PlotFile> ... </PlotFile>                                          *
// *         <EnergyLevel>1</EnergyLevel>   (0 default)                          *
// *         <PlotLabel>Plot Label!</PlotLabel>   (optional)                     *
// *         <Labels>label 1 | label 2</Labels>    (optional)                    *
// *         <SymbolColors> black blue </SymbolColor>  (default: blue)           *
// *         <SymbolTypes> circle square </SymbolTypes>  (default: circle)       *
// *         <QualityThreshold>0.05</QualityThreshold>  (0.1 default)            *
// *         <MaxErrorToPlot> ...</MaxErrorToPlot> (optional)                    *
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
  std::string plotlabel = "";
  std::vector<std::string> labels = {""};
  std::vector<std::string> symbolcolors = {"blue"};
  std::vector<std::string> symboltypes = {"circle"};
  double quality_threshold = 0.1;
  double maxerror = 0.;
  MCObsInfo chosen_fit;

  TminFitPlotInfo(XMLHandler& xmlin);

  TminFitPlotInfo(const std::string& in_plotfile) : plotfile(in_plotfile) {}
};


// *******************************************************************************
// *                                                                             *
// *        <DispersionFitPlotInfo>                                              *
// *           <PlotFile> ... </PlotFile>                                        *
// *           <ParticleName>pion</ParticleName>   (optional)                    *
// *           <SymbolColor> ... </SymbolColor>                                  *
// *           <SymbolType> ... </SymbolType>                                    *
// *           <Goodness>qual</Goodness>  "qual" or "chisq"                      *
// *        </DispersionFitPlotInfo>                                             *
// *                                                                             *
// *******************************************************************************

struct DispersionFitPlotInfo
{
  std::string plotfile;
  std::string particle_name = "";
  std::string symbolcolor = "blue";
  std::string symboltype = "circle";
  std::string goodness = "chisq";

  DispersionFitPlotInfo(XMLHandler& xmlin);

  DispersionFitPlotInfo(const std::string& in_plotfile) : plotfile(in_plotfile) {}
};

#endif

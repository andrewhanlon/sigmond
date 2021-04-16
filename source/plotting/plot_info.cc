#include "plot_info.h"

using namespace std;

FitEffEnergyPlotInfo::FitEffEnergyPlotInfo(XMLHandler& xmlin)
{
  xmlreadifchild(xmlin, "PlotFile", plotfile);
  plotfile = tidyString(plotfile);
  xmlreadifchild(xmlin, "CorrName", corrname);
  xmlreadifchild(xmlin, "TimeStep", timestep);
  xmlreadifchild(xmlin, "SymbolColor", symbolcolor);
  xmlreadifchild(xmlin, "SymbolType", symboltype);
  xmlreadifchild(xmlin, "MaxErrorToPlot", maxerror);
  xmlreadifchild(xmlin, "Goodness", goodness);
  if (xml_child_tag_count(xmlin, "ReferenceEnergy") == 1) {
    XMLHandler xmlref(xmlin, "ReferenceEnergy");
    string name;
    int index = 0;
    xmlreadchild(xmlref, "Name", name);
    xmlreadifchild(xmlref, "IDIndex", index);
    ref_energy = MCObsInfo(name, index);
  }
}

// *******************************************************************************

DataFitRatioPlotInfo::DataFitRatioPlotInfo(XMLHandler& xmlin)
{
  xmlreadifchild(xmlin, "PlotFile", plotfile);
  plotfile = tidyString(plotfile);
  xmlreadifchild(xmlin, "CorrName", corrname);
  xmlreadifchild(xmlin, "SymbolColor", symbolcolor);
  xmlreadifchild(xmlin, "SymbolType", symboltype);
  xmlreadifchild(xmlin, "MaxErrorToPlot", maxerror);
  xmlreadifchild(xmlin, "Goodness", goodness);
}

// *******************************************************************************

TminFitPlotInfo::TminFitPlotInfo(XMLHandler& xmlin)
{
  xmlreadifchild(xmlin, "PlotFile", plotfile);
  plotfile = tidyString(plotfile);
  xmlreadifchild(xmlin, "EnergyLevel", energy_level);
  xmlreadifchild(xmlin, "CorrName", corrname);
  xmlreadifchild(xmlin, "SymbolType", symboltype);
  xmlreadifchild(xmlin, "GoodFitSymbolColor", goodfit_symbolcolor);
  xmlreadifchild(xmlin, "BadFitSymbolColor", badfit_symbolcolor);
  correlatedfit_hollow = (xml_child_tag_count(xmlin, "CorrelatedFitSymbolHollow") > 0);
  uncorrelatedfit_hollow = (xml_child_tag_count(xmlin, "UncorrelatedFitSymbolHollow") > 0);
  xmlreadifchild(xmlin, "QualityThreshold", quality_threshold);
  xmlreadifchild(xmlin, "CorrelatedThreshold", correlated_threshold);
  if (xmlin.count_among_children("SubtractShift") > 0) subtract_shift = true;
  if (xml_child_tag_count(xmlin, "ShiftInfo") > 0) {
    XMLHandler xmlshift(xmlin, "ShiftInfo");
    string name;
    int index = 0;
    xmlreadchild(xmlshift, "Name", name);
    xmlreadifchild(xmlshift, "IDIndex", index);
    chosen_fit = MCObsInfo(name, index);
  }
  if (xml_child_tag_count(xmlin, "ChosenFitInfo") > 0) {
    XMLHandler xmlref(xmlin, "ChosenFitInfo");
    string name;
    int index = 0;
    xmlreadchild(xmlref, "Name", name);
    xmlreadifchild(xmlref, "IDIndex", index);
    chosen_fit = MCObsInfo(name, index);
  }
}

// *******************************************************************************

AnisotropyFitPlotInfo::AnisotropyFitPlotInfo(XMLHandler& xmlin)
{
  xmlreadifchild(xmlin, "PlotFile", plotfile);
  plotfile = tidyString(plotfile);
  xmlreadifchild(xmlin, "ParticleName", particle_name);
  xmlreadifchild(xmlin, "SymbolColor", symbolcolor);
  xmlreadifchild(xmlin, "SymbolType", symboltype);
  xmlreadifchild(xmlin, "Goodness", goodness);
}

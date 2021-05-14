#include <algorithm>
#include "plot_info.h"

using namespace std;

EffEnergyWithFitPlotInfo::EffEnergyWithFitPlotInfo(XMLHandler& xmlin)
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
  xmlreadifchild(xmlin, "PlotLabel", plotlabel);
  if (xml_child_tag_count(xmlin, "Labels")) {
    XMLHandler xml_labels(xmlin, "Labels");
    labels = ArgsHandler::split(xml_labels.get_text_content(), '|');
    for_each(labels.begin(), labels.end(), [](string &s){ tidyString(s); });
  }
  if (xml_child_tag_count(xmlin, "SymbolColors")) {
    XMLHandler xml_colors(xmlin, "SymbolColors");
    symbolcolors = ArgsHandler::split(xml_colors.get_text_content(), ' ');
    for_each(symbolcolors.begin(), symbolcolors.end(), [](string &s){ tidyString(s); });
  }
  if (xml_child_tag_count(xmlin, "SymbolTypes")) {
    XMLHandler xml_types(xmlin, "SymbolTypes");
    symboltypes = ArgsHandler::split(xml_types.get_text_content(), ' ');
    for_each(symboltypes.begin(), symboltypes.end(), [](string &s){ tidyString(s); });
  }
  xmlreadifchild(xmlin, "MaxErrorToPlot", maxerror);
}

// *******************************************************************************

TminFitPlotInfo::TminFitPlotInfo(XMLHandler& xmlin)
{
  xmlreadifchild(xmlin, "PlotFile", plotfile);
  plotfile = tidyString(plotfile);
  xmlreadifchild(xmlin, "EnergyLevel", energy_level);
  xmlreadifchild(xmlin, "PlotLabel", plotlabel);
  if (xml_child_tag_count(xmlin, "Labels")) {
    XMLHandler xml_labels(xmlin, "Labels");
    labels = ArgsHandler::split(xml_labels.get_text_content(), '|');
    for_each(labels.begin(), labels.end(), [](string &s){ tidyString(s); });
  }
  if (xml_child_tag_count(xmlin, "SymbolColors")) {
    XMLHandler xml_colors(xmlin, "SymbolColors");
    symbolcolors = ArgsHandler::split(xml_colors.get_text_content(), ' ');
    for_each(symbolcolors.begin(), symbolcolors.end(), [](string &s){ tidyString(s); });
  }
  if (xml_child_tag_count(xmlin, "SymbolTypes")) {
    XMLHandler xml_types(xmlin, "SymbolTypes");
    symboltypes = ArgsHandler::split(xml_types.get_text_content(), ' ');
    for_each(symboltypes.begin(), symboltypes.end(), [](string &s){ tidyString(s); });
  }
  xmlreadifchild(xmlin, "QualityThreshold", quality_threshold);
  xmlreadifchild(xmlin, "MaxErrorToPlot", maxerror);
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

DispersionFitPlotInfo::DispersionFitPlotInfo(XMLHandler& xmlin)
{
  xmlreadifchild(xmlin, "PlotFile", plotfile);
  plotfile = tidyString(plotfile);
  xmlreadifchild(xmlin, "ParticleName", particle_name);
  xmlreadifchild(xmlin, "SymbolColor", symbolcolor);
  xmlreadifchild(xmlin, "SymbolType", symboltype);
  xmlreadifchild(xmlin, "Goodness", goodness);
}

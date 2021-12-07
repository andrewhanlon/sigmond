#include "scalar_defs.h"
#include "ensemble_info.h"
#include "bins_info.h"
#include "bootstrapper.h"
#include "sampling_info.h"
#include "obs_get_handler.h"
#include "xml_handler.h"
#include "mcobs_info.h"
#include "mcobs_handler.h"
#include "operator_info.h"
#include "gen_irrep_operator_info.h"
#include "correlator_info.h"
#include "correlator_matrix_info.h"
#include "matrix.h"
#include "filelist_info.h"
#include "corr_data_handler.h"
#include "vev_data_handler.h"
#include "bins_handler.h"
#include "samplings_handler.h"
#include "mc_estimate.h"
#include "task_utils.h"
#include "momenta.h"
#include "minimizer.h"
#include "xml_handler.h"
#include "io_map.h"
#include "chisq_fit.h"
#include "chisq_base.h"
#include "chisq_3ptcorr.h"
#include "chisq_tcorr.h"
#include "plot_info.h"
#include "create_plots.h"

#include <vector>
#include <set>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

using namespace std;

// Note: Notice the 'xml' bindings make use of the ElementTree module in python.
//       I could not figure out a way to globally import this module for use
//       in each 'xml' binding (without resulting in seg faults). 
//       So, the import occurs every time the 'xml' bindings are called.
//       This doesn't seem ideal. Any solutions?

set<uint> getCorrelatorAvailableTimeSeps(MCObsHandler *moh, const CorrelatorInfo& corr, bool hermitian, ComplexArg arg)
{
  set<uint> tseps;
  getCorrelatorAvailableTimes(moh, tseps, corr, hermitian, arg);
  return tseps;
}

map<uint,set<uint> > getCorrelatorAvailableTimes(MCObsHandler *moh, const CorrelatorInfo& corr, bool hermitian, ComplexArg arg)
{
  map<uint,set<uint> > tavail;
  getCorrelatorAvailableTimes(moh, tavail, corr, hermitian, arg);
  return tavail;
}

map<uint,set<uint> > getCorrelatorRatioAvailableTimes(MCObsHandler *moh, const CorrelatorInfo& corr_3pt,
                            const OperatorInfo& two_pt_src_op, const OperatorInfo& two_pt_snk_op,
                            bool hermitian, ComplexArg arg)
{
  map<uint,set<uint> > tavail;
  getCorrelatorRatioAvailableTimes(moh, tavail, corr_3pt, two_pt_src_op, two_pt_snk_op, hermitian, arg);
  return tavail;
}

map<double,MCEstimate> getEffectiveEnergy(MCObsHandler *moh, const CorrelatorInfo& corr, 
                  bool hermitian, bool subtract_vev, ComplexArg arg, SamplingMode mode, uint step, 
                  uint efftype, double subtract_const)
{
  map<double,MCEstimate> results;
  getEffectiveEnergy(moh, corr, hermitian, subtract_vev, arg, mode, step, efftype, results, subtract_const);
  return results;
}

map<double,MCEstimate> getCorrelatorEstimates(MCObsHandler *moh, const CorrelatorInfo& corr,
                  bool hermitian, bool subtract_vev, ComplexArg arg, SamplingMode mode)
{
  map<double,MCEstimate> results;
  getCorrelatorEstimates(moh, corr, hermitian, subtract_vev, arg, mode, results);
  return results;
}

typedef pair<map<pair<double,double>,MCEstimate>,map<pair<double,double>,MCEstimate> > ratio_vals;
ratio_vals getCorrelatorRatioEstimates(MCObsHandler *moh, const CorrelatorInfo& corr_3pt, 
                            bool hermitian, bool subtract_vev, SamplingMode mode)
{
  map<pair<double,double>,MCEstimate> re_results;
  map<pair<double,double>,MCEstimate> im_results;
  getCorrelatorRatioEstimates(moh, corr_3pt, hermitian, subtract_vev, mode, re_results, im_results);
  return make_pair(re_results, im_results);
}

enum FileType {
  Correlator,
  VEV,
  Bins,
  Samplings,
  SinglePivot_CN,
  SinglePivot_RN,
};


FileType getFileID(const string& filename)
{
  ifstream fin(filename.c_str(), ios::binary);
  if (!fin)
    throw(invalid_argument("Could not find file '" + filename + "'"));

  char idstring[33];
  if (!fin.read(idstring,33))
    throw(invalid_argument("Invalid file '" + filename + "'"));

  string ID(&idstring[1],32);
  ID = tidyString(ID);

  if (ID=="Laph--CorrelatorFile")
    return Correlator;
  else if (ID=="Laph--VEVFile")
    return VEV;
  else if (ID=="Sigmond--BinsFile")
    return Bins;
  else if (ID=="Sigmond--SamplingsFile")
    return Samplings;
  else if (ID=="Sigmond--SinglePivotFile-CN")
    return SinglePivot_CN;
  else if (ID=="Sigmond--SinglePivotFile-RN")
    return SinglePivot_RN;
  else
    throw(invalid_argument("Invalid file '" + filename + "'"));
}

set<OperatorInfo> getOperatorBasis(const string& pivot_filename)
{
  XMLHandler xmlp;
  FileType file_id = getFileID(pivot_filename);
  if (file_id == SinglePivot_CN) {
    IOMap<UIntKey,Array<complex<double> > > iom;
    iom.openReadOnly(pivot_filename,"Sigmond--SinglePivotFile-CN");
    xmlp.set_from_string(iom.getHeader());}
  else if (file_id == SinglePivot_RN) {
    IOMap<UIntKey,Array<double> > iom;
    iom.openReadOnly(pivot_filename,"Sigmond--SinglePivotFile-RN");
    xmlp.set_from_string(iom.getHeader());}
  else
    throw(invalid_argument("File not a pivot file '" + pivot_filename + "'"));

  
  list<string> tagnames;
  tagnames.push_back("Operator");
  tagnames.push_back("OperatorString");
  tagnames.push_back("BLOperator");
  tagnames.push_back("BLOperatorString");
  tagnames.push_back("GIOperator");
  tagnames.push_back("GIOperatorString");
  xmlp.seek_child("MatrixDefinition");
  xmlp.seek_child("CorrelatorMatrixInfo");
  list<XMLHandler> opxml=xmlp.find_among_children(tagnames);
  
  set<OperatorInfo> opset;
  for (list<XMLHandler>::iterator ot=opxml.begin();ot!=opxml.end();++ot)
    opset.insert(OperatorInfo(*ot));


  return opset;
}

using RtcVector = vector<RealTemporalCorrelatorFit, allocator<RealTemporalCorrelatorFit> >;

PYBIND11_MODULE(sigmond, m) {

  m.doc() = "pybind11 wrapper for sigmond";

  // Functions
  m.def("createEffEnergyWithFitPlot", (void (*)(EffEnergyWithFitPlotInfo, RealTemporalCorrelatorFit&, FitResult&, MCObsHandler*, XMLHandler&)) &createEffEnergyWithFitPlot);
  m.def("createThreePointCorrelatorPlot", (void (*)(ThreePointCorrelatorPlotInfo, const CorrelatorInfo&, const OperatorInfo&, const OperatorInfo&, bool, MCObsHandler*, XMLHandler&)) &createThreePointCorrelatorPlot);
  m.def("createThreePointCorrelatorWithFitPlot", (void (*)(ThreePointCorrelatorWithFitPlotInfo, RealThreePointCorrelatorFit&, FitResult&, MCObsHandler*, XMLHandler&)) &createThreePointCorrelatorWithFitPlot);
  m.def("createDataFitRatioPlot", (void (*)(DataFitRatioPlotInfo, vector<RealTemporalCorrelatorFit>&, MCObsHandler*, XMLHandler&)) &createDataFitRatioPlot);
  m.def("createTminPlot", (void (*)(TminFitPlotInfo, vector<vector<RealTemporalCorrelatorFit> >&, vector<vector<FitResult> >&, MCObsHandler*, XMLHandler&)) &createTminPlot);
  m.def("createDispersionFitPlot", (void (*)(DispersionFitPlotInfo, DispersionFit&, FitResult&, MCObsHandler*, XMLHandler&, bool)) &createDispersionFitPlot);

  m.def("doRealTemporalCorrelatorFit", (FitResult (*)(RealTemporalCorrelatorFit&, ChiSquareMinimizerInfo&, vector<double>, bool, XMLHandler&)) &doChiSquareFitting);
  m.def("doRealThreePointCorrelatorFit", (FitResult (*)(RealThreePointCorrelatorFit&, ChiSquareMinimizerInfo&, vector<double>, bool, XMLHandler&)) &doChiSquareFitting);
  m.def("doDispersionFit", (FitResult (*)(DispersionFit&, ChiSquareMinimizerInfo&, vector<double>, bool, XMLHandler&)) &doChiSquareFitting);

  m.def("getCorrelatorAvailableTimeSeps", (set<uint> (*)(MCObsHandler*, const CorrelatorInfo&, bool, ComplexArg)) &getCorrelatorAvailableTimeSeps);
  m.def("getCorrelatorAvailableTimes", (map<uint,set<uint> > (*)(MCObsHandler*, const CorrelatorInfo&, bool, ComplexArg)) &getCorrelatorAvailableTimes);
  m.def("getCorrelatorRatioAvailableTimes", (map<uint,set<uint> > (*)(MCObsHandler*, const CorrelatorInfo&, const OperatorInfo&, const OperatorInfo&, bool, ComplexArg)) &getCorrelatorRatioAvailableTimes);
  m.def("getCorrelatorEstimates", (map<double,MCEstimate> (*)(MCObsHandler*, const CorrelatorInfo&, bool, bool, ComplexArg, SamplingMode)) &getCorrelatorEstimates);
  m.def("getEffectiveEnergy", (map<double,MCEstimate> (*)(MCObsHandler*, const CorrelatorInfo&, bool, bool, ComplexArg, SamplingMode, uint, uint, double)) &getEffectiveEnergy);
  m.def("getCorrelatorRatioEstimates", (ratio_vals (*)(MCObsHandler*, const CorrelatorInfo&, bool, bool, SamplingMode)) &getCorrelatorRatioEstimates);

  m.def("doRatioBySamplings", (void (*)(MCObsHandler&, const MCObsInfo&, const MCObsInfo&, const MCObsInfo&)) &doRatioBySamplings);
  m.def("doLogBySamplings", (void (*)(MCObsHandler&, const MCObsInfo&, const MCObsInfo&)) &doLogBySamplings);
  m.def("doBoostBySamplings", (void (*)(MCObsHandler&, const MCObsInfo&, double, const MCObsInfo&)) &doBoostBySamplings);
  m.def("doLinearSuperpositionBySamplings", (void (*)(MCObsHandler&, vector<MCObsInfo>&, vector<double>&, const MCObsInfo&)) &doLinearSuperpositionBySamplings);
  m.def("doSquareRootBySamplings", (void (*)(MCObsHandler&, const MCObsInfo&, const MCObsInfo& obs_out)) &doSquareRootBySamplings);

  m.def("getFileID", (FileType (*)(const string&)) &getFileID);
  m.def("getOperatorBasis", (set<OperatorInfo> (*)(const string&)) &getOperatorBasis);

  py::enum_<FileType>(m, "FileType")
    .value("Correlator", FileType::Correlator)
    .value("VEV", FileType::VEV)
    .value("Bins", FileType::Bins)
    .value("Samplings", FileType::Samplings)
    .value("SinglePivot_CN", FileType::SinglePivot_CN)
    .value("SinglePivot_RN", FileType::SinglePivot_RN);
  
  // Info classes
  py::class_<MCEnsembleInfo>(m, "MCEnsembleInfo")
    .def(py::init<const string &>())
    .def(py::init<const string &, const string &>())
    .def(py::init<const string &, uint, uint, uint, uint, uint, uint>())
    .def("getId", &MCEnsembleInfo::getId)
    .def("getLatticeTimeExtent", &MCEnsembleInfo::getLatticeTimeExtent)
    .def("getLatticeXExtent", &MCEnsembleInfo::getLatticeXExtent)
    .def("getLatticeYExtent", &MCEnsembleInfo::getLatticeYExtent)
    .def("getLatticeZExtent", &MCEnsembleInfo::getLatticeZExtent)
    .def("xml", [](const MCEnsembleInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.str()); })
    .def("__str__", &MCEnsembleInfo::getId)
    .def("__repr__", &MCEnsembleInfo::str)
    .def(py::self == py::self)
    .def(py::self != py::self);

  py::class_<MCBinsInfo>(m, "MCBinsInfo")
    .def(py::init<const MCEnsembleInfo &>())
    .def("setRebin", &MCBinsInfo::setRebin)
    .def("addOmission", &MCBinsInfo::addOmission)
    .def("addOmissions", &MCBinsInfo::addOmissions)
    .def("clearOmissions", &MCBinsInfo::clearOmissions)
    .def("getRebinFactor", &MCBinsInfo::getRebinFactor)
    .def("getOmissions", &MCBinsInfo::getOmissions)
    .def("getNumberOfBins", &MCBinsInfo::getNumberOfBins)
    .def("getMCEnsembleInfo", &MCBinsInfo::getMCEnsembleInfo)
    .def("getLatticeTimeExtent", &MCBinsInfo::getLatticeTimeExtent)
    .def("getLatticeXExtent", &MCBinsInfo::getLatticeXExtent)
    .def("getLatticeYExtent", &MCBinsInfo::getLatticeYExtent)
    .def("getLatticeZExtent", &MCBinsInfo::getLatticeZExtent)
    .def("xml", [](const MCBinsInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.str()); })
    .def("__str__", [](const MCBinsInfo &a) { return a.output(2); })
    .def("__repr__", &MCBinsInfo::str)
    .def(py::self == py::self)
    .def(py::self != py::self);

  // Note: Is there really no better way to initialize an enum with a string other
  //       than to use this 'create' static method?
  py::enum_<SamplingMode>(m, "SamplingMode")
    .value("Jackknife", SamplingMode::Jackknife)
    .value("Bootstrap", SamplingMode::Bootstrap)
    .def_static("create", [](string s) {
        if (s == "jackknife")
          return SamplingMode::Jackknife;
        else if (s == "bootstrap")
          return SamplingMode::Bootstrap;
        throw(invalid_argument("Bad SamplingMode '" + s + "'")); })
    .def("__str__", [](const SamplingMode &a) {
        switch(a) {
          case Jackknife : return "Jackknife";
          case Bootstrap : return "Bootstrap";
          default        : throw(invalid_argument("Bad SamplingMode"));
        }}, py::prepend());

  py::class_<Bootstrapper>(m, "Bootstrapper")
    .def(py::init<uint, uint,  unsigned long, uint, bool>());

  py::class_<MCSamplingInfo>(m, "MCSamplingInfo")
    .def(py::init<>())
    .def(py::init<uint, unsigned long, uint>())
    .def("isJackknifeMode", &MCSamplingInfo::isJackknifeMode)
    .def("isBootstrapMode", &MCSamplingInfo::isBootstrapMode)
    .def("setToJackknifeMode", &MCSamplingInfo::setToJackknifeMode)
    .def("setToBootstrapMode", &MCSamplingInfo::setToBootstrapMode)
    .def("getSamplingMode", &MCSamplingInfo::getSamplingMode)
    .def("getNumberOfReSamplings", &MCSamplingInfo::getNumberOfReSamplings)
    .def("getRNGSeed", &MCSamplingInfo::getRNGSeed)
    .def("getSkipValue", &MCSamplingInfo::getSkipValue)
    .def("xml", [](const MCSamplingInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.str()); })
    .def("__str__", [](const MCSamplingInfo &a) { return a.output(2); })
    .def("__repr__", &MCSamplingInfo::str)
    .def(py::self == py::self)
    .def(py::self != py::self);

  // Observables
  // Note: The way the hash functions are implemented aren't exactly ideal.
  //       For example, in the case of the 'BasicLapHOperatorInfo' it would be nice if
  //       one could use the 'icode' member for the '__repr__' method and hash that
  //       icode for the '__hash__' method. But, that isn't so straightforward.
  //       I think these hash implementations will suffice for now.
  py::enum_<ComplexArg>(m, "ComplexArg")
    .value("RealPart", ComplexArg::RealPart)
    .value("ImaginaryPart", ComplexArg::ImaginaryPart)
    .def("__str__", [](const ComplexArg &a) {
        switch(a) {
          case RealPart      : return "RealPart";
          case ImaginaryPart : return "ImaginaryPart";
          default            : throw(invalid_argument("Bad ComplexArg"));
        }}, py::prepend());

  py::enum_<OperatorInfo::OpKind>(m, "OpKind")
    .value("BasicLapH", OperatorInfo::OpKind::BasicLapH)
    .value("GenIrrep", OperatorInfo::OpKind::GenIrrep);

  py::class_<MCObsInfo>(m, "MCObsInfo")
    .def(py::init<const OperatorInfo &, ComplexArg>())
    .def(py::init<const OperatorInfo &, OperatorInfo &, int, bool, ComplexArg, bool>())
    .def(py::init<const OperatorInfo &, OperatorInfo &, OperatorInfo &, int, int, bool, ComplexArg, bool>())
    .def(py::init<const CorrelatorAtTimeInfo &, ComplexArg>())
    .def(py::init<const CorrelatorInfo &, int, bool, ComplexArg, bool>())
    .def(py::init<const string &, uint, bool, ComplexArg>())
    .def(py::init<const string &, uint>())
    .def("isVacuum", &MCObsInfo::isVacuum)
    .def("isVEV", &MCObsInfo::isVEV)
    .def("isCorrelatorAtTime", &MCObsInfo::isCorrelatorAtTime)
    .def("isHermitianCorrelatorAtTime", &MCObsInfo::isHermitianCorrelatorAtTime)
    .def("isRealPart", &MCObsInfo::isRealPart)
    .def("isImaginaryPart", &MCObsInfo::isImaginaryPart)
    .def("isSimple", &MCObsInfo::isSimple)
    .def("isNonSimple", &MCObsInfo::isNonSimple)
    .def("isPrimary", &MCObsInfo::isPrimary)
    .def("isSecondary", &MCObsInfo::isSecondary)
    .def("isBasicLapH", &MCObsInfo::isBasicLapH)
    .def("isGenIrrep", &MCObsInfo::isGenIrrep)
    .def("getVEVInfo", (OperatorInfo (MCObsInfo::*)() const) &MCObsInfo::getVEVInfo)
    .def("getCorrelatorAtTimeInfo", (CorrelatorAtTimeInfo (MCObsInfo::*)() const) &MCObsInfo::getCorrelatorAtTimeInfo)
    .def("getCorrelatorSourceInfo", &MCObsInfo::getCorrelatorSourceInfo)
    .def("getCorrelatorSinkInfo", &MCObsInfo::getCorrelatorSinkInfo)
    .def("getCorrelatorTimeIndex", &MCObsInfo::getCorrelatorTimeIndex)
    .def("getCorrelatorInfo", (CorrelatorInfo (MCObsInfo::*)() const) &MCObsInfo::getCorrelatorInfo)
    .def("getObsName", &MCObsInfo::getObsName)
    .def("getObsIndex", &MCObsInfo::getObsIndex)
    .def("xml", [](const MCObsInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("__str__", [](const MCObsInfo &a) { return a.output(false, 2); })
    .def("__repr__", &MCObsInfo::str)
    .def("__hash__", [](const MCObsInfo &a) { return hash<string>{}(a.str()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<OperatorInfo>(m, "OperatorInfo")
    .def(py::init<const string &, OperatorInfo::OpKind>())
    .def(py::init<const BasicLapHOperatorInfo &>())
    .def(py::init<const GenIrrepOperatorInfo &>())
    .def("isBasicLapH", &OperatorInfo::isBasicLapH)
    .def("isGenIrrep", &OperatorInfo::isGenIrrep)
    .def("getBasicLapH", &OperatorInfo::getBasicLapH)
    .def("getGenIrrep", &OperatorInfo::getGenIrrep)
    .def("isBackwards", &OperatorInfo::isBackwards)
    .def("setBackwards", &OperatorInfo::setBackwards)
    .def("setForwards", &OperatorInfo::setForwards)
    .def("long_xml", [](const OperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("xml", [](const OperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(false)); })
    .def("op_str", &OperatorInfo::short_output)
    .def("__str__", &OperatorInfo::short_output)
    .def("__repr__", &OperatorInfo::short_output)
    .def("__hash__", [](const OperatorInfo &a) { return hash<string>{}(a.short_output()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<BasicLapHOperatorInfo>(m, "BasicLapHOperatorInfo")
    .def(py::init<const string &>())
    .def("getNumberOfHadrons", &BasicLapHOperatorInfo::getNumberOfHadrons)
    .def("isGlueball", (bool (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::isGlueball)
    .def("isMeson", (bool (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::isMeson)
    .def("isBaryon", (bool (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::isBaryon)
    .def("isTetraquark", (bool (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::isTetraquark)
    .def("isMesonMeson", &BasicLapHOperatorInfo::isMesonMeson)
    .def("isMesonBaryon", &BasicLapHOperatorInfo::isMesonBaryon)
    .def("isBackwards", &BasicLapHOperatorInfo::isBackwards)
    .def("setBackwards", &BasicLapHOperatorInfo::setBackwards)
    .def("setForwards", &BasicLapHOperatorInfo::setForwards)
    .def("getMomentum", (Momentum (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getMomentum)
    .def("getXMomentum", (int (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getXMomentum)
    .def("getYMomentum", (int (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getYMomentum)
    .def("getZMomentum", (int (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getZMomentum)
    .def("getLGIrrep", (string (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getLGIrrep)
    .def("getLGClebschGordonIdNum", &BasicLapHOperatorInfo::getLGClebschGordonIdNum)
    .def("getLGIrrepRow", &BasicLapHOperatorInfo::getLGIrrepRow)
    .def("getIsospin", &BasicLapHOperatorInfo::getIsospin)
    .def("getIsospinClebschGordonIdNum", &BasicLapHOperatorInfo::getIsospinClebschGordonIdNum)
    .def("isBackwards", &BasicLapHOperatorInfo::isBackwards)
    .def("getFlavor", (string (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getFlavor)
    .def("getFlavorCode", &BasicLapHOperatorInfo::getFlavorCode)
    .def("getStrangeness", (int (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getStrangeness)
    .def("getTetraquarkColorType", (int (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getTetraquarkColorType)
    .def("getHadronFlavor", (string (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::getFlavor)
    .def("getHadronStrangeness", (int (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::getStrangeness)
    .def("isHadronGlueball", (bool (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::isGlueball)
    .def("isHadronMeson", (bool (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::isMeson)
    .def("isHadronBaryon", (bool (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::isBaryon)
    .def("isHadronFermion", &BasicLapHOperatorInfo::isFermion)
    .def("isHadronBoson", &BasicLapHOperatorInfo::isBoson)
    .def("getHadronLGIrrep", (string (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::getLGIrrep)
    .def("getHadronSpatialType", &BasicLapHOperatorInfo::getSpatialType)
    .def("getHadronSpatialIdNumber", &BasicLapHOperatorInfo::getSpatialIdNumber)
    .def("getHadronDisplacementLength", &BasicLapHOperatorInfo::getDisplacementLength)
    .def("getHadronMomentum", (Momentum (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::getMomentum)
    .def("getHadronXMomentum", (int (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::getXMomentum)
    .def("getHadronYMomentum", (int (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::getYMomentum)
    .def("getHadronZMomentum", (int (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::getZMomentum)
    .def("long_xml", [](const BasicLapHOperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("xml", [](const BasicLapHOperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(false)); })
    .def("op_str", &BasicLapHOperatorInfo::short_output)
    .def("__str__", &BasicLapHOperatorInfo::short_output)
    .def("__repr__", &BasicLapHOperatorInfo::short_output)
    .def("__hash__", [](const BasicLapHOperatorInfo &a) { return hash<string>{}(a.short_output()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<GenIrrepOperatorInfo>(m, "GenIrrepOperatorInfo")
    .def(py::init<const string &>())
    .def("getMomentum", &GenIrrepOperatorInfo::getMomentum)
    .def("getXMomentum", &GenIrrepOperatorInfo::getXMomentum)
    .def("getYMomentum", &GenIrrepOperatorInfo::getYMomentum)
    .def("getZMomentum", &GenIrrepOperatorInfo::getZMomentum)
    .def("isReferenceMomentum", &GenIrrepOperatorInfo::isReferenceMomentum)
    .def("getLGIrrep", &GenIrrepOperatorInfo::getLGIrrep)
    .def("getLGIrrepRow", &GenIrrepOperatorInfo::getLGIrrepRow)
    .def("isSU3flavor", &GenIrrepOperatorInfo::isSU3flavor)
    .def("getFlavor", &GenIrrepOperatorInfo::getFlavor)
    .def("getIDName", &GenIrrepOperatorInfo::getIDName)
    .def("getIDIndex", &GenIrrepOperatorInfo::getIDIndex)
    .def("resetIDIndex", (GenIrrepOperatorInfo& (GenIrrepOperatorInfo::*)(uint)) &GenIrrepOperatorInfo::resetIDIndex)
    .def("long_xml", [](const GenIrrepOperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("xml", [](const GenIrrepOperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(false)); })
    .def("op_str", &GenIrrepOperatorInfo::short_output)
    .def("__str__", &GenIrrepOperatorInfo::short_output)
    .def("__repr__", &GenIrrepOperatorInfo::short_output)
    .def("__hash__", [](const GenIrrepOperatorInfo &a) { return hash<string>{}(a.short_output()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<CorrelatorInfo>(m, "CorrelatorInfo")
    .def(py::init<const CorrelatorInfo &>())
    .def(py::init<const OperatorInfo &, const OperatorInfo &>())
    .def(py::init<const OperatorInfo &, const OperatorInfo &, const OperatorInfo &>())
    .def("getSource", &CorrelatorInfo::getSource)
    .def("getSink", &CorrelatorInfo::getSink)
    .def("isSinkSourceSame", &CorrelatorInfo::isSinkSourceSame)
    .def("isBackwards", &CorrelatorInfo::isBackwards)
    .def("setBackwards", &CorrelatorInfo::setBackwards)
    .def("setForwards", &CorrelatorInfo::setForwards)
    .def("xml", [](const CorrelatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("corr_str", [](const CorrelatorInfo &a) {
        return ("snk: " + a.getSink().short_output() + " src: " + a.getSource().short_output()); })
    .def("__str__", [](const CorrelatorInfo &a) { return a.output(false, 2); })
    .def("__repr__", [](const CorrelatorInfo &a) {
        return ("snk_" + a.getSink().short_output() + "-src_" + a.getSource().short_output()); })
    .def("__hash__", [](const CorrelatorInfo &a) { return hash<string>{}(a.str()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<CorrelatorAtTimeInfo>(m, "CorrelatorAtTimeInfo")
    .def(py::init<const OperatorInfo &, const OperatorInfo &, int, bool, bool>())
    .def(py::init<const OperatorInfo &, const OperatorInfo &, const OperatorInfo &, int, int, bool, bool>())
    .def(py::init<const CorrelatorInfo &, int, bool, bool>())
    .def(py::init<const CorrelatorInfo &, int, int, bool, bool>())
    .def ("resetTimeSeparation", &CorrelatorAtTimeInfo::resetTimeSeparation)
    .def ("resetTimeInsertion", &CorrelatorAtTimeInfo::resetTimeInsertion)
    .def("isBackwards", &CorrelatorAtTimeInfo::isBackwards)
    .def("setBackwards", &CorrelatorAtTimeInfo::setBackwards)
    .def("setForwards", &CorrelatorAtTimeInfo::setForwards)
    .def("xml", [](const CorrelatorAtTimeInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("__str__", [](const CorrelatorAtTimeInfo &a) { return a.output(false, 2); })
    .def("__repr__", &CorrelatorAtTimeInfo::str)
    .def("__hash__", [](const CorrelatorAtTimeInfo &a) { return hash<string>{}(a.str()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<CorrelatorMatrixInfo>(m, "CorrelatorMatrixInfo")
    .def(py::init<const set<OperatorInfo> &, bool, bool>())
    .def("long_xml", [](const CorrelatorMatrixInfo &a) {
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("xml", [](const CorrelatorMatrixInfo &a) {
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(false)); })
    .def("__str__", [](const CorrelatorMatrixInfo &a) { return a.output(false, 2); })
    .def("__repr__", &CorrelatorMatrixInfo::str)
    .def("__hash__", [](const CorrelatorMatrixInfo &a) { return hash<string>{}(a.str()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<Momentum>(m, "Momentum")
    .def(py::init<int, int, int>())
    .def_readwrite("x", &Momentum::x)
    .def_readwrite("y", &Momentum::y)
    .def_readwrite("z", &Momentum::z)
    .def("getPsq", (double (Momentum::*)(int, int, int) const) &Momentum::getPsq)
    .def("__str__", &Momentum::getMomentumString)
    .def("__repr__", &Momentum::getMomentumString)
    .def("__hash__", [](const Momentum &a) { return hash<string>{}(a.getMomentumString()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<ChiSquareMinimizerInfo>(m, "ChiSquareMinimizerInfo")
    .def(py::init<>())
    .def(py::init<char, double, double, int, char>())
    .def("xml", [](const ChiSquareMinimizerInfo &a) {
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output()); });

  py::class_<FitResult>(m, "FitResult")
    .def_readwrite("chisq_dof", &FitResult::chisq_dof)
    .def_readwrite("quality", &FitResult::quality);

  py::class_<RealTemporalCorrelatorFit>(m, "RealTemporalCorrelatorFit")
    .def(py::init<MCObsHandler&, OperatorInfo, bool, string, map<string,MCObsInfo>, uint, uint, double>())
    .def("getTmin", &RealTemporalCorrelatorFit::getTmin)
    .def("getTmax", &RealTemporalCorrelatorFit::getTmax)
    .def("getEffMassType", &RealTemporalCorrelatorFit::getEffMassType)
    .def("getFitParamInfos", &RealTemporalCorrelatorFit::getFitParamInfos)
    .def("setEnergyInfo", &RealTemporalCorrelatorFit::setEnergyInfo)
    .def("getObsMean", &RealTemporalCorrelatorFit::getObsMean)
    .def("getCovarianceMatrix", &RealTemporalCorrelatorFit::getCovarianceMatrix)
    .def("evalModelPoint", &RealTemporalCorrelatorFit::evalModelPoint)
    .def("addPriors", &RealTemporalCorrelatorFit::addPriors);

  py::class_<RealThreePointCorrelatorFit>(m, "RealThreePointCorrelatorFit")
    .def(py::init<MCObsHandler&, CorrelatorInfo, const OperatorInfo&, const OperatorInfo&, CorrelatorInfo, bool, ComplexArg, string, map<string,MCObsInfo>, map<uint,set<uint> > >())
    .def(py::init<MCObsHandler&, CorrelatorInfo, CorrelatorInfo, bool, ComplexArg, string, map<string,MCObsInfo>, map<uint,set<uint> > >())
    .def("getFitParamInfos", &RealThreePointCorrelatorFit::getFitParamInfos)
    .def("getCovarianceMatrix", &RealThreePointCorrelatorFit::getCovarianceMatrix)
    .def("evalModelPoint", &RealThreePointCorrelatorFit::evalModelPoint)
    .def("addPriors", &RealThreePointCorrelatorFit::addPriors);

  py::class_<DispersionFit>(m, "DispersionFit")
    .def(py::init<MCObsHandler&, const string&, map<MCObsInfo,double>, map<string,MCObsInfo> >())
    .def("getFitParamInfos", &DispersionFit::getFitParamInfos)
    .def("evalModelPoint", &DispersionFit::evalModelPoint)
    .def("addPriors", &DispersionFit::addPriors);

  py::class_<Prior>(m, "Prior")
    .def(py::init<double, double>())
    .def(py::init<MCObsInfo&, MCObsHandler&>());

  py::class_<EffEnergyWithFitPlotInfo>(m, "EffEnergyWithFitPlotInfo")
    .def(py::init<const string&>())
    .def_readwrite("plotfile", &EffEnergyWithFitPlotInfo::plotfile)
    .def_readwrite("corrname", &EffEnergyWithFitPlotInfo::corrname)
    .def_readwrite("timestep", &EffEnergyWithFitPlotInfo::timestep)
    .def_readwrite("symbolcolor", &EffEnergyWithFitPlotInfo::symbolcolor)
    .def_readwrite("symboltype", &EffEnergyWithFitPlotInfo::symboltype)
    .def_readwrite("maxerror", &EffEnergyWithFitPlotInfo::maxerror)
    .def_readwrite("goodness", &EffEnergyWithFitPlotInfo::goodness)
    .def_readwrite("ref_energy", &EffEnergyWithFitPlotInfo::ref_energy);

  py::class_<ThreePointCorrelatorPlotInfo>(m, "ThreePointCorrelatorPlotInfo")
    .def(py::init<const string&, const vector<uint>&>())
    .def_readwrite("plotfile", &ThreePointCorrelatorPlotInfo::plotfile)
    .def_readwrite("plotlabel", &ThreePointCorrelatorPlotInfo::plotlabel)
    .def_readwrite("time_seps", &ThreePointCorrelatorPlotInfo::time_seps)
    .def_readwrite("complex_arg", &ThreePointCorrelatorPlotInfo::complex_arg)
    .def_readwrite("labels", &ThreePointCorrelatorPlotInfo::labels)
    .def_readwrite("symbolcolors", &ThreePointCorrelatorPlotInfo::symbolcolors)
    .def_readwrite("symboltypes", &ThreePointCorrelatorPlotInfo::symboltypes);

  py::class_<ThreePointCorrelatorWithFitPlotInfo>(m, "ThreePointCorrelatorWithFitPlotInfo")
    .def(py::init<const string&, const vector<uint>&>())
    .def_readwrite("plotfile", &ThreePointCorrelatorWithFitPlotInfo::plotfile)
    .def_readwrite("plotlabel", &ThreePointCorrelatorWithFitPlotInfo::plotlabel)
    .def_readwrite("time_seps", &ThreePointCorrelatorWithFitPlotInfo::time_seps)
    .def_readwrite("goodness", &ThreePointCorrelatorWithFitPlotInfo::goodness)
    .def_readwrite("labels", &ThreePointCorrelatorWithFitPlotInfo::labels)
    .def_readwrite("symbolcolors", &ThreePointCorrelatorWithFitPlotInfo::symbolcolors)
    .def_readwrite("symboltypes", &ThreePointCorrelatorWithFitPlotInfo::symboltypes);

  py::class_<DataFitRatioPlotInfo>(m, "DataFitRatioPlotInfo")
    .def(py::init<const string&>())
    .def_readwrite("plotfile", &DataFitRatioPlotInfo::plotfile)
    .def_readwrite("plotlabel", &DataFitRatioPlotInfo::plotlabel)
    .def_readwrite("labels", &DataFitRatioPlotInfo::labels)
    .def_readwrite("symbolcolors", &DataFitRatioPlotInfo::symbolcolors)
    .def_readwrite("symboltypes", &DataFitRatioPlotInfo::symboltypes)
    .def_readwrite("maxerror", &DataFitRatioPlotInfo::maxerror);

  py::class_<TminFitPlotInfo>(m, "TminFitPlotInfo")
    .def(py::init<const string&>())
    .def_readwrite("plotfile", &TminFitPlotInfo::plotfile)
    .def_readwrite("energy_level", &TminFitPlotInfo::energy_level)
    .def_readwrite("plotlabel", &TminFitPlotInfo::plotlabel)
    .def_readwrite("labels", &TminFitPlotInfo::labels)
    .def_readwrite("symbolcolors", &TminFitPlotInfo::symbolcolors)
    .def_readwrite("symboltypes", &TminFitPlotInfo::symboltypes)
    .def_readwrite("quality_threshold", &TminFitPlotInfo::quality_threshold)
    .def_readwrite("maxerror", &TminFitPlotInfo::maxerror)
    .def_readwrite("chosen_fit", &TminFitPlotInfo::chosen_fit);

  py::class_<DispersionFitPlotInfo>(m, "DispersionFitPlotInfo")
    .def(py::init<const string&>())
    .def_readwrite("plotfile", &DispersionFitPlotInfo::plotfile)
    .def_readwrite("particle_name", &DispersionFitPlotInfo::particle_name)
    .def_readwrite("symbolcolor", &DispersionFitPlotInfo::symbolcolor)
    .def_readwrite("symboltype", &DispersionFitPlotInfo::symboltype)
    .def_readwrite("goodness", &DispersionFitPlotInfo::goodness);

  // Data Handlers
  py::class_<RVector>(m, "RVector")
    .def(py::init<>())
    .def(py::init<const vector<double> &>())
    .def("array", &RVector::c_vector);

  py::class_<MCEstimate>(m, "MCEstimate")
    .def(py::init<>())
    .def("getFullEstimate", &MCEstimate::getFullEstimate)
    .def("getAverageEstimate", &MCEstimate::getAverageEstimate)
    .def("getSymmetricError", &MCEstimate::getSymmetricError)
    .def("getRelativeError", &MCEstimate::getRelativeError)
    .def("isStatisticallyZero", [](const MCEstimate &a) {
        return (abs(a.getFullEstimate()) < a.getSymmetricError()); });

  py::class_<XMLHandler>(m, "XMLHandler")
    .def(py::init<>())
    .def(py::init<const string &>())
    .def(py::init<const string &, const string &>())
    .def("output", &XMLHandler::output)
    .def("set_from_string", &XMLHandler::set_from_string);

  py::class_<FileListInfo>(m, "FileListInfo")
    .def(py::init<const string &, int, int, bool>())
    .def("getFileStub", &FileListInfo::getFileStub)
    .def("getMaxFileNumber", &FileListInfo::getMaxFileNumber)
    .def("getMinFileNumber", &FileListInfo::getMinFileNumber)
    .def("isModeOverwrite", &FileListInfo::isModeOverwrite)
    .def("xml", [](const FileListInfo &a) {
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.str()); })
    .def("__str__", [](const FileListInfo &a) {
        return (a.getFileStub() + ".[" + to_string(a.getMinFileNumber()) + "," + to_string(a.getMaxFileNumber()) + "]"); })
    .def("__repr__", &FileListInfo::str)
    .def("__hash__", [](const FileListInfo &a) { return hash<string>{}(a.str()); })
    .def(py::self == py::self);

  py::class_<MCObsGetHandler>(m, "MCObsGetHandler")
    .def(py::init<XMLHandler&, const MCBinsInfo&, const MCSamplingInfo&>())
    .def(py::init<const MCBinsInfo&, const MCSamplingInfo&>())
    .def(py::init<const MCBinsInfo&, const MCSamplingInfo&, bool>())
    .def(py::init<const MCBinsInfo&, const MCSamplingInfo&,
                  list<FileListInfo>&, list<FileListInfo>&,
                  set<string>&, set<string>&,
                  bool>());

  py::class_<MCObsHandler>(m, "MCObsHandler")
    .def(py::init<MCObsGetHandler &, bool>())
    .def("getBins", (const RVector & (MCObsHandler::*)(const MCObsInfo &) ) &MCObsHandler::getBins)
    .def("getBin", (double (MCObsHandler::*)(const MCObsInfo &, int) ) &MCObsHandler::getBin)
    .def("getEstimate", (MCEstimate (MCObsHandler::*)(const MCObsInfo &) ) &MCObsHandler::getEstimate)
    .def("setSamplingBegin", &MCObsHandler::setSamplingBegin)
    .def("isSamplingEnd", &MCObsHandler::isSamplingEnd)
    .def("setSamplingNext", &MCObsHandler::setSamplingNext)
    .def("queryFullAndSamplings", (bool (MCObsHandler::*)(const MCObsInfo&) ) &MCObsHandler::queryFullAndSamplings)
    .def("getFullAndSamplingValues", (const RVector& (MCObsHandler::*)(const MCObsInfo&, SamplingMode) ) &MCObsHandler::getFullAndSamplingValues)
    .def("queryBins", (bool (MCObsHandler::*)(const MCObsInfo&) ) &MCObsHandler::queryBins)
    .def("getCurrentSamplingValue", &MCObsHandler::getCurrentSamplingValue)
    .def("putCurrentSamplingValue", &MCObsHandler::putCurrentSamplingValue)
    .def("writeSamplingValuesToFile", &MCObsHandler::writeSamplingValuesToFile)
    .def("getLatticeTimeExtent", &MCObsHandler::getLatticeTimeExtent)
    .def("getLatticeXExtent", &MCObsHandler::getLatticeXExtent)
    .def("getLatticeYExtent", &MCObsHandler::getLatticeYExtent)
    .def("getLatticeZExtent", &MCObsHandler::getLatticeZExtent)
    .def("putBins", &MCObsHandler::putBins)
    .def("writeBinsToFile", &MCObsHandler::writeBinsToFile)
    .def("clearData", &MCObsHandler::clearData)
    .def("clearSamplings", &MCObsHandler::clearSamplings)
    .def("eraseData", &MCObsHandler::eraseData)
    .def("eraseSamplings", &MCObsHandler::eraseSamplings);

  py::enum_<WriteMode>(m, "WriteMode")
    .value("Protect", WriteMode::Protect)
    .value("Update", WriteMode::Update)
    .value("Overwrite", WriteMode::Overwrite)
    .def_static("create", [](string s) {
        if (s == "protect")
          return WriteMode::Protect;
        else if (s == "update")
          return WriteMode::Update;
        else if (s == "overwrite")
          return WriteMode::Overwrite;
        throw(invalid_argument("Bad WriteMode")); })
    .def("__str__", [](const WriteMode &a) {
        switch(a) {
          case Protect   : return "protect";
          case Update    : return "update";
          case Overwrite : return "overwrite";
          default        : throw(invalid_argument("Bad WriteMode"));
        }}, py::prepend());

  py::class_<LaphEnv::BLCorrelatorDataHandler>(m, "BLCorrelatorDataHandler")
    .def(py::init<const list<FileListInfo> &, const set<CorrelatorInfo> &,
                  const set<CorrelatorInfo> &, const MCEnsembleInfo *, bool>())
    .def(py::init<const list<FileListInfo> &, const set<CorrelatorInfo> &,
                  const set<CorrelatorInfo> &, const MCEnsembleInfo *>())
    .def("getCorrelatorSet", &LaphEnv::BLCorrelatorDataHandler::getCorrelatorSet)
    .def("getFileName", &LaphEnv::BLCorrelatorDataHandler::getFileName)
    .def("getKeys", &LaphEnv::BLCorrelatorDataHandler::getKeys)
    .def("getOrderedKeys", &LaphEnv::BLCorrelatorDataHandler::getOrderedKeys)
    .def("getTimeSepRange", &LaphEnv::BLCorrelatorDataHandler::getTimeSepRange);

  py::class_<LaphEnv::BLCorrelatorDataHandler::RecordKey>(m, "BLCorrelatorRecordKey")
    .def(py::init<int, int>())
    .def("getTimeIndex", &LaphEnv::BLCorrelatorDataHandler::RecordKey::getTimeIndex)
    .def("getConfigSerialIndex", &LaphEnv::BLCorrelatorDataHandler::RecordKey::getConfigSerialIndex);

  py::class_<LaphEnv::BLVEVDataHandler>(m, "BLVEVDataHandler")
    .def(py::init<const list<FileListInfo> &, const set<OperatorInfo> &,
                  const MCEnsembleInfo *, bool>())
    .def(py::init<const list<FileListInfo> &, const set<OperatorInfo> &,
                  const MCEnsembleInfo *>())
    .def("getOperatorSet", &LaphEnv::BLVEVDataHandler::getOperatorSet)
    .def("getFileName", &LaphEnv::BLVEVDataHandler::getFileName)
    .def("getKeys", &LaphEnv::BLVEVDataHandler::getKeys);

  py::class_<BinsGetHandler>(m, "BinsGetHandler")
    .def(py::init<const MCBinsInfo &, const set<string> &, bool>())
    .def(py::init<const MCBinsInfo &, const set<string> &>())
    .def("getFileNames", &BinsGetHandler::getFileNames)
    .def("getKeys", &BinsGetHandler::getKeys);

  py::class_<BinsPutHandler>(m, "BinsPutHandler")
    .def(py::init<const MCBinsInfo &, const string &, WriteMode, bool>())
    .def(py::init<const MCBinsInfo &, const string &, WriteMode>())
    .def(py::init<const MCBinsInfo &, const string &>())
    .def("putData", &BinsPutHandler::putData)
    .def("close", &BinsPutHandler::close);

  py::class_<SamplingsGetHandler>(m, "SamplingsGetHandler")
    .def(py::init<const MCBinsInfo &, const MCSamplingInfo &,
                  const set<string> &, bool>())
    .def(py::init<const MCBinsInfo &, const MCSamplingInfo &,
                  const set<string> &>())
    .def("getFileNames", &SamplingsGetHandler::getFileNames)
    .def("getKeys", &SamplingsGetHandler::getKeys);
}

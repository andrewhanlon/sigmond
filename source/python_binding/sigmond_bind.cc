#include <scalar_defs.h>
#include <ensemble_info.h>
#include <bins_info.h>
#include <bootstrapper.h>
#include <sampling_info.h>
#include <obs_get_handler.h>
#include <xml_handler.h>
#include <mcobs_info.h>
#include <mcobs_handler.h>
#include <operator_info.h>
#include <gen_irrep_operator_info.h>
#include <correlator_info.h>
#include <correlator_matrix_info.h>
#include <matrix.h>
#include <filelist_info.h>
#include <corr_data_handler.h>
#include <vev_data_handler.h>
#include <bins_handler.h>
#include <samplings_handler.h>
#include <mc_estimate.h>
#include <task_utils.h>
#include <momenta.h>
#include <minimizer.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

// Note: Notice the 'xml' bindings make use of the ElementTree module in python.
//       I could not figure out a way to globally import this module for use
//       in each 'xml' binding (without resulting in seg faults). 
//       So, the import occurs every time the 'xml' bindings are called.
//       This doesn't seem ideal. Any solutions?
  
PYBIND11_MODULE(sigmondbind, m) {

  m.doc() = "pybind11 wrapper for sigmond";

  // Functions
  m.def("getEffectiveEnergy", (std::map<int,MCEstimate> (*)(MCObsHandler*, const CorrelatorInfo&, bool, bool, bool, ComplexArg, SamplingMode, uint, uint, double)) &getEffectiveEnergy);
  m.def("getCorrelatorEstimates", (std::map<int,MCEstimate> (*)(MCObsHandler*, const CorrelatorInfo&, bool, bool, bool, ComplexArg, SamplingMode)) &getCorrelatorEstimates);

  // Info classes
  py::class_<MCEnsembleInfo>(m, "MCEnsembleInfo")
    .def(py::init<const std::string &>())
    .def(py::init<const std::string &, uint, uint, uint, uint, uint, uint>())
    .def("getId", &MCEnsembleInfo::getId)
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
    .def("getMCEnsembleInfo", &MCBinsInfo::getMCEnsembleInfo)
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
    .def_static("create", [](std::string s) {
        if (s == "Jackknife")
          return SamplingMode::Jackknife;
        else if (s == "Bootstrap")
          return SamplingMode::Bootstrap;
        throw(std::invalid_argument("Bad SamplingMode")); })
    .def("__str__", [](const SamplingMode &a) {
        switch(a) {
          case Jackknife : return "Jackknife";
          case Bootstrap : return "Bootstrap";
          default        : throw(std::invalid_argument("Bad SamplingMode"));
        }});

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
          default            : throw(std::invalid_argument("Bad ComplexArg"));
        }});

  py::enum_<OperatorInfo::OpKind>(m, "OpKind")
    .value("BasicLapH", OperatorInfo::OpKind::BasicLapH)
    .value("GenIrrep", OperatorInfo::OpKind::GenIrrep);

  py::class_<MCObsInfo>(m, "MCObsInfo")
    .def(py::init<const OperatorInfo &, ComplexArg, bool>())
    .def(py::init<const OperatorInfo &, OperatorInfo &, int, bool, ComplexArg, bool, bool>())
    .def(py::init<const CorrelatorAtTimeInfo &, ComplexArg>())
    .def(py::init<const CorrelatorInfo &, int, bool, ComplexArg, bool, bool>())
    .def(py::init<const std::string &, uint>())
    .def("isVacuum", &MCObsInfo::isVacuum)
    .def("isVEV", &MCObsInfo::isVEV)
    .def("isCorrelatorAtTime", &MCObsInfo::isCorrelatorAtTime)
    .def("isHermitianCorrelatorAtTime", &MCObsInfo::isHermitianCorrelatorAtTime)
    .def("isReweightingFactor", &MCObsInfo::isReweightingFactor)
    .def("isRealPart", &MCObsInfo::isRealPart)
    .def("isImaginaryPart", &MCObsInfo::isImaginaryPart)
    .def("isSimple", &MCObsInfo::isSimple)
    .def("isNonSimple", &MCObsInfo::isNonSimple)
    .def("isPrimary", &MCObsInfo::isPrimary)
    .def("isSecondary", &MCObsInfo::isSecondary)
    .def("isBasicLapH", &MCObsInfo::isBasicLapH)
    .def("isGenIrrep", &MCObsInfo::isGenIrrep)
    .def("isVEVsubtractedCorrelatorAtTime", &MCObsInfo::isVEVsubtractedCorrelatorAtTime)
    .def("isReweightedCorrelatorAtTime", &MCObsInfo::isReweightedCorrelatorAtTime)
    .def("isReweightedVEV", &MCObsInfo::isReweightedVEV)
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
    .def("__hash__", [](const MCObsInfo &a) { return std::hash<std::string>{}(a.str()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<OperatorInfo>(m, "OperatorInfo")
    .def(py::init<const std::string &, OperatorInfo::OpKind>())
    .def(py::init<const BasicLapHOperatorInfo &>())
    .def(py::init<const GenIrrepOperatorInfo &>())
    .def("isBasicLapH", &OperatorInfo::isBasicLapH)
    .def("isGenIrrep", &OperatorInfo::isGenIrrep)
    .def("getBasicLapH", &OperatorInfo::getBasicLapH)
    .def("getGenIrrep", &OperatorInfo::getGenIrrep)
    .def("long_xml", [](const OperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("xml", [](const OperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(false)); })
    .def("op_str", &OperatorInfo::short_output)
    .def("__str__", &OperatorInfo::short_output)
    .def("__repr__", &OperatorInfo::short_output)
    .def("__hash__", [](const OperatorInfo &a) { return std::hash<std::string>{}(a.short_output()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<BasicLapHOperatorInfo>(m, "BasicLapHOperatorInfo")
    .def(py::init<const std::string &>())
    .def("getNumberOfHadrons", &BasicLapHOperatorInfo::getNumberOfHadrons)
    .def("isGlueball", (bool (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::isGlueball)
    .def("isMeson", (bool (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::isMeson)
    .def("isBaryon", (bool (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::isBaryon)
    .def("isTetraquark", (bool (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::isTetraquark)
    .def("isMesonMeson", &BasicLapHOperatorInfo::isMesonMeson)
    .def("isMesonBaryon", &BasicLapHOperatorInfo::isMesonBaryon)
    .def("getMomentum", (Momentum (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getMomentum)
    .def("getXMomentum", (int (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getXMomentum)
    .def("getYMomentum", (int (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getYMomentum)
    .def("getZMomentum", (int (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getZMomentum)
    .def("getLGIrrep", (std::string (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getLGIrrep)
    .def("getLGClebschGordonIdNum", &BasicLapHOperatorInfo::getLGClebschGordonIdNum)
    .def("getLGIrrepRow", &BasicLapHOperatorInfo::getLGIrrepRow)
    .def("getIsospin", &BasicLapHOperatorInfo::getIsospin)
    .def("getIsospinClebschGordonIdNum", &BasicLapHOperatorInfo::getIsospinClebschGordonIdNum)
    .def("getFlavor", (std::string (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getFlavor)
    .def("getFlavorCode", &BasicLapHOperatorInfo::getFlavorCode)
    .def("getStrangeness", (int (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getStrangeness)
    .def("getTetraquarkColorType", (int (BasicLapHOperatorInfo::*)() const) &BasicLapHOperatorInfo::getTetraquarkColorType)
    .def("getHadronFlavor", (std::string (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::getFlavor)
    .def("getHadronStrangeness", (int (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::getStrangeness)
    .def("isHadronGlueball", (bool (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::isGlueball)
    .def("isHadronMeson", (bool (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::isMeson)
    .def("isHadronBaryon", (bool (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::isBaryon)
    .def("isHadronFermion", &BasicLapHOperatorInfo::isFermion)
    .def("isHadronBoson", &BasicLapHOperatorInfo::isBoson)
    .def("getHadronLGIrrep", (std::string (BasicLapHOperatorInfo::*)(uint) const) &BasicLapHOperatorInfo::getLGIrrep)
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
    .def("__hash__", [](const BasicLapHOperatorInfo &a) { return std::hash<std::string>{}(a.short_output()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<GenIrrepOperatorInfo>(m, "GenIrrepOperatorInfo")
    .def(py::init<const std::string &>())
    .def("getMomentum", &GenIrrepOperatorInfo::getMomentum)
    .def("getXMomentum", &GenIrrepOperatorInfo::getXMomentum)
    .def("getYMomentum", &GenIrrepOperatorInfo::getYMomentum)
    .def("getZMomentum", &GenIrrepOperatorInfo::getZMomentum)
    .def("getMomentumSquared", &GenIrrepOperatorInfo::getMomentumSquared)
    .def("hasDefiniteMomentum", &GenIrrepOperatorInfo::hasDefiniteMomentum)
    .def("getLGIrrep", &GenIrrepOperatorInfo::getLGIrrep)
    .def("getLGIrrepRow", &GenIrrepOperatorInfo::getLGIrrepRow)
    .def("getIsospin", &GenIrrepOperatorInfo::getIsospin)
    .def("getStrangeness", &GenIrrepOperatorInfo::getStrangeness)
    .def("getIDName", &GenIrrepOperatorInfo::getIDName)
    .def("getIDIndex", &GenIrrepOperatorInfo::getIDIndex)
    .def("long_xml", [](const GenIrrepOperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("xml", [](const GenIrrepOperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(false)); })
    .def("op_str", &GenIrrepOperatorInfo::short_output)
    .def("__str__", &GenIrrepOperatorInfo::short_output)
    .def("__repr__", &GenIrrepOperatorInfo::short_output)
    .def("__hash__", [](const GenIrrepOperatorInfo &a) { return std::hash<std::string>{}(a.short_output()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<CorrelatorInfo>(m, "CorrelatorInfo")
    .def(py::init<const OperatorInfo &, const OperatorInfo &>())
    .def("getSource", &CorrelatorInfo::getSource)
    .def("getSink", &CorrelatorInfo::getSink)
    .def("isSinkSourceSame", &CorrelatorInfo::isSinkSourceSame)
    .def("xml", [](const CorrelatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("corr_str", [](const CorrelatorInfo &a) {
        return ("snk: " + a.getSink().short_output() + " src: " + a.getSource().short_output()); })
    .def("__str__", [](const CorrelatorInfo &a) { return a.output(false, 2); })
    .def("__repr__", [](const CorrelatorInfo &a) {
        return ("snk_" + a.getSink().short_output() + "-src_" + a.getSource().short_output()); })
    .def("__hash__", [](const CorrelatorInfo &a) { return std::hash<std::string>{}(a.str()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<CorrelatorAtTimeInfo>(m, "CorrelatorAtTimeInfo")
    .def(py::init<const OperatorInfo &, const OperatorInfo &, int, bool, bool, bool>())
    .def(py::init<const CorrelatorInfo &, int, bool, bool, bool>())
    .def ("resetTimeSeparation", &CorrelatorAtTimeInfo::resetTimeSeparation)
    .def("xml", [](const CorrelatorAtTimeInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("__str__", [](const CorrelatorAtTimeInfo &a) { return a.output(false, 2); })
    .def("__repr__", &CorrelatorAtTimeInfo::str)
    .def("__hash__", [](const CorrelatorAtTimeInfo &a) { return std::hash<std::string>{}(a.str()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<CorrelatorMatrixInfo>(m, "CorrelatorMatrixInfo")
    .def(py::init<const std::set<OperatorInfo> &, bool, bool, bool>())
    .def("long_xml", [](const CorrelatorMatrixInfo &a) {
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
    .def("xml", [](const CorrelatorMatrixInfo &a) {
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(false)); })
    .def("__str__", [](const CorrelatorMatrixInfo &a) { return a.output(false, 2); })
    .def("__repr__", &CorrelatorMatrixInfo::str)
    .def("__hash__", [](const CorrelatorMatrixInfo &a) { return std::hash<std::string>{}(a.str()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<Momentum>(m, "Momentum")
    .def(py::init<int, int, int>())
    .def_readwrite("x", &Momentum::x)
    .def_readwrite("y", &Momentum::y)
    .def_readwrite("z", &Momentum::z)
    .def("__str__", &Momentum::getMomentumString)
    .def("__repr__", &Momentum::getMomentumString)
    .def("__hash__", [](const Momentum &a) { return std::hash<std::string>{}(a.getMomentumString()); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<ChiSquareMinimizerInfo>(m, "MinimizerInfo")
    .def(py::init<>())
    .def(py::init<char, double, double, int, char>())
    .def("xml", [](const ChiSquareMinimizerInfo &a) {
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output()); });

  // Data Handlers
  py::class_<RVector>(m, "RVector")
    .def(py::init<>())
    .def(py::init<const std::vector<double> &>())
    .def("array", &RVector::c_vector);

  py::class_<MCEstimate>(m, "MCEstimate")
    .def(py::init<>())
    .def("getFullEstimate", &MCEstimate::getFullEstimate)
    .def("getAverageEstimate", &MCEstimate::getAverageEstimate)
    .def("getSymmetricError", &MCEstimate::getSymmetricError)
    .def("getRelativeError", [](const MCEstimate &a) {
        if (a.getFullEstimate() == 0.0)
          return a.getSymmetricError();

        return (a.getSymmetricError() / std::abs(a.getFullEstimate())); })
    .def("isStatisticallyZero", [](const MCEstimate &a) {
        return (std::abs(a.getFullEstimate()) < a.getSymmetricError()); });

  py::class_<XMLHandler>(m, "XMLHandler")
    .def(py::init<>())
    .def(py::init<const std::string &, const std::string &>())
    .def("set_from_string", &XMLHandler::set_from_string);

  py::class_<FileListInfo>(m, "FileListInfo")
    .def(py::init<const std::string &, int, int, bool>())
    .def("getFileStub", &FileListInfo::getFileStub)
    .def("getMaxFileNumber", &FileListInfo::getMaxFileNumber)
    .def("getMinFileNumber", &FileListInfo::getMinFileNumber)
    .def("isModeOverwrite", &FileListInfo::isModeOverwrite)
    .def("xml", [](const FileListInfo &a) {
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.str()); })
    .def("__str__", [](const FileListInfo &a) {
        return (a.getFileStub() + ".[" + std::to_string(a.getMinFileNumber()) + "," + std::to_string(a.getMaxFileNumber()) + "]"); })
    .def("__repr__", &FileListInfo::str)
    .def("__hash__", [](const FileListInfo &a) { return std::hash<std::string>{}(a.str()); })
    .def(py::self == py::self);

  py::class_<MCObsGetHandler>(m, "MCObsGetHandler")
    .def(py::init<XMLHandler &, const MCBinsInfo &, const MCSamplingInfo &>());

  py::class_<MCObsHandler>(m, "MCObsHandler")
    .def(py::init<MCObsGetHandler &, bool>())
    .def("getBins", (const RVector & (MCObsHandler::*)(const MCObsInfo &) ) &MCObsHandler::getBins)
    .def("getEstimate", (MCEstimate (MCObsHandler::*)(const MCObsInfo &) ) &MCObsHandler::getEstimate)
    .def("setSamplingBegin", &MCObsHandler::setSamplingBegin)
    .def("isSamplingEnd", &MCObsHandler::isSamplingEnd)
    .def("setSamplingNext", &MCObsHandler::setSamplingNext)
    .def("putCurrentSamplingValue", &MCObsHandler::putCurrentSamplingValue)
    .def("writeSamplingValuesToFile", &MCObsHandler::writeSamplingValuesToFile)
    .def("putBins", &MCObsHandler::putBins)
    .def("writeBinsToFile", &MCObsHandler::writeBinsToFile)
    .def("clearData", &MCObsHandler::clearData)
    .def("clearSamplings", &MCObsHandler::clearSamplings)
    .def("eraseData", &MCObsHandler::eraseData)
    .def("eraseSamplings", &MCObsHandler::eraseSamplings);

  py::class_<LaphEnv::BLCorrelatorDataHandler>(m, "BLCorrelatorDataHandler")
    .def(py::init<const std::list<FileListInfo> &, const std::set<CorrelatorInfo> &,
                  const std::set<CorrelatorInfo> &, const MCEnsembleInfo *, bool>())
    .def(py::init<const std::list<FileListInfo> &, const std::set<CorrelatorInfo> &,
                  const std::set<CorrelatorInfo> &, const MCEnsembleInfo *>())
    .def("getCorrelatorSet", &LaphEnv::BLCorrelatorDataHandler::getCorrelatorSet)
    .def("getFileName", &LaphEnv::BLCorrelatorDataHandler::getFileName)
    .def("getKeys", &LaphEnv::BLCorrelatorDataHandler::getKeys)
    .def("getOrderedKeys", &LaphEnv::BLCorrelatorDataHandler::getOrderedKeys);

  py::class_<LaphEnv::BLCorrelatorDataHandler::RecordKey>(m, "BLCorrelatorRecordKey")
    .def(py::init<int, int>())
    .def("getTimeIndex", &LaphEnv::BLCorrelatorDataHandler::RecordKey::getTimeIndex)
    .def("getConfigSerialIndex", &LaphEnv::BLCorrelatorDataHandler::RecordKey::getConfigSerialIndex);

  py::class_<LaphEnv::BLVEVDataHandler>(m, "BLVEVDataHandler")
    .def(py::init<const std::list<FileListInfo> &, const std::set<OperatorInfo> &,
                  const MCEnsembleInfo *, bool>())
    .def(py::init<const std::list<FileListInfo> &, const std::set<OperatorInfo> &,
                  const MCEnsembleInfo *>())
    .def("getOperatorSet", &LaphEnv::BLVEVDataHandler::getOperatorSet)
    .def("getFileName", &LaphEnv::BLVEVDataHandler::getFileName)
    .def("getKeys", &LaphEnv::BLVEVDataHandler::getKeys);

  py::class_<BinsGetHandler>(m, "BinsGetHandler")
    .def(py::init<const MCBinsInfo &, const std::set<std::string> &, bool>())
    .def(py::init<const MCBinsInfo &, const std::set<std::string> &>())
    .def("getFileNames", &BinsGetHandler::getFileNames)
    .def("getKeys", &BinsGetHandler::getKeys);

  py::class_<SamplingsGetHandler>(m, "SamplingsGetHandler")
    .def(py::init<const MCBinsInfo &, const MCSamplingInfo &,
                  const std::set<std::string> &, bool>())
    .def(py::init<const MCBinsInfo &, const MCSamplingInfo &,
                  const std::set<std::string> &>())
    .def("getFileNames", &SamplingsGetHandler::getFileNames)
    .def("getKeys", &SamplingsGetHandler::getKeys);
}

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
#include <matrix.h>
#include <filelist_info.h>
#include <corr_data_handler.h>
#include <vev_data_handler.h>
#include <bins_handler.h>
#include <samplings_handler.h>
#include <mc_estimate.h>
#include <task_utils.h>
#include <momenta.h>

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

  // Info classes
  py::class_<MCEnsembleInfo>(m, "MCEnsembleInfo")
    .def(py::init<const std::string &>())
    .def(py::init<const std::string &, uint, uint, uint, uint, uint, uint>())
    .def("xml", [](const MCEnsembleInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.str()); })
    .def("__str__", [](const MCEnsembleInfo &a) { return a.output(2); })
    .def("__repr__", &MCEnsembleInfo::str)
    .def(py::self == py::self)
    .def(py::self != py::self);

  py::class_<MCBinsInfo>(m, "MCBinsInfo")
    .def(py::init<const MCEnsembleInfo &>())
    .def("setRebin", &MCBinsInfo::setRebin)
    .def("addOmission", &MCBinsInfo::addOmission)
    .def("addOmissions", &MCBinsInfo::addOmissions)
    .def("clearOmissions", &MCBinsInfo::clearOmissions)
    .def("xml", [](const MCBinsInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.str()); })
    .def("__str__", [](const MCBinsInfo &a) { return a.output(2); })
    .def("__repr__", &MCBinsInfo::str)
    .def(py::self == py::self)
    .def(py::self != py::self);

  py::enum_<SamplingMode>(m, "SamplingMode")
    .value("Jackknife", SamplingMode::Jackknife)
    .value("Bootstrap", SamplingMode::Bootstrap);

  py::class_<Bootstrapper>(m, "Bootstrapper")
    .def(py::init<uint, uint,  unsigned long, uint, bool>());

  py::class_<MCSamplingInfo>(m, "MCSamplingInfo")
    .def(py::init<>())
    .def(py::init<uint, unsigned long, uint>())
    .def("setToJackknifeMode", &MCSamplingInfo::setToJackknifeMode)
    .def("setToBootstrapMode", &MCSamplingInfo::setToBootstrapMode)
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
    .value("ImaginaryPart", ComplexArg::ImaginaryPart);

  py::enum_<OperatorInfo::OpKind>(m, "OpKind")
    .value("BasicLapH", OperatorInfo::OpKind::BasicLapH)
    .value("GenIrrep", OperatorInfo::OpKind::GenIrrep);

  py::class_<MCObsInfo>(m, "MCObsInfo")
    .def(py::init<const OperatorInfo &, ComplexArg, bool>())
    .def(py::init<const OperatorInfo &, OperatorInfo &, int, bool, ComplexArg, bool, bool>())
    .def(py::init<const CorrelatorAtTimeInfo &, ComplexArg>())
    .def(py::init<const CorrelatorInfo &, int, bool, ComplexArg, bool, bool>())
    .def(py::init<const std::string &, uint, bool, ComplexArg>())
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
    .def("xml", [](const OperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
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
    .def("xml", [](const BasicLapHOperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
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
    .def("xml", [](const GenIrrepOperatorInfo &a) { 
        py::module ET = py::module::import("xml.etree.ElementTree");
        return ET.attr("fromstring")(a.output(true)); })
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
    .def("__str__", [](const CorrelatorInfo &a) { return a.output(false, 2); })
    .def("__repr__", &CorrelatorInfo::str)
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

  // Data Handlers
  py::class_<RVector>(m, "RVector")
    .def(py::init<>())
    .def(py::init<const std::vector<double> &>())
    .def("array", &RVector::c_vector);

  py::class_<MCEstimate>(m, "MCEstimate")
    .def(py::init<>())
    .def("getFullEstimate", &MCEstimate::getFullEstimate)
    .def("getAverageEstimate", &MCEstimate::getAverageEstimate)
    .def("getSymmetricError", &MCEstimate::getSymmetricError);

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
    .def("writeBinsToFile", &MCObsHandler::writeBinsToFile);

  py::class_<LaphEnv::BLCorrelatorDataHandler>(m, "BLCorrelatorDataHandler")
    .def(py::init<const std::list<FileListInfo> &, const std::set<CorrelatorInfo> &,
                  const std::set<CorrelatorInfo> &, const MCEnsembleInfo *, bool>())
    .def(py::init<const std::list<FileListInfo> &, const std::set<CorrelatorInfo> &,
                  const std::set<CorrelatorInfo> &, const MCEnsembleInfo *>())
    .def("getCorrelatorSet", &LaphEnv::BLCorrelatorDataHandler::getCorrelatorSet)
    .def("getFileName", &LaphEnv::BLCorrelatorDataHandler::getFileName)
    .def("getKeys", &LaphEnv::BLCorrelatorDataHandler::getKeys);

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
    .def("gitFileName", &LaphEnv::BLVEVDataHandler::getFileName)
    .def("getKeys", &LaphEnv::BLVEVDataHandler::getKeys);

  py::class_<BinsGetHandler>(m, "BinsGetHandler")
    .def(py::init<const MCBinsInfo &, const std::set<std::string> &, bool>())
    .def(py::init<const MCBinsInfo &, const std::set<std::string> &>())
    .def("getKeys", &BinsGetHandler::getKeys);

  py::class_<SamplingsGetHandler>(m, "SamplingsGetHandler")
    .def(py::init<const MCBinsInfo &, const MCSamplingInfo &,
                  const std::set<std::string> &, bool>())
    .def(py::init<const MCBinsInfo &, const MCSamplingInfo &,
                  const std::set<std::string> &>())
    .def("getKeys", &SamplingsGetHandler::getKeys);
}

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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(sigmondbind, m) {
  m.doc() = "pybind11 wrapper for sigmond";

  py::class_<RVector>(m, "RVector")
    .def(py::init<>())
    .def(py::init<const std::vector<double> &>());

  py::class_<MCEnsembleInfo>(m, "MCEnsembleInfo")
    .def(py::init<const std::string &>())
    .def(py::init<const std::string &, uint, uint, uint, uint, uint, uint>())
    .def("output", &MCEnsembleInfo::str);

  py::class_<MCBinsInfo>(m, "MCBinsInfo")
    .def(py::init<const MCEnsembleInfo &>())
    .def("setRebin", &MCBinsInfo::setRebin)
    .def("addOmission", &MCBinsInfo::addOmission)
    .def("addOmissions", &MCBinsInfo::addOmissions)
    .def("clearOmissions", &MCBinsInfo::clearOmissions)
    .def("output", &MCBinsInfo::str);

  py::class_<Bootstrapper>(m, "Bootstrapper")
    .def(py::init<uint, uint,  unsigned long, uint, bool>());

  py::class_<MCSamplingInfo>(m, "MCSamplingInfo")
    .def(py::init<>())
    .def(py::init<uint, unsigned long, uint>())
    .def("setToJackknifeMode", &MCSamplingInfo::setToJackknifeMode)
    .def("setToBootstrapMode", &MCSamplingInfo::setToBootstrapMode);

  py::class_<MCObsGetHandler>(m, "MCObsGetHandler")
    .def(py::init<XMLHandler &, const MCBinsInfo &, const MCSamplingInfo &>());

  py::class_<XMLHandler>(m, "XMLHandler")
    .def(py::init<const std::string &, const std::string &>());

  py::enum_<OperatorInfo::OpKind>(m, "OpKind")
    .value("BasicLapH", OperatorInfo::OpKind::BasicLapH)
    .value("GenIrrep", OperatorInfo::OpKind::GenIrrep);

  py::class_<OperatorInfo>(m, "OperatorInfo")
    .def(py::init<const std::string &, OperatorInfo::OpKind>())
    .def(py::init<const GenIrrepOperatorInfo &>());

  py::class_<GenIrrepOperatorInfo>(m, "GenIrrepOperatorInfo")
    .def(py::init<const std::string &>());

  py::class_<CorrelatorInfo>(m, "CorrelatorInfo")
    .def(py::init<const OperatorInfo &, const OperatorInfo &>())
    .def("str", &CorrelatorInfo::str);

  py::class_<CorrelatorAtTimeInfo>(m, "CorrelatorAtTimeInfo")
    .def(py::init<const OperatorInfo &, const OperatorInfo &, int, bool, bool, bool>())
    .def(py::init<const CorrelatorInfo &, int, bool, bool, bool>());

  py::enum_<ComplexArg>(m, "ComplexArg")
    .value("RealPart", ComplexArg::RealPart)
    .value("ImaginaryPart", ComplexArg::ImaginaryPart);

  py::class_<MCObsInfo>(m, "MCObsInfo")
    .def(py::init<const OperatorInfo &, ComplexArg, bool>())
    .def(py::init<const OperatorInfo &, OperatorInfo &, int, bool, ComplexArg, bool, bool>())
    .def(py::init<const CorrelatorAtTimeInfo &, ComplexArg>())
    .def(py::init<const CorrelatorInfo &, int, bool, ComplexArg, bool, bool>())
    .def(py::init<const std::string &, uint, bool, ComplexArg>());

  py::class_<MCObsHandler>(m, "MCObsHandler")
    .def(py::init<MCObsGetHandler &, bool>())
    .def("setSamplingBegin", &MCObsHandler::setSamplingBegin)
    .def("isSamplingEnd", &MCObsHandler::isSamplingEnd)
    .def("setSamplingNext", &MCObsHandler::setSamplingNext)
    .def("putCurrentSamplingValue", &MCObsHandler::putCurrentSamplingValue)
    .def("writeSamplingValuesToFile", &MCObsHandler::writeSamplingValuesToFile)
    .def("putBins", &MCObsHandler::putBins)
    .def("writeBinsToFile", &MCObsHandler::writeBinsToFile);

  py::class_<FileListInfo>(m, "FileListInfo")
    .def(py::init<const std::string &, int, int, bool>());

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
    .def("getConfigSerialIndex", &LaphEnv::BLCorrelatorDataHandler::RecordKey::getTimeIndex);

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

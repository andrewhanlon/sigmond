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
    .def(py::init<const OperatorInfo &, const OperatorInfo &>());

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
}

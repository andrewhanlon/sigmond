#include <scalar_defs.h>
#include <ensemble_info.h>
#include <bins_info.h>
#include <sampling_info.h>
#include <obs_get_handler.h>
#include <xml_handler.h>
#include <mcobs_info.h>
#include <mcobs_handler.h>
#include <operator_info.h>
#include <gen_irrep_operator_info.h>
#include <correlator_info.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(sigmondbind, m) {
  m.doc() = "pybind11 wrapper for sigmond";

  py::class_<MCEnsembleInfo>(m, "MCEnsembleInfo")
    .def(py::init<const std::string &, uint, uint, uint, uint, uint, uint>())
    .def("output", &MCEnsembleInfo::str);

  py::class_<MCBinsInfo>(m, "MCBinsInfo")
    .def(py::init<const MCEnsembleInfo &>())
    .def("setRebin", &MCBinsInfo::setRebin)
    .def("addOmission", &MCBinsInfo::addOmission)
    .def("addOmissions", &MCBinsInfo::addOmissions)
    .def("clearOmissions", &MCBinsInfo::clearOmissions)
    .def("output", &MCBinsInfo::str);

  py::class_<MCSamplingInfo>(m, "MCSamplingInfo")
    .def(py::init<>());

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
    .def(py::init<const OperatorInfo &, const OperatorInfo &, int, bool, bool>())
    .def(py::init<const CorrelatorInfo &, int, bool, bool>());

  py::enum_<ComplexArg>(m, "ComplexArg")
    .value("RealPart", ComplexArg::RealPart)
    .value("ImaginaryParty", ComplexArg::ImaginaryPart);

  py::class_<MCObsInfo>(m, "MCObsInfo")
    .def(py::init<const OperatorInfo &, ComplexArg>())
    .def(py::init<const OperatorInfo &, OperatorInfo &, int, bool, ComplexArg, bool>())
    .def(py::init<const CorrelatorAtTimeInfo &, ComplexArg>())
    .def(py::init<const CorrelatorInfo &, int, bool, ComplexArg, bool>())
    .def(py::init<const std::string &, uint, bool, ComplexArg>());

  py::class_<MCObsHandler>(m, "MCObsHandler")
    .def(py::init<MCObsGetHandler &, bool>())
    .def("putBins", &MCObsHandler::putBins)
    .def("writeBinsToFile", &MCObsHandler::writeBinsToFile);
}

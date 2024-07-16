#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "GaussianMixture.h"
#include "LindeBuzoGray.h"

namespace py = pybind11;

PYBIND11_MODULE(cstrlbg, m) {
    // Wrapper for RArg
    py::class_<RArg>(m, "RArg")
        .def(py::init<>())
        .def_readwrite("indexCmp", &RArg::indexCmp)
        .def_readwrite("indexCls", &RArg::indexCls)
        .def_readwrite("maxVar", &RArg::maxVar)
        .def_readwrite("Mean", &RArg::Mean)
        .def_readwrite("Var", &RArg::Var);

    // Wrapper for RModelStage
    py::class_<RModelStage>(m, "RModelStage")
        .def(py::init<>())
        .def_readwrite("Flag", &RModelStage::Flag)
        .def_readwrite("Cls", &RModelStage::Cls)
        .def_readwrite("indexCmp", &RModelStage::indexCmp)
        .def_readwrite("indexCls", &RModelStage::indexCls)
        .def_readwrite("maxVar", &RModelStage::maxVar)
        .def_readwrite("lnGaussPrb", &RModelStage::lnGaussPrb)
        .def_readwrite("Mean", &RModelStage::Mean)
        .def_readwrite("Var", &RModelStage::Var)
        .def_readwrite("lnVar", &RModelStage::lnVar);

    // Wrapper for TLindeBuzoGray
    py::class_<TLindeBuzoGray>(m, "TLindeBuzoGray")
        .def(py::init<>())
        .def("EvklDistance", &TLindeBuzoGray::EvklDistance)
        .def("RandomInitPoint", &TLindeBuzoGray::RandomInitPoint)
        .def("LindeBuzoGray_N", &TLindeBuzoGray::LindeBuzoGray_N)
        .def("K_MeanGeneral", &TLindeBuzoGray::K_MeanGeneral);

    // Wrapper for TGaussianMixture
    py::class_<TGaussianMixture>(m, "TGaussianMixture")
        .def(py::init<>())
        .def("EM_Alg", &TGaussianMixture::EM_Alg)
        .def("Likelihood", &TGaussianMixture::Likelihood);
}
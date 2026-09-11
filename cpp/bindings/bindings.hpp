#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <traffic_assignment/solve.hpp>

#include <cstdint>

namespace py = pybind11;
using namespace traffic_assignment;

using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;
using IndexArray = py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>;

template <typename T>
py::array_t<double> ToArray(const std::vector<T>& values) {
  py::array_t<double> result(static_cast<py::ssize_t>(values.size()));
  auto buffer = result.mutable_unchecked<1>();
  for (std::size_t i = 0; i < values.size(); ++i) buffer(i) = static_cast<double>(values[i]);
  return result;
}

void BindNetwork(py::module_& m);
void BindOptions(py::module_& m);
void BindResults(py::module_& m);
void BindSolvers(py::module_& m);

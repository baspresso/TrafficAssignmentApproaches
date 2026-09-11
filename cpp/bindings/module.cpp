#include "bindings.hpp"

PYBIND11_MODULE(_core, m) {
  m.doc() = "Python adapters for the traffic_assignment C++ library";
  BindNetwork(m);
  BindOptions(m);
  BindResults(m);
  BindSolvers(m);
}

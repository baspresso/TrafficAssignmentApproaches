#include "detail/tap/data/Link.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
  try {
    for (const long double power : {1.0L, 4.0L, 2.5L}) {
      TrafficAssignment::Link<long double> link(0, 1, 10, 0, 2, 1, power, 0, 0, 1);
      link.flow = 5;
      for (const long double capacity : {10.0L, 20.0L, 10.0L}) {
        link.capacity = capacity;
        const auto expected = 2 * power * std::pow(link.flow / capacity, power - 1) / capacity;
        for (int repeat = 0; repeat < 2; ++repeat) {
          if (std::abs(link.DelayDer() - expected) > 1e-15L) {
            throw std::runtime_error("derivative must follow capacity changes at fixed flow");
          }
        }
        // A trial-flow evaluation must not replace the current-flow cache.
        const auto trial = 2 * power * std::pow(7 / capacity, power - 1) / capacity;
        if (std::abs(link.DelayDer(7) - trial) > 1e-15L ||
            std::abs(link.DelayDer() - expected) > 1e-15L) {
          throw std::runtime_error("trial-flow evaluation corrupts cached derivative");
        }
      }
    }
    std::cout << "Capacity-dependent derivative checks passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

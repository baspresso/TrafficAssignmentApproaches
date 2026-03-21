#ifndef TAPAS_SHIFT_METHOD_FACTORY_H
#define TAPAS_SHIFT_METHOD_FACTORY_H

#include <memory>
#include <string>
#include <map>
#include <functional>
#include "components/TapasShiftMethod.h"
#include "components/TapasNewtonStepShiftMethod.h"
#include "components/TapasLineSearchShiftMethod.h"
#include "components/TapasAdvancedGradientDescentShiftMethod.h"
#include "components/TapasBisectionShiftMethod.h"

namespace TrafficAssignment {

#define REGISTER_TAPAS_SHIFT_METHOD(name, class_name) \
  Register(name, [](const std::vector<Link<T>>& links, T threshold) { \
    return std::make_unique<class_name<T>>(links, threshold); \
  })

/**
 * @brief Singleton factory for TAPAS PAS shift methods.
 *
 * Registers: "NewtonStep" (TapasNewtonStepShiftMethod), "LineSearch" (TapasLineSearchShiftMethod).
 * Note: TapasAdvancedGradientDescentShiftMethod is included but not registered.
 */
template <typename T>
class TapasShiftMethodFactory {
public:
  using Creator = std::function<std::unique_ptr<TapasShiftMethod<T>>(const std::vector<Link<T>>&, T)>;

  TapasShiftMethodFactory(const TapasShiftMethodFactory&) = delete;
  TapasShiftMethodFactory& operator=(const TapasShiftMethodFactory&) = delete;
    
  static TapasShiftMethodFactory& GetInstance() {
    static TapasShiftMethodFactory<T> instance;
    return instance;
  }
    
  void RegisterAll() {
    REGISTER_TAPAS_SHIFT_METHOD("NewtonStep", TapasNewtonStepShiftMethod);
    REGISTER_TAPAS_SHIFT_METHOD("LineSearch", TapasLineSearchShiftMethod);
    REGISTER_TAPAS_SHIFT_METHOD("Bisection", TapasBisectionShiftMethod);
  }
    
  void Register(const std::string& name, Creator creator) {
    registry_[name] = creator;
  }
    
  std::unique_ptr<TapasShiftMethod<T>> Create(const std::string& name, const std::vector<Link<T>>& links, T computation_threshold = T(1e-10)) {
    auto it = registry_.find(name);
    if (it != registry_.end()) {
      return it->second(links, computation_threshold);
    }
    throw std::runtime_error("Unknown shift method: " + name);
  }
    
  std::vector<std::string> GetRegisteredNames() const {
    std::vector<std::string> names;
    for (const auto& pair : registry_) {
      names.push_back(pair.first);
    }
    return names;
  }
    
private:

  TapasShiftMethodFactory() {
    RegisterAll();
  }
    
  std::map<std::string, Creator> registry_;
};

}

#endif // TAPAS_SHIFT_METHOD_FACTORY_H
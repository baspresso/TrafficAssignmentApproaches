#ifndef ROUTE_BASED_SHIFT_METHOD_FACTORY_H
#define ROUTE_BASED_SHIFT_METHOD_FACTORY_H

#include <memory>
#include <string>
#include <map>
#include <functional>
#include "components/RouteBasedShiftMethod.h"
//#include "components/RouteBasedKrylatov2023ShiftMethod.h"
#include "components/RouteBasedNewtonStepShiftMethod.h"

namespace TrafficAssignment {

#define REGISTER_SHIFT_METHOD(name, class_name) \
  Register(name, [](Network<T>& network) { \
    return std::make_unique<class_name<T>>(network); \
  })

template <typename T>
class RouteBasedShiftMethodFactory {
public:
  using Creator = std::function<std::unique_ptr<RouteBasedShiftMethod<T>>(Network<T>&)>;

  RouteBasedShiftMethodFactory(const RouteBasedShiftMethodFactory&) = delete;
  RouteBasedShiftMethodFactory& operator=(const RouteBasedShiftMethodFactory&) = delete;
    
  static RouteBasedShiftMethodFactory& GetInstance() {
    static RouteBasedShiftMethodFactory<T> instance;
    return instance;
  }
    
  void RegisterAll() {
    //REGISTER_SHIFT_METHOD("Krylatov2023", RouteBasedKrylatov2023ShiftMethod);
    REGISTER_SHIFT_METHOD("NewtonStep", RouteBasedNewtonStepShiftMethod);
  }
    
  void Register(const std::string& name, Creator creator) {
    registry_[name] = creator;
  }
    
  std::unique_ptr<RouteBasedShiftMethod<T>> Create(const std::string& name, Network<T>& network) {
    auto it = registry_.find(name);
    if (it != registry_.end()) {
      return it->second(network);
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

  RouteBasedShiftMethodFactory() {
    RegisterAll();
  }
    
  std::map<std::string, Creator> registry_;
};

}

#endif // ROUTE_BASED_SHIFT_METHOD_FACTORY_H
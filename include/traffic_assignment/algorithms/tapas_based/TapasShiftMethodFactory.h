// TapasShiftMethodRegistry.h
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

namespace TrafficAssignment {

#define REGISTER_SHIFT_METHOD(name, class_name) \
    Register(name, [](const std::vector<Link<T>>& links) { \
        return std::make_unique<class_name<T>>(links); \
    })

template <typename T>
class TapasShiftMethodFactory {
public:
    using Creator = std::function<std::unique_ptr<TapasShiftMethod<T>>(const std::vector<Link<T>>&)>;

    TapasShiftMethodFactory(const TapasShiftMethodFactory&) = delete;
    TapasShiftMethodFactory& operator=(const TapasShiftMethodFactory&) = delete;
    
    static TapasShiftMethodFactory& GetInstance() {
        static TapasShiftMethodFactory<T> instance;
        return instance;
    }
    
    void RegisterAll() {
        REGISTER_SHIFT_METHOD("NewtonStep", TapasNewtonStepShiftMethod);
        REGISTER_SHIFT_METHOD("LineSearch", TapasLineSearchShiftMethod);
        //REGISTER_SHIFT_METHOD("AdvancedGradientDescent", TapasAdvancedGradientDescentShiftMethod);       
    }
    
    void Register(const std::string& name, Creator creator) {
        registry_[name] = creator;
    }
    
    std::unique_ptr<TapasShiftMethod<T>> Create(const std::string& name, const std::vector<Link<T>>& links) {
        auto it = registry_.find(name);
        if (it != registry_.end()) {
            return it->second(links);
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

} // namespace TrafficAssignment

#endif // TAPAS_SHIFT_METHOD_FACTORY_H
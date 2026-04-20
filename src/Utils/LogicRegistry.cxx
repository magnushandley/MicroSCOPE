// LogicRegistry.cxx

#include "Utils/LogicRegistry.hxx"
#include <stdexcept>

LogicRegistry& LogicRegistry::Instance() {
    static LogicRegistry instance;
    return instance;
}

void LogicRegistry::Register(const std::string& name, OperationFunc func) {
    fRegistry[name] = func;
}

LogicRegistry::RNode LogicRegistry::Apply(RNode df, const LogicConfig& cfg) const
{

    auto it = fRegistry.find(cfg.type);
    if (it == fRegistry.end()) {
        throw std::runtime_error("Unknown logic type: " + cfg.type);
    }
    
    return it->second(df, cfg);

}
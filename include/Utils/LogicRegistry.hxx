// LogicRegistry.hxx

#pragma once

#include "Utils/LogicConfig.hxx"
#include <ROOT/RDataFrame.hxx>
#include <functional>
#include <unordered_map>

class LogicRegistry {
public:
    using RNode = ROOT::RDF::RNode;
    using OperationFunc = std::function<RNode(RNode, const LogicConfig&)>;

    static LogicRegistry& Instance();

    void Register(const std::string& name, OperationFunc func);

    RNode Apply(RNode df, const LogicConfig& cfg) const;

private:
    std::unordered_map<std::string, OperationFunc> fRegistry;
};
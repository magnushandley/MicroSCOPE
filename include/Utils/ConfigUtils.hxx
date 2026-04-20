// ConfigUtils.hxx

#pragma once
#include "Utils/LogicConfig.hxx"
#include <TEnv.h>
#include <sstream>

inline std::vector<LogicConfig> ParseLogicConfigs(const TEnv& cfg, const std::string& prefix)
{
    std::vector<LogicConfig> out;

    std::stringstream ss(cfg.GetValue((prefix + ".Operations").c_str(), ""));
    std::string name;

    while (ss >> name) {
        LogicConfig c;
        c.name = name;

        std::string base = prefix + ".Operation." + name;

        c.type   = cfg.GetValue((base + ".Type").c_str(), "");
        c.output = cfg.GetValue((base + ".Output").c_str(), "");

        std::stringstream ssInputs(cfg.GetValue((base + ".Inputs").c_str(), ""));
        std::string input;
        while (ssInputs >> input) {
            c.inputs.push_back(input);
        }

        c.expression = cfg.GetValue((base + ".Expression").c_str(), "");

        out.push_back(c);
    }

    return out;
}
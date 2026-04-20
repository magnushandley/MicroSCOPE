// LogicConfig.hxx

#pragma once
#include <string>
#include <vector>

struct LogicConfig {
    std::string name;                 // op1
    std::string type;                 // first_or_default_float
    std::string output;               // trk_score_v_first
    std::vector<std::string> inputs;  // [trk_score_v]
    std::string expression;           // optional
};
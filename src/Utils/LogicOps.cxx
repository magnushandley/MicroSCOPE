// LogicOps.cxx

#include "Utils/LogicRegistry.hxx"
#include <algorithm>

//--------------------------------------------
// First element or default
//--------------------------------------------
static auto first_or_default_float = [](auto df, const LogicConfig& cfg) {
    return df.Define(cfg.output,
        [](const std::vector<float>& v) {
            return v.empty() ? -9999.0f : v[0];
        },
        cfg.inputs
    );
};

//--------------------------------------------
// Sum int vector
//--------------------------------------------
static auto sum_int_vector = [](auto df, const LogicConfig& cfg) {
    return df.Define(cfg.output,
        [](const std::vector<int>& v) {
            int sum = 0;
            for (auto x : v) sum += x;
            return sum;
        },
        cfg.inputs
    );
};

//--------------------------------------------
// Value at max index (e.g. for maximum energy slice)
//--------------------------------------------
static auto value_at_max_index_float_int = [](auto df, const LogicConfig& cfg) {
    return df.Define(cfg.output,
        [](const std::vector<float>& var, const std::vector<int>& E) {
            if (var.empty() || E.empty()) return -9999.0f;
            int idx = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
            return var[idx];
        },
        cfg.inputs
    );
};

static auto value_at_max_index_int_int = [](auto df, const LogicConfig& cfg) {
    return df.Define(cfg.output,
        [](const std::vector<int>& var, const std::vector<int>& E) {
            if (var.empty() || E.empty()) return -9999;
            int idx = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
            return var[idx];
        },
        cfg.inputs
    );
};

static auto min_of_two_vectors_or_default_float = [](auto df, const LogicConfig& cfg) {
    float sentinel = -9999.0f;
    return df.Define(cfg.output,
        [sentinel](const std::vector<float>& v1, const std::vector<float>& v2) {
                if (v1.empty() || v2.empty()) return sentinel;

            float min1 = *std::min_element(v1.begin(), v1.end());
            float min2 = *std::min_element(v2.begin(), v2.end());

            return std::min(min1, min2);
        },
        cfg.inputs
    );
};

static auto max_of_two_vectors_or_default_float = [](auto df, const LogicConfig& cfg) {
    float sentinel = 9999.0f;
    return df.Define(cfg.output,
        [sentinel](const std::vector<float>& v1, const std::vector<float>& v2) {
            if (v1.empty() || v2.empty()) return sentinel;

            float max1 = *std::max_element(v1.begin(), v1.end());
            float max2 = *std::max_element(v2.begin(), v2.end());

            return std::max(max1, max2);
        },
        cfg.inputs
    );
};

static auto filter_expression = [](auto df, const LogicConfig& cfg) {
    std::cout << "Applying filter with expression: " << cfg.expression << std::endl;
    return df.Filter(cfg.expression);
};

static auto define_expression = [](auto df, const LogicConfig& cfg) {
    if (cfg.output.empty()) {
        throw std::runtime_error("define_expression requires Output for operation: " + cfg.name);
    }
    if (cfg.expression.empty()) {
        throw std::runtime_error("define_expression requires Expression for operation: " + cfg.name);
    }

    return df.Define(cfg.output, cfg.expression);

};

//--------------------------------------------
// Registration
//--------------------------------------------
namespace {
bool registered = [](){
    auto& reg = LogicRegistry::Instance();
    reg.Register("first_or_default_float", first_or_default_float);
    reg.Register("sum_int_vector", sum_int_vector);
    reg.Register("value_at_max_index_float_int", value_at_max_index_float_int);
    reg.Register("value_at_max_index_int_int", value_at_max_index_int_int);
    reg.Register("min_of_two_vectors_or_default_float", min_of_two_vectors_or_default_float);
    reg.Register("max_of_two_vectors_or_default_float", max_of_two_vectors_or_default_float);
    reg.Register("filter_expression", filter_expression);
    reg.Register("define_expression", define_expression);
    return true;
}();
}
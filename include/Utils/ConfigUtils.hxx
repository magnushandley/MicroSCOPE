// ConfigUtils.hxx

#pragma once
#include "Utils/LogicConfig.hxx"
#include <TEnv.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

struct PlotConfig {
    std::string name;
    std::string column;
    std::string xTitle;
    std::string yTitle = "Count";
    int         nBins = 0;
    double      xMin = 0.0;
    double      xMax = 0.0;
    std::string outputName;
    std::string histNamePrefix;
    bool        logY = false;
    std::string valueMode = "direct";
    std::string weightColumn;
    bool        enableSystematics = true;
    std::optional<double> showLowerCut;
    std::optional<double> showUpperCut;
};

enum class SampleType {
    BeamOff,
    Overlay,
    Dirt,
    Signal,
    Data,
    DetectorVariationCV,
    DetectorVariation
};

inline std::string ToLowerCopy(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

inline bool StringContainsCaseInsensitive(const std::string& value, const std::string& token)
{
    return ToLowerCopy(value).find(ToLowerCopy(token)) != std::string::npos;
}

inline std::string SampleTypeName(SampleType type)
{
    switch (type) {
    case SampleType::BeamOff: return "beamoff";
    case SampleType::Overlay: return "overlay";
    case SampleType::Dirt:    return "dirt";
    case SampleType::Signal:  return "signal";
    case SampleType::Data:    return "data";
    case SampleType::DetectorVariationCV: return "detvarcv";
    case SampleType::DetectorVariation:   return "detvar";
    }
    return "unknown";
}

inline SampleType ParseSampleType(const std::string& token, const std::string& key)
{
    if (token == "beamoff") return SampleType::BeamOff;
    if (token == "overlay") return SampleType::Overlay;
    if (token == "dirt")    return SampleType::Dirt;
    if (token == "signal")  return SampleType::Signal;
    if (token == "data")    return SampleType::Data;
    if (token == "detvarcv") return SampleType::DetectorVariationCV;
    if (token == "detvar")   return SampleType::DetectorVariation;

    throw std::runtime_error("[Config] Invalid sample type in " + key + ": " + token
                             + " (allowed: beamoff, overlay, dirt, signal, data, detvarcv, detvar)");
}

inline std::vector<SampleType> ParseSampleTypes(const std::string& typesString,
                                                const std::string& key)
{
    std::vector<SampleType> types;
    std::stringstream ssTypes{typesString};
    std::string typeToken;

    while (ssTypes >> typeToken) {
        if (!typeToken.empty() && typeToken.back() == ',') {
            typeToken.pop_back();
        }
        if (!typeToken.empty()) {
            types.push_back(ParseSampleType(typeToken, key));
        }
    }

    return types;
}

inline SampleType InferPlotSampleTypeFromLabel(const std::string& label)
{
    if (StringContainsCaseInsensitive(label, "signal")) return SampleType::Signal;
    if (StringContainsCaseInsensitive(label, "data"))   return SampleType::Data;
    if (StringContainsCaseInsensitive(label, "overlay")) return SampleType::Overlay;
    if (StringContainsCaseInsensitive(label, "dirt")) return SampleType::Dirt;
    return SampleType::BeamOff;
}

inline bool IsOverlaySample(SampleType type) { return type == SampleType::Overlay; }
inline bool IsDirtSample(SampleType type) { return type == SampleType::Dirt; }
inline bool IsSignalSample(SampleType type) { return type == SampleType::Signal; }
inline bool IsDataSample(SampleType type) { return type == SampleType::Data; }
inline bool IsDetectorVariationCVSample(SampleType type) { return type == SampleType::DetectorVariationCV; }
inline bool IsDetectorVariationSample(SampleType type) { return type == SampleType::DetectorVariation; }
inline bool IsDetectorVariationInputSample(SampleType type)
{
    return IsDetectorVariationCVSample(type) || IsDetectorVariationSample(type);
}
inline bool IsPlottableSample(SampleType type) { return !IsDetectorVariationInputSample(type); }

inline bool ConfigHasKey(const TEnv& cfg, const std::string& key)
{
    return cfg.Defined(key.c_str()) != 0;
}

inline std::string RequireConfigValue(const TEnv& cfg, const std::string& key)
{
    if (!ConfigHasKey(cfg, key)) {
        throw std::runtime_error("[Config] Missing required config key: " + key);
    }
    return cfg.GetValue(key.c_str(), "");
}

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

inline std::vector<PlotConfig> ParsePlotConfigs(const TEnv& cfg, const std::string& prefix)
{
    std::vector<PlotConfig> out;

    std::stringstream ss(cfg.GetValue((prefix + ".Plots").c_str(), ""));
    std::string name;
    std::unordered_set<std::string> seenNames;

    while (ss >> name) {
        if (!seenNames.insert(name).second) {
            throw std::runtime_error("[Config] Duplicate plot name in " + prefix + ".Plots: " + name);
        }

        PlotConfig plot;
        plot.name = name;

        const std::string base = prefix + ".Plot." + name;
        plot.column = RequireConfigValue(cfg, base + ".Column");
        plot.xTitle = cfg.GetValue((base + ".XTitle").c_str(), plot.column.c_str());
        plot.yTitle = cfg.GetValue((base + ".YTitle").c_str(), "Count");

        const std::string nBinsKey = base + ".NBins";
        const std::string xMinKey = base + ".XMin";
        const std::string xMaxKey = base + ".XMax";

        plot.nBins = std::stoi(RequireConfigValue(cfg, nBinsKey));
        plot.xMin = std::stod(RequireConfigValue(cfg, xMinKey));
        plot.xMax = std::stod(RequireConfigValue(cfg, xMaxKey));

        plot.outputName = cfg.GetValue((base + ".OutputName").c_str(),
                                       ("preselection_full_hist_" + name).c_str());
        plot.histNamePrefix = cfg.GetValue((base + ".HistNamePrefix").c_str(),
                                           ("preselection_hist_" + name + "_").c_str());
        plot.logY = cfg.GetValue((base + ".LogY").c_str(), false);
        plot.weightColumn = cfg.GetValue((base + ".WeightColumn").c_str(), "");
        plot.enableSystematics = cfg.GetValue((base + ".EnableSystematics").c_str(), true);

        const auto parseOptionalCut = [&](const std::string& setting) -> std::optional<double> {
            const std::string key = base + "." + setting;
            if (!ConfigHasKey(cfg, key)) {
                return std::nullopt;
            }

            const std::string rawValue = RequireConfigValue(cfg, key);
            std::stringstream valueStream(rawValue);
            double value = 0.0;
            std::string trailing;
            if (!(valueStream >> value) || (valueStream >> trailing) || !std::isfinite(value)) {
                throw std::runtime_error("[Config] Invalid numeric value for " + key
                                         + ": " + rawValue);
            }
            return value;
        };

        plot.showLowerCut = parseOptionalCut("ShowLowerCut");
        plot.showUpperCut = parseOptionalCut("ShowUpperCut");

        plot.valueMode = cfg.GetValue((base + ".ValueMode").c_str(), "direct");
        std::transform(plot.valueMode.begin(), plot.valueMode.end(), plot.valueMode.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

        if (plot.valueMode != "direct" && plot.valueMode != "first_element") {
            throw std::runtime_error("[Config] Invalid ValueMode for " + base
                                     + ".ValueMode: " + plot.valueMode
                                     + " (allowed: direct, first_element)");
        }

        if (plot.nBins <= 0) {
            throw std::runtime_error("[Config] Invalid NBins for " + base + ".NBins: must be > 0");
        }

        if (!(plot.xMax > plot.xMin)) {
            throw std::runtime_error("[Config] Invalid range for " + base
                                     + ": XMax must be greater than XMin");
        }

        const auto validateCutRange = [&](const std::optional<double>& cut,
                                          const std::string& setting) {
            if (cut && (*cut < plot.xMin || *cut > plot.xMax)) {
                throw std::runtime_error("[Config] Invalid " + base + "." + setting
                                         + ": cut position must be within [XMin, XMax]");
            }
        };
        validateCutRange(plot.showLowerCut, "ShowLowerCut");
        validateCutRange(plot.showUpperCut, "ShowUpperCut");

        if (plot.showLowerCut && plot.showUpperCut
            && *plot.showLowerCut > *plot.showUpperCut) {
            throw std::runtime_error("[Config] Invalid cut range for " + base
                                     + ": ShowLowerCut must not exceed ShowUpperCut");
        }

        out.push_back(std::move(plot));
    }

    return out;
}

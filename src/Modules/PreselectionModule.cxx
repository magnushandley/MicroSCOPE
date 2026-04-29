#include "Modules/PreselectionModule.hxx"
#include "Utils/Plotter.hxx"
#include "Utils/SystematicsUtil.hxx"
#include "Utils/TimingUtils.hxx"

#include <TEnv.h>
#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <algorithm>
#include <cmath>
#include <sstream>

using namespace Analysis;

namespace {

struct PlotRuntime {
    PlotConfig        config;
    std::vector<TH1D> histograms;
    std::vector<TH1D> systVarianceHists;
};

std::string DefaultPlotWeightColumn(SampleType sampleType)
{
    return (IsOverlaySample(sampleType) || IsDirtSample(sampleType)) ? "weight_cv" : "";
}

TH1D MakeEmptyVarianceHist(const TH1D& nominalHist, const std::string& suffix)
{
    TH1D emptyVar = nominalHist;
    emptyVar.SetName((std::string(nominalHist.GetName()) + suffix).c_str());
    emptyVar.Reset("ICES");
    return emptyVar;
}

} // namespace

//------------------------------------------------------------------------------
PreselectionModule::PreselectionModule(const TEnv& cfg)
    : Module(cfg)
    , fTreeName     (cfg.GetValue("Preselection.TreeName","nuselection/NeutrinoSelectionFilter"  ))
    , fRunLabel     (cfg.GetValue("Global.RunLabel","run_x") )
    , fMakePlots   (cfg.GetValue("Preselection.MakePlots", false))
{
    
    //----------------------------------------------------------------------
    // Parse the comma‑separated list of cuts, this is currently messy but was needed to ensure that the proper format of 
    // cutstring was passed to RDF::Filter
    //----------------------------------------------------------------------

    const std::string cutsString = cfg.GetValue("Preselection.Cuts", "");
    std::stringstream ssCuts{cutsString};
    std::string cutToken;

    while (std::getline(ssCuts, cutToken, ',')) {
        // trim leading/trailing whitespace (but keep quotes if present)
        cutToken.erase(0, cutToken.find_first_not_of(" \t\n\r"));
        cutToken.erase(cutToken.find_last_not_of(" \t\n\r") + 1);

        // Remove optional surrounding quotes so the expression is valid for RDataFrame::Filter
        if (cutToken.size() >= 2 &&
           ((cutToken.front() == '"'  && cutToken.back() == '"') ||
            (cutToken.front() == '\'' && cutToken.back() == '\'')))
        {
            cutToken = cutToken.substr(1, cutToken.size() - 2);
        }

        if (!cutToken.empty())
            cuts.push_back(cutToken);
    }
    
    if (cuts.empty()) {
        throw std::runtime_error("[Preselection] No cuts specified!");
    }

    std::stringstream ssKeep{cfg.GetValue("Preselection.Keep", "")};
    std::string keepItem;
    while (ssKeep >> keepItem) {
        if (keepItem.back()==',') keepItem.pop_back();
        fVarsToKeep.push_back(keepItem); 
    }

    std::stringstream ssInput{cfg.GetValue("Preselection.InputFiles", "")};
    std::string inputItem;
    while (ssInput >> inputItem) {
        if (inputItem.back()==',') inputItem.pop_back();
        fInputFiles.push_back(inputItem);
    }

    std::stringstream ssOut{cfg.GetValue("Preselection.Outputs", "")};
    std::string outItem;
    while (ssOut >> outItem) {
        if (outItem.back()==',') outItem.pop_back();
        fOutFiles.push_back(outItem);
    }

    const std::string labelsString = cfg.GetValue("Preselection.SampleLabels", "");
    std::stringstream ssLabels{labelsString};
    std::string labelToken;
    while (std::getline(ssLabels, labelToken, ',')) {
        // trim leading/trailing whitespace
        labelToken.erase(0, labelToken.find_first_not_of(" \t\n\r"));
        labelToken.erase(labelToken.find_last_not_of(" \t\n\r") + 1);

        // Remove optional surrounding quotes so labels like "Run 3 data" are kept intact
        if (labelToken.size() >= 2 &&
           ((labelToken.front() == '"'  && labelToken.back() == '"') ||
            (labelToken.front() == '\'' && labelToken.back() == '\'')))
        {
            labelToken = labelToken.substr(1, labelToken.size() - 2);
        }

        if (!labelToken.empty()) {
            fSampleLabels.push_back(labelToken);
        }
    }

    std::stringstream ssWeights{cfg.GetValue("Preselection.SampleWeights", "")};
    double weight;
    while (ssWeights >> weight) {
        fSampleWeights.push_back(weight);
    }

    const std::string sampleTypesKey = "Preselection.SampleTypes";
    fSampleTypes = ParseSampleTypes(RequireConfigValue(cfg, sampleTypesKey), sampleTypesKey);

    if (fVarsToKeep.empty()) {
        throw std::runtime_error("[Preselection] No variables to keep specified!");
    }

    const std::size_t nSamples = fInputFiles.size();
    if (nSamples == 0) {
        throw std::runtime_error("[Preselection] No input files specified!");
    }
    if (fOutFiles.size() != nSamples) {
        throw std::runtime_error("[Preselection] Outputs count (" + std::to_string(fOutFiles.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }
    if (fSampleLabels.size() != nSamples) {
        throw std::runtime_error("[Preselection] SampleLabels count (" + std::to_string(fSampleLabels.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }
    if (fSampleTypes.size() != nSamples) {
        throw std::runtime_error("[Preselection] SampleTypes count (" + std::to_string(fSampleTypes.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }
    if (fSampleWeights.size() != nSamples) {
        throw std::runtime_error("[Preselection] SampleWeights count (" + std::to_string(fSampleWeights.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }

    if (fMakePlots) {
        fPlotConfigs = ParsePlotConfigs(cfg, "Preselection");
        if (fPlotConfigs.empty()) {
            throw std::runtime_error("[Preselection] MakePlots is enabled but Preselection.Plots is empty.");
        }
    }
}

std::vector<std::unique_ptr<ROOT::RDataFrame>>
PreselectionModule::BuildDataFrames(const std::vector<std::string>& files,
                                     const std::string& treeName) const
{
    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfVec;
    for (const auto& fname : files) {
        auto file = TFile::Open(fname.c_str());
        if (!file || file->IsZombie()) {
            throw std::runtime_error("[Preselection] Cannot open file: " + fname);
        }
        auto tree = file->Get<TTree>(treeName.c_str());
        if (!tree) {
            throw std::runtime_error("[Preselection] Cannot find tree: " + treeName);
        }
        auto RDF = std::make_unique<ROOT::RDataFrame>(*tree);
        dfVec.push_back(std::move(RDF));
    }
    if (dfVec.empty())
        throw std::runtime_error("[Preselection] No data frames created!");

    return dfVec;
}

Long64_t PreselectionModule::EntryCount() const
{
    
    if (dfVec.size() == 0) {
        throw std::runtime_error("[Preselection] DataFrames not initialised!");
    }

    Long64_t totalEntries = 0;
    for (const auto& df : dfVec) {
        totalEntries += df->Count().GetValue();
    }
    return totalEntries;
}

void PreselectionModule::Initialise()
{
    dfVec = BuildDataFrames(fInputFiles, fTreeName);

     // One RNode per sample, initially pointing at the un-filtered DataFrame
     // This again feels messy because we initialise a set of pointers to rdataframes, but then RDF::Filter returns RNodes
     // and we have to keep track of those instead.
     
    std::vector<ROOT::RDF::RNode> nodes;
    nodes.reserve(dfVec.size());
    for (auto &dfPtr : dfVec) nodes.emplace_back(*dfPtr);

    // Apply every cut in sequence
    for (const auto &cut : cuts) {
        std::cout << "\n[Preselection] Cut: " << cut << '\n';

        for (std::size_t i = 0; i < nodes.size(); ++i) {
            auto before = nodes[i].Count().GetValue();
            std::cout << "    " << fSampleLabels[i] << " before: " << before << '\n';

            nodes[i] = nodes[i].Filter(cut);      
            auto after  = nodes[i].Count().GetValue();
            std::cout << "    " << fSampleLabels[i] << " after : " << after  << '\n';
        }
    }

    ROOT::RDF::RSnapshotOptions opt;
    opt.fMode = "RECREATE";
    opt.fCompressionAlgorithm = ROOT::kZLIB;
    opt.fCompressionLevel     = 4;

    // I manually add branches at the slimmer with different levels of weighting: for the tutorial, I disable systematics
    // but am leaving these here for reference
    SystematicsConfig systConfig;
    systConfig.genieMultisimBranch = "weightsGenie";
    systConfig.genieCVWeightBranch = "weight_cv_untuned";
    systConfig.genieGlobalCVWeightBranch = "weight_cv";
    systConfig.ppfxMultisimBranch = "weightsPPFX";
    systConfig.ppfxCVWeightBranch = "weight_cv_noppfx";
    systConfig.ppfxGlobalCVWeightBranch = "weight_cv";
    systConfig.reintMultisimBranch = "weightsReint";
    systConfig.reintCVWeightBranch = "weight_cv";
    systConfig.reintGlobalCVWeightBranch = "weight_cv";

    std::vector<PlotRuntime> plots;
    plots.reserve(fPlotConfigs.size());
    for (const auto& plotConfig : fPlotConfigs) {
        plots.push_back({plotConfig, {}, {}});
    }

    auto fillHistogramsForSample = [&](ROOT::RDF::RNode& node,
                                       const std::string& sampleLabel,
                                       SampleType sampleType) {
        for (auto& plot : plots) {
            const bool useFirstElement = plot.config.valueMode == "first_element";
            const std::string weightCol = plot.config.weightColumn.empty()
                ? DefaultPlotWeightColumn(sampleType)
                : plot.config.weightColumn;

            if (!weightCol.empty()) {
                std::cout << "    Applying " << weightCol
                          << " to plot " << plot.config.name
                          << " for sample: " << sampleLabel << "\n";
            }

            plot.histograms.push_back(
                Plotter::CreateTH1DFromRNode(
                    node,
                    plot.config.histNamePrefix + sampleLabel,
                    plot.config.column,
                    plot.config.xTitle,
                    plot.config.yTitle,
                    plot.config.nBins,
                    plot.config.xMin,
                    plot.config.xMax,
                    useFirstElement,
                    false,
                    weightCol));

            TH1D& nominalHist = plot.histograms.back();
            if (!plot.config.enableSystematics) {
                plot.systVarianceHists.push_back(MakeEmptyVarianceHist(nominalHist, "_systDisabledVar"));
            } else if (IsOverlaySample(sampleType)) {
                std::cout << "    Computing systematic variance histograms for sample: " << sampleLabel
                          << " and plot: " << plot.config.name << "\n";
                SystematicsUtil sysUtil;
                TH1D overlaySystHist = sysUtil.RunAllMultisimSystematics(
                    nominalHist,
                    node,
                    plot.config.column,
                    systConfig);
                plot.systVarianceHists.push_back(std::move(overlaySystHist));
            } else if (IsDirtSample(sampleType)) {
                TH1D dirtVar = nominalHist;
                dirtVar.SetName((std::string(nominalHist.GetName()) + "_dirtNormVar").c_str());
                dirtVar.Reset("ICES");

                for (int bin = 1; bin <= dirtVar.GetNbinsX(); ++bin) {
                    const double content = nominalHist.GetBinContent(bin);
                    const double variance = (0.75 * content) * (0.75 * content);
                    dirtVar.SetBinContent(bin, variance);
                    dirtVar.SetBinError(bin, 0.0);
                }

                plot.systVarianceHists.push_back(std::move(dirtVar));
            } else {
                plot.systVarianceHists.push_back(MakeEmptyVarianceHist(nominalHist, "_emptyVar"));
            }
        }
    };

    auto plotAllHistograms = [&]() {
        for (auto& plot : plots) {
            TH1D totalVarianceHist(
                plot.config.outputName.c_str(),
                (plot.config.outputName + " Total Variance").c_str(),
                plot.config.nBins,
                plot.config.xMin,
                plot.config.xMax);

            if (plot.systVarianceHists.size() != fSampleWeights.size()) {
                std::cerr << "[Preselection] ERROR: systVarianceHists size (" << plot.systVarianceHists.size()
                          << ") != fSampleWeights size (" << fSampleWeights.size()
                          << ") for plot " << plot.config.outputName << "\n";
            }

            const size_t n = std::min(plot.systVarianceHists.size(), fSampleWeights.size());
            for (size_t i = 0; i < n; ++i) {
                TH1D varHist = plot.systVarianceHists[i];
                const double weight = fSampleWeights[i];
                varHist.Scale(weight * weight);
                totalVarianceHist.Add(&varHist);
            }

            std::cout << "Total variance for plot " << plot.config.outputName << ":\n";
            for (int bin = 1; bin <= totalVarianceHist.GetNbinsX(); ++bin) {
                std::cout << "  Bin " << bin << ": " << totalVarianceHist.GetBinContent(bin) << "\n";
            }

            Plotter::FullDataMCSignalPlot(
                plot.histograms,
                fSampleLabels,
                fSampleTypes,
                plot.config.outputName,
                plot.config.logY,
                fSampleWeights,
                0.7,
                1.3,
                &totalVarianceHist);
        }
    };

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        std::cout << "\n[Preselection] Writing output for sample: " << fSampleLabels[i] << '\n';
        std::cout << "    to file: " << fOutFiles[i] << '\n';
        nodes[i].Snapshot(fTreeName, fOutFiles[i], fVarsToKeep, opt);

        if (!fMakePlots) continue;

        fillHistogramsForSample(nodes[i], fSampleLabels[i], fSampleTypes[i]);
    }

    if (!fMakePlots) return;
    plotAllHistograms();
}

void PreselectionModule::Finalise()
{
    // Nothing to do here
}

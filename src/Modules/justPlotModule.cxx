#include "Modules/justPlotModule.hxx"
#include "Utils/Plotter.hxx"

#include <TEnv.h>
#include <TFile.h>
#include <TString.h>
#include <algorithm>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <ROOT/RDataFrame.hxx>
#include <TTree.h>
#include <TH1D.h>

using namespace Analysis;

//------------------------------------------------------------------------------
PlotterModule::PlotterModule(const TEnv& cfg)
    : Module(cfg)
    , fTreeName   (cfg.GetValue("Plotter.TreeName", "nuselection/NeutrinoSelectionFilter"))
{
    
    std::stringstream ssInput{cfg.GetValue("Plotter.InputFiles", "")};
    std::string inputItem;
    while (ssInput >> inputItem) {
        if (inputItem.back()==',') inputItem.pop_back();
        std::cout << "Adding input file: " << inputItem << " (n=" << fInputFiles.size() << ")" << std::endl;
        fInputFiles.push_back(inputItem);
    }

    std::stringstream ssLabels{cfg.GetValue("Plotter.SampleLabels", "")};
    std::string label;
    while (ssLabels >> label) {
        if (label.back()==',') label.pop_back();
        fSampleLabels.push_back(label);
    }

    std::stringstream ssWeights{cfg.GetValue("Plotter.SampleWeights", "")};
    double weight;
    while (ssWeights >> weight) {
        fSampleWeights.push_back(weight);
    }
}

//------------------------------------------------------------------------------
std::vector<std::unique_ptr<ROOT::RDataFrame>>
PlotterModule::BuildDataFrames(const std::vector<std::string>& files,
                               const std::string& treeName) const
{
    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfVec;
    for (const auto& fname : files) {
        auto file = TFile::Open(fname.c_str());
        if (!file || file->IsZombie()) {
            throw std::runtime_error("[Plotter] Cannot open file: " + fname);
        }
        auto tree = file->Get<TTree>(treeName.c_str());
        if (!tree) {
            throw std::runtime_error("[Plotter] Cannot find tree: " + treeName);
        }
        auto RDF = std::make_unique<ROOT::RDataFrame>(*tree);
        dfVec.push_back(std::move(RDF));
    }
    if (dfVec.empty())
        throw std::runtime_error("[Plotter] No data frames created!");

    return dfVec;
}



Long64_t PlotterModule::EntryCount() const
{
    return 1; // Dummy, nothing per-event
}

//------------------------------------------------------------------------------
void PlotterModule::Initialise()
{
    auto dfVec = BuildDataFrames(fInputFiles, fTreeName);

    std::vector<ROOT::RDF::RNode> nodes;
    nodes.reserve(dfVec.size());
    for (auto &dfPtr : dfVec) nodes.emplace_back(*dfPtr);

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        auto before = nodes[i].Count().GetValue();
        std::cout << "    " << fSampleLabels[i] << " before: " << before << '\n';

        nodes[i] = nodes[i].Define("logit_bdt",
            [](float score) {
                const float eps = 1e-6f;
                const float s = std::min(std::max(score, eps), 1.0f - eps);
                return std::log(s / (1.0f - s));
            },
            {"bdt_score"});
    }

    std::vector<TH1D> bdtScoreVec;
    for (size_t i = 0; i < nodes.size(); ++i)
        bdtScoreVec.push_back(
            Plotter::CreateTH1DFromRNode(
                nodes[i],
                ("logit_bdt_" + fSampleLabels[i]).c_str(),
                "logit_bdt", 
                "Logit BDT Score",
                "Count",
                11, -5.0, 6.0)); 

    Plotter::FullDataMCSignalPlot(bdtScoreVec,
                        fSampleLabels,
                        "bdt_score_full_hist",
                        false, // logy
                        fSampleWeights);
    
}

//------------------------------------------------------------------------------
void PlotterModule::Finalise()
{
    // Nothing to do – Snapshot already wrote the slimmed tree.
}
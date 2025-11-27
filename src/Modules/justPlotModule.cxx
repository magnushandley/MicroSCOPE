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

// Helper function to save TH1D histograms to a ROOT file with proper weights
// and statistical errors.
void PlotterModule::SaveHistograms(const std::vector<TH1D>& hists,
                    const std::vector<std::string>& labels,
                    const std::vector<double>& weights,
                    const std::string& fileName)
{
    if (hists.size() != labels.size()) {
        throw std::runtime_error("[Plotter] Number of histograms and labels do not match in SaveHistograms");
    }

    TFile outFile(fileName.c_str(), "RECREATE");
    if (outFile.IsZombie()) {
        throw std::runtime_error("[Plotter] Cannot create output file: " + fileName);
    }

    for (std::size_t i = 0; i < hists.size(); ++i) {
        // Make a local copy so we do not modify the original histogram vector
        TH1D h = hists[i];

        // Ensure we have Sumw2 so that errors are stored and scaled correctly
        h.Sumw2();

        // Apply the sample weight if provided
        if (i < weights.size()) {
            h.Scale(weights[i]);
        }

        // Name the histogram according to the sample label
        h.SetName(labels[i].c_str());

        // Detach from any existing directory and write to the output file
        h.SetDirectory(&outFile);
        h.Write();
    }

    outFile.Close();
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
        
        //tmp for bdt testing
        float testFraction;
        if (fSampleLabels[i].find("data") != std::string::npos) {
            //data sample
            testFraction = 1.0;
            std::cout << "data sample found" << std::endl;
            std::cout << "testFraction: " << testFraction << std::endl;
        } else {
            //mc sample
            testFraction = 0.4;
            std::cout << "mc sample found" << std::endl;
            std::cout << "testFraction: " << testFraction << std::endl;
        }
        int nEntries = nodes[i].Count().GetValue();
        int nTest = static_cast<Long64_t>(std::round(nEntries * testFraction));
        nodes[i] = nodes[i].Define("logit_bdt",
            [](float score) {
                //const float eps = 1e-9f;
                //const float s = std::min(std::max(score, eps), 1.0f - eps);
                float s = (score + 1.0f) / 2.0f; //rescale from [-1,1] to [0,1]
                return std::log(s / (1.0f - s));
            },
            {"bdt_score"})
            .Range((nEntries-nTest), nEntries); //take only test sample
        std::cout << "    " << fSampleLabels[i] << " after: " << nodes[i].Count().GetValue() << '\n';
    }

    std::vector<TH1D> bdtScoreVec;
    for (size_t i = 0; i < nodes.size(); ++i)
        bdtScoreVec.push_back(
            Plotter::CreateTH1DFromRNode(
                nodes[i],
                ("logit_bdt_score_" + fSampleLabels[i]).c_str(),
                "logit_bdt", 
                "Logit BDT Score",
                "Count",
                10, -5.0, 5.0f)); 

    Plotter::FullDataMCSignalPlot(bdtScoreVec,
                        fSampleLabels,
                        "bdt_score_full_hist_tmva",
                        false, // logy
                        fSampleWeights);

    // Save the BDT score histograms to a ROOT file with proper weights and errors
    SaveHistograms(bdtScoreVec, fSampleLabels, fSampleWeights, "bdt_score_histograms.root");
}

//------------------------------------------------------------------------------
void PlotterModule::Finalise()
{
    // Nothing to do – Snapshot already wrote the slimmed tree.
}

//------------------------------------------------------------------------------

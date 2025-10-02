#include "Modules/BDTTrainModule.hxx"
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
#include <ROOT/RDFHelpers.hxx>
#include <TTree.h>
#include <TH1D.h>
#include <TChain.h>
#include <cstdio>

using namespace Analysis;

//------------------------------------------------------------------------------
BDTTrainModule::BDTTrainModule(const TEnv& cfg)
    : Module(cfg)
    , fTreeName   (cfg.GetValue("BDTTrainModule.TreeName", "nuselection/NeutrinoSelectionFilter"))
    , fTrainFraction(cfg.GetValue("BDTTrainModule.TrainFraction", 0.8f))
{

    std::stringstream ssInput{cfg.GetValue("BDTTrainModule.InputFiles", "")};
    std::string inputItem;
    while (ssInput >> inputItem) {
        if (inputItem.back()==',') inputItem.pop_back();
        std::cout << "Adding input file: " << inputItem << " (n=" << fInputFiles.size() << ")" << std::endl;
        fInputFiles.push_back(inputItem);
    }

    std::stringstream ssLabels{cfg.GetValue("BDTTrainModule.SampleLabels", "")};
    std::string label;
    while (ssLabels >> label) {
        if (label.back()==',') label.pop_back();
        fSampleLabels.push_back(label);
    }

    std::stringstream ssKeep{cfg.GetValue("BDTTrainModule.Keep", "")};
    std::string keepItem;
    while (ssKeep >> keepItem) {
        if (keepItem.back()==',') keepItem.pop_back();
        fVarsToKeep.push_back(keepItem);
    }

    std::stringstream ssTrainVars{cfg.GetValue("BDTTrainModule.TrainVariables", "")};
    std::string var;
    while (ssTrainVars >> var) {
        if (var.back()==',') var.pop_back();
        fTrainVars.push_back(var);
    }

    std::stringstream ssWeights{cfg.GetValue("BDTTrainModule.SampleWeights", "")};
    double weight;
    while (ssWeights >> weight) {
        fSampleWeights.push_back(weight);
    }
}

//------------------------------------------------------------------------------
std::vector<std::unique_ptr<ROOT::RDataFrame>>
BDTTrainModule::BuildDataFrames(const std::vector<std::string>& files,
                               const std::string& treeName) const
{
    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfVec;
    for (const auto& fname : files) {
        auto file = TFile::Open(fname.c_str());
        if (!file || file->IsZombie()) {
            throw std::runtime_error("[BDTTrainModule] Cannot open file: " + fname);
        }
        auto tree = file->Get<TTree>(treeName.c_str());
        if (!tree) {
            throw std::runtime_error("[BDTTrainModule] Cannot find tree: " + treeName);
        }
        auto RDF = std::make_unique<ROOT::RDataFrame>(*tree);
        dfVec.push_back(std::move(RDF));
    }
    if (dfVec.empty())
        throw std::runtime_error("[BDTTrainModule] No data frames created!");

    return dfVec;
}

std::vector<std::string>
BDTTrainModule::BuildTestTrainSamples(std::vector<ROOT::RDF::RNode> dfs,
                                      const std::vector<std::string> sampleLabels,
                                      const std::vector<double> sampleWeights,
                                      float trainFraction) const
{
    // Build FOUR output files, each aggregating events coming
    // from multiple input dataframes:

    if (dfs.size() != sampleLabels.size() || dfs.size() != sampleWeights.size()) {
        throw std::runtime_error("[BDTTrainModule] Mismatch in sizes of input vectors");
    }

    const std::string train_signal_File = "bdt_train_signal.root";
    const std::string train_bkg_File    = "bdt_train_bkg.root";
    const std::string test_signal_File  = "bdt_test_signal.root";
    const std::string test_bkg_File     = "bdt_test_bkg.root";

    // 1) Add a per-sample weight branch to each input dataframe and
    //    split them into signal/background groups.
    std::vector<ROOT::RDF::RNode> sigNodes;
    std::vector<ROOT::RDF::RNode> bkgNodes;
    sigNodes.reserve(dfs.size());
    bkgNodes.reserve(dfs.size());

    for (size_t i = 0; i < dfs.size(); ++i) {
        const std::string &label = sampleLabels[i];
        const bool isSignal = (label.find("signal") != std::string::npos);
        const double weight = sampleWeights[i];

        auto withWeight = dfs[i].Define("sample_weight", [weight]() { return weight; });
        if (isSignal) sigNodes.emplace_back(withWeight);
        else          bkgNodes.emplace_back(withWeight);
    }

    // Split each node individually, snapshot to temp files, then merge per-class.
    std::vector<std::string> tmp_train_sig, tmp_test_sig, tmp_train_bkg, tmp_test_bkg;
    tmp_train_sig.reserve(sigNodes.size());
    tmp_test_sig.reserve(sigNodes.size());
    tmp_train_bkg.reserve(bkgNodes.size());
    tmp_test_bkg.reserve(bkgNodes.size());

    // Keep weights in the snapshots
    std::vector<std::string> colsToKeep = fVarsToKeep;
    if (std::find(colsToKeep.begin(), colsToKeep.end(), std::string("sample_weight")) == colsToKeep.end())
        colsToKeep.push_back("sample_weight");

    auto createTempFiles = [&](const std::vector<ROOT::RDF::RNode>& nodes,
                            std::vector<std::string>& tmpTrain,
                            std::vector<std::string>& tmpTest,
                            const char* tag) {
        for (size_t i = 0; i < nodes.size(); ++i) {
            auto &n = const_cast<ROOT::RDF::RNode&>(nodes[i]);
            const auto nEntries = n.Count().GetValue();
            const auto nTrain  = static_cast<Long64_t>(std::round(nEntries * trainFraction));
            const auto nTest = nEntries - nTrain;

            std::cout << "[BDTTrainModule] Sample " << tag << " " << i
                      << ": total entries = " << nEntries
                      << ", train = " << nTrain
                      << ", test = " << nTest << std::endl;

            auto trainNode = n.Range(0, nTrain);
            auto testNode  = n.Range(nTrain, nEntries);

            ROOT::RDF::RSnapshotOptions opts;
            opts.fMode = "RECREATE"; // one snapshot per tmp file

            const std::string tmpTrainName = std::string("/tmp/") + std::to_string(reinterpret_cast<uintptr_t>(this)) +
                                             std::string("_") + tag + std::string("_train_") + std::to_string(i) + ".root";
            const std::string tmpTestName  = std::string("/tmp/") + std::to_string(reinterpret_cast<uintptr_t>(this)) +
                                             std::string("_") + tag + std::string("_test_")  + std::to_string(i) + ".root";

            trainNode.Snapshot("tree", tmpTrainName, colsToKeep, opts);
            testNode.Snapshot ("tree", tmpTestName,  colsToKeep, opts);

            tmpTrain.emplace_back(tmpTrainName);
            tmpTest.emplace_back(tmpTestName);
        }
    };

    createTempFiles(sigNodes, tmp_train_sig, tmp_test_sig, "sig");
    createTempFiles(bkgNodes, tmp_train_bkg, tmp_test_bkg, "bkg");

    auto mergeFiles = [](const std::vector<std::string>& inputs,
                         const std::string& output) {
        TChain chain("tree");
        for (auto &f : inputs) chain.Add(f.c_str());
        TFile out(output.c_str(), "RECREATE");
        auto *merged = chain.CloneTree(-1, "fast");
        merged->SetName("tree");
        merged->Write();
        out.Close();
    };

    auto cleanupFiles = [](const std::vector<std::string>& files) {
        for (const auto &f : files) {
            if (std::remove(f.c_str()) != 0) {
                // best-effort cleanup; ignore errors but you could log if desired
            }
        }
    };

    if (!tmp_train_sig.empty()) { mergeFiles(tmp_train_sig, train_signal_File); cleanupFiles(tmp_train_sig); }
    if (!tmp_test_sig.empty())  { mergeFiles(tmp_test_sig,  test_signal_File);  cleanupFiles(tmp_test_sig); }
    if (!tmp_train_bkg.empty()) { mergeFiles(tmp_train_bkg, train_bkg_File);    cleanupFiles(tmp_train_bkg); }
    if (!tmp_test_bkg.empty())  { mergeFiles(tmp_test_bkg,  test_bkg_File);     cleanupFiles(tmp_test_bkg); }

    return {train_signal_File, train_bkg_File, test_signal_File, test_bkg_File};
}



Long64_t BDTTrainModule::EntryCount() const
{
    return 1; // Dummy, nothing per-event
}

//------------------------------------------------------------------------------
void BDTTrainModule::Initialise()
{
    auto dfVec = BuildDataFrames(fInputFiles, fTreeName);

    std::vector<ROOT::RDF::RNode> nodes;
    nodes.reserve(dfVec.size());
    for (auto &dfPtr : dfVec) nodes.emplace_back(*dfPtr);

    std::cout << "Training fraction: " << fTrainFraction << std::endl;
    // Create training and testing samples, return the filenames
    std::vector<std::string> sampleVec = BuildTestTrainSamples(nodes, fSampleLabels, fSampleWeights, fTrainFraction);
    std::cout << "[BDTTrainModule] Created training and testing samples:\n";

    std::string train_signal_File = sampleVec[0];
    std::string train_bkg_File = sampleVec[1];
    std::string test_signal_File = sampleVec[2];
    std::string test_bkg_File = sampleVec[3];

    
    
}

//------------------------------------------------------------------------------
void BDTTrainModule::Finalise()
{
    // Nothing to do – Snapshot already wrote the slimmed tree.
}
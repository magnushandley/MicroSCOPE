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

#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Tools.h>

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

    std::stringstream ssTrainVars{cfg.GetValue("BDTTrainModule.TrainVars", "")};
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

std::string BDTTrainModule::BuildMethodString(int nTrees,
                                      int maxDepth,
                                      double learningRate,
                                      double minNodeSize,
                                    int nCuts) const
{
    std::ostringstream method;
    method << "!H:!V:NTrees=" << nTrees
           << ":MinNodeSize=" << minNodeSize << "%"
           << ":MaxDepth=" << maxDepth
           << ":BoostType=Grad"
           << ":Shrinkage=" << learningRate
           << ":nCuts=" << nCuts;
    return method.str();
}
std::string BDTTrainModule::FindOptimalCut(const std::string& trainSignalFile,
                              const std::string& trainBkgFile,
                              const std::string& testSignalFile,
                              const std::string& testBkgFile,
                              const std::vector<int> NTreesVec,
                              const std::vector<int> MaxDepthVec,
                              const std::vector<double> LearningRateVec,
                              const std::vector<double> MinNodeSizeVec,
                            const std::vector<int> nCutsVec)
{
    //Takes a range of possible hyperparameters and finds the optimal set based on test sample performance,
    // returning a string suitable for TMVA::Factory::BookMethod 
    std::string bestMethodString;

    // Loop over all combinations of hyperparameters
    double bestScore = -1.0;
    std::string bestMethod;
    for (const auto& nTrees : NTreesVec) {
        for (const auto& maxDepth : MaxDepthVec) {
            for (const auto& learningRate : LearningRateVec) {
                for (const auto& minNodeSize : MinNodeSizeVec) {
                    for (const auto& nCuts : nCutsVec) {
                        std::string methodString = BuildMethodString(nTrees, maxDepth, learningRate, minNodeSize, nCuts);
                        // Train BDT with these hyperparameters
                        double score = TrainBDT(trainSignalFile, trainBkgFile, testSignalFile, testBkgFile, methodString);
                        // Evaluate performance on test set, e.g. via AUC or significance
                        if (score > bestScore) {
                            bestScore = score;
                            bestMethodString = methodString;
                            std::cout << "[BDTTrainModule] New best score: " << bestScore
                                      << " with method: " << bestMethodString << std::endl;
                        }

                        std::cout << "[BDTTrainModule] Tested nTrees=" << nTrees
                                  << ", maxDepth=" << maxDepth
                                  << ", learningRate=" << learningRate
                                  << ", minNodeSize=" << minNodeSize
                                  << ", nCuts=" << nCuts
                                  << " => score: " << score << std::endl;
                        
                        //Check for value on the boundary of the hyperparameter ranges - may need to expand search space
                        
                        if (nTrees == NTreesVec.front()) {
                            std::cout << "[BDTTrainModule] Best nTrees at lower boundary: " << nTrees << std::endl;
                        }
                        if (nTrees == NTreesVec.back()) {
                            std::cout << "[BDTTrainModule] Best nTrees at upper boundary: " << nTrees << std::endl;
                        }
                        if (maxDepth == MaxDepthVec.front()) {
                            std::cout << "[BDTTrainModule] Best maxDepth at lower boundary: " << maxDepth << std::endl;
                        }
                        if (maxDepth == MaxDepthVec.back()) {
                            std::cout << "[BDTTrainModule] Best maxDepth at upper boundary: " << maxDepth << std::endl;
                        }
                        if (learningRate == LearningRateVec.front()) {
                            std::cout << "[BDTTrainModule] Best learningRate at lower boundary: " << learningRate << std::endl;
                        }
                        if (learningRate == LearningRateVec.back()) {
                            std::cout << "[BDTTrainModule] Best learningRate at upper boundary: " << learningRate << std::endl;
                        }
                        if (minNodeSize == MinNodeSizeVec.front()) {
                            std::cout << "[BDTTrainModule] Best minNodeSize at lower boundary: " << minNodeSize << std::endl;
                        }
                        if (minNodeSize == MinNodeSizeVec.back()) {
                            std::cout << "[BDTTrainModule] Best minNodeSize at upper boundary: " << minNodeSize << std::endl;
                        }
                        if (nCuts == nCutsVec.front()) {
                            std::cout << "[BDTTrainModule] Best nCuts at lower boundary: " << nCuts << std::endl;
                        }
                        if (nCuts == nCutsVec.back()) {
                            std::cout << "[BDTTrainModule] Best nCuts at upper boundary: " << nCuts << std::endl;
                        }   
                    }
                }
            }
        }
    }

    std::cout << "[BDTTrainModule] Optimal method string: " << bestMethodString << std::endl;
    return bestMethodString;

}

double BDTTrainModule::TrainBDT(const std::string& trainSignalFile,
                              const std::string& trainBkgFile,
                              const std::string& testSignalFile,
                              const std::string& testBkgFile,
                              std::string methodString) const
{
    //Trains the BDT and returns a figure of merit (e.g. ROC) on the test set.
    // Build chains so we can control the split with SplitMode=Block
    TChain sigChain("tree");
    TChain bkgChain("tree");

    // Add in the order: TRAIN first, then TEST
    sigChain.Add(trainSignalFile.c_str());
    sigChain.Add(testSignalFile.c_str());

    bkgChain.Add(trainBkgFile.c_str());
    bkgChain.Add(testBkgFile.c_str());

    // Count how many entries belong to the training portion
    Long64_t nTrainSig = 0, nTrainBkg = 0;
    {
        std::unique_ptr<TFile> fSig{TFile::Open(trainSignalFile.c_str())};
        if (!fSig || fSig->IsZombie()) throw std::runtime_error("[BDTTrainModule] Cannot open " + trainSignalFile);
        TTree* tSig = nullptr; fSig->GetObject("tree", tSig);
        if (!tSig) throw std::runtime_error("[BDTTrainModule] 'tree' not found in " + trainSignalFile);
        nTrainSig = tSig->GetEntries();
        for (auto& v : fTrainVars) if (!tSig->GetBranch(v.c_str()))
        std::cerr << "[BDTTrainModule] SIGNAL MISSING branch: " << v << "\n";
    }
    {
        std::unique_ptr<TFile> fBkg{TFile::Open(trainBkgFile.c_str())};
        if (!fBkg || fBkg->IsZombie()) throw std::runtime_error("[BDTTrainModule] Cannot open " + trainBkgFile);
        TTree* tBkg = nullptr; fBkg->GetObject("tree", tBkg);
        if (!tBkg) throw std::runtime_error("[BDTTrainModule] 'tree' not found in " + trainBkgFile);
        nTrainBkg = tBkg->GetEntries();
        for (auto& v : fTrainVars) if (!tBkg->GetBranch(v.c_str()))
        std::cerr << "[BDTTrainModule] BACKGROUND MISSING branch: " << v << "\n";
    }

    // Determine whether a per-event weight exists
    bool hasSampleWeight = false;
    {
        std::unique_ptr<TFile> fProbe{TFile::Open(trainSignalFile.c_str())};
        TTree* tProbe = nullptr; fProbe->GetObject("tree", tProbe);
        if (tProbe && tProbe->GetBranch("sample_weight")) hasSampleWeight = true;
    }

    // TMVA setup
    TMVA::Tools::Instance();
    std::unique_ptr<TFile> outFile{TFile::Open("tmva_training_output.root", "RECREATE")};

    TMVA::Factory factory("TMVAClassification", outFile.get(),
                          "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");
    TMVA::DataLoader loader("dataset");

    // Register training variables from fTrainVars
    std::cout << "[BDTTrainModule] Registering training variables:\n";
    for (const auto& var : fTrainVars) {
        std::cout << "[BDTTrainModule] Adding training variable: " << var << std::endl;
        loader.AddVariable(var.c_str(), 'F');
    }

    // Add the signal and background trees
    loader.AddSignalTree(&sigChain, 1.0);
    loader.AddBackgroundTree(&bkgChain, 1.0);

    // Optional event weights
    if (hasSampleWeight) {
        loader.SetSignalWeightExpression("sample_weight");
        loader.SetBackgroundWeightExpression("sample_weight");
    }

    // Use Block split so the first N entries (our TRAIN files) are used for training
    std::ostringstream prep;
    prep << "nTrain_Signal=" << nTrainSig
         << ":nTrain_Background=" << nTrainBkg
         << ":SplitMode=Block:NormMode=None:!V";

    loader.PrepareTrainingAndTestTree("", "", prep.str().c_str());

    // Book a simple BDTG. You can later move these options into the TEnv cfg.
    //"!H:!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=3:BoostType=Grad:"
                       //"Shrinkage=0.1:nCuts=20"
    factory.BookMethod(&loader, TMVA::Types::kBDT, "BDTG",
                         methodString.c_str());

    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();
    // Retrieve a figure of merit on the test set
    double fom = factory.GetROCIntegral(&loader, "BDTG", /*iClass=*/0, TMVA::Types::kTesting);
    return fom;

    // XML weights will be in: dataset/weights/TMVAClassification_BDTG.weights.xml
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
    
    //"!H:!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=3:BoostType=Grad:"
                       //"Shrinkage=0.1:nCuts=20"
    //Optimise hyperparameters
    std::vector<int> NTreesRange = {150, 200, 250};
    std::vector<int> MaxDepthRange = {2, 3, 4};
    std::vector<double> LearningRateRange = {0.05, 0.1, 1.5};
    std::vector<double> MinNodeSizeRange = {1.5, 2.5, 3.5};
    std::vector<int> nCutsRange = {10, 20, 30};
    std::string methodString = FindOptimalCut(train_signal_File, train_bkg_File, test_signal_File, test_bkg_File,
                                              NTreesRange, MaxDepthRange, LearningRateRange, MinNodeSizeRange, nCutsRange);
    // Train the BDT using the selected variables (fTrainVars) and the prepared samples
    double fom = TrainBDT(train_signal_File, train_bkg_File, test_signal_File, test_bkg_File, methodString);

    std::cout << "[BDTTrainModule] TMVA training complete. Weights XML written under dataset/weights/." << std::endl;
}

//------------------------------------------------------------------------------
void BDTTrainModule::Finalise()
{
    // Nothing to do – Snapshot already wrote the slimmed tree.
}
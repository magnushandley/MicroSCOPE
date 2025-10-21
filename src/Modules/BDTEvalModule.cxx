#include "Modules/BDTEvalModule.hxx"

#include <TMVA/Reader.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TDirectory.h>

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <cstdio>

using namespace Analysis;

//---------------------------------------------
static std::vector<std::string> split_ws_or_commas(const std::string& raw) {
    std::vector<std::string> out;
    std::stringstream ss(raw);
    std::string tok;
    while (ss >> tok) {
        if (!tok.empty() && tok.back()==',') tok.pop_back();
        if (!tok.empty()) out.push_back(tok);
    }
    return out;
}

//---------------------------------------------
BDTEvalModule::BDTEvalModule(const TEnv& cfg)
: Module(cfg)
, fTreeName   (cfg.GetValue("BDTEvalModule.TreeName", "nuselection/NeutrinoSelectionFilter"))
, fWeightsXML (cfg.GetValue("BDTEvalModule.WeightsXML", "dataset/weights/TMVAClassification_BDTG.weights.xml"))
, fMethodName (cfg.GetValue("BDTEvalModule.MethodName", "BDTG"))
, fOutputTag  (cfg.GetValue("BDTEvalModule.OutputTag",  "_bdt"))
{
    // Input files: allow spaces and/or commas
    fInputFiles = split_ws_or_commas(cfg.GetValue("BDTEvalModule.InputFiles", ""));
    for (const auto& file : fInputFiles) {
        std::cout << "[BDTEvalModule] Input file: " << file << "\n";
    }
    // Variables to evaluate (must match training names!)
    fEvalVars   = split_ws_or_commas(cfg.GetValue("BDTEvalModule.EvalVars", ""));

    if (fInputFiles.empty())
        throw std::runtime_error("[BDTEvalModule] No input files provided (BDTEvalModule.InputFiles).");

    if (fEvalVars.empty())
        throw std::runtime_error("[BDTEvalModule] No EvalVariables provided — must match training variables.");
}

//---------------------------------------------
std::vector<std::string> BDTEvalModule::TokeniseCSV(const std::string& s) {
    return split_ws_or_commas(s);
}

//---------------------------------------------
Long64_t BDTEvalModule::EntryCount() const { return 1; }

//---------------------------------------------
void BDTEvalModule::ProcessOneFile(const std::string& inPath) const
{
    // open input
    std::cout << "Loop!" << std::endl;
    std::cout << "[BDTEvalModule] Opening input file: " << inPath << "\n";
    std::unique_ptr<TFile> inFile{TFile::Open(inPath.c_str(), "READ")};
    if (!inFile || inFile->IsZombie())
        throw std::runtime_error("[BDTEvalModule] Cannot open input file: " + inPath);

    // fetch tree (allow for directory path in fTreeName)
    TTree* inTree = nullptr;
    {
        TObject* obj = inFile->Get(fTreeName.c_str());
        if (!obj) {
            throw std::runtime_error("[BDTEvalModule] Cannot find tree: " + fTreeName
                                      + " in file " + inPath);
        }
        inTree = dynamic_cast<TTree*>(obj);
        if (!inTree) {
            throw std::runtime_error("[BDTEvalModule] Object at " + fTreeName + " is not a TTree.");
        }
    }

    // check variables exist and cache TLeaf*
    std::vector<TLeaf*> leaves;
    leaves.reserve(fEvalVars.size());
    for (const auto& v : fEvalVars) {
        // Try both leaf and (branch->leaf) name resolution
        TLeaf* leaf = inTree->GetLeaf(v.c_str());
        if (!leaf) {
            // Some TTrees store leaves as "branch.leaf". Try to find by scanning all leaves.
            bool found = false;
            TObjArray* leafs = inTree->GetListOfLeaves();
            for (int i = 0; i < leafs->GetEntries(); ++i) {
                auto* L = static_cast<TLeaf*>(leafs->At(i));
                if (std::string(L->GetName()) == v) { leaf = L; found = true; break; }
            }
            if (!found) {
                std::ostringstream msg;
                msg << "[BDTEvalModule] Missing variable leaf '" << v
                    << "' in tree '" << fTreeName << "' (file " << inPath << ").";
                throw std::runtime_error(msg.str());
            }
        }
        leaves.push_back(leaf);
    }

    // prepare TMVA::Reader with float buffers (kept in scope)
    TMVA::Reader reader("!Color:!Silent");
    std::vector<float> varBuf(fEvalVars.size(), 0.f);
    for (size_t i = 0; i < fEvalVars.size(); ++i) {
        reader.AddVariable(fEvalVars[i].c_str(), &varBuf[i]);
    }

    // book MVA
    reader.BookMVA(fMethodName.c_str(), fWeightsXML.c_str());

    // output file name
    std::string outPath = inPath;
    auto pos = outPath.find_last_of('.');
    if (pos == std::string::npos) outPath += fOutputTag + ".root";
    else                          outPath.insert(pos, fOutputTag); // e.g. input.root -> input_bdt.root

    std::unique_ptr<TFile> outFile{TFile::Open(outPath.c_str(), "RECREATE")};
    if (!outFile || outFile->IsZombie())
        throw std::runtime_error("[BDTEvalModule] Cannot create output file: " + outPath);

    // clone tree structure and add bdt_score branch
    outFile->cd();
    std::unique_ptr<TTree> outTree{inTree->CloneTree(0)};
    float bdt_score = -999.f;
    TBranch* brScore = outTree->Branch("bdt_score", &bdt_score, "bdt_score/F");

    // main loop
    const Long64_t nEntries = inTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        inTree->GetEntry(i);

        // populate var buffers from leaves (ROOT returns double; cast to float)
        for (size_t j = 0; j < leaves.size(); ++j) {
            // GetValue(0) handles both scalar and first element of arrays; if your
            // variables are arrays, change to GetValue(index) as needed.
            varBuf[j] = static_cast<float>(leaves[j]->GetValue(0));
        }

        bdt_score = reader.EvaluateMVA(fMethodName.c_str());
        brScore->Fill();
        outTree->Fill(); // copies all original branches + our new one
    }

    outTree->Write("", TObject::kOverwrite);
    // Destroy the tree before closing the file; otherwise ROOT may segfault
    // when the tree tries to detach from a closed TDirectory/TFile.
    outTree.reset(nullptr);

    outFile->Write();
    outFile->Close();

    std::cout << "[BDTEvalModule] Wrote: " << outPath
              << "  (entries: " << nEntries << ")\n";
}

//---------------------------------------------
void BDTEvalModule::Initialise() {
    std::cout << "[BDTEvalModule] Using weights xml: " << fWeightsXML << "\n";
    std::cout << "[BDTEvalModule] Method: " << fMethodName << "\n";
    std::cout << "[BDTEvalModule] Tree: " << fTreeName << "\n";
    std::cout << "[BDTEvalModule] Variables (" << fEvalVars.size() << "): ";
    for (auto& v : fEvalVars) std::cout << v << " ";
    std::cout << "\n";

    for (const auto& f : fInputFiles) {
        std::cout << "[BDTEvalModule] Will loop over input file: " << f << "\n";
    }
    for (const auto& f : fInputFiles) {
        std::cout << "[BDTEvalModule] Processing: " << f << "\n";
        ProcessOneFile(f);
    }
    std::cout << "[BDTEvalModule] Done.\n";
}

//---------------------------------------------
void BDTEvalModule::Finalise() {
    // Nothing to do
}
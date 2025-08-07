#include "Modules/SlimmerModule.hxx"
#include "Utils/Plotter.hxx"

#include <TEnv.h>
#include <TFile.h>
#include <TString.h>
#include <algorithm>
#include <sstream>

using namespace Analysis;

//------------------------------------------------------------------------------
SlimmerModule::SlimmerModule(const TEnv& cfg)
    : Module(cfg)
    , fInputFiles   (cfg.GetValue("Slimmer.Input",   ""        ))
    , fTreeName     (cfg.GetValue("Slimmer.TreeName","nuselection/NeutrinoSelectionFilter"  ))
    , fOutFile      (cfg.GetValue("Slimmer.Output", "slimmed.root"))
    , fRunLabel     (cfg.GetValue("Global.RunLabel","run_x") )
{
    
    std::stringstream ss{cfg.GetValue("Slimmer.Keep", "")};
    std::string item;
    while (ss >> item) {
        if (item.back()==',') item.pop_back();
        fVarsToKeep.push_back(item); 
    }
}

//------------------------------------------------------------------------------
std::unique_ptr<TChain>
SlimmerModule::BuildInputChain(const std::string& files,
                               const std::string& treeName) const
{
    // The functionality in this module is different from the others in that here we chain several input files together
    // into a single TChain, which is then passed to RDataFrame. This is because we want to create a single output file
    // with a single tree, rather than one output file per input file. Other modules (e.g. Preselection) create one output file per input file,
    // internally using a vector of RDFs.

    auto chain = std::make_unique<TChain>(treeName.c_str());
    std::stringstream ss{files};
    std::string fname;
    while (std::getline(ss, fname, ',')) {
        if (!fname.empty()) chain->Add(fname.c_str());
    }
    if (chain->GetEntries()==0)
        throw std::runtime_error("[Slimmer] Input chain is empty!");
    return chain;
}

Long64_t SlimmerModule::EntryCount() const
{
    if (!fChain) {
        throw std::runtime_error("[Slimmer] Chain not initialized!");
    }

    return fChain->GetEntries();
}

//------------------------------------------------------------------------------
void SlimmerModule::Initialise()
{
    std::cout << "[Slimmer] Initialising with input files: " << fInputFiles
              << " and output file: " << fOutFile << "\n";
    fChain = BuildInputChain(fInputFiles, fTreeName);

    // RDataFrame takes ownership of the TChain pointer
    fRDF = std::make_unique<ROOT::RDataFrame>(*fChain);

    //----------------------------------------------------------------------
    // 1.  Define derived variables
    //----------------------------------------------------------------------
    auto& df = *fRDF;

    // Fiducial variables to assess containment (taken from HNL analysis). The whole mess with big and small
    // values is because in ext files the trk_sce_start_x_v vectors can be empty if there is no neutrino slice.
    // In overlay this doesn't happen, but need to for data/ext files I think, so I set min/max to values outside the fiducial volume.

    using VecF = const std::vector<float>&;

    auto df1 = df
    .Define("min_x",
        [](VecF a, VecF b) {
            const float small = -9999;
            float aMin = a.empty() ? small : *std::min_element(a.begin(), a.end());
            float bMin = b.empty() ? small : *std::min_element(b.begin(), b.end());
            return std::min(aMin, bMin);
        },
        {"trk_sce_start_x_v", "trk_sce_end_x_v"})
    .Define("max_x",
            [](VecF a, VecF b) {
                const float big = 9999;
                float aMax = a.empty() ? big : *std::max_element(a.begin(), a.end());
                float bMax = b.empty() ? big : *std::max_element(b.begin(), b.end());
                return std::max(aMax, bMax);
            },
            {"trk_sce_start_x_v", "trk_sce_end_x_v"})
    .Define("min_y",
            [](VecF a, VecF b) {
                const float small = -9999;
                float aMin = a.empty() ? small : *std::min_element(a.begin(), a.end());
                float bMin = b.empty() ? small : *std::min_element(b.begin(), b.end());
                return std::min(aMin, bMin);
            },
            {"trk_sce_start_y_v", "trk_sce_end_y_v"})
    .Define("max_y",
            [](VecF a, VecF b) {
                const float big = 9999;
                float aMax = a.empty() ? big : *std::max_element(a.begin(), a.end());
                float bMax = b.empty() ? big : *std::max_element(b.begin(), b.end());
                return std::max(aMax, bMax);
            },
            {"trk_sce_start_y_v", "trk_sce_end_y_v"})
    .Define("min_z",
            [](VecF a, VecF b) {
                const float small = -9999;
                float aMin = a.empty() ? small : *std::min_element(a.begin(), a.end());
                float bMin = b.empty() ? small : *std::min_element(b.begin(), b.end());
                return std::min(aMin, bMin);
            },
            {"trk_sce_start_z_v", "trk_sce_end_z_v"})
    .Define("max_z",
            [](VecF a, VecF b) {
                const float big = 9999;
                float aMax = a.empty() ? big : *std::max_element(a.begin(), a.end());
                float bMax = b.empty() ? big : *std::max_element(b.begin(), b.end());
                return std::max(aMax, bMax);
            },
            {"trk_sce_start_z_v", "trk_sce_end_z_v"});

    ROOT::RDF::RSnapshotOptions opt;
    opt.fMode = "RECREATE";
    opt.fCompressionAlgorithm = ROOT::kZLIB;
    opt.fCompressionLevel     = 4;

    df1.Snapshot(fTreeName, fOutFile, fVarsToKeep, opt);

    Plotter::SaveHist(
        df1.Histo1D({"sub_hist", ";run_number;Count", 50, 0, 600}, "sub").GetPtr(),
        "slimmer_"+fRunLabel+"_run_histogram" , "prelim");
}

//------------------------------------------------------------------------------
void SlimmerModule::Finalise()
{
    // Nothing to do – Snapshot already wrote the slimmed tree.
}
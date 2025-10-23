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
    , fTreeName     (cfg.GetValue("Slimmer.TreeName","nuselection/NeutrinoSelectionFilter"  ))
    , fRunLabel     (cfg.GetValue("Global.RunLabel","run_x") )
{

    std::stringstream ssInput{cfg.GetValue("Slimmer.InputFiles", "")};
    std::string inputItem;
    while (ssInput >> inputItem) {
        if (inputItem.back()==',') inputItem.pop_back();
        fInputFiles.push_back(inputItem);
    }

    std::stringstream ssOutput{cfg.GetValue("Slimmer.OutputFiles", "")};
    std::string outputItem;
    while (ssOutput >> outputItem) {
        if (outputItem.back()==',') outputItem.pop_back();
        fOutputFiles.push_back(outputItem);
    }
    
    std::stringstream ss{cfg.GetValue("Slimmer.Keep", "")};
    std::string item;
    while (ss >> item) {
        if (item.back()==',') item.pop_back();
        fVarsToKeep.push_back(item); 
    }
}

std::vector<std::unique_ptr<ROOT::RDataFrame>>
SlimmerModule::BuildDataFrames(const std::vector<std::string>& files,
                                     const std::string& treeName) const
{
    // Create a vector of unique pointers to RDataFrames, one per input file.

    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfVec;
    for (const auto& fname : files) {
        std::cout << "Creating RDF for file: " << fname << std::endl;
        auto file = TFile::Open(fname.c_str());
        if (!file || file->IsZombie()) {
            throw std::runtime_error("[Slimmer] Cannot open file: " + fname);
        }
        auto tree = file->Get<TTree>(treeName.c_str());
        if (!tree) {
            throw std::runtime_error("[Slimmer] Cannot find tree: " + treeName);
        }
        auto RDF = std::make_unique<ROOT::RDataFrame>(*tree);
        dfVec.push_back(std::move(RDF));
    }
    if (dfVec.empty())
        throw std::runtime_error("[Slimmer] No data frames created!");

    return dfVec;
}

Long64_t SlimmerModule::EntryCount() const
{
    
    if (dfVec.size() == 0) {
        throw std::runtime_error("[Slimmer] DataFrames not initialised!");
    }

    Long64_t totalEntries = 0;
    for (const auto& df : dfVec) {
        totalEntries += df->Count().GetValue();
    }
    return totalEntries;
}

//------------------------------------------------------------------------------
void SlimmerModule::Initialise()
{
    //std::cout << "[Slimmer] Initialising with input files: " << fInputFiles
              //<< " and output file: " << fOutFile << "\n";
    //fChain = BuildInputChain(fInputFiles, fTreeName);

    std::cout << "[Slimmer] Initialising with input files: \n";
    for (const auto& file : fInputFiles) {
        std::cout << "  " << file << "\n";
    }
    dfVec = BuildDataFrames(fInputFiles, fTreeName);

    //Silly games to convert vector of unique_ptr<RDataFrame> to vector of RNodes
    std::vector<ROOT::RDF::RNode> nodes;
    nodes.reserve(dfVec.size());
    for (auto &dfPtr : dfVec) nodes.emplace_back(*dfPtr);

    // RDataFrame takes ownership of the TChain pointer
    //fRDF = std::make_unique<ROOT::RDataFrame>(*fChain);

    //----------------------------------------------------------------------
    // 1.  Define derived variables
    //----------------------------------------------------------------------
    int fileIndex = 0;
    for (auto df : nodes) {
        std::cout << "[Slimmer] Number of entries in input file: " << df.Count().GetValue() << '\n';
        std::string fOutFile = fOutputFiles[fileIndex];
        std::cout << "[Slimmer] Will write slimmed tree to: " << fOutFile << '\n';


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
                {"trk_sce_start_z_v", "trk_sce_end_z_v"})
        .Filter("swtrig==1"); // keep only events passing the software trigger

        ROOT::RDF::RSnapshotOptions opt;
        opt.fMode = "RECREATE";
        opt.fCompressionAlgorithm = ROOT::kZLIB;
        opt.fCompressionLevel     = 4;

        df1.Snapshot(fTreeName, fOutFile, fVarsToKeep, opt);

        Plotter::SaveHist(
            df1.Histo1D({"sub_hist", ";run_number;Count", 50, 0, 600}, "sub").GetPtr(),
            "slimmer_"+fRunLabel+"_run_histogram" , "prelim");

        fileIndex++;
    }
}

//------------------------------------------------------------------------------
void SlimmerModule::Finalise()
{
    // Nothing to do – Snapshot already wrote the slimmed tree.
}
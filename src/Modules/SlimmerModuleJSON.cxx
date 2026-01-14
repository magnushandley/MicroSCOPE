#include "Modules/SlimmerModuleJSON.hxx"
#include "Utils/Plotter.hxx"

#include <TFile.h>
#include <TString.h>
#include <algorithm>
#include <sstream>

using namespace Analysis;

//------------------------------------------------------------------------------
SlimmerModuleJSON::SlimmerModuleJSON(const nlohmann::json& cfgJson) : Module(cfgJson) 
{
    // Nothing to do here, all in Initialise()
}

//------------------------------------------------------------------------------
Long64_t SlimmerModuleJSON::EntryCount() const
{
    bool found_vec = false;
    for (const auto& pair : fRdfMap) {
        if (pair.second.size() > 0) {
            found_vec = true;
            break;
        }
    }
    if (!found_vec) {
        throw std::runtime_error("[Slimmer] DataFrames not initialised!");
    }

    Long64_t totalEntries = 0;
    for (const auto& pair : fRdfMap) {
        for (const auto& df : pair.second) {
            totalEntries += df->Count().GetValue();
        }
    }
    return totalEntries;
}

//------------------------------------------------------------------------------
void SlimmerModuleJSON::Initialise()
{

    std::cout << "[Slimmer] Reading configuration from JSON file.\n";
    // Load configuration from JSON
    auto globalSettings = fCfgJson["GlobalSettings"];
    auto variables = fCfgJson["Variables"];

    std::cout << "[Slimmer] General categories read.\n";

    auto tempFiles = globalSettings["InputFiles"];
    for (const auto& file : tempFiles) {
        fInputFiles.push_back(file.get<std::string>());
    }
    auto recoAlgs = globalSettings["RecoAlgs"];
    fOutputFiles = globalSettings["OutputFile"];
    fRunLabel = globalSettings["RunLabel"];


    std::cout << "[Slimmer] Initialising with input files: \n";
    for (const auto& file : fInputFiles) {
        std::cout << "  " << file << "\n";
    }

    ROOT::RDF::RSnapshotOptions opt;
    opt.fMode = "UPDATE";
    opt.fCompressionAlgorithm = ROOT::kZLIB;
    opt.fCompressionLevel = 4;

    //Silly games to convert vector of unique_ptr<RDataFrame> to vector of RNodes
    std::map<std::string, std::vector<ROOT::RDF::RNode>> nodesMap;

    //Prepare dataframes and list of variables to keep
    for (const auto& alg : recoAlgs) {
        fVarsToKeep.clear();
        if (variables.contains(alg)){
            for (const auto& varCat : variables[alg]) {
                for (const auto& var : varCat) {
                    std::cout << "[Slimmer] Will keep variable: " << var << "\n";
                    fVarsToKeep.push_back(var.get<std::string>());
                }
            }
        }
        if(alg == "pandora")
            fVarsToKeep.insert(fVarsToKeep.end(), { "min_x", "max_x", "min_y", "max_y", "min_z", "max_z"});;



        fRdfMap[alg] = BuildDataFrames(fInputFiles, treeMap[alg]);

        //Silly games to convert vector of unique_ptr<RDataFrame> to vector of RNodes
        std::vector<ROOT::RDF::RNode> nodes;
        nodes.reserve(fRdfMap[alg].size());
        for (auto &dfPtr : fRdfMap[alg]) nodes.emplace_back(*dfPtr);

        int fileIndex = 0;
        for (auto df : nodes) {
            std::cout << "[Slimmer] Number of entries in input file: " << df.Count().GetValue() << '\n';
            std::string fOutFile = fOutputFiles[fileIndex];
            std::cout << "[Slimmer] Will write slimmed tree to: " << fOutFile << '\n';
            if(alg == "pandora"){
                std::cout << "[Slimmer] Define variables for pandora tree \n";
                using VecF = const std::vector<float>&;
                auto df1 = df.Define("min_x", [](VecF a, VecF b) {
                            const float small = -9999;
                            float aMin = a.empty() ? small : *std::min_element(a.begin(), a.end());
                            float bMin = b.empty() ? small : *std::min_element(b.begin(), b.end());
                            return std::min(aMin, bMin);
                        }, {"trk_sce_start_x_v", "trk_sce_end_x_v"})
                        .Define("max_x", [](VecF a, VecF b) {
                            const float big = 9999;
                            float aMax = a.empty() ? big : *std::max_element(a.begin(), a.end());
                            float bMax = b.empty() ? big : *std::max_element(b.begin(), b.end());
                            return std::max(aMax, bMax);
                        }, {"trk_sce_start_x_v", "trk_sce_end_x_v"})
                        .Define("min_y",
                            [](VecF a, VecF b) {
                                const float small = -9999;
                                float aMin = a.empty() ? small : *std::min_element(a.begin(), a.end());
                                float bMin = b.empty() ? small : *std::min_element(b.begin(), b.end());
                                return std::min(aMin, bMin);
                            },{"trk_sce_start_y_v", "trk_sce_end_y_v"})
                        .Define("max_y",
                            [](VecF a, VecF b) {
                                const float big = 9999;
                                float aMax = a.empty() ? big : *std::max_element(a.begin(), a.end());
                                float bMax = b.empty() ? big : *std::max_element(b.begin(), b.end());
                                return std::max(aMax, bMax);
                            }, {"trk_sce_start_y_v", "trk_sce_end_y_v"})
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

                df1.Snapshot(treeMap[alg], fOutFile, fVarsToKeep, opt);
            }
            else df.Snapshot(treeMap[alg], fOutFile, fVarsToKeep, opt);
            ++fileIndex;
            
        }

    }    
}


//------------------------------------------------------------------------------
std::vector<std::unique_ptr<ROOT::RDataFrame>>
SlimmerModuleJSON::BuildDataFrames(const std::vector<std::string>& files,
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
        else {
            std::cout << "[Slimmer] Found tree: " << treeName << " with "
                      << tree->GetEntries() << " entries" << std::endl;
        }   
        auto RDF = std::make_unique<ROOT::RDataFrame>(*tree);
        dfVec.push_back(std::move(RDF));
    }
    if (dfVec.empty())
        throw std::runtime_error("[Slimmer] No data frames created!");
    else
        std::cout << "[Slimmer] Created " << dfVec.size() << " data frames." << std::endl;

    return dfVec;
}
//------------------------------------------------------------------------------
void SlimmerModuleJSON::Finalise()
{
    // Nothing to do – Snapshot already wrote the slimmed tree.
}
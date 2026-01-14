#include "Modules/PreselectionModuleJSON.hxx"
#include "Utils/Plotter.hxx"

#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <algorithm>
#include <sstream>

using namespace Analysis;

//------------------------------------------------------------------------------
PreselectionModuleJSON::PreselectionModuleJSON(const nlohmann::json& cfgJson)
: Module(cfgJson) {}

//------------------------------------------------------------------------------
std::vector<std::unique_ptr<ROOT::RDataFrame>>
PreselectionModuleJSON::BuildDataFrames(const std::vector<std::string>& files,
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

//------------------------------------------------------------------------------
Long64_t PreselectionModuleJSON::EntryCount() const
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

//------------------------------------------------------------------------------
void PreselectionModuleJSON::Initialise()
{
    ROOT::RDF::RSnapshotOptions opt;
    opt.fMode = "UPDATE";
    opt.fCompressionAlgorithm = ROOT::kZLIB;
    opt.fCompressionLevel = 4;

    // Load configuration from JSON
    std::cout << "[Preselection] Start\n";
    auto globalSettings = fCfgJson["GlobalSettings"];
    auto allCuts = fCfgJson["Cuts"];
   
    auto tempFiles = globalSettings["InputFiles"];
    for (const auto& file : tempFiles) {
        fInputFiles.push_back(file.get<std::string>());
    }
    tempFiles = globalSettings["OutputFiles"];
    for (const auto& file : tempFiles) {
        fOutFiles.push_back(file.get<std::string>());
    }

    std::cout << "[Preselection] General Categories read\n";

    auto recoAlgs = globalSettings["RecoAlgs"];
    fSampleLabels = globalSettings["SampleLabels"];

    std::cout << "[Preselection] Start algrorithm loop\n";

    for (const auto& alg : recoAlgs) {
        cuts.clear();

        fRdfMap[alg] = BuildDataFrames(fInputFiles, treeMap[alg]);

        auto tempCuts = allCuts[alg];
        for (const auto& cut : tempCuts) {
            cuts.push_back(cut.get<std::string>());
        }

        std::vector<ROOT::RDF::RNode> nodes;
        nodes.reserve(fRdfMap[alg].size());
        for (auto &dfPtr : fRdfMap[alg]) nodes.emplace_back(*dfPtr);

        int fileIndex = 0;
        for (auto df : nodes) {
            std::cout << "[Preselection] Number of entries in input file: " << df.Count().GetValue() << '\n';
            std::string fOutFile = fOutFiles[fileIndex];
            std::cout << "[Preselection] Will write pre-selected tree to: " << fOutFile << '\n';

            for (const auto &cut : cuts) {
                std::cout << "\n[Preselection] Cut: " << cut << '\n';
                df = df.Filter(cut, cut);   
            }


            //Print cut reports
            auto cutReport = df.Report();
            cutReport->Print();  
            std::cout << "\n[Preselection] Writing output for sample: " << fSampleLabels[fileIndex] << '\n';
            std::cout << "    to file: " << fOutFile << '\n';
            df.Snapshot(treeMap[alg], fOutFile, ".*", opt);
            ++fileIndex;
           
        }  
    }
}

//------------------------------------------------------------------------------
void PreselectionModuleJSON::Finalise()
{
    // Nothing to do here
}
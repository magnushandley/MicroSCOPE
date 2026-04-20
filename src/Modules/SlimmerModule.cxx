#include "Modules/SlimmerModule.hxx"
#include "Utils/Plotter.hxx"
#include "Utils/TimingUtils.hxx"
#include "Utils/ConfigUtils.hxx"
#include "Utils/LogicRegistry.hxx"

#include <TEnv.h>
#include <TFile.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <algorithm>
#include <sstream>
#include <TRandom.h>
#include <unordered_map>
#include <cmath>

using namespace Analysis;

//------------------------------------------------------------------------------
SlimmerModule::SlimmerModule(const TEnv& cfg)
    : Module(cfg)
    , fTreeName     (cfg.GetValue("Slimmer.TreeName","nuselection/NeutrinoSelectionFilter"))
    , fRunLabel     (cfg.GetValue("Global.RunLabel","run_x") )
    , fMakePlots   (cfg.GetValue("Slimmer.MakePlots", false))
    , fConfig(const_cast<TEnv&>(cfg)) // Store a reference to the config for use in logic operations
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

    std::stringstream ssLabels{cfg.GetValue("Slimmer.SampleLabels", "")};
    std::string label;
    while (ssLabels >> label) {
        if (label.back()==',') label.pop_back();
        fSampleLabels.push_back(label);
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

    std::cout << "[Slimmer] Initialising with input files: \n";
    for (const auto& file : fInputFiles) {
        std::cout << "  " << file << "\n";
    }
    dfVec = BuildDataFrames(fInputFiles, fTreeName);

    //Silly games to convert vector of unique_ptr<RDataFrame> to vector of RNodes
    std::vector<ROOT::RDF::RNode> nodes;
    nodes.reserve(dfVec.size());
    for (auto &dfPtr : dfVec) nodes.emplace_back(*dfPtr);

    //----------------------------------------------------------------------
    // 1.  Define derived variables
    //----------------------------------------------------------------------

    int fileIndex = 0;
    for (auto df : nodes) {
        const double beamSpillPeriod = fBeamSpillPeriod;
        std::cout << "[Slimmer] Number of entries in input file: " << df.Count().GetValue() << '\n';
        std::string fOutFile = fOutputFiles[fileIndex];
        std::cout << "[Slimmer] Will write slimmed tree to: " << fOutFile << '\n';

        auto logicConfigs = ParseLogicConfigs(fConfig, "Slimmer");
        std::cout << "[Slimmer] Parsed " << logicConfigs.size() << " logic operations from config.\n";

        ROOT::RDF::RNode df_int = df; // Intermediate RNode to apply logic operations to, before snapshotting

        //Create derived varibles based on config.
        //for (const auto& op : logicConfigs) {
        //    df_int = LogicRegistry::Instance().Apply(df_int, op);
        //}

        //For the tutorial, apply hard coded operations - in practice, these should be defined in the config and applied via the registry as above.
        //-----------------------------------------------------------------------
        //Add logic here:
        // e.g. df_int = df_int.Define("new_variable", "expression or lambda function");
        //

        ROOT::RDF::RNode dfOut = df_int;

        // Could apply file-specific logic separate to global filters and definitions.
        

        //----------------------------------------------------------------------
        // 2.  Snapshot only the variables we want to keep
        //----------------------------------------------------------------------
        
        ROOT::RDF::RSnapshotOptions opt;
        opt.fMode = "RECREATE";
        opt.fCompressionAlgorithm = ROOT::kZLIB;
        opt.fCompressionLevel     = 4;

        std::cout << "[Slimmer] Writing slimmed tree to file: " << fOutFile << '\n';
        std::cout << "[Slimmer] Variables to keep in slimmed tree:\n";
        
        for (const auto& var : fVarsToKeep) {
            std::cout << "  " << var << "\n";
        }

        dfOut.Snapshot(fTreeName, fOutFile, fVarsToKeep, opt);

        if (fMakePlots) {
            std::cout << "Creating plots for file: " << fOutFile << std::endl;
            // Could add some quick diagnostic plots here.
        }

        fileIndex++;
    }
}

//------------------------------------------------------------------------------
void SlimmerModule::Finalise()
{
    // Nothing to do – Snapshot already wrote the slimmed tree.
}

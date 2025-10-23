#include "Modules/PreselectionModule.hxx"
#include "Utils/Plotter.hxx"

#include <TEnv.h>
#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <algorithm>
#include <sstream>

using namespace Analysis;

//------------------------------------------------------------------------------
PreselectionModule::PreselectionModule(const TEnv& cfg)
    : Module(cfg)
    , fTreeName     (cfg.GetValue("Preselection.TreeName","nuselection/NeutrinoSelectionFilter"  ))
    , fRunLabel     (cfg.GetValue("Global.RunLabel","run_x") )
{
    
        // --------------------------------------------------------------------
    // Parse the comma‑separated list of cuts, this is currently messy but was needed to ensure that the proper format of 
    // cutstring was passed to RDF::Filter
    // --------------------------------------------------------------------

    const std::string cutsString = cfg.GetValue("Preselection.Cuts", "");
    std::stringstream ssCuts{cutsString};
    std::string cutToken;

    while (std::getline(ssCuts, cutToken, ',')) {
        // trim leading/trailing whitespace (but keep quotes if present)
        cutToken.erase(0, cutToken.find_first_not_of(" \t\n\r"));
        cutToken.erase(cutToken.find_last_not_of(" \t\n\r") + 1);

        // Remove optional surrounding quotes so the expression is valid for RDataFrame::Filter
        if (cutToken.size() >= 2 &&
           ((cutToken.front() == '"'  && cutToken.back() == '"') ||
            (cutToken.front() == '\'' && cutToken.back() == '\'')))
        {
            cutToken = cutToken.substr(1, cutToken.size() - 2);
        }

        if (!cutToken.empty())
            cuts.push_back(cutToken);

    }
    
    if (cuts.empty()) {
        throw std::runtime_error("[Preselection] No cuts specified!");
    }

    std::stringstream ssKeep{cfg.GetValue("Preselection.Keep", "")};
    std::string keepItem;
    while (ssKeep >> keepItem) {
        if (keepItem.back()==',') keepItem.pop_back();
        fVarsToKeep.push_back(keepItem); 
    }

    std::stringstream ssInput{cfg.GetValue("Preselection.InputFiles", "")};
    std::string inputItem;
    while (ssInput >> inputItem) {
        if (inputItem.back()==',') inputItem.pop_back();
        fInputFiles.push_back(inputItem);
    }

    std::stringstream ssOut{cfg.GetValue("Preselection.Outputs", "")};
    std::string outItem;
    while (ssOut >> outItem) {
        if (outItem.back()==',') outItem.pop_back();
        fOutFiles.push_back(outItem);
    }

    std::stringstream ssLabels{cfg.GetValue("Preselection.SampleLabels", "")};
    std::string label;
    while (ssLabels >> label) {
        if (label.back()==',') label.pop_back();
        fSampleLabels.push_back(label);
    }

    std::stringstream ssWeights{cfg.GetValue("Preselection.SampleWeights", "")};
    double weight;
    while (ssWeights >> weight) {
        fSampleWeights.push_back(weight);
    }

    if (fVarsToKeep.empty()) {
        throw std::runtime_error("[Preselection] No variables to keep specified!");
    }
}

std::vector<std::unique_ptr<ROOT::RDataFrame>>
PreselectionModule::BuildDataFrames(const std::vector<std::string>& files,
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

Long64_t PreselectionModule::EntryCount() const
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

void PreselectionModule::Initialise()
{
    dfVec = BuildDataFrames(fInputFiles, fTreeName);

     // One RNode per sample, initially pointing at the un-filtered DataFrame
     // This again feels messy because we initialise a set of pointers to rdataframes, but then RDF::Filter returns RNodes
     // and we have to keep track of those instead.
     
    std::vector<ROOT::RDF::RNode> nodes;
    nodes.reserve(dfVec.size());
    for (auto &dfPtr : dfVec) nodes.emplace_back(*dfPtr);

    // Apply every cut in sequence
    for (const auto &cut : cuts) {
        std::cout << "\n[Preselection] Cut: " << cut << '\n';

        for (std::size_t i = 0; i < nodes.size(); ++i) {
            auto before = nodes[i].Count().GetValue();
            std::cout << "    " << fSampleLabels[i] << " before: " << before << '\n';

            nodes[i] = nodes[i].Filter(cut);      
            auto after  = nodes[i].Count().GetValue();
            std::cout << "    " << fSampleLabels[i] << " after : " << after  << '\n';
        }
    }

    ROOT::RDF::RSnapshotOptions opt;
    opt.fMode = "RECREATE";
    opt.fCompressionAlgorithm = ROOT::kZLIB;
    opt.fCompressionLevel     = 4;

    // Now we can write each filtered RNode to a new TTree in the output file

    std::vector<TH1D> preSelectednpfpsVec;
    std::vector<TH1D> preSelectedNuE2Vec;
    std::vector<TH1D> preSelectedFlashMatchScoreVec;
    std::vector<TH1D> preSelectedTopologicalScoreVec;
    std::vector<TH1D> preSelectedShrPhivVec;
    std::vector<TH1D> preSelectedShrFitPzFracVec;
    std::vector<TH1D> preSelectedShrFitThetaVec;

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        std::cout << "\n[Preselection] Writing output for sample: " << fSampleLabels[i] << '\n';
        std::cout << "    to file: " << fOutFiles[i] << '\n';
        nodes[i].Snapshot(fTreeName, fOutFiles[i], fVarsToKeep, opt);


        //ROOT::RDF::TH1DModel NuE2Model(
        //    ("pre_hist_" + fSampleLabels[i]).c_str(),
        //    ";Neutrino Energy [MeV];Count",
        //    20, 0.0, 500.0);

        //ROOT::RDF::TH1DModel FlashMatchScoreModel(
        //    ("pre_hist_" + fSampleLabels[i]).c_str(),
        //    ";Flash Match Score;Count",
        //    20, 0.0, 15.0);

        //ROOT::RDF::TH1DModel TopologicalScoreModel(
        //    ("pre_hist_" + fSampleLabels[i]).c_str(),
        //    ";Topological Score;Count",
        //    30, 0.0, 1.0);

        //preSelectedHistVec.push_back(
            //nodes[i].Histo1D(enModel, "NeutrinoEnergy2").GetPtr());
        //TH1D histNuE2 = nodes[i].Histo1D(NuE2Model, "NeutrinoEnergy2").GetValue();
        //histNuE2.SetName(("preselection_hist_" + fSampleLabels[i]).c_str());
        //preSelectedNuE2Vec.push_back(histNuE2);
        //TH1D histFlashMatchScore = nodes[i].Histo1D(FlashMatchScoreModel, "nu_flashmatch_score").GetValue();
        //histFlashMatchScore.SetName(("preselection_hist_" + fSampleLabels[i]).c_str());
        //preSelectedFlashMatchScoreVec.push_back(histFlashMatchScore);
        //TH1D histTopoScore = nodes[i].Histo1D(TopologicalScoreModel, "topological_score").GetValue();

        preSelectednpfpsVec.push_back(
            Plotter::CreateTH1DFromRNode(
                nodes[i],
                ("preselection_hist_npfps_" + fSampleLabels[i]).c_str(),
                "n_pfps", 
                "Number of PFParticles",
                "Count",
                5, 0.5, 5.5));

        preSelectedNuE2Vec.push_back(
            Plotter::CreateTH1DFromRNode(
                nodes[i],
                ("preselection_hist_NeutrinoEnergy2_" + fSampleLabels[i]).c_str(),
                "NeutrinoEnergy2",
                "Neutrino Energy [MeV]",
                "Count",
                20, 0.0, 500.0));

        preSelectedFlashMatchScoreVec.push_back(
            Plotter::CreateTH1DFromRNode(
                nodes[i],
                ("preselection_hist_FlashMatchScore_" + fSampleLabels[i]).c_str(),
                "nu_flashmatch_score", 
                "Flash Match Score",
                "Count",
                20, 0.0, 15.0));

        preSelectedTopologicalScoreVec.push_back(
            Plotter::CreateTH1DFromRNode(
                nodes[i],
                ("preselection_hist_TopologicalScore_" + fSampleLabels[i]).c_str(),
                "topological_score", 
                "Topological Score",
                "Count",
                30, 0.0, 1.0));

        preSelectedShrPhivVec.push_back(
            Plotter::CreateTH1DFromRNode(
                nodes[i],
                ("preselection_hist_ShrPhiv_" + fSampleLabels[i]).c_str(),
                "shr_phi_v", 
                "Shr Phi [rad]",
                "Count",
                20, -3.14, 3.14,
                true)); // remove vector duplicates by taking first element only

        preSelectedShrFitPzFracVec.push_back(
            Plotter::CreateTH1DFromRNode(
                nodes[i],
                ("preselection_hist_ShrFitPzFrac_" + fSampleLabels[i]).c_str(),
                "shr_pz_v",
                "Shr Fit Pz Frac",
                "Count",
                20, -1.0, 1.0,
                true)); // remove vector duplicates by taking first element only

        preSelectedShrFitThetaVec.push_back(
            Plotter::CreateTH1DFromRNode(
                nodes[i],
                ("preselection_hist_ShrFitTheta_" + fSampleLabels[i]).c_str(),
                "shr_theta_v",
                "Shr Fit Theta [rad]",
                "Count",
                20, 0.0, 3.14,
                true)); // remove vector duplicates by taking first element only

    }

    // Example of creating a stacked histogram

    Plotter::FullDataMCSignalPlot(preSelectednpfpsVec,
                        fSampleLabels,
                        "preselection_full_hist_npfps",
                        false, // logy
                        fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedNuE2Vec,
                         fSampleLabels,
                         "preselection_full_hist_NeutrinoEnergy2",
                         false, // logy
                         fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedFlashMatchScoreVec,
                         fSampleLabels,
                         "preselection_full_hist_FlashMatchScore",
                         false, // logy
                         fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedShrPhivVec,
                            fSampleLabels,
                            "preselection_full_hist_ShrPhiv",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedTopologicalScoreVec,
                         fSampleLabels,
                         "preselection_full_hist_TopologicalScore",
                         false, // logy
                         fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedShrFitPzFracVec,
                            fSampleLabels,
                            "preselection_full_hist_ShrFitPzFrac",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedShrFitThetaVec,
                            fSampleLabels,
                            "preselection_full_hist_ShrFitTheta",
                            false, // logy
                            fSampleWeights);
}

void PreselectionModule::Finalise()
{
    // Nothing to do here
}
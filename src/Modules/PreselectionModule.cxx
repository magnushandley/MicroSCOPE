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
    , fMakePlots   (cfg.GetValue("Preselection.MakePlots", false))
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
    std::vector<TH1D> preSelectednTracksVec;
    std::vector<TH1D> preSelectedNuE2Vec;
    std::vector<TH1D> preSelectedSliceCaloE2Vec;
    std::vector<TH1D> preSelectedShrdedxmaxVec;
    std::vector<TH1D> preSelectedShrETotVec;
    std::vector<TH1D> preSelectedTrkETotVec;
    std::vector<TH1D> preSelectedMaxTrkEVec;
    std::vector<TH1D> preSelectedShrClusDir2Vec;
    std::vector<TH1D> preSelectedFlashMatchScoreVec;
    std::vector<TH1D> preSelectedTopologicalScoreVec;
    std::vector<TH1D> preSelectedShrPhivVec;
    std::vector<TH1D> preSelectedShrFitPzFracVec;
    std::vector<TH1D> preSelectedShrFitThetaVec;
    std::vector<TH1D> pfng2shravrg_maxEVec;
    std::vector<TH1D> preSelectedpfng2mipfrac_maxEVec;
    std::vector<TH1D> preSelectedpfng2hipfrac_maxEVec;
    std::vector<TH1D> preSelectedpfng2bkgfrac_maxEVec;
    std::vector<TH1D> preSelectedMergedTimeVec;
    std::vector<TH1D> preSelectedInteractionTimeVec;
    std::vector<TH1D> preSelectedUPlaneHitsVec;
    std::vector<TH1D> preSelectedVPlaneHitsVec;
    std::vector<TH1D> preSelectedYPlaneHitsVec;
    std::vector<TH1D> preSelectedShrFitThetaMaxEVec;
    std::vector<TH1D> preSelectedShrFitPzFracMaxEVec;
    std::vector<TH1D> preSelectedShrPhivMaxEVec;
    std::vector<TH1D> preSelectedTrkFitThetaMaxEVec;
    std::vector<TH1D> preSelectedTrkFitPzFracMaxEVec;
    std::vector<TH1D> preSelectedTrkPhivMaxEVec;
    std::vector<TH1D> preSelectedpi0MassYVec;
    std::vector<TH1D> preselectedFlashTimeVec;

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
        
        if (fMakePlots) {
            preSelectednpfpsVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_npfps_" + fSampleLabels[i]).c_str(),
                    "n_pfps", 
                    "Number of PFParticles",
                    "Count",
                    5, 0.5, 5.5));

            preSelectednTracksVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_nTracks_" + fSampleLabels[i]).c_str(),
                    "n_tracks", 
                    "Number of Tracks",
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

            preSelectedSliceCaloE2Vec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_SliceCaloE2_" + fSampleLabels[i]).c_str(),
                    "SliceCaloEnergy2",
                    "Slice Calorimetric Energy [MeV]",
                    "Count",
                    20, 0.0, 500.0));

            preSelectedShrdedxmaxVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_Shrdedxmax_" + fSampleLabels[i]).c_str(),
                    "shr_tkfit_dedx_max",
                    "Max Shower dE/dx [MeV/cm]",
                    "Count",
                    20, 0.0, 10.0));

            preSelectedShrETotVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_ShrETot_" + fSampleLabels[i]).c_str(),
                    "shr_energy_tot",
                    "Total Shower Energy [MeV]",
                    "Count",
                    20, 0.0, 0.25));

            preSelectedTrkETotVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_TrkETot_" + fSampleLabels[i]).c_str(),
                    "trk_energy_tot",
                    "Total Track Energy [MeV]",
                    "Count",
                    20, 0.0, 0.75));

            preSelectedMaxTrkEVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_MaxTrkE_" + fSampleLabels[i]).c_str(),
                    "trk_energy",
                    "Max Track Energy [MeV]",
                    "Count",
                    20, 0.0, 500.0));

            preSelectedShrClusDir2Vec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_ShrClusDir2_" + fSampleLabels[i]).c_str(),
                    "shrclusdir2",
                    "Average Shower Cluster Direction [degrees]",
                    "Count",
                    20, 0, 360.0));

            preSelectedFlashMatchScoreVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_FlashMatchScore_" + fSampleLabels[i]).c_str(),
                    "nu_flashmatch_score", 
                    "Flash Match Score",
                    "Count",
                    30, 0.0, 30.0));

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

            pfng2shravrg_maxEVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_pfng2shravrg_maxE_" + fSampleLabels[i]).c_str(),
                    "pfng2shravrg_maxE",
                    "PF NG2 Shr Average Score",
                    "Count",
                    20, 0.0, 1.0));

            preSelectedpfng2mipfrac_maxEVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_pfng2mipfrac_maxE_" + fSampleLabels[i]).c_str(),
                    "pfng2mipfrac_maxE",
                    "PFP NG2 MIP Fraction",
                    "Count",
                    20, 0.0, 1.0));

            preSelectedpfng2hipfrac_maxEVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_pfng2hipfrac_maxE_" + fSampleLabels[i]).c_str(),
                    "pfng2hipfrac_maxE",
                    "PFP NG2 HIP Fraction",
                    "Count",
                    20, 0.0, 1.0));

            preSelectedpfng2bkgfrac_maxEVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_pfng2bkgfrac_maxE_" + fSampleLabels[i]).c_str(),
                    "pfng2bkgfrac_maxE",
                    "PFP NG2 Background Fraction",
                    "Count",
                    20, 0.0, 1.0));

            preSelectedMergedTimeVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_MergedTime_" + fSampleLabels[i]).c_str(),
                    "interaction_time_merged",
                    "Merged Interaction Time [ns]",
                    "Count",
                    20, 0.0, 18.831));

            preSelectedInteractionTimeVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_InteractionTime_" + fSampleLabels[i]).c_str(),
                    "interaction_time_abs",
                    "Absolute Interaction Time [ns]",
                    "Count",
                    500, 0.0, 20000.0));

            preSelectedUPlaneHitsVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_UPlaneHits_" + fSampleLabels[i]).c_str(),
                    "pfnplanehits_U",
                    "Number of U Plane Hits",
                    "Count",
                    30, 0.0, 300.0));

            preSelectedVPlaneHitsVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_VPlaneHits_" + fSampleLabels[i]).c_str(),
                    "pfnplanehits_V",
                    "Number of V Plane Hits",
                    "Count",
                    30, 0.0, 300.0));

            preSelectedYPlaneHitsVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],   
                    ("preselection_hist_YPlaneHits_" + fSampleLabels[i]).c_str(),
                    "pfnplanehits_Y",
                    "Number of Y Plane Hits",
                    "Count",
                    30, 0.0, 300.0));

            preSelectedShrFitThetaMaxEVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_ShrFitThetaMaxE_" + fSampleLabels[i]).c_str(),
                    "shr_theta_v_maxE",
                    "Shr Fit Theta (max E object) [rad]",
                    "Count",
                    20, 0.0, 3.14));

            preSelectedShrFitPzFracMaxEVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_ShrFitPzFracMaxE_" + fSampleLabels[i]).c_str(),
                    "shr_pz_v_maxE",
                    "Shr Fit Pz Frac (max E object)",
                    "Count",
                    20, -1.0, 1.0));

            preSelectedShrPhivMaxEVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_ShrPhivMaxE_" + fSampleLabels[i]).c_str(),
                    "shr_phi_v_maxE",
                    "Shr Phi (max E object) [rad]",
                    "Count",
                    20, -3.14, 3.14));

            preSelectedTrkFitThetaMaxEVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_TrkFitThetaMaxE_" + fSampleLabels[i]).c_str(),
                    "trk_theta_v_maxE",
                    "Track Fit Theta (max E object) [rad]",
                    "Count",
                    20, 0.0, 3.14));

            preSelectedTrkFitPzFracMaxEVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_TrkFitPzFracMaxE_" + fSampleLabels[i]).c_str(),
                    "trk_dir_z_v_maxE",
                    "Track Fit Pz Frac (max E object)",
                    "Count",
                    20, -1.0, 1.0));

            preSelectedTrkPhivMaxEVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_TrkPhivMaxE_" + fSampleLabels[i]).c_str(),
                    "trk_phi_v_maxE",
                    "Track Phi (max E object) [rad]",
                    "Count",
                    20, -3.14, 3.14));

            preSelectedpi0MassYVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_pi0MassY_" + fSampleLabels[i]).c_str(),
                    "pi0_mass_Y",
                    "Reconstructed #pi^{0} Mass [MeV]",
                    "Count",
                    30, 0.0, 200.0));

            preselectedFlashTimeVec.push_back(
                Plotter::CreateTH1DFromRNode(
                    nodes[i],
                    ("preselection_hist_FlashTime_" + fSampleLabels[i]).c_str(),
                    "flash_time_flash_matching",
                    "Flash Time [ns]",
                    "Count",
                    30, 6.5, 16.5));

        }

    }

    // Example of creating a stacked histogram
    if (!fMakePlots) return;
    Plotter::FullDataMCSignalPlot(preSelectednpfpsVec,
                        fSampleLabels,
                        "preselection_full_hist_npfps",
                        false, // logy
                        fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectednTracksVec,
                         fSampleLabels,
                         "preselection_full_hist_nTracks",
                         false, // logy
                         fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedNuE2Vec,
                         fSampleLabels,
                         "preselection_full_hist_NeutrinoEnergy2",
                         false, // logy
                         fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedSliceCaloE2Vec,
                         fSampleLabels,
                         "preselection_full_hist_SliceCaloE2",
                         false, // logy
                         fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedShrdedxmaxVec,
                         fSampleLabels,
                         "preselection_full_hist_Shrdedxmax",
                         false, // logy
                         fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedShrETotVec,
                         fSampleLabels,
                         "preselection_full_hist_ShrETot",
                         false, // logy
                         fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedTrkETotVec,
                         fSampleLabels,
                         "preselection_full_hist_TrkETot",
                         false, // logy
                         fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedMaxTrkEVec,
                         fSampleLabels,
                         "preselection_full_hist_MaxTrkE",
                         false, // logy
                         fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedShrClusDir2Vec,
                         fSampleLabels,
                         "preselection_full_hist_ShrClusDir2",
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

    Plotter::FullDataMCSignalPlot(pfng2shravrg_maxEVec,
                            fSampleLabels,
                            "preselection_full_hist_pfng2shravrg_maxE",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedpfng2mipfrac_maxEVec,
                            fSampleLabels,
                            "preselection_full_hist_pfng2mipfrac_maxE",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedpfng2hipfrac_maxEVec,
                            fSampleLabels,
                            "preselection_full_hist_pfng2hipfrac_maxE",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedpfng2bkgfrac_maxEVec,
                            fSampleLabels,
                            "preselection_full_hist_pfng2bkgfrac_maxE",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedMergedTimeVec,
                            fSampleLabels,
                            "preselection_full_hist_MergedTime",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedInteractionTimeVec,
                            fSampleLabels,
                            "preselection_full_hist_InteractionTime",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedUPlaneHitsVec,
                            fSampleLabels,
                            "preselection_full_hist_UPlaneHits",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedVPlaneHitsVec,
                            fSampleLabels,
                            "preselection_full_hist_VPlaneHits",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedYPlaneHitsVec,
                            fSampleLabels,
                            "preselection_full_hist_YPlaneHits",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedShrFitThetaMaxEVec,
                            fSampleLabels,
                            "preselection_full_hist_ShrFitThetaMaxE",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedShrFitPzFracMaxEVec,
                            fSampleLabels,
                            "preselection_full_hist_ShrFitPzFracMaxE",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedTrkPhivMaxEVec,
                            fSampleLabels,
                            "preselection_full_hist_TrkPhivMaxE",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedTrkFitThetaMaxEVec,
                            fSampleLabels,
                            "preselection_full_hist_TrkFitThetaMaxE",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedTrkFitPzFracMaxEVec,
                            fSampleLabels,
                            "preselection_full_hist_TrkFitPzFracMaxE",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preSelectedpi0MassYVec,
                            fSampleLabels,
                            "preselection_full_hist_pi0MassY",
                            false, // logy
                            fSampleWeights);

    Plotter::FullDataMCSignalPlot(preselectedFlashTimeVec,
                            fSampleLabels,
                            "preselection_full_hist_FlashTime",
                            false, // logy
                            fSampleWeights);

    
}

void PreselectionModule::Finalise()
{
    // Nothing to do here
}
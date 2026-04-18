#include "Modules/PreselectionModule.hxx"
#include "Utils/Plotter.hxx"
#include "Utils/SystematicsUtil.hxx"
#include "Utils/TimingUtils.hxx"

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

    const std::string labelsString = cfg.GetValue("Preselection.SampleLabels", "");
    std::stringstream ssLabels{labelsString};
    std::string labelToken;
    while (std::getline(ssLabels, labelToken, ',')) {
        // trim leading/trailing whitespace
        labelToken.erase(0, labelToken.find_first_not_of(" \t\n\r"));
        labelToken.erase(labelToken.find_last_not_of(" \t\n\r") + 1);

        // Remove optional surrounding quotes so labels like "Run 3 data" are kept intact
        if (labelToken.size() >= 2 &&
           ((labelToken.front() == '"'  && labelToken.back() == '"') ||
            (labelToken.front() == '\'' && labelToken.back() == '\'')))
        {
            labelToken = labelToken.substr(1, labelToken.size() - 2);
        }

        if (!labelToken.empty()) {
            fSampleLabels.push_back(labelToken);
        }
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

    // Test of whether any further timing alignment is needed. Use the timing utility function from TimingUtils to find the
    // mean timing offsets. Should be zero if no additional shift happens.
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        //Create vector of times from the dataframe
        std::vector<double> times = nodes[i].Take<double>("interaction_time_merged").GetValue();

        auto [A, mu, sigma, C, muError] = TimingUtils::WrappedGaussianFit(times, /*period=*/18.831, /*K=*/3);
        std::cout << "[Preselection] Timing fit results for sample " << fSampleLabels[i] << ":\n";
        std::cout << "  A     = " << A << "\n";
        std::cout << "  mu    = " << mu << "\n";
        std::cout << "  sigma = " << sigma << "\n";
        std::cout << "  C     = " << C << "\n";
        std::cout << "  muError = " << muError << "\n";

        //From this, there is a differenc of 0.5915 ns. For diagnostics, add another branch to the
        //RNode, adding this as a correction to the overlay sample only
        if (fSampleLabels[i].find("overlay") != std::string::npos) {
            nodes[i] = nodes[i].Redefine("interaction_time_merged",
                [](double t) {
                    double corrected_time = t + 0.5915; // Apply the timing correction
                    double remerged_time = std::fmod(corrected_time, 18.831); // Wrap around using the spill period
                    if (remerged_time < 0) remerged_time += 18.831; // Ensure non-negative
                    return remerged_time;
                },
                {"interaction_time_merged"}
            );
            std::cout << "[Preselection] Applied timing correction of " << mu << " ns to sample " << fSampleLabels[i] << ".\n";
        }
    }

    std::vector<TH1D> preSelectednpfpsVec;
    std::vector<TH1D> preSelectednTracksVec;
    std::vector<TH1D> preSelectedNuE2Vec;
    std::vector<TH1D> preSelectedNuE2HigherRangeVec;
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
    std::vector<TH1D> preSelectedUPlaneHitsVec;
    std::vector<TH1D> preSelectedVPlaneHitsVec;
    std::vector<TH1D> preSelectedYPlaneHitsVec;
    std::vector<TH1D> preSelectedShrFitThetaMaxEVec;
    std::vector<TH1D> preSelectedShrFitPzFracMaxEVec;
    std::vector<TH1D> preSelectedShrPhivMaxEVec;
    std::vector<TH1D> preSelectedTrkFitThetaMaxEVec;
    std::vector<TH1D> preSelectedTrkFitPzFracMaxEVec;
    std::vector<TH1D> preSelectedTrkPhivMaxEVec;
    std::vector<TH1D> preSelectedMergedTimeVec;
    std::vector<TH1D> preSelectedFlashTimeVec;
    std::vector<TH1D> preSelectedNG2ShrAvrgMaxEVec;
    std::vector<TH1D> preSelectedNG2ShrAvrgMaxEHighScoresVec;
    std::vector<TH1D> preSelectedpi0MassYVec;

    struct HistSpec {
        std::vector<TH1D>* vec;
        std::string        histNamePrefix;   // e.g. "preselection_hist_npfps_"
        std::string        colName;
        std::string        xTitle;
        std::string        yTitle;
        int                nBins;
        double             xMin;
        double             xMax;
        bool               isAngle;
        bool               removeVectorDuplicates;
        std::string        fullPlotName;     // e.g. "preselection_full_hist_npfps"
        bool               logY;
        std::vector<TH1D>  systVarianceHists; //Histograms with systematic variances per bin
    };

    SystematicsConfig systConfig;
    systConfig.genieMultisimBranch = "weightsGenie";
    systConfig.genieCVWeightBranch = "weight_cv_untuned";
    systConfig.genieGlobalCVWeightBranch = "weight_cv";
    systConfig.ppfxMultisimBranch = "weightsPPFX";
    systConfig.ppfxCVWeightBranch = "weight_cv_noppfx";
    systConfig.ppfxGlobalCVWeightBranch = "weight_cv";
    systConfig.reintMultisimBranch = "weightsReint";
    systConfig.reintCVWeightBranch = "weight_cv";
    systConfig.reintGlobalCVWeightBranch = "weight_cv";

    // Define all histogram specifications in one place
    std::vector<HistSpec> histSpecs;
    histSpecs.reserve(32);
    //histSpecs.push_back({&preSelectednpfpsVec,            "preselection_hist_npfps_",             "n_pfps",             "Number of PFParticles",                 "Count",  5,  0.5,  5.5,  false, false, "preselection_full_hist_npfps",             false});
    //histSpecs.push_back({&preSelectednTracksVec,          "preselection_hist_nTracks_",           "n_tracks",           "Number of Tracks",                      "Count",  5,  0.5,  5.5,  false, false, "preselection_full_hist_nTracks",           false});
    histSpecs.push_back({&preSelectedNuE2Vec,             "preselection_hist_NeutrinoEnergy2_",   "NeutrinoEnergy2",     "Reconstructed Neutrino Energy [MeV]",                 "Count", 20,  0.0,  500.0,false, false, "preselection_full_hist_NeutrinoEnergy2",   false});
    histSpecs.push_back({&preSelectedNuE2HigherRangeVec, "preselection_hist_NuE2HigherRange_",   "NeutrinoEnergy2",     "Reconstructed Neutrino Energy [MeV]",                 "Count", 30,  0.0,  1500.0,false, false, "preselection_full_hist_NuE2HigherRange",   false});
    //histSpecs.push_back({&preSelectedSliceCaloE2Vec,      "preselection_hist_SliceCaloE2_",       "SliceCaloEnergy2",    "Slice Calorimetric Energy [MeV]",        "Count", 20,  0.0,  500.0,false, false, "preselection_full_hist_SliceCaloE2",       false});
    //histSpecs.push_back({&preSelectedShrdedxmaxVec,       "preselection_hist_Shrdedxmax_",        "shr_tkfit_dedx_max",  "Max Shower dE/dx [MeV/cm]",             "Count", 20,  0.0,  10.0, false, false, "preselection_full_hist_Shrdedxmax",        false});
    //histSpecs.push_back({&preSelectedShrETotVec,          "preselection_hist_ShrETot_",           "shr_energy_tot",      "Total Shower Energy [MeV]",              "Count", 20,  0.0,  0.25, false, false, "preselection_full_hist_ShrETot",           false});
    //histSpecs.push_back({&preSelectedTrkETotVec,          "preselection_hist_TrkETot_",           "trk_energy_tot",      "Total Track Energy [MeV]",               "Count", 20,  0.0,  0.75, false, false, "preselection_full_hist_TrkETot",           false});
    //histSpecs.push_back({&preSelectedMaxTrkEVec,          "preselection_hist_MaxTrkE_",           "trk_energy",          "Max Track Energy [MeV]",                 "Count", 20,  0.0,  500.0,false, false, "preselection_full_hist_MaxTrkE",           false});
    //histSpecs.push_back({&preSelectedShrClusDir2Vec,      "preselection_hist_ShrClusDir2_",       "shrclusdir2",         "Average Shower Cluster Direction [degrees]","Count",20,  0.0,  360.0,false, false, "preselection_full_hist_ShrClusDir2",       false});
    //histSpecs.push_back({&preSelectedFlashMatchScoreVec,  "preselection_hist_FlashMatchScore_",   "nu_flashmatch_score", "Flash Match Score",                      "Count", 30,  0.0,  30.0, false, false, "preselection_full_hist_FlashMatchScore",   false});
    //histSpecs.push_back({&preSelectedTopologicalScoreVec, "preselection_hist_TopologicalScore_",  "topological_score",   "Topological Score",                      "Count", 30,  0.0,  1.0,  false, false, "preselection_full_hist_TopologicalScore",  false});
    //histSpecs.push_back({&preSelectedShrPhivVec,          "preselection_hist_ShrPhiv_",           "shr_phi_v",           "Shr Phi [rad]",                          "Count", 20, -3.14, 3.14, true,  false, "preselection_full_hist_ShrPhiv",           false});
    //histSpecs.push_back({&preSelectedShrFitPzFracVec,     "preselection_hist_ShrFitPzFrac_",      "shr_pz_v",            "Shr Fit Pz Frac",                        "Count", 20, -1.0,  1.0,  true,  false, "preselection_full_hist_ShrFitPzFrac",      false});
    //histSpecs.push_back({&preSelectedShrFitThetaVec,      "preselection_hist_ShrFitTheta_",       "shr_theta_v",         "Shr Fit Theta [rad]",                    "Count", 20,  0.0,  3.14, true,  false, "preselection_full_hist_ShrFitTheta",       false});
    //histSpecs.push_back({&preSelectedUPlaneHitsVec,       "preselection_hist_UPlaneHits_",        "pfnplanehits_U",      "Number of U Plane Hits",                 "Count", 30,  0.0,  300.0,false, false, "preselection_full_hist_UPlaneHits",        false});
    //histSpecs.push_back({&preSelectedVPlaneHitsVec,       "preselection_hist_VPlaneHits_",        "pfnplanehits_V",      "Number of V Plane Hits",                 "Count", 30,  0.0,  300.0,false, false, "preselection_full_hist_VPlaneHits",        false});
    //histSpecs.push_back({&preSelectedYPlaneHitsVec,       "preselection_hist_YPlaneHits_",        "pfnplanehits_Y",      "Number of Y Plane Hits",                 "Count", 30,  0.0,  300.0,false, false, "preselection_full_hist_YPlaneHits",        false});
    histSpecs.push_back({&preSelectedShrFitThetaMaxEVec,  "preselection_hist_ShrFitThetaMaxE_",   "shr_theta_v_maxE",    "Shr Fit Theta (max E object) [rad]",     "Count", 20,  0.0,  3.14, false, false, "preselection_full_hist_ShrFitThetaMaxE",   false});
    histSpecs.push_back({&preSelectedShrFitPzFracMaxEVec, "preselection_hist_ShrFitPzFracMaxE_",  "shr_pz_v_maxE",       "Momentum Fraction in Forward Direction",         "Count", 20, -1.0,  1.0,  false, false, "preselection_full_hist_ShrFitPzFracMaxE",  false});
    histSpecs.push_back({&preSelectedShrPhivMaxEVec,      "preselection_hist_ShrPhivMaxE_",       "shr_phi_v_maxE",      "Shr Phi (max E object) [rad]",           "Count", 20, -3.14, 3.14, false, false, "preselection_full_hist_ShrPhivMaxE",       false});
    //histSpecs.push_back({&preSelectedTrkFitThetaMaxEVec,  "preselection_hist_TrkFitThetaMaxE_",   "trk_theta_v_maxE",    "Track Fit Theta (max E object) [rad]",   "Count", 20,  0.0,  3.14, false, false, "preselection_full_hist_TrkFitThetaMaxE",   false});
    //histSpecs.push_back({&preSelectedTrkFitPzFracMaxEVec, "preselection_hist_TrkFitPzFracMaxE_",  "trk_dir_z_v_maxE",    "Track Fit Pz Frac (max E object)",       "Count", 20, -1.0,  1.0,  false, false, "preselection_full_hist_TrkFitPzFracMaxE",  false});
    //histSpecs.push_back({&preSelectedTrkPhivMaxEVec,      "preselection_hist_TrkPhivMaxE_",       "trk_phi_v_maxE",      "Track Phi (max E object) [rad]",         "Count", 20, -3.14, 3.14, false, false, "preselection_full_hist_TrkPhivMaxE",       false});
    //histSpecs.push_back({&preSelectedMergedTimeVec,       "preselection_hist_MergedTime_",        "interaction_time_merged","Merged Interaction Time [ns]",   "Count", 20,  0.0,  18.831,false, false, "preselection_full_hist_MergedTime",        false});
    //histSpecs.push_back({&preSelectedFlashTimeVec,        "preselection_hist_FlashTime_",         "flash_time_flash_matching", "Flash Match Time [ns]",           "Count", 30,  0.0,  20.0, false, false, "preselection_full_hist_FlashTime",         false});
    histSpecs.push_back({&preSelectedNG2ShrAvrgMaxEVec,   "preselection_hist_NG2ShrAvrgMaxE_",    "pfng2shravrg_maxE",  "NuGraph Average Shower Score", "Count", 20,  0, 1.0, false, false, "preselection_full_hist_NG2ShrAvrgMaxE",   false});
    histSpecs.push_back({&preSelectedNG2ShrAvrgMaxEHighScoresVec,   "preselection_hist_NG2ShrAvrgMaxEHighScores_",    "pfng2shravrg_maxE",  "NuGraph Average Shower Score", "Count", 20,  0.5, 1.0, false, false, "preselection_full_hist_NG2ShrAvrgMaxEHighScores",   false});
    //histSpecs.push_back({&preSelectedpi0MassYVec,         "preselection_hist_pi0MassY_",          "pi0_mass_Y",          "Reconstructed pi0 Mass (Y Plane) [MeV]",             "Count", 20,  0.0,  250.0,false, false, "preselection_full_hist_pi0MassY",          false});

    auto fillHistogramsForSample = [&](ROOT::RDF::RNode &node,
                                       const std::string &sampleLabel,
                                       const std::string &weightCol) {
        for (auto &spec : histSpecs) {
            spec.vec->push_back(
                Plotter::CreateTH1DFromRNode(
                    node,
                    (spec.histNamePrefix + sampleLabel).c_str(),
                    spec.colName,
                    spec.xTitle,
                    spec.yTitle,
                    spec.nBins, spec.xMin, spec.xMax,
                    spec.isAngle,
                    spec.removeVectorDuplicates,
                    weightCol));
            // Fill systematic variance histograms - special treatment per sample: overlay gets all, dirt gets 75% normalisation, others get none
            TH1D &nominalHist = spec.vec->back();
            if (sampleLabel == "Run 4b in-cryo nu (overlay)") {
                std::cout << "    Computing systematic variance histograms for sample: " << sampleLabel << "\n";
                SystematicsUtil sysUtil;
                TH1D overlaySystHist = sysUtil.RunAllMultisimSystematics(
                    nominalHist,
                    node,
                    spec.colName,
                    systConfig);
                spec.systVarianceHists.push_back(overlaySystHist);
            } else if (sampleLabel == "Run 4b out-of-cryo nu (dirt)") {
                // Create per-bin variance histogram for a 75% normalisation uncertainty:
                //   var = (0.75 * yield)^2
                TH1D dirtVar = nominalHist;
                dirtVar.SetName((std::string(nominalHist.GetName()) + "_dirtNormVar").c_str());
                dirtVar.Reset("ICES");

                for (int bin = 1; bin <= dirtVar.GetNbinsX(); ++bin) {
                    const double content  = nominalHist.GetBinContent(bin);
                    const double variance = (0.75 * content) * (0.75 * content);
                    dirtVar.SetBinContent(bin, variance);
                    dirtVar.SetBinError(bin, 0.0);
                }

                spec.systVarianceHists.push_back(dirtVar);
            } else {
                // No systematic variance for other samples
                TH1D emptyVar = nominalHist;
                emptyVar.SetName((std::string(nominalHist.GetName()) + "_emptyVar").c_str());
                emptyVar.Reset("ICES");
                spec.systVarianceHists.push_back(emptyVar);
            }
        }
    };

    auto plotAllHistograms = [&]() {
        for (auto &spec : histSpecs) {
            //Compute weighted total variance histograms for each sample
            TH1D totalVarianceHist(spec.fullPlotName.c_str(), (spec.fullPlotName + " Total Variance").c_str(), spec.nBins, spec.xMin, spec.xMax);
            if (spec.systVarianceHists.size() != fSampleWeights.size()) {
                std::cerr << "[Preselection] ERROR: systVarianceHists size (" << spec.systVarianceHists.size()
                          << ") != fSampleWeights size (" << fSampleWeights.size()
                          << ") for plot " << spec.fullPlotName << "\n";
            }
            const size_t n = std::min(spec.systVarianceHists.size(), fSampleWeights.size());
            for (size_t i = 0; i < n; ++i) {
                TH1D varHist = spec.systVarianceHists[i];
                double weight = fSampleWeights[i];
                varHist.Scale(weight * weight); // Scale variance
                totalVarianceHist.Add(&varHist);
            }
            // Debug: Print total variance
            std::cout << "Total variance for plot " << spec.fullPlotName << ":\n";
            for (int bin = 1; bin <= totalVarianceHist.GetNbinsX(); ++bin) {
                std::cout << "  Bin " << bin << ": " << totalVarianceHist.GetBinContent(bin) << "\n";
            }
            Plotter::FullDataMCSignalPlot(
                *spec.vec,
                fSampleLabels,
                spec.fullPlotName,
                spec.logY,
                fSampleWeights,
                0.7, 1.3,
                &totalVarianceHist);
        }
    };

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        std::cout << "\n[Preselection] Writing output for sample: " << fSampleLabels[i] << '\n';
        std::cout << "    to file: " << fOutFiles[i] << '\n';
        nodes[i].Snapshot(fTreeName, fOutFiles[i], fVarsToKeep, opt);

        if (!fMakePlots) continue;

        // Weighting policy for plots
        std::string weightCol;
        if (fSampleLabels[i] == "Run 4b in-cryo nu (overlay)" || fSampleLabels[i] == "Run 4b out-of-cryo nu (dirt)") {
            // Apply CV weight to MC samples that require it
            std::cout << "    Applying weight_cv to MC histograms for sample: " << fSampleLabels[i] << "\n";
            weightCol = "weight_cv";
        } else {
            weightCol = "";
        }

        // Fill all plot histograms in one go
        fillHistogramsForSample(nodes[i], fSampleLabels[i], weightCol);
    }

    // Example of creating stacked histograms
    if (!fMakePlots) return;
    plotAllHistograms();
}

void PreselectionModule::Finalise()
{
    // Nothing to do here
}
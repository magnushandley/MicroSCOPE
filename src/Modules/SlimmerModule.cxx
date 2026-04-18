#include "Modules/SlimmerModule.hxx"
#include "Utils/Plotter.hxx"
#include "Utils/TimingUtils.hxx"

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
    , fTreeName     (cfg.GetValue("Slimmer.TreeName","nuselection/NeutrinoSelectionFilter"  ))
    , fRunLabel     (cfg.GetValue("Global.RunLabel","run_x") )
    , fMakePlots   (cfg.GetValue("Slimmer.MakePlots", false))
    , fBeamSpillPeriod (cfg.GetValue("Global.BeamSpillPeriod", 18.831))
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


        // Fiducial variables to assess containment (taken from HNL analysis). The whole mess with big and small
        // values is because in ext files the trk_sce_start_x_v vectors can be empty if there is no neutrino slice.
        // In overlay this doesn't happen, but need to for data/ext files I think, so I set min/max to values outside the fiducial volume.

        using VecF = const std::vector<float>&;
        using VecI = const std::vector<int>&;

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
	  .Define("trk_score_v_first",
            [](VecF v) {
                return v.empty() ? -9999.0f : v[0];
            },
            {"trk_score_v"})
        .Define("shr_theta_v_first",
            [](VecF v) {
                return v.empty() ? -9999.0f : v[0];
            },
            {"shr_theta_v"})
        .Define("shr_px_v_first",
            [](VecF v) {
                return v.empty() ? -9999.0f : v[0];
            },
            {"shr_px_v"})
        .Define("trk_end_x_v_first",
            [](VecF v) {
                return v.empty() ? -9999.0f : v[0];
            },
            {"trk_end_x_v"})
        .Define("shr_phi_v_first",
            [](VecF v) {
                return v.empty() ? -9999.0f : v[0];
            },
            {"shr_phi_v"})
        .Define("shr_pz_v_first",
            [](VecF v) {
                return v.empty() ? -9999.0f : v[0];
            },
            {"shr_pz_v"})
        .Define("trk_theta_v_first",
            [](VecF v) {
                return v.empty() ? -9999.0f : v[0];
            },
            {"trk_theta_v"})
        .Define("trk_phi_v_first",
            [](VecF v) {
                return v.empty() ? -9999.0f : v[0];
            },
            {"trk_phi_v"})
        .Define("trk_dir_z_v_first",
            [](VecF v) {
                return v.empty() ? -9999.0f : v[0];
            },
            {"trk_dir_z_v"})
        .Define("trk_calo_energy_u_v_first",
            [](VecF v) {
                return v.empty() ? -9999.0f : v[0];
            },
            {"trk_calo_energy_u_v"})
        .Define("pfnplanehits_U_sum",
            [](VecI v) {
                int sum = 0;
                for (const auto& val : v) sum += val;
                return sum;
            },
            {"pfnplanehits_U"})
        .Define("pfnplanehits_V_sum",
            [](VecI v) {
                int sum = 0;
                for (const auto& val : v) sum += val;
                return sum;
            },
            {"pfnplanehits_V"})
        .Define("pfnplanehits_Y_sum",
            [](VecI v) {
                int sum = 0;
                for (const auto& val : v) sum += val;
                return sum;
            },
            {"pfnplanehits_Y"})
        .Define("trk_score_v_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"trk_score_v", "pfnplanehits_Y"})
        .Define("shr_theta_v_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"shr_theta_v", "pfnplanehits_Y"})
        .Define("shr_px_v_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"shr_px_v", "pfnplanehits_Y"})
        .Define("trk_end_x_v_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"trk_end_x_v", "pfnplanehits_Y"})
        .Define("shr_phi_v_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"shr_phi_v", "pfnplanehits_Y"})
        .Define("shr_pz_v_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"shr_pz_v", "pfnplanehits_Y"})
        .Define("trk_theta_v_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"trk_theta_v", "pfnplanehits_Y"})
        .Define("trk_phi_v_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"trk_phi_v", "pfnplanehits_Y"})
        .Define("trk_dir_z_v_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"trk_dir_z_v", "pfnplanehits_Y"})
        .Define("trk_calo_energy_u_v_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"trk_calo_energy_u_v", "pfnplanehits_Y"})
        .Define("pfnplanehits_U_maxE",
            [](VecI var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -1 : var[maxEIndex];
            },
            {"pfnplanehits_U", "pfnplanehits_Y"})
        .Define("pfnplanehits_V_maxE",
            [](VecI var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -1 : var[maxEIndex];
            },
            {"pfnplanehits_V", "pfnplanehits_Y"})
        .Define("pfnplanehits_Y_maxE",
            [](VecI var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -1 : var[maxEIndex];
            },
            {"pfnplanehits_Y", "pfnplanehits_Y"})
        .Define("pfng2shravrg_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"pfng2shravrg", "pfnplanehits_Y"})
        .Define("pfng2mipfrac_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"pfng2mipfrac", "pfnplanehits_Y"})
        .Define("pfng2hipfrac_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"pfng2hipfrac", "pfnplanehits_Y"})
        .Define("pfng2bkgfrac_maxE",
            [](VecF var, VecI E){
                int maxEIndex = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
                return var.empty() ? -9999.0f : var[maxEIndex];
            },
            {"pfng2bkgfrac", "pfnplanehits_Y"});
        
        //.Filter("(par_decay_vz > 70000 && par_decay_pz < 0.01)");

        // Append timing-related columns onto this node and snapshot from it
        ROOT::RDF::RNode dfOut = df1;
        

        // Special handling for timing variables
        // Ext files have garbage times, so just add a random time within the spill window
        if (fSampleLabels[fileIndex].find("beamoff") != std::string::npos) {
            std::cout << "[Slimmer] Adding random timing offsets for beam-off file.\n";
            dfOut = df1.Define(
                "interaction_time_merged",
                [beamSpillPeriod](float /*time*/) {
                    const double random_offset = gRandom->Uniform(0.0, beamSpillPeriod);
                    //std::cout << "[Slimmer] Assigned random timing offset for beam-off event: " << random_offset << " ns\n";
                    return random_offset;
                },
                {"interaction_time_abs"}
            )
            .Redefine("Med_TT3", []() {
                return 99999.0f;   // or whatever placeholder you want
            });
        }
        else if (fSampleLabels[fileIndex].find("data") != std::string::npos || fSampleLabels[fileIndex].find("overlay") != std::string::npos || fSampleLabels[fileIndex].find("dirt") != std::string::npos) {
            std::cout << "[Slimmer] Fitting timing offsets for data file.\n";
            auto run_numbers = df1.Take<int>("run").GetValue();
            auto times_f     = df1.Take<float>("interaction_time_abs").GetValue();
            std::cout << "[Slimmer] Fitting timing offsets for " << run_numbers.size() << " events.\n";
            std::vector<double> times;
            times.reserve(times_f.size());
            for (const float t : times_f) times.push_back(static_cast<double>(t));
            std::cout << "[Slimmer] Creating run offset map.\n";
            std::unordered_map<int, std::pair<double, double>> runOffsetMap =
                Analysis::TimingUtils::CreateRunOffsetMap(
                    run_numbers,
                    times,
                    beamSpillPeriod,
                    1,  // K
                    50  // runWindowSize
                );
            std::cout << "[Slimmer] Created run offset map with " << runOffsetMap.size() << " entries.\n";
            auto runOffsetMapPtr = std::make_shared<const std::unordered_map<int, std::pair<double, double>>>(std::move(runOffsetMap));
            for (const auto& [run, offset] : *runOffsetMapPtr) {
                std::cout << "[Slimmer] Run " << run << " has timing offset: " << offset.first << " ns (uncertainty: " << offset.second << " ns)\n";
            }

            // Scatter plot of run vs timing offset
            std::vector<double> runs;
            std::vector<double> offsets;
            std::vector<double> uncertainties;

            runs.reserve(runOffsetMapPtr->size());
            offsets.reserve(runOffsetMapPtr->size());
            uncertainties.reserve(runOffsetMapPtr->size());

            for (const auto& [run, offset] : *runOffsetMapPtr) {
                runs.push_back(static_cast<double>(run));
                offsets.push_back(offset.first);
                uncertainties.push_back(offset.second);
            }

            // Make scatter plot
            auto c = std::make_unique<TCanvas>("c_run_offset", "Run vs timing offset", 800, 600);

            std::vector<double> xerr(runs.size(), 0.0);
            auto g = std::make_unique<TGraphErrors>(
                static_cast<int>(runs.size()),
                runs.data(),
                offsets.data(),
                xerr.data(),
                uncertainties.data()
            );

            g->SetTitle("Timing offset vs run;Run number;Timing offset [ns]");
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.9);
            g->SetMarkerColor(kBlue+1);

            g->Draw("AP");

            c->SaveAs(("timing_offset_vs_run_" + fSampleLabels[fileIndex] + ".pdf").c_str());


            std::cout << "[Slimmer] Applying run-by-run timing offsets for data file.\n";
            dfOut = df1
                .Define(
                    "interaction_time_merged",
                    [beamSpillPeriod, runOffsetMapPtr](float time, int run) {
                        // Wrap into spill period
                        const auto it = runOffsetMapPtr->find(run);
                        const double offset = (it == runOffsetMapPtr->end()) ? 0.0 : it->second.first;
                        //std::cout << "[Slimmer] Run " << run << " applying offset: " << offset << " ns\n";
                        double time_corrected = time - offset;
                        double wrappedTime = std::fmod(time_corrected, beamSpillPeriod);
                        if (wrappedTime < 0) wrappedTime += beamSpillPeriod;
                        return wrappedTime;
                    },
                    {"interaction_time_abs", "run"}
                );
        }
        else {
            // For non-data, just wrap into spill period
            std::cout << "[Slimmer] Wrapping timing for non-data or ext file.\n";
            dfOut = df1.Define(
                "interaction_time_merged",
                [beamSpillPeriod](float time) {
                    double wrappedTime = std::fmod(static_cast<double>(time), beamSpillPeriod);
                    if (wrappedTime < 0) wrappedTime += beamSpillPeriod;
                    return wrappedTime;
                },
                {"interaction_time_abs"}
            );
        }

        std::cout << "[Slimmer] Added interaction_time_merged variable.\n";

        if (fSampleLabels[fileIndex].find("overlay") != std::string::npos || fSampleLabels[fileIndex].find("dirt") != std::string::npos) {
            // Temporary way to add in central value weights for overlay and dirt, otherwise set to 1.0, should do this by label in future
            dfOut = dfOut.Define("weight_cv",
                [](float w1, float w2, int npi0) {
                    float safeWeight1 = (w1 > 0.0f && !std::isnan(w1) && !std::isinf(w1) && w1 < 100) ? w1 : 1.0f;
                    float safeWeight2 = (w2 > 0.0f && !std::isnan(w2) && !std::isinf(w2) && w2 < 100) ? w2 : 1.0f;
                    if (npi0 > 0) {
                        safeWeight2 = safeWeight2 * 0.759; // Scaling on files with pi0's, as done by David
                    }
                    return safeWeight1 * safeWeight2;
                }, {"weightSplineTimesTune", "ppfx_cv", "npi0"})
                .Define("weight_cv_untuned",
                [](float w1, float w2, int npi0) {
                    float safeWeight1 = (w1 > 0.0f && !std::isnan(w1) && !std::isinf(w1) && w1 < 100) ? w1 : 1.0f;
                    float safeWeight2 = (w2 > 0.0f && !std::isnan(w2) && !std::isinf(w2) && w2 < 100) ? w2 : 1.0f;
                    if (npi0 > 0) {
                        safeWeight2 = safeWeight2 * 0.759; // Scaling on files with pi0's, as done by David
                    }
                    return safeWeight1 * safeWeight2;
                }, {"weightSpline", "ppfx_cv", "npi0"})
                .Define("weight_cv_nosplineortune",
                [](float w, int npi0) {
                    float safeWeight = (w > 0.0f && !std::isnan(w) && !std::isinf(w) && w < 100) ? w : 1.0f;
                    if (npi0 > 0) {
                        safeWeight = safeWeight * 0.759; // Scaling on files with pi0's, as done by David
                    }
                    return safeWeight;
                }, {"ppfx_cv", "npi0"})
                .Define("weight_cv_noppfx",
                [](float w1, int npi0) {
                    float safeWeight = (w1 > 0.0f && !std::isnan(w1) && !std::isinf(w1) && w1 < 100) ? w1 : 1.0f;
                    if (npi0 > 0) {
                        safeWeight = safeWeight * 0.759; // Scaling on files with pi0's, as done by David
                    }
                    return safeWeight;
                }, {"weightSplineTimesTune", "npi0"});
        }
        else {
            dfOut = dfOut.Define("weight_cv",
                []() {
                    return 1.0f;
                });
        }
        
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
            Plotter::SaveHist(
                dfOut.Histo1D({"sub_hist", ";run_number;Count", 50, 0, 600}, "sub").GetPtr(),
                "slimmer_"+fRunLabel+"_run_histogram" , "prelim");
            Plotter::SaveHist(
                dfOut.Histo1D({"par_decay_vz_hist", ";Kaon Decay Vertex Z [cm];Count", 50, 0, 75000}, "par_decay_vz").GetPtr(),
                "slimmer_"+fRunLabel+"_kaon_decay_vz_histogram" , "prelim");
            Plotter::SaveHist(
                dfOut.Histo1D({"par_decay_pz_hist", ";Kaon Decay Pz [GeV];Count", 50, 0, 12}, "par_decay_pz").GetPtr(),
                "slimmer_"+fRunLabel+"_kaon_decay_pz_histogram" , "prelim");
            Plotter::SaveHist(
                dfOut.Histo1D({"reco_minus_true_t_hist", ";Reconstructed - True Time [ns];Count", 100, 4000, 4200}, "reco_minus_true_time").GetPtr(),
                "slimmer_"+fRunLabel+"_reco_minus_true_t_histogram" , "prelim");
            Plotter::SaveHist(
                dfOut.Histo1D({"pmt_time_hist", ";PMT Time [ns];Count", 100, 0, 20}, "pmt_time").GetPtr(),
                "slimmer_"+fRunLabel+"_pmt_time_histogram" , "prelim");
            std::cout << "[Slimmer] Checking if mc_interaction_time column exists for plotting...\n";
            Plotter::SaveHist(
                dfOut.Histo1D({"mc_interaction_time_hist", ";True Interaction Time [ns];Count", 1000, -10000, 20000}, "mc_interaction_time").GetPtr(),
                "slimmer_"+fRunLabel+"_mc_interaction_time_histogram" , "prelim");
        }

        fileIndex++;
    }
}

//------------------------------------------------------------------------------
void SlimmerModule::Finalise()
{
    // Nothing to do – Snapshot already wrote the slimmed tree.
}

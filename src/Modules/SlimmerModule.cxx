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
    , fTreeName     (cfg.GetValue("Slimmer.TreeName","nuselection/NeutrinoSelectionFilter"  ))
    , fRunLabel     (cfg.GetValue("Global.RunLabel","run_x") )
    , fMakePlots   (cfg.GetValue("Slimmer.MakePlots", false))
    , fBeamSpillPeriod (cfg.GetValue("Global.BeamSpillPeriod", 18.831))
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

        ROOT::RDF::RNode df1 = df;

        //Create derived varibles based on config.
        for (const auto& op : logicConfigs) {
            df1 = LogicRegistry::Instance().Apply(df1, op);
        }

        // Append timing-related columns onto this node and snapshot from it
        ROOT::RDF::RNode dfOut = df1;
        

        // Special handling for timing variables
        // Ext files have garbage times, so just add a random time within the spill window
        //if (fSampleLabels[fileIndex].find("beamoff") != std::string::npos) {
        //    std::cout << "[Slimmer] Adding random timing offsets for beam-off file.\n";
        //    dfOut = df1.Define(
        //        "interaction_time_merged",
        //        [beamSpillPeriod](float) {
        //            const double random_offset = gRandom->Uniform(0.0, beamSpillPeriod);
        //            //std::cout << "[Slimmer] Assigned random timing offset for beam-off event: " << random_offset << " ns\n";
        //            return random_offset;
        //        },
        //        {"interaction_time_abs"}
        //    )
        //    .Redefine("Med_TT3", []() {
        //        return 99999.0f;   // arbitrary placeholder
        //    });
        //}

        //If not ext file, apply offsets.
        if (fSampleLabels[fileIndex].find("beamoff") == std::string::npos && fSampleLabels[fileIndex].find("signal") == std::string::npos) {
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

        //Adding info required for systematics

        if (fSampleLabels[fileIndex].find("overlay") != std::string::npos || fSampleLabels[fileIndex].find("dirt") != std::string::npos || fSampleLabels[fileIndex].find("detvar") != std::string::npos || fSampleLabels[fileIndex].find("signal") != std::string::npos) {
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
                })
                .Define("weight_cv_untuned",
                []() {
                    return 1.0f;
                })
                .Define("weight_cv_nosplineortune",
                []() {
                    return 1.0f;
                })
                .Define("weight_cv_noppfx",
                []() {
                    return 1.0f;
                })
                .Define("ppfx_cv",
                []() {
                    return 1.0f;
                })
                .Define("weightsFlux",
                []() {
                    return 1.0f;
                })
                .Define("weightsGenie",
                []() {
                    return 1.0f;
                })
                .Define("weightsReint",
                []() {
                    return 1.0f;
                })
                .Define("weightsPPFX",
                []() {
                    return 1.0f;
                })
                .Define("weightSplineTimesTune",
                []() {
                    return 1.0f;
                })
                .Define("weightSpline",
                []() {
                    return 1.0f;
                })
                .Define("weightTune",
                []() {
                    return 1.0f;
                });
        }

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

        //----------------------------------------------------------------------
        // 3.  Example debugging plots (you can see me debugging some ns timing variables here...)
        //----------------------------------------------------------------------

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

#include "Modules/HistFitModule.hxx"
#include "Utils/Plotter.hxx"

#include <TEnv.h>
#include <TFile.h>
#include <TString.h>
#include <algorithm>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TTree.h>
#include <TH1D.h>
#include <TChain.h>
#include <cstdio>

#include <RooStats/HistFactory/MakeModelAndMeasurementsFast.h>
#include <RooStats/HistFactory/Measurement.h>

#include <RooStats/AsymptoticCalculator.h>
#include <RooStats/HypoTestInverter.h>
#include <RooStats/HypoTestInverterResult.h>
#include <RooStats/ProfileLikelihoodTestStat.h>
#include <RooStats/HypoTestInverterPlot.h> 

#include "TRandom3.h"


using namespace Analysis;

HistFitModule::HistFitModule(const TEnv& cfg)
    : Module(cfg)
    , fTreeName   (cfg.GetValue("HistFitModule.TreeName", "nuselection/NeutrinoSelectionFilter"))
    , fDataPOT    (cfg.GetValue("HistFitModule.DataPOT", 1.0e20))
    , fSignalPOT  (cfg.GetValue("HistFitModule.SignalPOT", 1.0e20)) // We only need data and signal POT because the backgrounds are scaled to data POT anyway
    , fSimulatedSignalU2(cfg.GetValue("HistFitModule.SimulatedSignalU2", 1.0e-4))
{

    std::stringstream ssInput{cfg.GetValue("HistFitModule.InputFiles", "")};
    std::string inputItem;
    while (ssInput >> inputItem) {
        if (inputItem.back()==',') inputItem.pop_back();
        fInputFiles.push_back(inputItem);
    }

    std::stringstream ssLabels{cfg.GetValue("HistFitModule.SampleLabels", "")};
    std::string label;
    while (ssLabels >> label) {
        if (label.back()==',') label.pop_back();
        fSampleLabels.push_back(label);
    }

    std::stringstream ssWeights{cfg.GetValue("HistFitModule.SampleWeights", "")};
    double weight;
    while (ssWeights >> weight) {
        fSampleWeights.push_back(weight);
    }

    // If provided, then use only the subset of events corresponding to the test samples of the bdt training
    std::stringstream ssTestFractions{cfg.GetValue("HistFitModule.TestFractions", "")};
    double testFraction;
    while (ssTestFractions >> testFraction) {
        fTestFractions.push_back(testFraction);
    }
}


//------------------------------------------------------------------------------
std::vector<ROOT::RDF::RNode>
HistFitModule::BuildDataFrames(const std::vector<std::string>& files,
                               const std::string& treeName,
                               const std::vector<double>& testFractions) const
{
    //std::vector<std::unique_ptr<ROOT::RNode>> dfVec;

    //Because we need to filter, the return type is an RNode
    std::vector<ROOT::RDF::RNode> nodes;
    std::cout << "[HistFitModule] Building DataFrames for " << files.size() << " input files\n";
    nodes.reserve(files.size());
    std::cout << "[HistFitModule] Nodes size reserved: " << nodes.size() << std::endl;

    int it = 0;
    for (const auto& fname : files) {
        auto file = TFile::Open(fname.c_str());
        if (!file || file->IsZombie()) {
            throw std::runtime_error("[HistFitModule] Cannot open file: " + fname);
        }
        auto tree = file->Get<TTree>(treeName.c_str());
        if (!tree) {
            throw std::runtime_error("[HistFitModule] Cannot find tree: " + treeName);
        }
        auto RDF = std::make_unique<ROOT::RDataFrame>(*tree);
        
        // Different from all other modules: apply test sample fraction if provided
        // Default is 1.0 if not provided
        double testFraction = 1.0;
        if (testFractions.size() != 0) {
            testFraction = testFractions[it];
        }
        int nEvents = RDF->Count().GetValue();
        int nKeep = static_cast<int>(nEvents * std::abs(testFraction));
        // Depending on whether we want to keep the first or last fraction, the sign of testFraction is different
        //ROOT::RDF::RNode filteredRDF;
        if (testFraction > 0)
            nodes.push_back(RDF->Range(0, nKeep)); // Take first nKeep events
        else
            nodes.push_back(RDF->Range(nEvents - nKeep, nEvents)); // Take last nKeep events

        ++it;
    }
    if (nodes.empty())
        throw std::runtime_error("[HistFitModule] No data frames created!");

    return nodes;
}

void HistFitModule::SaveHistograms(const std::vector<TH1D>& hists,
                    const std::vector<std::string>& labels,
                    const std::vector<double>& weights,
                    const std::string& fileName)
{
    if (hists.size() != labels.size()) {
        throw std::runtime_error("[Plotter] Number of histograms and labels do not match in SaveHistograms");
    }

    TFile outFile(fileName.c_str(), "RECREATE");
    if (outFile.IsZombie()) {
        throw std::runtime_error("[Plotter] Cannot create output file: " + fileName);
    }

    for (std::size_t i = 0; i < hists.size(); ++i) {
        // Make a local copy so we do not modify the original histogram vector
        TH1D h = hists[i];

        // Ensure we have Sumw2 so that errors are stored and scaled correctly
        h.Sumw2();

        // Apply the sample weight if provided
        if (i < weights.size()) {
            h.Scale(weights[i]);
        }

        // Name the histogram according to the sample label
        h.SetName(labels[i].c_str());

        // Detach from any existing directory and write to the output file
        h.SetDirectory(&outFile);
        h.Write();
    }

    outFile.Close();
}

 std::unique_ptr<RooWorkspace> HistFitModule::BuildModelWorkspace(
    std::vector<TH1D>& histVec,
    const std::vector<std::string>& labels,
    const std::string& inputFile) const
{
    //std::string InputFile = "/Users/magnus/Documents/PhD/MicroSCOPE/build/run/bdt_score_histograms.root";
   // in case the file is not found

    bool bfile = gSystem->AccessPathName(inputFile.c_str());
    if (bfile) {
        std::cout << "Input file is not found - run prepareHistFactory script " << std::endl;
        gROOT->ProcessLine(".! prepareHistFactory .");
        bfile = gSystem->AccessPathName(inputFile.c_str());
        if (bfile) {
            std::cout << "Still no " << inputFile << ", giving up.\n";
            exit(1);
        }
    }

    // Build a RooStats HistFactory model from the input histograms, creating s and s+b models
    // Convention: last histogram in histVec is data, second last is signal, all others are backgrounds
    std::cout << "[HistFitModule] Building model workspace from histograms" << std::endl;
    RooStats::HistFactory::Measurement meas("meas", "meas");
    std::cout << "[HistFitModule] Measurement created" << std::endl;
    meas.SetOutputFilePrefix("./results/example_UsingC"); // Set this from config later
    meas.SetPOI("SigXsecOverSim");
    //meas.AddConstantParam("alpha_syst1");
    //meas.AddConstantParam("Lumi");
    std::cout << "[HistFitModule] Measurement configured" << std::endl;
    
    meas.SetLumi(1.0);
    //meas.SetLumiRelErr(0.3); // Investigate impact of this
    //meas.SetBinHigh(2); // Investigate impact of this

    // Create a channel

    RooStats::HistFactory::Channel chan("channel1");
    //TH1D& hData = histVec.back();
    chan.SetData(labels.back(), inputFile);
    std::cout << "[HistFitModule] Data sample set: " << labels.back() << std::endl;
    chan.SetStatErrorConfig(0.02, "Poisson"); // Investigate impact of this


    // Now, create some samples

    // Create the signal sample
    //TH1D& hsig = histVec[histVec.size() - 2];
    //RooStats::HistFactory::Sample signal(hsig.GetName());
    RooStats::HistFactory::Sample signal(labels[histVec.size() - 2], labels[histVec.size() - 2], inputFile);
    std::cout << "[HistFitModule] Signal sample created" << std::endl;
    std::cout << "[HistFitModule] Signal sample label: " << labels[histVec.size() - 2] << std::endl;
    //signal.SetHisto(&hsig);
    signal.AddNormFactor("SigXsecOverSim", 1, 0, 0.01); 
    chan.AddSample(signal);

    std::cout << "[HistFitModule] Signal sample added" << std::endl;

    // Background 1
    //TH1D& h1 = histVec[0];
    RooStats::HistFactory::Sample background1(labels[0], labels[0], inputFile);
    std::cout << "[HistFitModule] Background 1 sample created" << std::endl;
    std::cout << "[HistFitModule] Background 1 label: " << labels[0] << std::endl;
    //background1.SetHisto(&h1);
    background1.ActivateStatError();
    //background1.AddOverallSys("syst2", 0.95, 1.05);
    chan.AddSample(background1);

    std::cout << "[HistFitModule] Background 1 sample added" << std::endl;

    // Background 2
    //TH1D& h2 = histVec[1];
    RooStats::HistFactory::Sample background2(labels[1], labels[1], inputFile);
    //background2.SetHisto(&h2);
    background2.ActivateStatError();
    //background2.AddOverallSys("syst3", 0.95, 1.05);
    chan.AddSample(background2);

    std::cout << "[HistFitModule] Background 2 sample added" << std::endl;

    // Background 3
    //TH1D& h3 = histVec[2];
    RooStats::HistFactory::Sample background3(labels[2], labels[2], inputFile);
    //background3.SetHisto(&h3);
    background3.ActivateStatError();
    //background3.AddOverallSys("syst3", 0.95, 1.05);
    chan.AddSample(background3);

    std::cout << "[HistFitModule] Background 3 sample added" << std::endl;

    // Done with this channel
    // Add it to the measurement:
    meas.AddChannel(chan);

    std::cout << "[HistFitModule] Channel added to measurement" << std::endl;

    // Collect the histograms from their files,
    // print some output,
    meas.CollectHistograms();
    meas.PrintTree();

    std::unique_ptr<RooWorkspace> ws{MakeModelAndMeasurementFast(meas)};
    std::cout << "[HistFitModule] Model workspace built" << std::endl;
    return ws;
}

RooStats::ModelConfig* HistFitModule::GetSPlusBModel(RooWorkspace* ws) const
{
    // Retrieve the signal plus background model from the workspace
    RooStats::ModelConfig* sbModel = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!sbModel) {
        throw std::runtime_error("[HistFitModule] Cannot retrieve ModelConfig from workspace");
    }
    return sbModel;
}

RooStats::ModelConfig* HistFitModule::GetBOnlyModel(RooWorkspace* ws) const
{
    RooStats::ModelConfig* sbModel = GetSPlusBModel(ws);
    RooStats::ModelConfig* bModel = new RooStats::ModelConfig(*sbModel);
    bModel->SetName("BModel");

    RooRealVar* poi = static_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
    poi->setVal(0.0);

    // snapshot only the POI at μ=0 for bModel
    RooArgSet poiSet(*poi);
    bModel->SetSnapshot(poiSet);

    return bModel;
}

double HistFitModule::CLsOutputToU2(double cls, double simulatedU2, double dataPOT, double signalPOT) const
{
    // Based on the output of the CLs hypothesis test, which simply gives a limit 
    // on the ratio of signal strength to that in the simulated signal sample, convert
    // this to a limit on U^2.
    double POTratio = signalPOT / dataPOT;
    double UsquaredLimit = sqrt(cls * POTratio) * simulatedU2;
    return UsquaredLimit;
}

Long64_t HistFitModule::EntryCount() const
{
    return 1; // Dummy, nothing per-event
}


//------------------------------------------------------------------------------

void HistFitModule::Initialise()
{
    RNodes = BuildDataFrames(fInputFiles, fTreeName, fTestFractions);

    // Create new sample weights based on the fact that we might be only using part of the samples

    std::vector<double> sampleWeights;
    for (size_t i = 0; i < RNodes.size(); ++i) {
        sampleWeights.push_back(fSampleWeights[i]/std::abs(fTestFractions[i]));
    }

    // Raw BDT scores are saved as "bdt_score", with a range of -1 to 1
    // We apply the logit transformation to spread out high bdt scores

    for (int i=0; i<RNodes.size(); i++){
        RNodes[i] = RNodes[i].Define("logit_bdt",
            [](float score) {
                float s = (score + 1.0f) / 2.0f; //rescale from [-1,1] to [0,1]
                return std::log(s / (1.0f - s));
            },
            {"bdt_score"});
    }

    //Want to initialise histograms based on the maximum value of overlay/data
    
    double maxBkg = 10.0; //initial high value
    for (int i=1; i < 5; ++i){
        auto maxInRNode = RNodes[i].Max("logit_bdt").GetValue();
        if (maxInRNode < maxBkg) maxBkg = maxInRNode;
    }

    std::vector<TH1D> bdtScoreVec;
    for (size_t i = 0; i < RNodes.size(); ++i)
        bdtScoreVec.push_back(
            Plotter::CreateTH1DFromRNode(
                RNodes[i],
                ("logit_bdt_score_" + fSampleLabels[i]).c_str(),
                "logit_bdt", 
                "Logit BDT Score",
                "Count",
                9.0, -5.0, 4.0,
                false, // removeVectorDuplicates
                true));   // createOverFlowBin

    // TEST OF CLS INFRASTRUCTURE - manually set data histogram to be the sum of bkgd histograms
    
    TRandom3 rng(0); //  
    
    double signalStrength = 0.0;

    TH1D fakeDataHist("fakeDataHist", "Fake Data Histogram", 10, -5.0, 5.0f);
    for (int i = 1; i < 11; ++i){
        double fakeDataBinEntry = 0.0;
        double fakeDataBinError = 0.0;
        double weightedBkgBin = 0.0;
        double SignalBin = 0.0;
        double weightedSignalBin = 0.0;
        for (int j = 0; j < 3; ++j){
            double rawBkgBin = 0.0;
            rawBkgBin = bdtScoreVec[j].GetBinContent(i);
            SignalBin = bdtScoreVec[bdtScoreVec.size() - 2].GetBinContent(i);
            SignalBin *= signalStrength; // Scale signal by some strength
            std::cout << "Raw bkg bin content for bin " << i << " of sample " << j << ": " << rawBkgBin << std::endl;
            weightedBkgBin = rawBkgBin * sampleWeights[j];
            fakeDataBinEntry += (weightedBkgBin + SignalBin);
        }
        fakeDataBinError = sqrt(fakeDataBinEntry); // Poisson errors
        //Random gaussian fluctuation to add

        double fluctuation = 1.0 *rng.Gaus(0.0, fakeDataBinError);
        fakeDataBinEntry += fluctuation;
        std::cout << "Bin " << i << ": fake data bin entry before fluctuation: " << fakeDataBinEntry - fluctuation << ", after fluctuation: " << fakeDataBinEntry << std::endl;
        if (i == 10) std::cout << "Fluctuation in last bin: " << fluctuation << std::endl;
        if (fakeDataBinEntry < 0.0) fakeDataBinEntry = 0.0; // No negative entries
        fakeDataHist.SetBinContent(i, fakeDataBinEntry);
        //fakeDataHist.SetBinError(i, fakeDataBinError);
        fakeDataHist.SetBinError(i, 0.002);
    }
    // Overwrite the data histogram in the vector
    //bdtScoreVec[bdtScoreVec.size() - 1] = fakeDataHist;

    // End of test code
             
    HistFitModule::SaveHistograms(bdtScoreVec, fSampleLabels, sampleWeights, "bdt_score_histograms_tmp_4.root");

    //std::vector<double> placeholderweights = {1.0, 1.0, 1.0, 1.0, 1.0};
    Plotter::FullDataMCSignalPlot(bdtScoreVec,
                        fSampleLabels,
                        "bdt_score_full_hist_tmva_histfitmodule",
                        false, // logy
                        sampleWeights);

    // Save histograms to temp file to match implementation in the RooFit examples. Could maybe
    // be done directly in memory but right now it's nice to verify you're passing in correctly
    // weighted histograms by saving them with weights applied, and this allows you to manually
    // inspect the saved histograms and errors too.
    


    std::unique_ptr<RooWorkspace> ws = BuildModelWorkspace(bdtScoreVec, fSampleLabels, "bdt_score_histograms_tmp_4.root");
    //RooStats::ModelConfig* sbModel = GetSPlusBModel(ws.get());
    //RooStats::ModelConfig* bModel = GetBOnlyModel(ws.get());
    

    // Observed data
    //RooAbsData *data = ws->data("obsData");
    //if (!data) throw std::runtime_error("[HistFitModule] Cannot retrieve observed data from workspace");

    //auto *poi = static_cast<RooRealVar*>(sbModel->GetParametersOfInterest()->first());
    //if (!poi) throw std::runtime_error("[HistFitModule] Cannot retrieve POI from ModelConfig");

    // Set up the hypothesis test calculator
    //RooStats::AsymptoticCalculator ac(*data, *bModel, *sbModel);
    //ac.SetOneSided(true);

    //RooStats::HypoTestInverter calc(ac);
    //calc.SetConfidenceLevel(0.95);
    //calc.UseCLs(true);
    //calc.SetVerbose(false);

    // Scan from 0 to, say, 0.07
    //const double muMin   = 0.0;
    //const double muMax   = 0.07;
    //const int    nPoints = 60;
    //calc.SetFixedScan(nPoints, muMin, muMax);
   
    //std::unique_ptr<RooStats::HypoTestInverterResult> htres{calc.GetInterval()};
        //if (!htres) {
        //    std::cerr << "ERROR: HypoTestInverterResult is null\n";
        //    return;
        //}

        //double upperLimit     = htres->UpperLimit();
        //double upperLimitErr  = htres->UpperLimitEstimatedError();

        //std::cout << "\n====================================================\n";
        //std::cout << " Asymptotic CLs 95% upper limit on "
        //        << poi->GetName() << " in [" << muMin << ", " << muMax << "] : "
        //        << upperLimit << " +/- " << upperLimitErr << "\n";
        //std::cout << "====================================================\n\n";

        // --- Optional: plot CLs vs SigXsecOverSim ---

        //TCanvas *c_cls = new TCanvas("c_cls", "CLs vs SigXsecOverSim", 800, 600);
        //RooStats::HypoTestInverterPlot *clsPlot =
        //    new RooStats::HypoTestInverterPlot("clsPlot", "CLs scan", htres.get());

        // "CLs" plots CLs(µ); you can also use "CLb2CLs", "CLb", "CLs+b"
        //clsPlot->Draw("CLs");
        //c_cls->SetLogy();
        //c_cls->Update();
        //c_cls->SaveAs("cls_plot_histfitmodule.png");

        // All the above is legacy code from the example, kept for reference for now.

    ws->Print();
    RooAbsData* data = ws->data("obsData");
    RooStats::ModelConfig* sbModel = (RooStats::ModelConfig*) ws->obj("ModelConfig");
    RooStats::ModelConfig* bModel = (RooStats::ModelConfig*) sbModel->Clone("BonlyModel");
    RooRealVar* poi = (RooRealVar*) bModel->GetParametersOfInterest()->first();
    poi->setVal(0);
    bModel->SetSnapshot(*poi);

    RooStats::AsymptoticCalculator  asympCalc(*data, *bModel, *sbModel);
    asympCalc.SetOneSided(true);

    //RooStats::FrequentistCalculator  freqCalc(*data, *bModel, *sbModel);

    RooStats::HypoTestInverter inverter(asympCalc);

    inverter.SetConfidenceLevel(0.95);
    inverter.UseCLs(true);  
    inverter.SetVerbose(false);
    inverter.SetFixedScan(60, 0.0, 0.01);
        
    RooStats::HypoTestInverterResult* result =  inverter.GetInterval();

    std::cout << 100*inverter.ConfidenceLevel() << "%  upper limit : " << result->UpperLimit() << std::endl;

    std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
    std::cout << " expected limit (median) " << result->GetExpectedUpperLimit(0) << std::endl;
    std::cout << " expected limit (-1 sig) " << result->GetExpectedUpperLimit(-1) << std::endl;
    std::cout << " expected limit (+1 sig) " << result->GetExpectedUpperLimit(1) << std::endl;
    std::cout << " expected limit (-2 sig) " << result->GetExpectedUpperLimit(-2) << std::endl;
    std::cout << " expected limit (+2 sig) " << result->GetExpectedUpperLimit(2) << std::endl;

    std::cout << "Converting to U^2 limits: " << std::endl;
    double clsU2 = CLsOutputToU2(result->UpperLimit(), fSimulatedSignalU2, fDataPOT, fSignalPOT);
    std::cout << " Observed U^2 limit: " << clsU2 << std::endl;

    TCanvas* c_limit = new TCanvas("c_limit", "HypoTestInverter Result", 800, 600);
    RooStats::HypoTestInverterPlot* plot = new RooStats::HypoTestInverterPlot("HTI_Result_Plot","HypoTest Scan Result",result);
    plot->Draw("CLb 2CL");  // plot also CLb and CLs+b
    c_limit->SetLogy();
    c_limit->Draw();
    c_limit->SaveAs("hypotestinverter_result_histfitmodule.png");
}

void HistFitModule::Finalise()
{
    // Nothing here
}
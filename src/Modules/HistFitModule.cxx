#include "Modules/HistFitModule.hxx"
#include "Utils/Plotter.hxx"

#include <TEnv.h>
#include <TFile.h>
#include <TString.h>
#include <algorithm>
#include <cctype>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
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

using namespace Analysis;

namespace {

std::string TrimCopy(const std::string& value)
{
    const auto first = value.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";

    const auto last = value.find_last_not_of(" \t\n\r");
    return value.substr(first, last - first + 1);
}

std::string StripOptionalQuotes(std::string value)
{
    if (value.size() >= 2 &&
        ((value.front() == '"'  && value.back() == '"') ||
         (value.front() == '\'' && value.back() == '\''))) {
        value = value.substr(1, value.size() - 2);
    }
    return value;
}

std::vector<std::string> ParseSampleLabels(const std::string& labelsString)
{
    std::vector<std::string> labels;

    if (labelsString.find(',') != std::string::npos) {
        std::stringstream ssLabels{labelsString};
        std::string labelToken;
        while (std::getline(ssLabels, labelToken, ',')) {
            labelToken = StripOptionalQuotes(TrimCopy(labelToken));
            if (!labelToken.empty()) {
                labels.push_back(labelToken);
            }
        }
    } else {
        std::stringstream ssLabels{labelsString};
        std::string labelToken;
        while (ssLabels >> labelToken) {
            labelToken = StripOptionalQuotes(TrimCopy(labelToken));
            if (!labelToken.empty()) {
                labels.push_back(labelToken);
            }
        }
    }

    return labels;
}

} // namespace

HistFitModule::HistFitModule(const TEnv& cfg)
    : Module(cfg)
    , fTreeName   (cfg.GetValue("HistFitModule.TreeName", "nuselection/NeutrinoSelectionFilter"))
    , fDataPOT    (cfg.GetValue("HistFitModule.DataPOT", 1.0e20))
    , fSignalPOT  (cfg.GetValue("HistFitModule.SignalPOT", 1.0e20)) // We only need data and signal POT because the backgrounds are scaled to data POT anyway
    , fSimulatedSignalU2(cfg.GetValue("HistFitModule.SimulatedSignalU2", 1.0e-4))
    , fBlindData  (cfg.GetValue("HistFitModule.BlindData", false))
    , fRateScaling(cfg.GetValue("HistFitModule.RateScaling", 1.0))
{

    std::stringstream ssInput{cfg.GetValue("HistFitModule.InputFiles", "")};
    std::string inputItem;
    while (ssInput >> inputItem) {
        if (inputItem.back()==',') inputItem.pop_back();
        fInputFiles.push_back(inputItem);
    }

    fSampleLabels = ParseSampleLabels(cfg.GetValue("HistFitModule.SampleLabels", ""));

    const std::string sampleTypesKey = "HistFitModule.SampleTypes";
    fSampleTypes = ParseSampleTypes(RequireConfigValue(cfg, sampleTypesKey), sampleTypesKey);

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

    const std::size_t nSamples = fInputFiles.size();
    if (nSamples == 0) {
        throw std::runtime_error("[HistFitModule] No input files specified.");
    }
    if (fSampleLabels.size() != nSamples) {
        throw std::runtime_error("[HistFitModule] SampleLabels count (" + std::to_string(fSampleLabels.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }
    if (fSampleTypes.size() != nSamples) {
        throw std::runtime_error("[HistFitModule] SampleTypes count (" + std::to_string(fSampleTypes.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }
    if (fSampleWeights.size() != nSamples) {
        throw std::runtime_error("[HistFitModule] SampleWeights count (" + std::to_string(fSampleWeights.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }
    if (!fTestFractions.empty() && fTestFractions.size() != nSamples) {
        throw std::runtime_error("[HistFitModule] TestFractions count (" + std::to_string(fTestFractions.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }
    if (fRateScaling <= 0.0) {
        throw std::runtime_error("[HistFitModule] RateScaling must be positive.");
    }

    for (std::size_t i = 0; i < fTestFractions.size(); ++i) {
        const double testFraction = fTestFractions[i];
        if (testFraction == 0.0) {
            throw std::runtime_error("[HistFitModule] TestFractions entry " + std::to_string(i)
                                     + " is zero; use a non-zero fraction or omit TestFractions.");
        }
        if (std::abs(testFraction) > 1.0) {
            throw std::runtime_error("[HistFitModule] TestFractions entry " + std::to_string(i)
                                     + " has magnitude greater than 1.");
        }
    }

    const auto nDataSamples = std::count_if(
        fSampleTypes.begin(),
        fSampleTypes.end(),
        [](SampleType type) { return IsDataSample(type); });
    const auto nSignalSamples = std::count_if(
        fSampleTypes.begin(),
        fSampleTypes.end(),
        [](SampleType type) { return IsSignalSample(type); });
    const auto nFitBackgroundSamples = std::count_if(
        fSampleTypes.begin(),
        fSampleTypes.end(),
        [this](SampleType type) { return IsFitBackground(type); });
    const auto nDetVarCVSamples = std::count_if(
        fSampleTypes.begin(),
        fSampleTypes.end(),
        [](SampleType type) { return IsDetectorVariationCVSample(type); });
    const auto nDetVarSamples = std::count_if(
        fSampleTypes.begin(),
        fSampleTypes.end(),
        [](SampleType type) { return IsDetectorVariationSample(type); });

    if (nDataSamples != 1) {
        throw std::runtime_error("[HistFitModule] SampleTypes must contain exactly one data sample.");
    }
    if (nSignalSamples == 0) {
        throw std::runtime_error("[HistFitModule] SampleTypes must contain at least one signal sample.");
    }
    if (nFitBackgroundSamples == 0) {
        throw std::runtime_error("[HistFitModule] SampleTypes must contain at least one beamoff, overlay, or dirt background sample.");
    }
    if (nDetVarSamples > 0 && nDetVarCVSamples == 0) {
        throw std::runtime_error("[HistFitModule] SampleTypes contains detvar samples but no detvarcv sample.");
    }
    if (nDetVarCVSamples > 1) {
        throw std::runtime_error("[HistFitModule] SampleTypes must contain at most one detvarcv sample.");
    }
    if (nDetVarCVSamples == 1 && nDetVarSamples == 0) {
        std::cout << "[HistFitModule] Warning: detvarcv sample configured without detvar samples; "
                  << "detector variation systematics will be skipped in stage 1.\n";
    }
}


//------------------------------------------------------------------------------
std::vector<ROOT::RDF::RNode>
HistFitModule::BuildDataFrames(const std::vector<std::string>& files,
                               const std::string& treeName) const
{
    std::vector<ROOT::RDF::RNode> nodes;
    std::cout << "[HistFitModule] Building DataFrames for " << files.size() << " input files\n";
    nodes.reserve(files.size());

    for (std::size_t i = 0; i < files.size(); ++i) {
        const auto& fname = files[i];
        auto file = TFile::Open(fname.c_str());
        if (!file || file->IsZombie()) {
            throw std::runtime_error("[HistFitModule] Cannot open file: " + fname);
        }
        auto tree = file->Get<TTree>(treeName.c_str());
        if (!tree) {
            throw std::runtime_error("[HistFitModule] Cannot find tree: " + treeName);
        }
        auto RDF = std::make_unique<ROOT::RDataFrame>(*tree);
        
        const double testFraction = EffectiveTestFraction(i);
        const auto nEvents = RDF->Count().GetValue();
        const auto nKeep = static_cast<ULong64_t>(nEvents * std::abs(testFraction));

        // Depending on whether we want to keep the first or last fraction, the sign of testFraction is different
        if (testFraction > 0) {
            nodes.push_back(RDF->Range(0, nKeep)); // Take first nKeep events
        } else {
            nodes.push_back(RDF->Range(nEvents - nKeep, nEvents)); // Take last nKeep events
        }
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
        throw std::runtime_error("[HistFitModule] Number of histograms and labels do not match in SaveHistograms");
    }
    if (weights.size() != hists.size()) {
        throw std::runtime_error("[HistFitModule] Number of histograms and weights do not match in SaveHistograms");
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

        h.Scale(weights[i]);

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
    const std::vector<std::string>& histNames,
    const std::vector<SampleType>& sampleTypes,
    const std::string& inputFile) const
{
    if (histVec.size() != histNames.size()) {
        throw std::runtime_error("[HistFitModule] Histogram count does not match histogram name count.");
    }
    if (histVec.size() != sampleTypes.size()) {
        throw std::runtime_error("[HistFitModule] Histogram count does not match sample type count.");
    }

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

    // Build a RooStats HistFactory model from the input histograms, creating s and s+b models.
    // Samples are selected by SampleTypes rather than by position in the input list.
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
    //meas.SetBinLow(16); // Investigate impact of this
    //int nBins = histVec[0].GetNbinsX();
    //meas.SetBinHigh(nBins); // Investigate impact of this

    // Create a channel

    RooStats::HistFactory::Channel chan("channel1");

    auto dataIt = std::find_if(
        sampleTypes.begin(),
        sampleTypes.end(),
        [](SampleType type) { return IsDataSample(type); });
    if (dataIt == sampleTypes.end()) {
        throw std::runtime_error("[HistFitModule] Cannot build model without a data histogram.");
    }

    const std::size_t dataIndex = static_cast<std::size_t>(std::distance(sampleTypes.begin(), dataIt));
    chan.SetData(histNames[dataIndex], inputFile);
    std::cout << "[HistFitModule] Data sample set: " << histNames[dataIndex] << std::endl;
    chan.SetStatErrorConfig(0.02, "Poisson"); // Investigate impact of this

    for (std::size_t i = 0; i < histVec.size(); ++i) {
        const SampleType sampleType = sampleTypes[i];
        if (IsDataSample(sampleType)) {
            continue;
        }
        if (IsDetectorVariationInputSample(sampleType)) {
            std::cout << "[HistFitModule] Skipping detector variation sample in stage-1 model: "
                      << histNames[i] << std::endl;
            continue;
        }

        RooStats::HistFactory::Sample sample(histNames[i], histNames[i], inputFile);
        if (IsSignalSample(sampleType)) {
            sample.AddNormFactor("SigXsecOverSim", 1, 0, 0.01);
            std::cout << "[HistFitModule] Added signal sample: " << histNames[i] << std::endl;
        } else if (IsFitBackground(sampleType)) {
            sample.ActivateStatError();
            std::cout << "[HistFitModule] Added background sample: " << histNames[i]
                      << " (" << SampleTypeName(sampleType) << ")" << std::endl;
        } else {
            throw std::runtime_error("[HistFitModule] Unsupported sample type in HistFactory model: "
                                     + SampleTypeName(sampleType));
        }

        chan.AddSample(sample);
    }

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

double HistFitModule::EffectiveTestFraction(std::size_t sampleIndex) const
{
    if (sampleIndex >= fInputFiles.size()) {
        throw std::runtime_error("[HistFitModule] EffectiveTestFraction sample index out of range.");
    }
    return fTestFractions.empty() ? 1.0 : fTestFractions[sampleIndex];
}

double HistFitModule::EffectiveSampleWeight(std::size_t sampleIndex) const
{
    if (sampleIndex >= fSampleWeights.size()) {
        throw std::runtime_error("[HistFitModule] EffectiveSampleWeight sample index out of range.");
    }
    return fSampleWeights[sampleIndex] / std::abs(EffectiveTestFraction(sampleIndex));
}

bool HistFitModule::IsFitBackground(SampleType type) const
{
    return IsOverlaySample(type) || IsDirtSample(type) || type == SampleType::BeamOff;
}

std::string HistFitModule::SanitiseHistName(const std::string& label) const
{
    std::string safeName;
    safeName.reserve(label.size());

    for (unsigned char ch : label) {
        safeName.push_back(std::isalnum(ch) || ch == '_' ? static_cast<char>(ch) : '_');
    }

    if (safeName.empty()) {
        safeName = "sample";
    }
    if (std::isdigit(static_cast<unsigned char>(safeName.front()))) {
        safeName = "h_" + safeName;
    }

    return safeName;
}

double HistFitModule::BasicSensitivityEstimate(const std::vector<TH1D>& bdtScoreVec,
    const std::vector<SampleType>& sampleTypes,
    const std::vector<double>& sampleWeights,
    double signalBinThreshold) const
{
    // Simple sensitivity estimate based on S/sqrt(B) in the high scoring region of the BDT score histogram
    if (bdtScoreVec.size() != sampleTypes.size() || bdtScoreVec.size() != sampleWeights.size()) {
        throw std::runtime_error("[HistFitModule] BasicSensitivityEstimate input vector sizes do not match.");
    }

    double signalCount = 0.0;
    double bkgCount = 0.0;

    for (size_t i = 0; i < bdtScoreVec.size(); ++i) {
        const TH1D& hist = bdtScoreVec[i];
        const SampleType sampleType = sampleTypes[i];
        const double weight = sampleWeights[i];

        // Sum entries above the threshold
        for (int bin = 1; bin <= hist.GetNbinsX(); ++bin) {
            double binCenter = hist.GetBinCenter(bin);
            if (binCenter >= signalBinThreshold) {
                double binContent = hist.GetBinContent(bin) * weight;
                if (IsSignalSample(sampleType)) {
                    signalCount += binContent;
                } else if (IsFitBackground(sampleType)) {
                    bkgCount += binContent;
                }
            }
        }
    }

    std::cout << "Basic sensitivity estimate: Signal count = " << signalCount
              << ", Background count = " << bkgCount << std::endl;
    if (bkgCount <= 0.0) {
        std::cout << "[HistFitModule] Background count is zero; returning sensitivity estimate of 0." << std::endl;
        return 0.0;
    }
    double sensitivity = signalCount / std::sqrt(bkgCount);
    return sensitivity;

}

//------------------------------------------------------------------------------

void HistFitModule::Initialise()
{
    RNodes = BuildDataFrames(fInputFiles, fTreeName);

    std::vector<std::string> histNames;
    histNames.reserve(fSampleLabels.size());
    std::unordered_map<std::string, int> histNameCounts;

    for (const auto& label : fSampleLabels) {
        const std::string baseName = SanitiseHistName(label);
        const int duplicateIndex = histNameCounts[baseName]++;
        histNames.push_back(duplicateIndex == 0 ? baseName : baseName + "_" + std::to_string(duplicateIndex));
    }

    // Raw BDT scores are saved as "bdt_score", with a range of -1 to 1
    // We apply the logit transformation to spread out high bdt scores

    std::vector<TH1D> bdtScoreVec;
    std::vector<std::string> fitLabels;
    std::vector<std::string> fitHistNames;
    std::vector<SampleType> fitSampleTypes;
    std::vector<double> fitSampleWeights;

    for (std::size_t i = 0; i < RNodes.size(); ++i) {
        if (IsDetectorVariationInputSample(fSampleTypes[i])) {
            std::cout << "[HistFitModule] Skipping detector variation input in stage-1 nominal fit: "
                      << fSampleLabels[i] << " (" << SampleTypeName(fSampleTypes[i]) << ")\n";
            continue;
        }

        RNodes[i] = RNodes[i].Define("logit_bdt",
            [](float score) {
                //float s = (score + 1.0f) / 2.0f; //rescale from [-1,1] to [0,1]
                return std::log(score / (1.0f - score));
            },
            {"bdt_score"});

        bdtScoreVec.push_back(
            Plotter::CreateTH1DFromRNode(
                RNodes[i],
                ("logit_bdt_score_" + histNames[i]).c_str(),
                "logit_bdt",
                "Logit BDT Score",
                "Count",
                10.0, -5.0, 5.0,
                false, // removeVectorDuplicates
                true));   // createOverFlowBin

        fitLabels.push_back(fSampleLabels[i]);
        fitHistNames.push_back(histNames[i]);
        fitSampleTypes.push_back(fSampleTypes[i]);
        fitSampleWeights.push_back(EffectiveSampleWeight(i));
    }

    if (bdtScoreVec.empty()) {
        throw std::runtime_error("[HistFitModule] No fit histograms were created.");
    }

    //Scale every element by rate scaling and every bin error by sqrt(rate scaling)
    std::vector<TH1D> bdtScoreVecRateScaled;
    for (size_t i = 0; i < bdtScoreVec.size(); ++i){
        TH1D hOriginal = bdtScoreVec[i];
        TH1D hScaled = hOriginal;
        for (int bin = 1; bin <= hOriginal.GetNbinsX(); ++bin){
            double originalBinContent = hOriginal.GetBinContent(bin);
            double originalBinError = hOriginal.GetBinError(bin);
            double scaledBinContent = originalBinContent * fRateScaling;
            double scaledBinError = originalBinError * sqrt(fRateScaling);
            hScaled.SetBinContent(bin, scaledBinContent);
            hScaled.SetBinError(bin, scaledBinError);
        }
        bdtScoreVecRateScaled.push_back(hScaled);
    }

    HistFitModule::SaveHistograms(bdtScoreVecRateScaled, fitHistNames, fitSampleWeights, "bdt_score_histograms_tmp_4.root");

    //std::vector<double> placeholderweights = {1.0, 1.0, 1.0, 1.0, 1.0};

    Plotter::BlindedMCSignalPlot(bdtScoreVecRateScaled,
                        fitLabels,
                        "bdt_score_blinded_hist_tmva_histfitmodule",
                        false, // logy
                        fitSampleWeights); // blinded data

    // Save histograms to temp file to match implementation in the RooFit examples. Could maybe
    // be done directly in memory but right now it's nice to verify you're passing in correctly
    // weighted histograms by saving them with weights applied, and this allows you to manually
    // inspect the saved histograms and errors too.


    double sensitivity = BasicSensitivityEstimate(bdtScoreVec,
        fitSampleTypes,
        fitSampleWeights,
        3.0); // signal bin threshold

    std::cout << "Estimated basic sensitivity (S/sqrt(B)) in logit BDT > 3.0 region: " << sensitivity << std::endl;

    std::unique_ptr<RooWorkspace> ws = BuildModelWorkspace(bdtScoreVecRateScaled, fitHistNames, fitSampleTypes, "bdt_score_histograms_tmp_4.root");

    std::cout << "HypoTestInverter starting..." << std::endl;

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

    if (!fBlindData){
        std::cout << 100*inverter.ConfidenceLevel() << "%  upper limit : " << result->UpperLimit() << std::endl;
    }
    else{
        std::cout << "Data is blinded, not showing observed limit." << std::endl;
    }

    std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
    std::cout << " expected limit (median) " << result->GetExpectedUpperLimit(0) << std::endl;
    std::cout << " expected limit (-1 sig) " << result->GetExpectedUpperLimit(-1) << std::endl;
    std::cout << " expected limit (+1 sig) " << result->GetExpectedUpperLimit(1) << std::endl;
    std::cout << " expected limit (-2 sig) " << result->GetExpectedUpperLimit(-2) << std::endl;
    std::cout << " expected limit (+2 sig) " << result->GetExpectedUpperLimit(2) << std::endl;

    if (!fBlindData){
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
    else{
        //Median expected limit plot only
        double expectedLimit = result->GetExpectedUpperLimit(0);        
        double expectedLimitMinus1Sigma = result->GetExpectedUpperLimit(-1);
        double expectedLimitPlus1Sigma = result->GetExpectedUpperLimit(1);
        double expectedLimitMinus2Sigma = result->GetExpectedUpperLimit(-2);
        double expectedLimitPlus2Sigma = result->GetExpectedUpperLimit(2);
        std::cout << "Converting expected limits to U^2: " << std::endl;
        std::cout << "Simulated POT: " << fSignalPOT << ", Data POT: " << fDataPOT << ", Simulated U^2: " << fSimulatedSignalU2 << std::endl;
        double clsU2 = CLsOutputToU2(expectedLimit, fSimulatedSignalU2, fDataPOT, fSignalPOT);
        double clsU2Minus1Sigma = CLsOutputToU2(expectedLimitMinus1Sigma, fSimulatedSignalU2, fDataPOT, fSignalPOT);
        double clsU2Plus1Sigma = CLsOutputToU2(expectedLimitPlus1Sigma, fSimulatedSignalU2, fDataPOT, fSignalPOT);
        double clsU2Minus2Sigma = CLsOutputToU2(expectedLimitMinus2Sigma, fSimulatedSignalU2, fDataPOT, fSignalPOT);
        double clsU2Plus2Sigma = CLsOutputToU2(expectedLimitPlus2Sigma, fSimulatedSignalU2, fDataPOT, fSignalPOT);
        std::cout << " Expected U^2 limit (median): " << clsU2 << std::endl;
        std::cout << " Expected U^2 limit (-1 sigma): " << clsU2Minus1Sigma << std::endl;
        std::cout << " Expected U^2 limit (+1 sigma): " << clsU2Plus1Sigma << std::endl;
        std::cout << " Expected U^2 limit (-2 sigma): " << clsU2Minus2Sigma << std::endl;
        std::cout << " Expected U^2 limit (+2 sigma): " << clsU2Plus2Sigma << std::endl;
    }
}

void HistFitModule::Finalise()
{
    // Nothing here
}

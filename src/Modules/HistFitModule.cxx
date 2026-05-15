#include "Modules/HistFitModule.hxx"
#include "Utils/Plotter.hxx"
#include "Utils/SystematicsUtil.hxx"

#include <TEnv.h>
#include <TFile.h>
#include <TString.h>
#include <algorithm>
#include <cctype>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TTree.h>
#include <TH1D.h>
#include <TChain.h>
#include <TObject.h>
#include <TVectorD.h>
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

std::string DefaultPlotWeightColumn(SampleType sampleType)
{
    return (IsOverlaySample(sampleType) || IsDirtSample(sampleType) || IsDetectorVariationInputSample(sampleType))
        ? "weight_cv"
        : "";
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
    , fPlotSystematicsDebug(cfg.GetValue("HistFitModule.PlotSystematicsDebug", false))
    , fBDTScoreMinX(cfg.GetValue("HistFitModule.BDTScoreMinX", -5.0))
    , fBDTScoreBinsBelowOverflow(cfg.GetValue("HistFitModule.BDTScoreBinsBelowOverflow", 9))
    , fBDTScoreOverflowBackgroundEvents(cfg.GetValue("HistFitModule.BDTScoreOverflowBackgroundEvents", 5.0))
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
    if (!std::isfinite(fBDTScoreMinX)) {
        throw std::runtime_error("[HistFitModule] BDTScoreMinX must be finite.");
    }
    if (fBDTScoreBinsBelowOverflow <= 0) {
        throw std::runtime_error("[HistFitModule] BDTScoreBinsBelowOverflow must be positive.");
    }
    if (!std::isfinite(fBDTScoreOverflowBackgroundEvents) || fBDTScoreOverflowBackgroundEvents <= 0.0) {
        throw std::runtime_error("[HistFitModule] BDTScoreOverflowBackgroundEvents must be finite and positive.");
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

std::vector<HistFitModule::HistoSysVariation>
HistFitModule::WriteOverlayHistoSysVariations(
    const TH1D& overlayHist,
    const std::string& overlayHistName,
    double sampleWeight,
    const TMatrixD& covariance,
    const std::string& inputFile) const
{
    const int nBins = overlayHist.GetNbinsX();
    if (covariance.GetNrows() != nBins || covariance.GetNcols() != nBins) {
        throw std::runtime_error("[HistFitModule] Overlay covariance dimensions ("
                                 + std::to_string(covariance.GetNrows()) + "x"
                                 + std::to_string(covariance.GetNcols())
                                 + ") do not match overlay histogram bins ("
                                 + std::to_string(nBins) + ").");
    }

    TH1D nominal = overlayHist;
    nominal.SetDirectory(nullptr);
    nominal.Scale(sampleWeight);

    SystematicsUtil sysUtil;
    auto [eigenVectors, eigenValues] = sysUtil.EigenDecomposition(covariance);

    double maxDiag = 0.0;
    for (int i = 0; i < covariance.GetNrows(); ++i) {
        maxDiag = std::max(maxDiag, std::abs(covariance(i, i)));
    }
    const double negativeEigenvalueTolerance = std::max(1e-9, 1e-10 * maxDiag);
    const double zeroShiftTolerance = 1e-12;

    TFile outFile(inputFile.c_str(), "UPDATE");
    if (outFile.IsZombie()) {
        throw std::runtime_error("[HistFitModule] Cannot update histogram file with HistoSys variations: "
                                 + inputFile);
    }

    std::vector<HistoSysVariation> variations;
    for (int eig = 0; eig < eigenValues.GetNrows(); ++eig) {
        double eigenValue = eigenValues(eig);
        if (!std::isfinite(eigenValue)) {
            std::cout << "[HistFitModule] Skipping non-finite overlay covariance eigenvalue "
                      << eig << ".\n";
            continue;
        }
        if (eigenValue < -negativeEigenvalueTolerance) {
            throw std::runtime_error("[HistFitModule] Overlay covariance has materially negative eigenvalue "
                                     + std::to_string(eig) + ": " + std::to_string(eigenValue));
        }
        if (eigenValue < 0.0) {
            eigenValue = 0.0;
        }

        const double shiftScale = std::sqrt(eigenValue);
        double maxAbsShift = 0.0;
        for (int bin = 0; bin < nBins; ++bin) {
            const double shift = eigenVectors(bin, eig) * shiftScale;
            if (!std::isfinite(shift)) {
                maxAbsShift = std::numeric_limits<double>::infinity();
                break;
            }
            maxAbsShift = std::max(maxAbsShift, std::abs(shift));
        }
        if (!std::isfinite(maxAbsShift) || maxAbsShift <= zeroShiftTolerance) {
            continue;
        }

        const std::string systName = "overlay_multisim_eig" + std::to_string(eig);
        const std::string lowHistName = overlayHistName + "_multisim_eig" + std::to_string(eig) + "_low";
        const std::string highHistName = overlayHistName + "_multisim_eig" + std::to_string(eig) + "_high";

        TH1D lowHist = nominal;
        TH1D highHist = nominal;
        lowHist.SetName(lowHistName.c_str());
        highHist.SetName(highHistName.c_str());
        lowHist.SetDirectory(nullptr);
        highHist.SetDirectory(nullptr);


        //Debug: work out mean of shift histogram - if negative, flip shiftScale to make "up" variation always above nominal
        double meanShift = 0.0;
        for (int bin = 1; bin <= nBins; ++bin) {
            const double shift = eigenVectors(bin - 1, eig) * shiftScale;
            meanShift += shift;
        }
        meanShift /= nBins;
        //double shiftSign = (meanShift >= 0.0) ? 1.0 : -1.0;
        double shiftSign = 1.0; // For now, do not flip sign of shift even if mean is negative, to avoid confusion with interpretation of eigenvectors

        for (int bin = 1; bin <= nBins; ++bin) {
            const double nominalBin = nominal.GetBinContent(bin);
            const double shift = eigenVectors(bin - 1, eig) * shiftScale * shiftSign;
            //Debug, set shift to 30% of nominal bin content
            //const double shift = 0.1 * nominalBin;
            lowHist.SetBinContent(bin, std::max(0.0, nominalBin - shift));
            highHist.SetBinContent(bin, std::max(0.0, nominalBin + shift));
        }

        //Debug - 
        

        outFile.cd();
        lowHist.Write("", TObject::kOverwrite);
        highHist.Write("", TObject::kOverwrite);
        variations.push_back({systName, lowHistName, highHistName});
    }

    outFile.Close();
    return variations;
}

 std::unique_ptr<RooWorkspace> HistFitModule::BuildModelWorkspace(
    std::vector<TH1D>& histVec,
    const std::vector<std::string>& histNames,
    const std::vector<SampleType>& sampleTypes,
    const std::vector<double>& sampleWeights,
    const std::optional<TMatrixD>& overlayMultisimCovariance,
    const std::string& inputFile) const
{
    if (histVec.size() != histNames.size()) {
        throw std::runtime_error("[HistFitModule] Histogram count does not match histogram name count.");
    }
    if (histVec.size() != sampleTypes.size()) {
        throw std::runtime_error("[HistFitModule] Histogram count does not match sample type count.");
    }
    if (histVec.size() != sampleWeights.size()) {
        throw std::runtime_error("[HistFitModule] Histogram count does not match sample weight count.");
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

    //Setting this to zero causes the minimisation to take forever and gives weird results, so set to a small non-zero value for now. Need to investigate further.
    meas.SetLumiRelErr(0.01);
    
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

    if (overlayMultisimCovariance) {
        const auto nOverlaySamples = std::count_if(
            sampleTypes.begin(),
            sampleTypes.end(),
            [](SampleType type) { return IsOverlaySample(type); });
        if (nOverlaySamples != 1) {
            throw std::runtime_error("[HistFitModule] Overlay covariance HistoSys construction requires exactly one overlay sample.");
        }
    }

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
            sample.AddNormFactor("SigXsecOverSim", 0.001, 0, 0.005);
            sample.ActivateStatError();
            sample.AddOverallSys("signal_norm_30pct", 0.7, 1.3);
            std::cout << "[HistFitModule] Added signal sample: " << histNames[i] << std::endl;
        } else if (IsFitBackground(sampleType)) {
            sample.ActivateStatError();
            std::cout << "[HistFitModule] Added background sample: " << histNames[i]
                      << " (" << SampleTypeName(sampleType) << ")" << std::endl;
            if (IsOverlaySample(sampleType) && overlayMultisimCovariance) {
                const auto variations = WriteOverlayHistoSysVariations(
                    histVec[i],
                    histNames[i],
                    sampleWeights[i],
                    *overlayMultisimCovariance,
                    inputFile);
                for (const auto& variation : variations) {
                    std::cout << "Dummy add" << std::endl;
                    sample.AddHistoSys(
                        variation.systName,
                        variation.lowHistName,
                        inputFile,
                        "",
                        variation.highHistName,
                        inputFile,
                        "");
                }
                std::cout << "[HistFitModule] Attached " << variations.size()
                          << " overlay multisim HistoSys variations to "
                          << histNames[i] << ".\n";
            }
        } else {
            throw std::runtime_error("[HistFitModule] Unsupported sample type in HistFactory model: "
                                     + SampleTypeName(sampleType));
        }

        //Diagnostic dump of HistoSys list
        if (IsOverlaySample(sampleType)) {
            std::cout << "[HistFitModule] Overlay sample HistoSys list:\n";
            for (const auto& hs : sample.GetHistoSysList()) {
                std::cout << "  " << hs.GetName()
                        << " low=" << hs.GetHistoNameLow()
                        << " high=" << hs.GetHistoNameHigh()
                        << " lowFile=" << hs.GetInputFileLow()
                        << " highFile=" << hs.GetInputFileHigh()
                        << "\n";
            }
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
    ws->Print("t");

    auto* mc = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (mc && mc->GetNuisanceParameters()) {
        std::cout << "[HistFitModule] Nuisance parameters:\n";
        mc->GetNuisanceParameters()->Print("v");
    }

    auto* nuis = mc ? mc->GetNuisanceParameters() : nullptr;
    if (nuis) {
        std::cout << "\n[HistFitModule] Nuisance parameter details:\n";

        TIterator* it = nuis->createIterator();
        TObject* obj = nullptr;
        while ((obj = it->Next())) {
            auto* var = dynamic_cast<RooRealVar*>(obj);
            if (!var) continue;

            std::cout << "  " << var->GetName()
                    << " value=" << var->getVal()
                    << " error=" << var->getError()
                    << " range=[" << var->getMin() << ", " << var->getMax() << "]"
                    << " constant=" << var->isConstant()
                    << "\n";
        }
        delete it;
    }

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

HistFitModule::DynamicBDTBinning
HistFitModule::ComputeDynamicBDTBinning(const std::vector<ROOT::RDF::RNode>& nodes) const
{
    std::vector<std::pair<double, double>> weightedBackgroundScores;

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (!IsFitBackground(fSampleTypes[i])) {
            continue;
        }

        const double globalScale = EffectiveSampleWeight(i) * fRateScaling;
        const std::string weightCol = DefaultPlotWeightColumn(fSampleTypes[i]);
        const std::string predictedWeightCol = "histfit_dynamic_bdt_weight_" + std::to_string(i);
        const std::string globalScaleString = std::to_string(globalScale);

        ROOT::RDF::RNode baseNode = nodes[i];
        ROOT::RDF::RNode weightedNode = weightCol.empty()
            ? baseNode.Define(predictedWeightCol, globalScaleString)
            : baseNode.Define(predictedWeightCol, "(" + weightCol + ") * (" + globalScaleString + ")");

        const auto scoreValues = weightedNode.Take<double>("logit_bdt").GetValue();
        const auto weightValues = weightedNode.Take<double>(predictedWeightCol).GetValue();
        if (scoreValues.size() != weightValues.size()) {
            throw std::runtime_error("[HistFitModule] Dynamic binning score/weight vector size mismatch.");
        }

        for (std::size_t j = 0; j < scoreValues.size(); ++j) {
            const double score = scoreValues[j];
            const double predictedWeight = weightValues[j];
            if (!std::isfinite(score) || !std::isfinite(predictedWeight) || predictedWeight <= 0.0) {
                continue;
            }
            if (score >= fBDTScoreMinX) {
                weightedBackgroundScores.push_back({score, predictedWeight});
            }
        }
    }

    if (weightedBackgroundScores.empty()) {
        throw std::runtime_error("[HistFitModule] Cannot compute dynamic BDT binning: no weighted background scores above BDTScoreMinX.");
    }

    std::sort(
        weightedBackgroundScores.begin(),
        weightedBackgroundScores.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });

    double cumulativeYield = 0.0;
    std::optional<double> overflowEdge;
    for (const auto& [score, weight] : weightedBackgroundScores) {
        cumulativeYield += weight;
        if (cumulativeYield > fBDTScoreOverflowBackgroundEvents) {
            overflowEdge = score;
            break;
        }
    }

    if (!overflowEdge) {
        throw std::runtime_error("[HistFitModule] Cannot compute dynamic BDT binning: target overflow background yield "
                                 + std::to_string(fBDTScoreOverflowBackgroundEvents)
                                 + " exceeds total predicted background yield above BDTScoreMinX.");
    }
    if (*overflowEdge <= fBDTScoreMinX) {
        throw std::runtime_error("[HistFitModule] Dynamic BDT overflow edge is not greater than BDTScoreMinX.");
    }

    double overflowYield = 0.0;
    for (const auto& [score, weight] : weightedBackgroundScores) {
        if (score >= *overflowEdge) {
            overflowYield += weight;
        } else {
            break;
        }
    }

    DynamicBDTBinning binning;
    binning.xMin = fBDTScoreMinX;
    binning.overflowEdge = *overflowEdge;
    binning.binsBelowOverflow = fBDTScoreBinsBelowOverflow;
    binning.binWidth = (binning.overflowEdge - binning.xMin) / static_cast<double>(binning.binsBelowOverflow);
    binning.xMax = binning.overflowEdge + binning.binWidth;
    binning.overflowBackgroundYield = overflowYield;
    return binning;
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
    fOverlayMultisimCovariance.reset();

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
    std::vector<std::size_t> fitSampleIndices;

    for (std::size_t i = 0; i < RNodes.size(); ++i) {
        if (IsDetectorVariationInputSample(fSampleTypes[i])) {
            std::cout << "[HistFitModule] Skipping detector variation input in stage-1 nominal fit: "
                      << fSampleLabels[i] << " (" << SampleTypeName(fSampleTypes[i]) << ")\n";
            continue;
        }

        RNodes[i] = RNodes[i].Define("logit_bdt",
            [](float score) -> double {
                //float s = (score + 1.0f) / 2.0f; //rescale from [-1,1] to [0,1]
                const double scoreDouble = static_cast<double>(score);
                return std::log(scoreDouble / (1.0 - scoreDouble));
            },
            {"bdt_score"});
        fitSampleIndices.push_back(i);
    }

    const DynamicBDTBinning bdtBinning = ComputeDynamicBDTBinning(RNodes);
    std::cout << "[HistFitModule] Dynamic BDT binning: min=" << bdtBinning.xMin
              << ", overflow edge=" << bdtBinning.overflowEdge
              << ", bins below overflow=" << bdtBinning.binsBelowOverflow
              << ", target overflow bkg=" << fBDTScoreOverflowBackgroundEvents
              << ", actual overflow bkg=" << bdtBinning.overflowBackgroundYield
              << ".\n";

    for (const std::size_t i : fitSampleIndices) {
        const double overflowEdge = bdtBinning.overflowEdge;
        const double overflowBinCenter = bdtBinning.overflowEdge + 0.5 * bdtBinning.binWidth;
        RNodes[i] = RNodes[i].Define(
            "logit_bdt_hist",
            [overflowEdge, overflowBinCenter](double score) {
                return score >= overflowEdge ? overflowBinCenter : score;
            },
            {"logit_bdt"});

        const std::string weightCol = DefaultPlotWeightColumn(fSampleTypes[i]);

        bdtScoreVec.push_back(
            Plotter::CreateTH1DFromRNode(
                RNodes[i],
                ("logit_bdt_score_" + histNames[i]).c_str(),
                "logit_bdt_hist",
                "Logit BDT Score",
                "Count",
                bdtBinning.binsBelowOverflow + 1,
                bdtBinning.xMin,
                bdtBinning.xMax,
                false, // removeVectorDuplicates
                false, // createOverFlowBin
                weightCol));

        fitLabels.push_back(fSampleLabels[i]);
        fitHistNames.push_back(histNames[i]);
        fitSampleTypes.push_back(fSampleTypes[i]);
        fitSampleWeights.push_back(EffectiveSampleWeight(i));

        if (IsOverlaySample(fSampleTypes[i])) {
            const double overlayGlobalScale = EffectiveSampleWeight(i) * fRateScaling;
            TH1D overlayNominalForCov = bdtScoreVec.back();
            overlayNominalForCov.SetDirectory(nullptr);
            overlayNominalForCov.Scale(overlayGlobalScale);

            SystematicsUtil sysUtil;
            std::vector<TH1D> genieUniverses = sysUtil.createMultiSimUniverses(
                overlayNominalForCov,
                RNodes[i],
                "logit_bdt_hist",
                systConfig.genieMultisimBranch,
                systConfig.genieCVWeightBranch,
                systConfig.genieGlobalCVWeightBranch,
                overlayGlobalScale);
            std::vector<TH1D> ppfxUniverses = sysUtil.createMultiSimUniverses(
                overlayNominalForCov,
                RNodes[i],
                "logit_bdt_hist",
                systConfig.ppfxMultisimBranch,
                systConfig.ppfxCVWeightBranch,
                systConfig.ppfxGlobalCVWeightBranch,
                overlayGlobalScale);
            std::vector<TH1D> reintUniverses = sysUtil.createMultiSimUniverses(
                overlayNominalForCov,
                RNodes[i],
                "logit_bdt_hist",
                systConfig.reintMultisimBranch,
                systConfig.reintCVWeightBranch,
                systConfig.reintGlobalCVWeightBranch,
                overlayGlobalScale);

            TMatrixD genieCov = sysUtil.covarianceMatrixFromMultisims(genieUniverses, overlayNominalForCov);
            TMatrixD ppfxCov = sysUtil.covarianceMatrixFromMultisims(ppfxUniverses, overlayNominalForCov);
            TMatrixD reintCov = sysUtil.covarianceMatrixFromMultisims(reintUniverses, overlayNominalForCov);

            fOverlayMultisimCovariance = genieCov + ppfxCov + reintCov;

            //Test case of zero covariance to verify that we can handle this without issue
            //fOverlayMultisimCovariance = TMatrixD(genieCov.GetNrows(), genieCov.GetNcols());
            //for (int row = 0; row < fOverlayMultisimCovariance->GetNrows(); ++row) {
            //    for (int col = 0; col < fOverlayMultisimCovariance->GetNcols(); ++col) {
            //        (*fOverlayMultisimCovariance)(row, col) = 0.01;
            //    }
            //}

            std::cout << "[HistFitModule] Built overlay multisim covariance for "
                      << fSampleLabels[i] << " with dimensions "
                      << fOverlayMultisimCovariance->GetNrows() << "x"
                      << fOverlayMultisimCovariance->GetNcols() << ".\n";

            if (fPlotSystematicsDebug) {
                sysUtil.PlotMatrix(*fOverlayMultisimCovariance, "histfit_overlay_multisim_cov");
                sysUtil.PlotFractionalCovarianceMatrix(
                    *fOverlayMultisimCovariance,
                    overlayNominalForCov,
                    "histfit_overlay_multisim_frac_cov");
            }
        }
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
                        fitSampleWeights);

    // Save histograms to temp file to match implementation in the RooFit examples. Could maybe
    // be done directly in memory but right now it's nice to verify you're passing in correctly
    // weighted histograms by saving them with weights applied, and this allows you to manually
    // inspect the saved histograms and errors too.


    double sensitivity = BasicSensitivityEstimate(bdtScoreVec,
        fitSampleTypes,
        fitSampleWeights,
        3.0); // signal bin threshold

    std::cout << "Estimated basic sensitivity (S/sqrt(B)) in logit BDT > 3.0 region: " << sensitivity << std::endl;

    std::unique_ptr<RooWorkspace> ws = BuildModelWorkspace(
        bdtScoreVecRateScaled,
        fitHistNames,
        fitSampleTypes,
        fitSampleWeights,
        fOverlayMultisimCovariance,
        "bdt_score_histograms_tmp_4.root");

    std::cout << "HypoTestInverter starting..." << std::endl;

    ws->Print();
    RooAbsData* data = ws->data("obsData");

    //Legacy code ----------------------
    //RooStats::ModelConfig* sbModel = (RooStats::ModelConfig*) ws->obj("ModelConfig");
    //RooStats::ModelConfig* bModel = (RooStats::ModelConfig*) sbModel->Clone("BonlyModel");
    //RooRealVar* poi = (RooRealVar*) bModel->GetParametersOfInterest()->first();
    //poi->setVal(0);
    //bModel->SetSnapshot(*poi);

    ////RooStats::AsymptoticCalculator  asympCalc(*data, *bModel, *sbModel);
    //RooStats::AsymptoticCalculator asympCalc(*data, *bModel, *sbModel, true);
    //asympCalc.SetOneSided(true);
    //End of legacy code-----------------------
    
    RooStats::ModelConfig* sbModel =
    static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));

    auto* poi =
        dynamic_cast<RooRealVar*>(sbModel->GetParametersOfInterest()->first());

    if (!poi) {
        throw std::runtime_error("POI is not a RooRealVar");
    }

    // Save nominal S+B POI snapshot before constructing B-only snapshot.
    poi->setVal(0.001);
    sbModel->SetSnapshot(RooArgSet(*poi));

    RooStats::ModelConfig* bModel =
        static_cast<RooStats::ModelConfig*>(sbModel->Clone("BonlyModel"));

    auto* poiB =
        dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());

    const double oldVal = poiB->getVal();
    poiB->setVal(0.0);
    bModel->SetSnapshot(RooArgSet(*poiB));
    poiB->setVal(oldVal);

    // For nominal expected limits:
    RooStats::AsymptoticCalculator asympCalc(*data, *bModel, *sbModel, true);
    asympCalc.SetOneSided(true);

    // Dignostic output
    //asympCalc.Initialize();
    //std::cout << "\nBest-fit parameters used by AsymptoticCalculator:\n";
    //asympCalc.GetBestFitParams().Print("v");
    //std::cout << "\nBest-fit POI:\n";
    //asympCalc.GetBestFitPoi().Print("v");

    //RooStats::FrequentistCalculator  freqCalc(*data, *bModel, *sbModel);

    RooStats::HypoTestInverter inverter(asympCalc);

    inverter.SetConfidenceLevel(0.95);
    inverter.UseCLs(true);  
    inverter.SetVerbose(true);
    inverter.SetFixedScan(60, 0.0, 0.005);
        
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

#ifndef HISTFIT_MODULE_HXX
#define HISTFIT_MODULE_HXX

#include "Framework/Module.hxx"
#include "Utils/Plotter.hxx"
#include "Utils/ConfigUtils.hxx"

#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <memory>
#include <optional>
#include <unordered_map>
#include <string>
#include <vector>
#include <TMatrixD.h>

#include <RooStats/HistFactory/MakeModelAndMeasurementsFast.h>
#include <RooStats/HistFactory/Measurement.h>

#include <RooStats/AsymptoticCalculator.h>
#include <RooStats/HypoTestInverter.h>
#include <RooStats/HypoTestInverterResult.h>
#include <RooStats/ProfileLikelihoodTestStat.h>
#include <RooStats/HypoTestInverterPlot.h> 

namespace Analysis {

class HistFitModule final : public Module {
public:
    explicit HistFitModule(const TEnv& cfg);

    // Module interface
    void Initialise() override;
    void Execute(Long64_t /*entry*/) override {}   // nothing per-event
    void Finalise()   override;
    Long64_t EntryCount() const override;

    std::string Name() const override { return "HistFit"; }

private:
    struct HistoSysVariation {
        std::string systName;
        std::string lowHistName;
        std::string highHistName;
    };
    struct DynamicBDTBinning {
        double xMin = 0.0;
        double overflowEdge = 0.0;
        double binWidth = 0.0;
        double xMax = 0.0;
        double overflowBackgroundYield = 0.0;
        int binsBelowOverflow = 0;
    };
    struct ChannelInput {
        std::string name;
        std::vector<std::size_t> sampleIndices;
    };
    struct ChannelFitInputs {
        std::string name;
        std::vector<TH1D> hists;
        std::vector<std::string> labels;
        std::vector<std::string> histNames;
        std::vector<SampleType> sampleTypes;
        std::vector<double> sampleWeights;
        std::optional<TMatrixD> overlayShapeCovariance;
        bool overlayShapeCovarianceIncludesDetVars = false;
    };

    // Helper: build the dataframe vector from a file list
    std::vector<ROOT::RDF::RNode> BuildDataFrames(const std::vector<std::string>& files,
                                            const std::string& treeName) const;

    void SaveHistograms(const std::vector<TH1D>& hists,
                    const std::vector<std::string>& labels,
                    const std::vector<double>& weights,
                    const std::string& fileName);

    std::unique_ptr<RooWorkspace> BuildModelWorkspace(
                                            std::vector<ChannelFitInputs>& channels,
                                            const std::string& inputFile) const;
    std::vector<HistoSysVariation> WriteOverlayHistoSysVariations(
                                            const TH1D& overlayHist,
                                            const std::string& overlayHistName,
                                            double sampleWeight,
                                            const TMatrixD& covariance,
                                            const std::string& systNamePrefix,
                                            const std::string& inputFile) const;
    std::string WriteOverlayShapeSysUncertainty(
                                            const TH1D& overlayHist,
                                            const std::string& overlayHistName,
                                            double sampleWeight,
                                            double relativeUncertainty,
                                            const std::string& inputFile) const;
    void PrintSourceFractionalUncertainties(
                                            const std::string& channelName,
                                            const std::string& sourceName,
                                            const TMatrixD& covariance,
                                            const TH1D& denominatorHist) const;
    TMatrixD FractionalCovarianceFromAbsolute(
                                            const TMatrixD& covariance,
                                            const TH1D& scaledNominalHist) const;
    TMatrixD AbsoluteCovarianceFromFractional(
                                            const TMatrixD& fractionalCovariance,
                                            const TH1D& scaledNominalHist) const;

    RooStats::ModelConfig* GetSPlusBModel(RooWorkspace* ws) const;
    RooStats::ModelConfig* GetBOnlyModel(RooWorkspace* ws) const;
    double CLsOutputToU2(double cls, double simulatedU2, double dataPOT, double signalPOT) const;
    double EffectiveTestFraction(std::size_t sampleIndex) const;
    double EffectiveSampleWeight(std::size_t sampleIndex) const;
    bool IsFitBackground(SampleType type) const;
    std::string SanitiseHistName(const std::string& label) const;
    std::vector<ChannelInput> BuildChannelInputs() const;
    void ValidateChannelInputs(const std::vector<ChannelInput>& channels) const;
    void ValidateDetVarCovarianceTransfers(const std::vector<ChannelInput>& channels) const;
    DynamicBDTBinning ComputeDynamicBDTBinning(
        const std::vector<ROOT::RDF::RNode>& nodes,
        const std::vector<std::size_t>& sampleIndices) const;
    double BasicSensitivityEstimate(const std::vector<TH1D>& bdtScoreVec,
        const std::vector<SampleType>& sampleTypes,
        const std::vector<double>& sampleWeights,
        double signalBinThreshold) const;

    /// Configuration
    std::vector<std::string>      fInputFiles;
    std::string        fTreeName;        ///< name of the input TTree
    bool fBlindData;     
    std::vector<std::string> fSampleLabels; ///< Labels for the samples, e.g. "data", "overlay", "signal"
    std::vector<SampleType> fSampleTypes; ///< Analysis role for each sample
    std::vector<std::string> fSampleChannels; ///< HistFactory channel label for each sample
    std::vector<double> fSampleWeights; ///< Weights for each sample to normalise to POT
    double fDataPOT;    ///< POT for the data sample
    double fSignalPOT;  ///< POT for the signal MC sample
    double fSimulatedSignalU2; ///< The U^2 value used in the generator
    std::vector<double> fTestFractions; ///< Fractions of events to keep for each sample (for BDT test samples)
    double fRateScaling; ///< Optional global rate scaling for histogram contents
    bool fPlotSystematicsDebug; ///< Whether to write covariance diagnostic plots
    double fBDTScoreMinX; ///< Lower edge for dynamic BDT-score histograms
    int fBDTScoreBinsBelowOverflow; ///< Number of BDT-score bins before the overflow-like bin
    double fBDTScoreOverflowBackgroundEvents; ///< Target predicted background yield in overflow-like bin
    bool fLegacySingleChannelMode; ///< True when SampleChannels is omitted and legacy names should be preserved
    std::unordered_map<std::string, std::string> fDetVarCovarianceTransfers; ///< target channel -> source channel
    std::string fDetVarCovarianceTransferMode; ///< Transfer mode for temporary detector covariance reuse


    /// Working objects
    std::vector<ROOT::RDF::RNode> RNodes; ///< DataFrames for each input file
};

} // namespace Analysis
#endif

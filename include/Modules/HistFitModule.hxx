#ifndef HISTFIT_MODULE_HXX
#define HISTFIT_MODULE_HXX

#include "Framework/Module.hxx"
#include "Utils/Plotter.hxx"

#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <memory>
#include <vector>

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

    // Helper: build the dataframe vector from a file list
    std::vector<ROOT::RDF::RNode> BuildDataFrames(const std::vector<std::string>& files,
                                            const std::string& treeName,
                                            const std::vector<double>& testFractions) const;

    void SaveHistograms(const std::vector<TH1D>& hists,
                    const std::vector<std::string>& labels,
                    const std::vector<double>& weights,
                    const std::string& fileName);

    std::unique_ptr<RooWorkspace> BuildModelWorkspace(
                                            std::vector<TH1D>& histVec,
                                            const std::vector<std::string>& labels,
                                            const std::string& inputFile) const;

    RooStats::ModelConfig* GetSPlusBModel(RooWorkspace* ws) const;
    RooStats::ModelConfig* GetBOnlyModel(RooWorkspace* ws) const;

    /// Configuration
    std::vector<std::string>      fInputFiles;
    std::string        fTreeName;        ///< name of the input TTree
    std::vector<std::string> fSampleLabels; ///< Labels for the samples, e.g. "data", "overlay", "signal"
    std::vector<double> fSampleWeights; ///< Weights for each sample to normalise to POT
    double fDataPOT;    ///< POT for the data sample
    double fSignalPOT;  ///< POT for the signal MC sample
    std::vector<double> fTestFractions; ///< Fractions of events to keep for each sample (for BDT test samples)


    /// Working objects
    std::vector<ROOT::RDF::RNode> RNodes; ///< DataFrames for each input file
};

} // namespace Analysis
#endif
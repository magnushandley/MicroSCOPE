#ifndef PRESELECTION_MODULE_HXX
#define PRESELECTION_MODULE_HXX

#include "Framework/Module.hxx"
#include "Utils/ConfigUtils.hxx"
#include "Utils/Plotter.hxx"

#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <memory>
#include <vector>

namespace Analysis {

class PreselectionModule final : public Module {
public:
    explicit PreselectionModule(const TEnv& cfg);

    // Module interface
    void Initialise() override;
    void Execute(Long64_t /*entry*/) override {}   // nothing per-event
    void Finalise()   override;
    Long64_t EntryCount() const override;

    std::string Name() const override { return "Preselection"; }

private:
    // Helper: build the input chain from a comma-separated list
    std::vector<std::unique_ptr<ROOT::RDataFrame>> BuildDataFrames(const std::vector<std::string>& files,
                                            const std::string& treeName) const;

    /// Configuration
    std::vector<std::string> fInputFiles;      ///< comma-separated list
    std::string        fTreeName;        ///< name of the input TTree
    std::vector<std::string>        fOutFiles;         ///< result file
    std::vector<std::string> cuts;       ///< List of cuts strings to apply
    std::vector<std::string> fSampleLabels; ///< Legend/log labels for each sample
    std::vector<SampleType> fSampleTypes; ///< Analysis role for each sample
    std::vector<double> fSampleWeights; ///< Weights for each sample to normalise to POT
    std::vector<std::string> fVarsToKeep;///< thin list, incl. derived vars
    std::string        fRunLabel;        ///< “numi_run4b”, …
    bool              fMakePlots;      ///< whether to create output plots
    std::vector<PlotConfig> fPlotConfigs; ///< Config-driven plot definitions

    /// Working objects
    std::unique_ptr<TChain>     fChain;
    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfVec; ///< DataFrames for each input file
};

} // namespace Analysis
#endif

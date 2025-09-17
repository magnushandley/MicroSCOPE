#ifndef PLOTTER_MODULE_HXX
#define PLOTTER_MODULE_HXX

#include "Framework/Module.hxx"
#include "Utils/Plotter.hxx"

#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <memory>
#include <vector>

namespace Analysis {

class PlotterModule final : public Module {
public:
    explicit PlotterModule(const TEnv& cfg);

    // Module interface
    void Initialise() override;
    void Execute(Long64_t /*entry*/) override {}   // nothing per-event
    void Finalise()   override;
    Long64_t EntryCount() const override;

    std::string Name() const override { return "Plotter"; }

private:
    // Helper: build the input chain from a comma-separated list
    std::vector<std::unique_ptr<ROOT::RDataFrame>> BuildDataFrames(const std::vector<std::string>& files,
                                            const std::string& treeName) const;

    /// Configuration
    std::vector<std::string>      fInputFiles;
    std::string        fTreeName;        ///< name of the input TTree
    std::vector<std::string> fSampleLabels; ///< Labels for the samples, e.g. "data", "overlay", "signal"
    std::vector<double> fSampleWeights; ///< Weights for each sample to normalise to POT
    std::vector<std::string> fVarsToKeep;///< thin list, incl. derived vars

    /// Working objects
    std::unique_ptr<TChain>     fChain;
    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfVec; ///< DataFrames for each input file
};

} // namespace Analysis
#endif
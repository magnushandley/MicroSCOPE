#ifndef SLIMMER_MODULE_HXX
#define SLIMMER_MODULE_HXX

#include "Framework/Module.hxx"
#include "Utils/Plotter.hxx"

#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <memory>
#include <vector>

namespace Analysis {

class SlimmerModule final : public Module {
public:
    explicit SlimmerModule(const TEnv& cfg);

    // Module interface
    void Initialise() override;
    void Execute(Long64_t /*entry*/) override {}   // nothing per-event
    void Finalise()   override;
    Long64_t EntryCount() const override;

    std::string Name() const override { return "Slimmer"; }

private:
    // Helper: build the input vector of RDataFrames
    std::vector<std::unique_ptr<ROOT::RDataFrame>> BuildDataFrames(const std::vector<std::string>& files,
                                            const std::string& treeName) const;

    /// Configuration
    std::vector<std::string> fInputFiles;      ///< comma-separated list
    std::string        fTreeName;        ///< name of the input TTree
    std::vector<std::string> fOutputFiles;         ///< result files
    std::vector<std::string> fVarsToKeep;///< thin list, incl. derived vars
    std::string        fRunLabel;        ///< “run_x”, …

    /// Working objects
    std::unique_ptr<TChain>     fChain;
    std::unique_ptr<ROOT::RDataFrame> fRDF;
    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfVec; ///< DataFrames for each input file
};

} // namespace Analysis
#endif
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
    // Helper: build the input chain from a comma-separated list
    std::unique_ptr<TChain> BuildInputChain(const std::string& files,
                                            const std::string& treeName) const;

    /// Configuration
    std::string        fInputFiles;      ///< comma-separated list
    std::string        fTreeName;        ///< name of the input TTree
    std::string        fOutFile;         ///< result file
    std::vector<std::string> fVarsToKeep;///< thin list, incl. derived vars
    std::string        fRunLabel;        ///< “data2023B”, “bnb_overlay”, …

    /// Working objects
    std::unique_ptr<TChain>     fChain;
    std::unique_ptr<ROOT::RDataFrame> fRDF;
};

} // namespace Analysis
#endif
#ifndef SLIMMER_MODULE_JSON_HXX
#define SLIMMER_MODULE_JSON_HXX

#include "Framework/Module.hxx"
#include "Utils/Plotter.hxx"

#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <memory>
#include <vector>

namespace Analysis {

class SlimmerModuleJSON final : public Module {
public:
    explicit SlimmerModuleJSON(const nlohmann::json& cfgJson);

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

    // Helper: build the input vector of RDataFrames
    std::vector<std::unique_ptr<ROOT::RDataFrame>> BuildDataFrames(const std::vector<std::string>& files,
                                            const std::string& treeName) const;                                            

    /// Configuration
    std::vector<std::string> fInputFiles; ///< comma-separated list
    std::string        fTreeNames;        ///< name of the input TTree
    std::vector<std::string>       fOutputFiles;         ///< result file
    std::vector<std::string> fVarsToKeep;///< thin list, incl. derived vars
    std::string        fRunLabel;        ///< “data2023B”, “bnb_overlay”, …

    /// Working objects
    std::unique_ptr<TChain>     fChain;
    std::vector<ROOT::RDF::RNode> fRDFs;

    std::map<std::string, std::vector<std::unique_ptr<ROOT::RDataFrame>>> fRdfMap;
    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfVec; ///< DataFrames for each input file


    std::map<std::string, std::string> treeMap = {
        {"pandora", "nuselection/NeutrinoSelectionFilter"},
        {"wirecell", "wcpselection/T_PFeval"},
        {"lantern", "lantern/EventTree"}
    };

};

} // namespace Analysis
#endif
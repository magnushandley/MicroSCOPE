#ifndef BDTTRAIN_MODULE_HXX
#define BDTTRAIN_MODULE_HXX
/*--------------------------------------------------------------------------*
 *  Helper functions to train and evaluate BDTs using TMVA
 *--------------------------------------------------------------------------*/
#include "Framework/Module.hxx"
#include "Utils/Plotter.hxx"

#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <memory>
#include <vector>

namespace Analysis {

class BDTTrainModule final : public Module {
public:
    explicit BDTTrainModule(const TEnv& cfg);

    // Module interface
    void Initialise() override;
    void Execute(Long64_t /*entry*/) override {}   // nothing per-event
    void Finalise()   override;
    Long64_t EntryCount() const override;

    std::string Name() const override { return "BDTTrainModule"; }

private:
    // Helper: build the input chain from a comma-separated list
    std::vector<std::unique_ptr<ROOT::RDataFrame>> BuildDataFrames(const std::vector<std::string>& files,
                                            const std::string& treeName) const;

     // Helper: Creates temporary snapshots to disk with training and testing samples, and weights added as new branch. 
     // returns the filenames of train_signal, train_bkg, test_signal, test_bkg files.                                           
    std::vector<std::string> BuildTestTrainSamples(std::vector<ROOT::RDF::RNode> dfs, std::vector<std::string> sampleLabels, std::vector<double> sampleWeights, float testFraction) const;

    /// Configuration
    std::vector<std::string>      fInputFiles;
    std::string        fTreeName;        ///< name of the input TTree
    std::vector<std::string> fSampleLabels; ///< Labels for the samples, e.g. "data", "overlay", "signal"
    std::vector<double> fSampleWeights; ///< Weights for each sample to normalise to POT
    std::vector<std::string> fVarsToKeep;///< thin list, incl. derived vars
    std::vector<std::string> fTrainVars;///< variables to use for training
    float fTrainFraction; ///< Fraction of events to use for training (rest for testing)

    /// Working objects
    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfVec; ///< DataFrames for each input file
    std::vector<ROOT::RDF::RNode> testNodes; ///< RNodes for each input file
    std::vector<ROOT::RDF::RNode> trainNodes; ///< RNodes for each input file
};

} // namespace Analysis
#endif
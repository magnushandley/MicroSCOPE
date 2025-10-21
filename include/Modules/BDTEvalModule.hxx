#ifndef BDTEVAL_MODULE_HXX
#define BDTEVAL_MODULE_HXX

#include "Framework/Module.hxx"
#include "Utils/Plotter.hxx"

#include <TEnv.h>
#include <string>
#include <vector>

namespace Analysis {

class BDTEvalModule final : public Module {
public:
    explicit BDTEvalModule(const TEnv& cfg);

    Long64_t EntryCount() const override;
    void Initialise() override;
    void Execute(Long64_t /*entry*/) override {}   // nothing per-event
    void Finalise() override;

    std::string Name() const override { return "BDTEvalModule"; }

private:
    // config
    std::string fTreeName;
    std::string fWeightsXML;
    std::string fMethodName;               // default "BDTG"
    std::string fOutputTag;                // appended to filenames, default "_bdt"
    std::vector<std::string> fInputFiles;
    std::vector<std::string> fEvalVars;    // must match training variable names

    // helpers
    static std::vector<std::string> TokeniseCSV(const std::string& s);
    void ProcessOneFile(const std::string& inPath) const;
};

} // namespace Analysis
#endif
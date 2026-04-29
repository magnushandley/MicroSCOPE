#ifndef PYTHON_MODULE_HXX
#define PYTHON_MODULE_HXX

#include "Framework/Module.hxx"

#include <TEnv.h>
#include <string>
#include <vector>

namespace Analysis {

class PythonModule final : public Module {
public:
    PythonModule(const TEnv& cfg, std::string configPath);

    void Initialise() override;
    void Execute(Long64_t /*entry*/) override {}
    void Finalise() override {}
    Long64_t EntryCount() const override { return 1; }

    std::string Name() const override { return "PythonModule"; }

private:
    std::string fConfigPath;
    std::string fPythonExecutable;
    std::string fScript;
    std::string fWorkingDirectory;
    bool fPassConfig = false;
    std::string fConfigArg;
    std::vector<std::string> fCommand;

    void ParseConfig();
    std::vector<std::string> BuildCommand() const;
    int RunCommand() const;

    static std::vector<std::string> SplitWhitespace(const std::string& value);
    static std::string RequireValue(const TEnv& cfg, const std::string& key);
    static std::string ToCliFlag(const std::string& name);
    static std::string FormatCommand(const std::vector<std::string>& command);
    static void StreamWithPrefix(int fd, const std::string& prefix);
};

} // namespace Analysis

#endif

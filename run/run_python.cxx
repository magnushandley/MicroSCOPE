#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Framework/Module.hxx"
#include "Framework/ModuleManager.hxx"
#include "Modules/PythonModule.hxx"

#include <TEnv.h>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: run_python <config.cfg>\n";
        return 1;
    }

    try {
        const std::string configPath = argv[1];
        TEnv cfg(configPath.c_str());

        std::vector<std::unique_ptr<Analysis::Module>> modules;
        modules.emplace_back(std::make_unique<Analysis::PythonModule>(cfg, configPath));

        Analysis::ModuleManager mgr(std::move(modules));
        mgr.Run();
    } catch (const std::exception& ex) {
        std::cerr << "[run_python] " << ex.what() << "\n";
        return 1;
    }
    return 0;
}

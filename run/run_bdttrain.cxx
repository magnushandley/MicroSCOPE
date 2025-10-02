#include <memory>
#include <vector>
#include "Framework/Module.hxx"
#include "Framework/ModuleManager.hxx"
#include "Modules/BDTTrainModule.hxx"
#include <TEnv.h>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc!=2) {
        std::cerr << "Usage: run_preselection <config.cfg>\n";
        return 1;
    }
    TEnv cfg(argv[1]);
    
    // Build the vector explicitly (was having problems before)

    std::vector<std::unique_ptr<Analysis::Module>> modules;
    modules.emplace_back(std::make_unique<Analysis::BDTTrainModule>(cfg));

    Analysis::ModuleManager mgr(std::move(modules));
    mgr.Run();
    return 0;
}
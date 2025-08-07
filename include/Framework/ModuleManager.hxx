#ifndef ANALYSIS_FRAMEWORK_MODULEMANAGER_HXX
#define ANALYSIS_FRAMEWORK_MODULEMANAGER_HXX
/*--------------------------------------------------------------------------*
 *  Orchestrates a list of Modules: Initialise → Event loop → Finalise
 *--------------------------------------------------------------------------*/



#include <TStopwatch.h>
#include <memory>
#include <vector>
#include "Framework/Module.hxx"

namespace Analysis {

class ModuleManager {
public:
    // Construct from a list of already created modules. 
    explicit ModuleManager(std::vector<std::unique_ptr<Module>> mods);

    // Add another module after construction.
    void Add(std::unique_ptr<Module> mod);

    // Run the full life-cycle for all registered modules.
    void Run();

private:
    // Determine how many events to loop over.
    Long64_t DetermineNEntries() const;

    std::vector<std::unique_ptr<Module>> fModules;
};

} // namespace Analysis
#endif
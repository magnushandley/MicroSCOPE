#include "Framework/ModuleManager.hxx"

#include <iostream>
#include <iomanip>

using namespace Analysis;

//----------------------------------------------------------------------------//
ModuleManager::ModuleManager(std::vector<std::unique_ptr<Module>> mods)
    : fModules(std::move(mods)) {}

//----------------------------------------------------------------------------//
void ModuleManager::Add(std::unique_ptr<Module> mod)
{
    fModules.emplace_back(std::move(mod));
}

//----------------------------------------------------------------------------//
Long64_t ModuleManager::DetermineNEntries() const
{
    for (const auto& m : fModules) {
        const Long64_t n = m->EntryCount();
        if (n > 0) return n;
    }
    throw std::runtime_error(
        "[ModuleManager] Could not determine number of events – "
        "not all modules returned a valid EntryCount().");
}

//----------------------------------------------------------------------------//
void ModuleManager::Run()
{
    if (fModules.empty())
        throw std::runtime_error("[ModuleManager] No modules registered!");

    //------------------------------------------------------------------//
    // 1.  Initialise
    //------------------------------------------------------------------//
    for (auto& m : fModules) {
        std::cout << "  ↳ Initialising " << m->Name() << " …\n";
        m->Initialise();

    // Have to call this after Initialise() to ensure all modules are ready - in the slimmer this requires
    // the TChain to be built first, in the others it requires available data frames.

        const Long64_t nEntries = DetermineNEntries();
        std::cout << "[ModuleManager] Will process " << nEntries
                << " entries.\n";

        //------------------------------------------------------------------//
        // 2.  Event loop - stopwatch only gives useful information if we actually loop over events in Execute()
        //------------------------------------------------------------------//

        TStopwatch sw;
        sw.Start();
        for (Long64_t i = 0; i < nEntries; ++i) {
            if (i%10000==0)
                std::cout << "\r[ModuleManager] " << std::setw(7) << i
                        << " / " << nEntries << std::flush;


            m->Execute(i);
        }
        sw.Stop();
        const Double_t cpuTime = sw.CpuTime();
        std::cout << "\r[ModuleManager] Finished loop in "
                << std::fixed << std::setprecision(3)
                << cpuTime << " s (" << (cpuTime/nEntries) * 1e3
                << " ms / evt)\n";

        //------------------------------------------------------------------//
        // 3.  Finalise
        //------------------------------------------------------------------//
        std::cout << "  ↳ Finalising " << m->Name() << " …\n";
        m->Finalise();
    }
}

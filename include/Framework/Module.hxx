#ifndef ANALYSIS_FRAMEWORK_MODULE_HXX
#define ANALYSIS_FRAMEWORK_MODULE_HXX
/*--------------------------------------------------------------------------*
 *  Abstract base class for every processing stage (Slimmer, Preselection, etc.)
 *--------------------------------------------------------------------------*/

#include <TEnv.h>
#include <Rtypes.h>          
#include <string>


namespace Analysis {

// All modules overwrite these general methods
// and can be run by the ModuleManager.
// The ModuleManager will call Initialize() once, then Execute() for each event, 
// then Finalize() once.

class Module {
public:
    explicit Module(const TEnv& cfg) : fCfg(&cfg) {}
    virtual ~Module() = default;

    virtual void Initialise() = 0;
    virtual void Execute(Long64_t entry) = 0;
    virtual void Finalise() = 0;

    
    virtual std::string Name() const = 0;

    virtual Long64_t EntryCount() const { return -1; }

    //Access to the configuration object
    const TEnv& Cfg() const { return *fCfg; }

protected:
    const TEnv* fCfg;   
};

} // namespace Analysis
#endif
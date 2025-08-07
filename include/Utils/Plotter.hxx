#ifndef ANALYSIS_UTILS_PLOTTER_HXX
#define ANALYSIS_UTILS_PLOTTER_HXX

/*--------------------------------------------------------------------------*
 *  Helper functions to draw and save histograms and graphs
 *--------------------------------------------------------------------------*/

#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <TH1D.h>

class TH1;
class THStack;
class TGraph;
class TCanvas;

namespace Analysis {

class Plotter {
public:
    // ------------------------------------------------------------------
    //  Single-histogram helpers
    // ------------------------------------------------------------------
    /** Draw a TH1 and write <basename>.png + .pdf into the current dir. */
    static void SaveHist(TH1* h,
                         const std::string& basename,
                         const std::string& style = "default");

    // ------------------------------------------------------------------
    //  Stacked histogram helpers
    /** Draw a stack of TH1s with labels and save as <basename>.png + .pdf. */

    static void StackedHist(std::vector<TH1D>& hists,
                        const std::vector<std::string>& labels,
                        const std::string& basename,
                        bool logy = false,
                        const std::vector<double> weights = {});

private:
    
    Plotter()  = default;
    ~Plotter() = default;

    // ----  internal helpers --------------------------------------
    static std::unique_ptr<TCanvas> MakeCanvas(const std::string& title);
    static void ApplyStyle(const std::string& style);
};

} // namespace Analysis
#endif 
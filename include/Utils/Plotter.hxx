#ifndef ANALYSIS_UTILS_PLOTTER_HXX
#define ANALYSIS_UTILS_PLOTTER_HXX

/*--------------------------------------------------------------------------*
 *  Helper functions to draw and save histograms and graphs
 *--------------------------------------------------------------------------*/

#include "Utils/ConfigUtils.hxx"

#include <string>
#include <vector>
#include <memory>
#include <TH1.h>
#include <TH1D.h>
#include <THStack.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFile.h>
#include <TSystem.h>
#include <iostream>
#include <TLine.h>
#include <ROOT/RDataFrame.hxx>

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
    
        static TH1D CreateTH1DFromRNode(
        ROOT::RDF::RNode node,
        const std::string& name,
        const std::string& varName,
        const std::string& xLabel,
        const std::string& yLabel,
        int nBins,
        double xMin,
        double xMax,
        bool removeVectorDuplicates = false,
        bool createOverFlowBin = false,
        std::string weightCol = "");

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

    static void FullDataMCSignalPlot(std::vector<TH1D>& hists,
                    const std::vector<std::string>& labels,
                    const std::string& basename,
                    bool logy,
                    const std::vector<double> weights,
                    double ratioYMin = 0.5,
                    double ratioYMax = 1.5,
                    const TH1D* bkgSysVarHist = nullptr,
                    const std::string& MicroBooNELabel = "NuMI 2.6e20 POT");

    static void FullDataMCSignalPlot(std::vector<TH1D>& hists,
                    const std::vector<std::string>& labels,
                    const std::vector<SampleType>& sampleTypes,
                    const std::string& basename,
                    bool logy,
                    const std::vector<double> weights,
                    double ratioYMin = 0.5,
                    double ratioYMax = 1.5,
                    const TH1D* bkgSysVarHist = nullptr,
                    const std::string& MicroBooNELabel = "NuMI 2.6e20 POT");

    static void BlindedMCSignalPlot(std::vector<TH1D>& rawHists,
                          const std::vector<std::string>& labels,
                          const std::string& basename,
                          bool logy = false,
                          const std::vector<double> weights = {},
                          const std::string& MicroBooNELabel = "MicroBooNE Simulation In Progress");

private:
    
    Plotter()  = default;
    ~Plotter() = default;

    // ----  internal helpers --------------------------------------
    static std::unique_ptr<TCanvas> MakeCanvas(const std::string& title,
                                               int width = 800,
                                               int height = 600);
    static void ApplyStyle(const std::string& style);
};

} // namespace Analysis
#endif 

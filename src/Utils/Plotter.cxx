#include "Utils/Plotter.hxx"

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

using namespace Analysis;

// ----------------------------------------------------------------------//
std::unique_ptr<TCanvas> Plotter::MakeCanvas(const std::string& title)
{
    auto c = std::make_unique<TCanvas>(title.c_str(), title.c_str(), 800, 600);
    c->SetTicks(1, 1);
    c->SetMargin(0.12, 0.02, 0.12, 0.08);
    return c;
}

// ----------------------------------------------------------------------//
void Plotter::ApplyStyle(const std::string& style)
{
    if (style == "prelim") {
        gStyle->SetOptStat(0);
        gStyle->SetLineWidth(2);
    }
    else if (style == "mdh_nice"){
        gStyle->SetTextFont(132);
        gStyle->SetLabelFont(132, "XYZ");   
        gStyle->SetTitleFont(132, "XYZ");   
        gStyle->SetTitleFont(132, "Title");
        gStyle->SetLegendFont(132);
        gStyle->SetLineWidth(2);
        gStyle->SetOptStat(0);
    }
    else {           // "default"
        gStyle->SetOptStat(0); 
    }
}

TH1D Plotter::CreateTH1DFromRNode(
    ROOT::RDF::RNode node,
    const std::string& name,
    const std::string& varName,
    const std::string& xLabel,
    const std::string& yLabel,
    int nBins,
    double xMin,
    double xMax,
    bool removeVectorDuplicates)
{
    std::string axisString = std::string(";") + xLabel + ";" + yLabel;
    ROOT::RDF::TH1DModel model(name.c_str(), axisString.c_str(), nBins, xMin, xMax);

    // If the variable is a vector, can remove duplicates by taking only the first element
    if (removeVectorDuplicates) {
        // Define a new column that extracts the first element of the vector,
        // if vector is empty, return default value -9999.0
        auto firstElementCol = node.Define((varName + "_first").c_str(),
            [](const ROOT::VecOps::RVec<float>& vec) {
                return vec.empty() ? -9999.0f : vec[0];
            }, {varName.c_str()});
        TH1D hist = firstElementCol.Histo1D(model, (varName + "_first").c_str()).GetValue();
        hist.SetDirectory(nullptr);   // decouple from any current file
        hist.SetName(name.c_str());
        return hist;
    }

    TH1D hist = node.Histo1D(model, varName).GetValue();
    hist.SetDirectory(nullptr);   // decouple from any current file
    hist.SetName(name.c_str());
    return hist;
}

// ----------------------------------------------------------------------//
void Plotter::SaveHist(TH1* h,
                       const std::string& basename,
                       const std::string& style)
{
    if (!h) return;
    ApplyStyle(style);

    auto c = MakeCanvas(basename);
    h->Draw("HIST");
    c->SaveAs((basename + ".png").c_str());
    c->SaveAs((basename + ".pdf").c_str());
}

// ----------------------------------------------------------------------//
void Plotter::StackedHist(std::vector<TH1D>& hists,
                      const std::vector<std::string>& labels,
                      const std::string& basename,
                      bool logy,
                      const std::vector<double> weights)
{
    // I am using TH1D objects rather than pointers because the pointers returned by RDataFrame are
    // smart pointers, which I think are deallocated before we can use them here. Had a bunch of seg faults...

    std::cout << "[Plotter] Creating stacked histogram: " << basename << std::endl;
    static const Int_t colors[] = {TColor::GetColor("#e69f00"),TColor::GetColor("#5664e9"),TColor::GetColor("#009e73"), kOrange, kViolet, kCyan, kMagenta, kYellow};
    
    if (hists.empty() || hists.size() != labels.size()) return;

    ApplyStyle("mdh_nice");
    auto c = MakeCanvas(basename);
    if (logy) c->SetLogy();

    std::string stackTitle = "Stacked Histogram: " + basename;
    THStack *hs = new THStack("hs", stackTitle.c_str());

    std::cout << "[Plotter] Adding " << hists.size() << " histograms to stack:\n";

    for (size_t i = 0; i < hists.size(); ++i) {
        auto& hist = hists[i];   
        std::cout << "  Adding hist number = " << i
                  << " name=" << hist.GetName()
                  << " entries=" << hist.GetEntries()
                  << std::endl;
        hist.SetFillColor(colors[i % (sizeof(colors)/sizeof(colors[0]))]);
        hist.SetLineColor(hist.GetFillColor());
        hist.SetFillStyle(1001);
        std::cout << " [Plotter] Applying weight for histogram " << i << ": ";
        if (i < weights.size() && weights[i] != 1.0) {
            hist.Scale(weights[i]);
            std::cout << "    Scaled by weight: " << weights[i] << std::endl;
        }
        hs->Add(&hist);
    }

    TList* histList = hs->GetHists();

    // Iterate through the list and print information about each histogram
    if (histList) {
        TIterator* nextHist = histList->MakeIterator();
        TObject* obj;
        std::cout << "Histograms in THStack '" << hs->GetName() << "':" << std::endl;
        while ((obj = nextHist->Next())) {
            TH1* hist = dynamic_cast<TH1*>(obj);
            if (hist) {
                std::cout << "  - Name: " << hist->GetName() << ", Title: " << hist->GetTitle() << std::endl;
            }
        }
        delete nextHist; // Clean up the iterator
    }

    

    // Draw the first hist to set axes
    //hists[0].Draw("HIST");
    //for (size_t i = 1; i < hists.size(); ++i)
    //    hists[i].Draw("HIST SAME");

    hs->Draw("HIST");
    hs->GetXaxis()->SetTitle(hists[0].GetXaxis()->GetTitle());
    hs->GetYaxis()->SetTitle(hists[0].GetYaxis()->GetTitle());
    // Legend
    auto leg = new TLegend(0.7, 0.7, 0.88, 0.88);
    for (size_t i = 0; i < hists.size(); ++i)
        leg->AddEntry(&hists[i], labels[i].c_str(), "l");   // pass pointer

    leg->Draw();

    c->SaveAs((basename + ".png").c_str());
    c->SaveAs((basename + ".pdf").c_str());
}

// ----------------------------------------------------------------------//
// ----------------------------------------------------------------------//
void Plotter::FullDataMCSignalPlot(std::vector<TH1D>& hists,
                      const std::vector<std::string>& labels,
                      const std::string& basename,
                      bool logy,
                      const std::vector<double> weights)
{
    // I am using TH1D objects rather than pointers because the pointers returned by RDataFrame are
    // smart pointers, which I think are deallocated before we can use them here. Had a bunch of seg faults...

    std::cout << "[Plotter] Creating full stacked histogram: " << basename << std::endl;
    static const Int_t colors[] = {TColor::GetColor("#e69f00"),TColor::GetColor("#5664e9"),TColor::GetColor("#009e73"), kOrange, kViolet, kCyan, kMagenta, kYellow};
    
    std::cout << "Number of entries in " << labels[0] << ": " << hists.back().GetEntries() << std::endl;
    std::cout << "Hists size: " << hists.size() << ", Labels size: " << labels.size() << std::endl;
    if (hists.empty() || hists.size() != labels.size()) return;
    std::cout << "Number of histograms: " << hists.size() << std::endl;

    ApplyStyle("prelim");
    auto c = MakeCanvas(basename);
    if (logy) c->SetLogy();

    std::string stackTitle = "Stacked Histogram: " + basename;
    THStack *hs = new THStack("hs", stackTitle.c_str());

    std::cout << "[Plotter] Adding " << hists.size() << " histograms to stack:\n";

    std::vector<TH1D*> signalHists;
    std::vector<TH1D*> dataHists;
    std::vector<TH1D*> bkgHists;

    //TCanvas *stackCanvas = new TCanvas("stackCanvas","MC Stack with Signal and Data",800,800);
    c->Divide(1,2);

    // Top pad: main plot                                                             
    c->cd(1);
    gPad->UseCurrentStyle();
    gPad->SetPad(0.0, 0.3, 1.0, 1.0);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.15);

    for (size_t i = 0; i < hists.size(); ++i) {
        auto& hist = hists[i];
        
        // Because the labels are used for plotting, I check for the substrings rather than whole things.
        auto lower = [](std::string s){ std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s; };
        const std::string llbl = lower(labels[i]);
        const bool isSignal = llbl.find("signal") != std::string::npos;
        const bool isData   = llbl.find("data")   != std::string::npos;

        if (isSignal) {
            // Style the signal and queue for overlay after the stack is drawn
            hist.SetLineColor(kRed);
            hist.SetLineWidth(2);
            hist.SetFillStyle(0);
            // Optionally scale signal here (currently unity)
            if (i < weights.size() && weights[i] != 1.0) {
                hist.Scale(weights[i]);
                std::cout << "    [Plotter] Signal scaled by weight: " << weights[i] << std::endl;
            }
            signalHists.push_back(&hist);
            std::cout << " [Plotter] Queued signal histogram " << i << " for overlay." << std::endl;
            std::cout << "Number of signal entries: " << hist.GetEntries() << std::endl;
        }
        else if (isData) {
            // Style the data and queue for overlay after the stack is drawn
            hist.Sumw2();
            hist.SetMarkerStyle(20);
            hist.SetMarkerSize(1);
            hist.SetLineColor(kBlack);
            hist.SetFillStyle(0);
            if (i < weights.size() && weights[i] != 1.0) {
                hist.Scale(weights[i]);
                std::cout << "    [Plotter] Data scaled by weight: " << weights[i] << std::endl;
            }
            dataHists.push_back(&hist);
            std::cout << " [Plotter] Queued data histogram " << i << " for overlay." << std::endl;
            std::cout << "Number of data entries: " << hist.GetEntries() << std::endl;
            std::cout << "Maximum bin content in data: " << hist.GetMaximum() << std::endl;
        }
        else {
            hist.Sumw2();
            std::cout << "  Adding hist number = " << i
                      << " name=" << hist.GetName()
                      << " entries=" << hist.GetEntries()
                      << std::endl;
            hist.SetFillColor(colors[i % (sizeof(colors)/sizeof(colors[0]))]);
            hist.SetLineColor(hist.GetFillColor());
            hist.SetFillStyle(1001);
            std::cout << " [Plotter] Applying weight for histogram " << i << ": ";
            if (i < weights.size() && weights[i] != 1.0) {
                hist.Scale(weights[i]);
                std::cout << "    Scaled by weight: " << weights[i] << std::endl;
            }
            hs->Add(&hist);
            bkgHists.push_back(&hist);
        }
    }

    TList* histList = hs->GetHists();

    // Iterate through the list and print information about each histogram
    if (histList) {
        TIterator* nextHist = histList->MakeIterator();
        TObject* obj;
        std::cout << "Histograms in THStack '" << hs->GetName() << "':" << std::endl;
        while ((obj = nextHist->Next())) {
            TH1* hist = dynamic_cast<TH1*>(obj);
            if (hist) {
                std::cout << "  - Name: " << hist->GetName() << ", Title: " << hist->GetTitle() << std::endl;
            }
        }
        delete nextHist; // Clean up the iterator
    }


    // Draw the stacked backgrounds first
    hs->Draw("HIST");

    // Build total background histogram and its statistical uncertainty band
    TH1D* hBkgTotal = nullptr;
    if (!bkgHists.empty()) {
        hBkgTotal = (TH1D*)bkgHists[0]->Clone((std::string(basename) + "_bkg_total").c_str());
        hBkgTotal->SetDirectory(nullptr);
        hBkgTotal->Sumw2();
        for (size_t ib = 1; ib < bkgHists.size(); ++ib) {
            hBkgTotal->Add(bkgHists[ib]);
        }
    }

    double yMax = hs->GetMaximum();
    for (auto* sh : signalHists) yMax = std::max(yMax, sh->GetMaximum());
    for (auto* dh : dataHists)   yMax = std::max(yMax, dh->GetMaximum());
    if (hBkgTotal) {
        // consider content+error to avoid clipping the band
        double yBandMax = 0.0;
        for (int b = 1; b <= hBkgTotal->GetNbinsX(); ++b) {
            yBandMax = std::max(yBandMax, hBkgTotal->GetBinContent(b) + hBkgTotal->GetBinError(b));
        }
        yMax = std::max(yMax, yBandMax);
    }
    hs->SetMaximum(1.2 * yMax);

    // Draw statistical uncertainty band for the total background
    if (hBkgTotal) {
        TH1D* hBkgBand = (TH1D*)hBkgTotal->Clone((std::string(basename) + "_bkg_band").c_str());
        hBkgBand->SetDirectory(nullptr);
        hBkgBand->SetFillColorAlpha(kGray+1, 0.5);
        //hBkgBand->SetFillStyle(3001); // hatched/transparent style
        hBkgBand->SetLineColor(kGray+2);
        hBkgBand->SetLineWidth(1);
        hBkgBand->SetMarkerSize(0);
        hBkgBand->Draw("E2 SAME"); // draw band as content ± error
    }

    // Now overlay signal(s) and data on top of the stack
    for (auto* sh : signalHists) sh->Draw("HIST SAME");
    for (auto* dh : dataHists)   dh->Draw("E SAME");

    // Axis titles come from the first histogram
    hs->GetXaxis()->SetTitle(hists[0].GetXaxis()->GetTitle());
    hs->GetYaxis()->SetTitle(hists[0].GetYaxis()->GetTitle());

    auto leg = new TLegend(0.7, 0.7, 0.88, 0.88);
    for (size_t i = 0; i < hists.size(); ++i) {
        auto lower = [](std::string s){ std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s; };
        const std::string llbl = lower(labels[i]);
        const bool isSignal = llbl.find("signal") != std::string::npos;
        const bool isData   = llbl.find("data")   != std::string::npos;
        const char* opt = isData ? "lep" : (isSignal ? "l" : "f");
        leg->AddEntry(&hists[i], labels[i].c_str(), opt);
    }
    if (hBkgTotal) {
        // Add a dummy clone for legend styling consistency
        TH1D* hBkgBandLegend = (TH1D*)hBkgTotal->Clone((std::string(basename) + "_bkg_band_legend").c_str());
        hBkgBandLegend->SetFillColor(kGray+1);
        hBkgBandLegend->SetFillStyle(3002);
        hBkgBandLegend->SetLineColor(kGray+2);
        hBkgBandLegend->SetMarkerSize(0);
        leg->AddEntry(hBkgBandLegend, "Bkg. stat. unc.", "f");
    }

    leg->Draw();

    // Bottom pad: ratio panel                                                        
    c->cd(2);
    gPad->UseCurrentStyle();
    gPad->SetPad(0.0, 0.0, 1.0, 0.3);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.3);
    gPad->SetLeftMargin(0.15);

    // Build and draw ratio histogram                                                 
    TH1D *ratio = (TH1D*)dataHists[0]->Clone("ratio");
    ratio->SetDirectory(0);
    ratio->SetTitle("");
    ratio->Sumw2();
    ratio->Divide(hBkgTotal);
    ratio->GetYaxis()->SetRangeUser(0.6,1.4);

    ratio->SetMarkerStyle(20);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->GetYaxis()->SetTitle("Data / MC");
    //ratio->GetYaxis()->SetTitleSize(0.08);                                          
    //ratio->GetYaxis()->SetLabelSize(0.06);                                          
    ratio->GetYaxis()->SetLabelOffset(0.0);
    ratio->GetYaxis()->SetTitleOffset(0.5);
    ratio->GetXaxis()->SetTitle(hists[0].GetXaxis()->GetTitle());
    ratio->GetXaxis()->SetTitleSize(0.08);
    ratio->GetXaxis()->SetLabelSize(0.08);
    ratio->GetYaxis()->SetTitleSize(0.08);
    ratio->GetYaxis()->SetLabelSize(0.08);

    // Draw MC statistical uncertainty band in ratio panel                            
    TH1D *ratioMC = (TH1D*)hBkgTotal->Clone("ratioMC");
    ratioMC->SetDirectory(0);
    ratioMC->SetTitle("");
    ratioMC->GetYaxis()->SetTitleOffset(0.5);
    ratioMC->GetXaxis()->SetTitle(hists[0].GetXaxis()->GetTitle());
    ratioMC->Sumw2();
    ratioMC->GetYaxis()->SetRangeUser(0.6,1.4);
    ratioMC->GetYaxis()->SetTitle("Data / MC");
    ratioMC->Divide(hBkgTotal);
    //ratioMC->SetFillStyle(3001);                    
    ratioMC->SetFillColorAlpha(kGray + 2, 0.3);     
    ratioMC->SetLineColor(kBlack);
    ratioMC->GetXaxis()->SetTitleSize(0.08);
    ratioMC->GetXaxis()->SetLabelSize(0.08);
    ratioMC->GetYaxis()->SetTitleSize(0.08);
    ratioMC->GetYaxis()->SetLabelSize(0.08);
    ratioMC->Draw("E2 SAME");

    // Draw ratio points on top                                                       
    ratio->Draw("E1 SAME");

    // Horizontal line at 1
    double xlow = ratio->GetXaxis()->GetXmin();
    double xhigh = ratio->GetXaxis()->GetXmax();                                                           
    TLine *unity = new TLine(xlow,1.0,xhigh,1.0);
    unity->SetLineStyle(2);
    unity->SetLineColor(kRed);
    unity->Draw();
    c->Update();

    c->SaveAs((basename + ".png").c_str());
    c->SaveAs((basename + ".pdf").c_str());
}
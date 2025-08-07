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
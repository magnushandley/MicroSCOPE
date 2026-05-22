// /Users/magnus/Documents/PhD/MicroSCOPE/src/Utils/SystematicsUtil.cxx
#include "Utils/SystematicsUtil.hxx"

#include <cmath>
#include <numeric>
#include <sstream>
#include <TH2D.h>
#include <TCanvas.h>
#include "ROOT/RVec.hxx"
#include <TLegend.h>
#include <stdexcept>
#include <string>

#include "TMath.h"

namespace Analysis
{
    SystematicsUtil::SystematicsUtil(std::string name)
    {
        (void)name;
    }


    std::vector<TH1D> SystematicsUtil::createMultiSimUniverses(
        const TH1D& nominalHist,
        const ROOT::RDF::RNode& rawDataFrame,
        const std::string& variableName,
        const std::string& multisimWeightColumn,
        const std::string& CVWeightColumn,
        const std::string& GlobalCVWeightColumn,
        const double& globalScaleFactor)
    {
        std::vector<TH1D> universes;
        ROOT::RDF::RNode dataFrame = rawDataFrame;

        // Ensure the variable we histogram is a scalar double column to avoid TTreeReader type mismatches
        // (e.g. int branches like n_pfps being read as float).
        std::string xCol = variableName;
        try {
            const std::string vType = dataFrame.GetColumnType(variableName);
            const std::string safeNameBase = variableName + "_asDouble";
            if (vType.find("RVec") != std::string::npos || vType.find("vector") != std::string::npos) {
                xCol = safeNameBase + "_first";
                dataFrame = dataFrame.Define(
                    xCol,
                    ("(" + variableName + ".size() > 0) ? (double)(" + variableName + "[0]) : -9999.0").c_str());
            } else {
                xCol = safeNameBase;
                dataFrame = dataFrame.Define(xCol, ("(double)(" + variableName + ")").c_str());
            }
        } catch (...) {
            xCol = variableName + "_asDouble";
            dataFrame = dataFrame.Define(xCol, ("(double)(" + variableName + ")").c_str());
        }

        //int nUniverses = dataFrame.Take(multisimWeightColumn)[0].GetSize(); 
        auto wTake = dataFrame.Take<ROOT::VecOps::RVec<unsigned short>>(multisimWeightColumn);
        const auto& allWeights = wTake.GetValue();     
        const auto& firstWeights = allWeights.at(0);   
        const int nUniverses = static_cast<int>(firstWeights.size());

        universes.reserve(nUniverses);
        double xMin = nominalHist.GetXaxis()->GetXmin();
        double xMax = nominalHist.GetXaxis()->GetXmax();
        int nBins = nominalHist.GetNbinsX();

        for (int u = 0; u < nUniverses; ++u)
        {
            universes.emplace_back(Form("%s_u%d", nominalHist.GetName(), u),
                                  nominalHist.GetTitle(), nBins, xMin, xMax);
            universes.back().Sumw2();
        }

        //Debug: create histogram which will contain the fractional difference between multisim mean and nominal weights
        TH1D fracDiffHist("fracDiffHist", "Fractional Difference between mean multisim weight and tune weight;Fractional Difference;Events", 200, -0.3, 0.4);
        TH1D rawMultiSimWeightHist("rawMultiSimWeightHist", "Raw Multisim Weights;Weight;Events", 500, -1000, 7000);
        TH1D rawMultiSimNormalisedbyCV("rawMultiSimNormalisedbyCV", "Raw Multisim Weights Normalised by CV tune weight;Normalised multisim weight;Events", 500, -1, 10);

        dataFrame.Foreach(
            [&](double x, const ROOT::VecOps::RVec<unsigned short>& mswvec, const float& cvWeight, const float& globalCVWeight, const float& splineWeight, const float& tuneWeight, const float& splineTimesTuneWeight, const float& ppfxWeight, const int& npi0)
            {
                const int nThis = static_cast<int>(mswvec.size());
                const int nFill = std::min(nUniverses, nThis);
                //const int nFill = 100; // For testing, only fill first 100 universes to speed up
                double mswSum = 0.0;
                int nSanitised = 0;
                
                for (int u = 0; u < nFill; ++u)
                {
                    double safeMSW = 0.0;
                    // Need to make mswvec[u] safe (could be zero, negative, or nan/inf). Need to be careful with type conversions because casting to
                    // unsigned short was causing loss of information.
                    rawMultiSimWeightHist.Fill(static_cast<double>(mswvec[u]));
                    if (mswvec[u] > 0 && mswvec[u] < 60000 && mswvec[u] != 1000){
                        rawMultiSimNormalisedbyCV.Fill(static_cast<double>(mswvec[u]) / (tuneWeight * 1000));
                        safeMSW = static_cast<double>(mswvec[u]) / 1000.0;
                    }
                    else {
                        safeMSW = 1.0;
                    }
                    //double safeMSW = (mswvec[u] > 0 && mswvec[u] < 60000.0 && mswvec[u] != 1000) ? static_cast<double>(mswvec[u])/1000.0 : 1000;
                    bool sanitised = (safeMSW == 1.0);
                    if (!sanitised) {
                        universes[u].Fill(x, static_cast<double>(safeMSW) * globalScaleFactor * cvWeight);
                        mswSum += safeMSW;
                    }
                    else {
                        universes[u].Fill(x, static_cast<double>(globalCVWeight) * globalScaleFactor);
                        nSanitised++;
                    }
                }
                //if (nSanitised > 10) {
                    //std::cout << "Warning: large number of sanitised universes" << std::endl;
                    //Print all weights for this event for debugging
                    //std::cout << "Event with x = " << x << " has " << nSanitised << " sanitised universes. Weights: ";
                    //for (int u = 0; u < nFill; ++u) {
                    //    std::cout << "Universe " << u << ": ";
                    //    std::cout << "MSW: " << mswvec[u];
                    //}
                    //std::cout << std::endl;
                //}
                double mswMean = mswSum / (nFill - nSanitised);
                double fracDiff = (mswMean - static_cast<double>(tuneWeight)) / static_cast<double>(tuneWeight);
                fracDiffHist.Fill(fracDiff);
                //if (fracDiff < -0.05 || fracDiff > 0.05) {
                //    std::cout << "Warning: large fractional difference between multisim mean and CV weight: ";
                //    std::cout <<"CV weight: " << cvWeight << " Mean multisim weight: " << mswMean << " Spline Weight: " << splineWeight << " Tune Weight: " << tuneWeight << " Spline Times Tune Weight: " << splineTimesTuneWeight << " PPFX Weight: " << ppfxWeight << " nPi0: " << npi0 << " Number of sanitised universes: " << nSanitised << std::endl;
                //}
            },
            {xCol, multisimWeightColumn, CVWeightColumn, GlobalCVWeightColumn, "weightSpline", "weightTune", "weightSplineTimesTune", "ppfx_cv", "npi0"}
        );
        //Plot fractional difference histogram to check sanity of multisim weights compared to CV weight
        //TCanvas cFrac("cFrac", "Fractional Difference", 800, 600);
        //fracDiffHist.Draw("HIST");
        //cFrac.Update();
        //cFrac.SaveAs("debug_multisim_fractional_difference.png");

        //Plot raw multisim weight distribution to check for any pathological weights
        //TCanvas cRaw("cRaw", "Raw Multisim Weights", 800, 600);
        //rawMultiSimWeightHist.Draw("HIST");
        //cRaw.Update();
        //cRaw.SaveAs("debug_multisim_raw_weights.png");

        //Plot raw multisim weight distribution normalised by CV weight
        //TCanvas cRawNorm("cRawNorm", "Raw Multisim Weights Normalised by CV", 800, 600);
        //rawMultiSimNormalisedbyCV.Draw("HIST");
        //cRawNorm.Update();
        //cRawNorm.SaveAs("debug_multisim_raw_weights_normalised_by_cv.png");

        //Debug: print out nominal histogram contents
        //std::cout << "Nominal histogram contents: ";
        //for (int b = 1; b <= nominalHist.GetNbinsX(); ++b)
        //{
        //    std::cout << nominalHist.GetBinContent(b) << " ";
        //}
        //std::cout << std::endl;
        //Debug: print out the first few universes
        //for (int u = 0; u < std::min(20, nUniverses); ++u)
        //{
        //    std::cout << "Universe " << u << " contents: ";
        //    for (int b = 1; b <= universes[u].GetNbinsX(); ++b)
        //    {
        //        std::cout << universes[u].GetBinContent(b) << " ";
        //    }
        //    std::cout << std::endl;
        //}
        //Debug: calculate the mean universe and fractional difference from nominal
        //Also plot each universe overlaid with nominal for a few universes to visually check sanity
        TCanvas c("c", "c", 800, 600);
        
        // Draw the first universe first to setup axes, but ensure title/labels are correct
        // Or better: Draw valid frame or nominal first. 
        // Let's modify logic to ensure axis labels apply.
        
        TH1D nomHistClone = nominalHist; // Make a clone to avoid modifying original
        nomHistClone.SetLineColor(kRed);
        nomHistClone.SetLineWidth(2);
        nomHistClone.SetFillStyle(0);
        // Set title and axis labels here
        nomHistClone.SetTitle("Multisim Universes Comparison;Logit BDT Score;Events");
        nomHistClone.SetStats(0); // Optional: remove stats box for cleaner plot

        //Debug, set y axis to log scale to better see differences in low-stat bins
        c.SetLogy();

        // Calculate mean universe first so we can add it to legend
        TH1D meanUniverse("mean_universe", "Mean of Universes", nBins, xMin, xMax);
        //meanUniverse.Sumw2();
        for (int u = 0; u < nUniverses; ++u)
        {
            meanUniverse.Add(&universes[u]);
        }
        meanUniverse.Scale(1.0 / nUniverses);
        
        // Draw Nominal First to set the axes range and labels
        nomHistClone.Draw("HIST"); 

        for (int u = 0; u < nUniverses; ++u)
        {
            universes[u].SetLineColor(kGray+1); // Lighter color for background lines often looks better
            universes[u].SetLineWidth(1);
            universes[u].SetFillStyle(0);
            universes[u].Draw("HIST SAME");
        }
        
        // Redraw Nominal on top
        nomHistClone.Draw("HIST SAME");

        // ... existing code for printing values ...
        
        // Plot mean universe
        meanUniverse.SetLineColor(kBlue);
        meanUniverse.SetLineWidth(2);
        meanUniverse.SetFillStyle(0);
        meanUniverse.Draw("HIST SAME");

        // Add Legend
        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(&nomHistClone, "Nominal CV", "l");
        leg->AddEntry(&meanUniverse, "Mean Universe", "l");
        if (nUniverses > 0) leg->AddEntry(&universes[0], "Multisim Universes", "l");
        leg->Draw();

        c.Update();
        c.SaveAs("debug_universes_overlay.png");

        // Clean up legend
        delete leg; 

        return universes;
    }


    using EventKey = uint64_t;
    using EventSet = std::unordered_set<EventKey>;

    // Fast hash map for run, sub, evt -> single number so we can do matching.
    // Caveat: only works for run < 2^22, subrun < 2^21, event < 2^21.
    // If we exceed these, then we will have bit collisions (bad).

    EventKey MakeKey(int run, int subrun, int event) {
        return (uint64_t(run)   << 42) |
            (uint64_t(subrun) << 21) |
            uint64_t(event);
    }

    EventSet BuildCVEventSet(ROOT::RDF::RNode cv_df,
                        const std::string& run_col = "run",
                        const std::string& subrun_col = "sub",
                        const std::string& event_col = "evt")
    {
        EventSet cv_events;
        cv_df.Foreach(
            [&cv_events](int run,
                        int subrun,
                        int event) {
                cv_events.insert(MakeKey(run, subrun, event));
            },
            {run_col, subrun_col, event_col}
        );
        return cv_events;
    }

    ROOT::RDF::RNode FilterToCVEvents(ROOT::RDF::RNode detvar_df,
                    std::shared_ptr<const EventSet> cv_events,
                    const std::string& run_col = "run",
                    const std::string& subrun_col = "sub",
                    const std::string& event_col = "evt")
    {
        return detvar_df.Filter(
            [cv_events](int run,
                        int subrun,
                        int event) {
                return cv_events->count(MakeKey(run, subrun, event)) > 0;
            },
            {run_col, subrun_col, event_col},
            "Keep only detector variation events present in CV"
        );
    }

    std::vector<TH1D> SystematicsUtil::createMatchedDetVarHists(
        const ROOT::RDF::RNode& CVNode,
        const std::vector<ROOT::RDF::RNode>& detVarNodes,
        const std::vector<std::string>& detVarNames,
        std::string branchName,
        const std::string& weightColumn,
        const TH1D& CVNominalHist)
        {

        using EventKey = uint64_t;
        using EventSet = std::unordered_set<EventKey>;

        auto cv_events = std::make_shared<EventSet>(BuildCVEventSet(CVNode, "run", "sub", "evt"));
        std::vector<TH1D> hists;
        hists.reserve(detVarNodes.size());

        auto hist_model = CVNominalHist; // Use CV nominal as model for binning, etc.

        for (size_t i = 0; i < detVarNodes.size(); ++i) {
            auto detvar_df = detVarNodes[i];
            auto nEntriesBefore = detvar_df.Count().GetValue();
            auto matched_df = FilterToCVEvents(detvar_df, cv_events,
                                "run", "sub", "evt");
            auto nEntriesAfter = matched_df.Count().GetValue();
            std::cout << "DetVar sample " << detVarNames[i] << ": "
                      << nEntriesBefore << " entries before matching, "
                      << nEntriesAfter << " entries after matching.\n";
            std::cout << "  Matched " << nEntriesAfter << " / " << nEntriesBefore
                      << " events (" << (100.0 * nEntriesAfter / nEntriesBefore) << "%) in detvar sample " << detVarNames[i] << "\n";

            ROOT::RDF::RResultPtr<TH1D> histPtr = weightColumn.empty()
                ? matched_df.Histo1D(hist_model, branchName)
                : matched_df.Histo1D(hist_model, branchName, weightColumn);

            TH1D hist = *histPtr;
            hist.SetDirectory(nullptr);

            if (i < detVarNames.size()) {
                hist.SetName(Form("%s_%s", CVNominalHist.GetName(), detVarNames[i].c_str()));
                hist.SetTitle(Form("%s %s", CVNominalHist.GetTitle(), detVarNames[i].c_str()));
            }

            hists.push_back(std::move(hist));
        }
        return hists;
    }



    // ------------------------------------------------------------------------
    // Compute covariance from a set of universe histograms.
    // ------------------------------------------------------------------------
    TMatrixD SystematicsUtil::covarianceMatrixFromMultisims(
        const std::vector<TH1D>& universeHists,
        const TH1D& nominalHist)
    {
        const int nUniverses = static_cast<int>(universeHists.size());
        const int nBins = nominalHist.GetNbinsX(); // excluding UF/OF

        TMatrixD cov(nBins, nBins);
        cov.Zero();

        if (nUniverses == 0)
        {
            std::cout << "[SystematicsUtil] Warning: zero universes provided for covariance calculation. Returning zero matrix." << std::endl;
            return cov;
        }
        // Compute covariance matrix
        for (int i = 0; i < nBins; ++i)
        {
            const double nom_i = nominalHist.GetBinContent(i + 1); // bins are 1-indexed
            for (int j = 0; j < nBins; ++j)
            {
                const double nom_j = nominalHist.GetBinContent(j + 1); // bins are 1-indexed
                double cov_ij = 0.0;
                for (int u = 0; u < nUniverses; ++u)
                {
                    const double uni_i = universeHists[u].GetBinContent(i + 1);
                    const double uni_j = universeHists[u].GetBinContent(j + 1);
                    cov_ij += (uni_i - nom_i) * (uni_j - nom_j);
                }
                cov_ij /= static_cast<double>(nUniverses); 
                
                double frac_cov_ij = (nom_i != 0.0 && nom_j != 0.0) ? cov_ij / (nom_i * nom_j) : 0.0;
                //std::cout << "Covariance[" << i << "," << j << "] = " << cov_ij << std::endl;
                //std::cout << " Fractional Covariance[" << i << "," << j << "] = "
                //          << frac_cov_ij << std::endl;
                cov(i, j) = cov_ij;
            }
        }
        //Print out fractional uncertainties (diagonal elements)
        std::cout << "Fractional uncertainties per bin:" << std::endl;
        for (int i = 0; i < nBins; ++i)
        {
            const double nom_i = nominalHist.GetBinContent(i + 1);
            const double var_i = cov(i, i);
            double frac_uncert_i = (nom_i != 0.0) ? std::sqrt(var_i) / nom_i : 0.0;
            std::cout << " Bin " << i << ": " << frac_uncert_i << std::endl;
        }
        return cov;
    }

    TMatrixD SystematicsUtil::covarianceMatrixFromDetVars(
        const TH1D& nominalHist,
        const std::vector<TH1D>& detVarHists,
        const std::vector<double>& detVarWeights,
        const double nomHistScaleFactor,
        const double uncertaintyCap)
        {
            // Create a covariance matrix from a set of set of detvars, but also have the option to apply a cap to the uncertainty to avoid pathological bins dominating fits. 
            
            const int nBins = nominalHist.GetNbinsX();
            const int nDetVars = static_cast<int>(detVarHists.size());

            TMatrixD cov(nBins, nBins);
            cov.Zero();

            for (int i = 0; i < nBins; ++i)
            {
                const double nom_i = nominalHist.GetBinContent(i + 1) * nomHistScaleFactor;
                for (int j = 0; j < nBins; ++j)
                {
                    const double nom_j = nominalHist.GetBinContent(j + 1) * nomHistScaleFactor;
                    double cov_ij = 0.0;
                    for (int v = 0; v < nDetVars; ++v)
                    {
                        const double weight = detVarWeights[v];
                        const double var_i = detVarHists[v].GetBinContent(i + 1) * weight;
                        const double var_j = detVarHists[v].GetBinContent(j + 1) * weight;
                        double delta_i = var_i - nom_i;
                        //double maxDelta_i = uncertaintyCap * std::abs(nom_i);

                        //if (std::abs(delta_i) > maxDelta_i) {
                        //    std::cout << "Warning: detvar variation in bin " << i << " exceeds cap. Nominal: " << nom_i << ", variation: " << var_i
                        //              << ", delta: " << delta_i << ", max allowed delta: " << maxDelta_i << std::endl;
                        //    delta_i = std::copysign(maxDelta_i, delta_i);
                        //}
                        double delta_j = var_j - nom_j;
                        //double maxDelta_j = uncertaintyCap * std::abs(nom_j);
                        //if (std::abs(delta_j) > maxDelta_j) {
                        //    std::cout << "Warning: detvar variation in bin " << j << " exceeds cap. Nominal: " << nom_j << ", variation: " << var_j
                        //              << ", delta: " << delta_j << ", max allowed delta: " << maxDelta_j << std::endl;
                        //    delta_j = std::copysign(maxDelta_j, delta_j);
                        //}
                        cov_ij += delta_i * delta_j;
                    }
                    cov(i, j) = cov_ij;
                    //double frac_cov_ij = (nom_i != 0.0 && nom_j != 0.0) ? cov_ij / (nom_i * nom_j) : 0.0;
                    //if (frac_cov_ij > uncertaintyCap*uncertaintyCap) {
                    //    std::cout << "Warning: fractional covariance between bin " << i << " and bin " << j << " is " << frac_cov_ij
                    //              << ", which exceeds the cap of " << uncertaintyCap << ". Capping covariance to " << uncertaintyCap * uncertaintyCap * nom_i * nom_j << std::endl;
                    //    cov(i, j) = uncertaintyCap * uncertaintyCap * nom_i * nom_j;
                    //    cov(j, i) = cov(i, j); // Ensure covariance matrix remains symmetric
                    //}
                    //else {
                    //    std::cout << "Detvar Covariance[" << i << "," << j << "] = " << cov_ij << std::endl;
                    //    std::cout << " Fractional Covariance[" << i << "," << j << "] = "
                    //              << frac_cov_ij << std::endl;
                    //}

                }
            }
            std::vector<double> scale(nBins, 1.0);
            for (int i = 0; i < nBins; ++i) {
                const double variance = cov(i, i);
                const double sigma = variance > 0.0 ? std::sqrt(variance) : 0.0;
                const double nom = nominalHist.GetBinContent(i + 1) * nomHistScaleFactor;
                const double capSigma = uncertaintyCap * std::abs(nom);

                if (sigma > 0.0 && sigma > capSigma) {
                    scale[i] = capSigma / sigma;
                }
            }

            // 3. Apply C' = D C D once.
            for (int i = 0; i < nBins; ++i) {
                for (int j = 0; j < nBins; ++j) {
                    cov(i, j) *= scale[i] * scale[j];
                }
            }
            return cov;
        }

    TMatrixD SystematicsUtil::BuildDetVarCovariance(
        const TH1D& nominalHist,
        const ROOT::RDF::RNode& cvNode,
        const std::vector<ROOT::RDF::RNode>& detVarNodes,
        const std::vector<std::string>& detVarNames,
        const std::string& variableName,
        const std::vector<double>& detVarGlobalWeights,
        double nomHistScaleFactor,
        const std::string& weightColumn)
    {
        if (detVarNodes.size() != detVarNames.size()) {
            throw std::runtime_error("[SystematicsUtil] Detector variation node/name count mismatch.");
        }
        if (detVarNodes.size() != detVarGlobalWeights.size()) {
            throw std::runtime_error("[SystematicsUtil] Detector variation node/weight count mismatch.");
        }

        const int nBins = nominalHist.GetNbinsX();
        TMatrixD zeroCov(nBins, nBins);
        zeroCov.Zero();

        if (detVarNodes.empty()) {
            std::cout << "[SystematicsUtil] Warning: zero detector variations provided. Returning zero matrix." << std::endl;
            return zeroCov;
        }

        auto cvEvents = std::make_shared<EventSet>(BuildCVEventSet(cvNode, "run", "sub", "evt"));
        std::vector<TH1D> detVarHists;
        std::vector<double> detVarWeights;
        detVarHists.reserve(detVarNodes.size());
        detVarWeights.reserve(detVarNodes.size());

        TH1D histModel = nominalHist;

        for (size_t i = 0; i < detVarNodes.size(); ++i) {
            ROOT::RDF::RNode detVarNode = detVarNodes[i];
            const auto nEntriesBefore = detVarNode.Count().GetValue();
            if (nEntriesBefore == 0) {
                throw std::runtime_error("[SystematicsUtil] Detector variation sample has zero events before matching: "
                                         + detVarNames[i]);
            }

            ROOT::RDF::RNode matchedNode = FilterToCVEvents(detVarNode, cvEvents, "run", "sub", "evt");
            const auto nEntriesAfter = matchedNode.Count().GetValue();
            if (nEntriesAfter == 0) {
                throw std::runtime_error("[SystematicsUtil] Detector variation sample has zero matched events: "
                                         + detVarNames[i]);
            }

            std::cout << "Detector Variation: " << detVarNames[i]
                      << " Initial events: " << nEntriesBefore
                      << " Final events: " << nEntriesAfter
                      << " Fraction kept: "
                      << (static_cast<double>(nEntriesAfter) / static_cast<double>(nEntriesBefore))
                      << std::endl;

            ROOT::RDF::RResultPtr<TH1D> histPtr = weightColumn.empty()
                ? matchedNode.Histo1D(histModel, variableName)
                : matchedNode.Histo1D(histModel, variableName, weightColumn);

            TH1D hist = *histPtr;
            hist.SetDirectory(nullptr);
            hist.SetName(Form("%s_%s", nominalHist.GetName(), detVarNames[i].c_str()));
            hist.SetTitle(Form("%s %s", nominalHist.GetTitle(), detVarNames[i].c_str()));
            detVarHists.push_back(std::move(hist));

            const double matchingScale = static_cast<double>(nEntriesBefore) / static_cast<double>(nEntriesAfter);
            const double totalScale = detVarGlobalWeights[i] * matchingScale;
            detVarWeights.push_back(totalScale);

            std::cout << "Detector Variation: " << detVarNames[i]
                      << " Base global scale: " << detVarGlobalWeights[i]
                      << " Matching scale: " << matchingScale
                      << " Total covariance scale: " << totalScale
                      << std::endl;
        }

        return covarianceMatrixFromDetVars(nominalHist, detVarHists, detVarWeights, nomHistScaleFactor);
    }


    TMatrixD SystematicsUtil::combineCovarianceMatrices(
        const std::vector<TMatrixD>& matrices)
    {
        if (matrices.empty())
        {
            std::cout << "[SystematicsUtil] Warning: no matrices provided for combination. Returning zero matrix." << std::endl;
            return TMatrixD();
        }

        const int nRows = matrices.front().GetNrows();
        const int nCols = matrices.front().GetNcols();

        TMatrixD total(nRows, nCols);
        total.Zero();

        for (size_t k = 0; k < matrices.size(); ++k)
        {
            const auto& m = matrices[k];
            if (m.GetNrows() != nRows || m.GetNcols() != nCols)
            {
                std::ostringstream ss;
                ss << "combineCovarianceMatrices: matrix " << k
                   << " has dimension (" << m.GetNrows() << "," << m.GetNcols()
                   << ") but expected (" << nRows << "," << nCols << ").";
                throw std::runtime_error(ss.str());
            }
            total += m;
        }

        return total;
    }

    void SystematicsUtil::PlotMatrix(const TMatrixD& matrix, const std::string& name)
    {
        int nRows = matrix.GetNrows();
        int nCols = matrix.GetNcols();

        TH2D hist(Form("%s_hist", name.c_str()), Form("Covariance Matrix: %s", name.c_str()),
                  nCols, 0.5, nCols + 0.5,
                  nRows, 0.5, nRows + 0.5);

        for (int i = 0; i < nRows; ++i)
        {
            for (int j = 0; j < nCols; ++j)
            {
                hist.SetBinContent(j + 1, i + 1, matrix(i, j));
            }
        }

        TCanvas c(Form("%s_canvas", name.c_str()), Form("Covariance Matrix Canvas: %s", name.c_str()), 800, 600);
        hist.GetXaxis()->SetTitle("Bin Index");
        hist.GetYaxis()->SetTitle("Bin Index");
        hist.Draw("COLZ");
        c.SaveAs(Form("%s.pdf", name.c_str()));
    }

    void SystematicsUtil::PlotFractionalCovarianceMatrix(
    const TMatrixD& cov,
    const TH1D& nominalHist,
    const std::string& outName) const
    {
        const int nRows = cov.GetNrows();
        const int nCols = cov.GetNcols();
        const int nBins = nominalHist.GetNbinsX(); 

        if (nRows != nCols)
        {
            std::ostringstream ss;
            ss << "PlotFractionalCovarianceMatrix: covariance matrix is not square ("
            << nRows << "x" << nCols << ")";
            throw std::runtime_error(ss.str());
        }

        if (nRows != nBins)
        {
            std::ostringstream ss;
            ss << "PlotFractionalCovarianceMatrix: matrix dimension (" << nRows
            << ") does not match nominal histogram bins (" << nBins << ")";
            throw std::runtime_error(ss.str());
        }

        TH2D hist(
            Form("%s_hist", outName.c_str()),
            Form("Fractional covariance matrix: %s;Bin index;Bin index", outName.c_str()),
            nCols, 0.5, nCols + 0.5,
            nRows, 0.5, nRows + 0.5);

        hist.SetStats(0);

        for (int i = 0; i < nRows; ++i)
        {
            const double nom_i = nominalHist.GetBinContent(i + 1);
            for (int j = 0; j < nCols; ++j)
            {
                const double nom_j = nominalHist.GetBinContent(j + 1);
                const double denom = nom_i * nom_j;

                // fractional cov_ij = cov_ij / (N_i * N_j)
                const double fracCov = (denom != 0.0) ? (cov(i, j) / denom) : 0.0;
                hist.SetBinContent(j + 1, i + 1, fracCov);
            }
        }

        TCanvas c(Form("%s_canvas", outName.c_str()), "Fractional Covariance Matrix", 900, 800);
        hist.Draw("COLZ");
        c.SaveAs(Form("%s.pdf", outName.c_str()));
        c.SaveAs(Form("%s.png", outName.c_str()));
    }



    std::pair<TMatrixD, TVectorD>
    SystematicsUtil::EigenDecomposition(const TMatrixD& cov)
    {
        TVectorD eigenValues;
        TMatrixD eigenVectors = cov.EigenVectors(eigenValues); // construct with correct size
        return {eigenVectors, eigenValues};
    }

    std::vector<TVectorD> SystematicsUtil::CreateNuisanceParams(const TMatrixD& cov)
    {
        // Perform eigen decomposition, then return a vector of nuisance parameters corresponding to the eigenvectors,
        // with uncertainties given by the square root of the eigenvalues.
        auto [eigenVectors, eigenValues] = EigenDecomposition(cov);
        std::vector<TVectorD> nuisanceParams;
        TVectorD vk(cov.GetNrows());
        for (int k = 0; k < eigenValues.GetNrows(); ++k) {
            for (int i = 0; i < cov.GetNrows(); ++i) {
                vk(i) = eigenVectors(i, k);
            }
            nuisanceParams.push_back(vk * std::sqrt(eigenValues(k)));
        }
        return nuisanceParams;
    }

    void SystematicsUtil::PlotNuisanceParams(
    const std::vector<TVectorD>& nuisanceParams,
    const TH1D& nominalHist,
    const std::string& outPrefix,
    const int maxOverlay) const
    {
        const int nPars = static_cast<int>(nuisanceParams.size());
        const int nBins = nominalHist.GetNbinsX(); // excluding UF/OF

        if (nPars == 0)
        {
        std::cout << "[SystematicsUtil] PlotNuisanceParams: no nuisance parameters provided. Nothing to plot.\n";
        return;
        }

        // Validate dimensions
        for (int k = 0; k < nPars; ++k)
        {
        if (nuisanceParams[k].GetNrows() != nBins)
        {
            std::ostringstream ss;
            ss << "PlotNuisanceParams: nuisance param " << k
            << " has length " << nuisanceParams[k].GetNrows()
            << " but nominal histogram has " << nBins << " bins.";
            throw std::runtime_error(ss.str());
        }
        }

        const double xMin = 0.5;
        const double xMax = nBins + 0.5;

        // ------------------------------------------------------------------
        // 2D heatmaps: parameter index vs bin index
        // ------------------------------------------------------------------
        TH2D hAbs(
        Form("%s_nuis_abs", outPrefix.c_str()),
        "Nuisance parameters (absolute);Bin index;Nuisance parameter index",
        nBins, xMin, xMax,
        nPars, 0.5, nPars + 0.5);

        TH2D hFrac(
        Form("%s_nuis_frac", outPrefix.c_str()),
        "Nuisance parameters (fractional vs nominal);Bin index;Nuisance parameter index",
        nBins, xMin, xMax,
        nPars, 0.5, nPars + 0.5);

        for (int k = 0; k < nPars; ++k)
        {
        for (int i = 0; i < nBins; ++i)
        {
            const double absVal = nuisanceParams[k](i);
            const double nom    = nominalHist.GetBinContent(i + 1);
            const double fracVal = (nom != 0.0) ? (absVal / nom) : 0.0;

            hAbs.SetBinContent(i + 1, k + 1, absVal);
            hFrac.SetBinContent(i + 1, k + 1, fracVal);
        }
        }

        {
        TCanvas cAbs(Form("%s_cAbs", outPrefix.c_str()), "Nuisance parameters absolute", 1100, 700);
        hAbs.SetStats(0);
        hAbs.Draw("COLZ");
        cAbs.SaveAs(Form("%s_nuisance_params_absolute.pdf", outPrefix.c_str()));
        cAbs.SaveAs(Form("%s_nuisance_params_absolute.png", outPrefix.c_str()));
        }

        {
        TCanvas cFrac(Form("%s_cFrac", outPrefix.c_str()), "Nuisance parameters fractional", 1100, 700);
        hFrac.SetStats(0);
        hFrac.Draw("COLZ");
        cFrac.SaveAs(Form("%s_nuisance_params_fractional.pdf", outPrefix.c_str()));
        cFrac.SaveAs(Form("%s_nuisance_params_fractional.png", outPrefix.c_str()));
        }

        // ------------------------------------------------------------------
        // Overlay plots (first N nuisance parameters) for quick inspection
        // ------------------------------------------------------------------
        const int nOverlay = std::min(maxOverlay, nPars);
        if (nOverlay <= 0) return;

        // Helper to compute overall min/max across the first nOverlay parameters
        auto computeMinMax = [&](const bool fractional) {
            double vmin = 0.0;
            double vmax = 0.0;
            bool first = true;

            for (int k = 0; k < nOverlay; ++k)
            {
                for (int i = 0; i < nBins; ++i)
                {
                    const double nom = nominalHist.GetBinContent(i + 1);
                    const double val = fractional ? ((nom != 0.0) ? (nuisanceParams[k](i) / nom) : 0.0)
                                                  : nuisanceParams[k](i);
                    if (first)
                    {
                        vmin = val;
                        vmax = val;
                        first = false;
                    }
                    else
                    {
                        vmin = std::min(vmin, val);
                        vmax = std::max(vmax, val);
                    }
                }
            }

            double span = vmax - vmin;
            if (span <= 0.0) span = (vmax != 0.0) ? std::abs(vmax) : 1.0;
            const double margin = 0.10 * span;
            return std::pair<double, double>(vmin - margin, vmax + margin);
        };

        // Absolute overlay
        {
        TCanvas cOvAbs(Form("%s_cOvAbs", outPrefix.c_str()), "Nuisance parameters overlay (absolute)", 1100, 700);

        TH1D frame(
            Form("%s_frame_abs", outPrefix.c_str()),
            "Nuisance parameters (absolute, overlay);Bin index;Absolute nuisance shift",
            nBins, xMin, xMax);
        frame.SetStats(0);
        frame.SetLineColor(0);
        frame.SetFillStyle(0);
        const auto [absMin, absMax] = computeMinMax(false);
        frame.SetMinimum(absMin);
        frame.SetMaximum(absMax);
        frame.Draw("AXIS");

        TLegend leg(0.70, 0.70, 0.93, 0.93);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);

        std::vector<TH1D> hists;
        hists.reserve(nOverlay);

        for (int k = 0; k < nOverlay; ++k)
        {
            hists.emplace_back(Form("%s_nuis_abs_k%d", outPrefix.c_str(), k), "", nBins, xMin, xMax);
            hists.back().SetStats(0);
            hists.back().SetLineWidth(2);
            hists.back().SetLineColor(k + 1);
            hists.back().SetFillStyle(0);

            for (int i = 0; i < nBins; ++i)
            hists.back().SetBinContent(i + 1, nuisanceParams[k](i));

            hists.back().Draw("HIST SAME");
            leg.AddEntry(&hists.back(), Form("k=%d", k), "l");
        }

        leg.Draw();
        cOvAbs.SaveAs(Form("%s_nuisance_params_absolute_overlay_first%d.pdf", outPrefix.c_str(), nOverlay));
        cOvAbs.SaveAs(Form("%s_nuisance_params_absolute_overlay_first%d.png", outPrefix.c_str(), nOverlay));
        }

        // Fractional overlay
        {
        TCanvas cOvFrac(Form("%s_cOvFrac", outPrefix.c_str()), "Nuisance parameters overlay (fractional)", 1100, 700);

        TH1D frame(
            Form("%s_frame_frac", outPrefix.c_str()),
            "Nuisance parameters (fractional, overlay);Bin index;Fractional nuisance shift",
            nBins, xMin, xMax);
        frame.SetStats(0);
        frame.SetLineColor(0);
        frame.SetFillStyle(0);
        const auto [fracMin, fracMax] = computeMinMax(true);
        frame.SetMinimum(fracMin);
        frame.SetMaximum(fracMax);
        frame.Draw("AXIS");

        TLegend leg(0.70, 0.70, 0.93, 0.93);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);

        std::vector<TH1D> hists;
        hists.reserve(nOverlay);

        for (int k = 0; k < nOverlay; ++k)
        {
            hists.emplace_back(Form("%s_nuis_frac_k%d", outPrefix.c_str(), k), "", nBins, xMin, xMax);
            hists.back().SetStats(0);
            hists.back().SetLineWidth(2);
            hists.back().SetLineColor(k + 1);
            hists.back().SetFillStyle(0);

            for (int i = 0; i < nBins; ++i)
            {
            const double nom = nominalHist.GetBinContent(i + 1);
            const double frac = (nom != 0.0) ? (nuisanceParams[k](i) / nom) : 0.0;
            hists.back().SetBinContent(i + 1, frac);
            }

            hists.back().Draw("HIST SAME");
            leg.AddEntry(&hists.back(), Form("k=%d", k), "l");
        }

        leg.Draw();
        cOvFrac.SaveAs(Form("%s_nuisance_params_fractional_overlay_first%d.pdf", outPrefix.c_str(), nOverlay));
        cOvFrac.SaveAs(Form("%s_nuisance_params_fractional_overlay_first%d.png", outPrefix.c_str(), nOverlay));
        }
    }

    TH1D SystematicsUtil::RunAllMultisimSystematics(
        const TH1D& nominalHist,
        const ROOT::RDF::RNode& rawDataFrame,
        const std::string& variableName,
        const SystematicsConfig& systConfig
    )
    {
        //Struct contains all the branch names for the genie, ppfx and reint
        std::string genieMultisimBranch = systConfig.genieMultisimBranch;
        std::string genieCVWeightBranch = systConfig.genieCVWeightBranch;
        std::string genieGlobalCVWeightBranch = systConfig.genieGlobalCVWeightBranch;
        std::string ppfxMultisimBranch = systConfig.ppfxMultisimBranch;
        std::string ppfxCVWeightBranch = systConfig.ppfxCVWeightBranch;
        std::string ppfxGlobalCVWeightBranch = systConfig.ppfxGlobalCVWeightBranch;
        std::string reintMultisimBranch = systConfig.reintMultisimBranch;
        std::string reintCVWeightBranch = systConfig.reintCVWeightBranch;
        std::string reintGlobalCVWeightBranch = systConfig.reintGlobalCVWeightBranch;

        // Create multisim universes for each systematic category

        std::vector<TH1D> genieUniverses = createMultiSimUniverses(nominalHist, rawDataFrame, variableName, genieMultisimBranch, genieCVWeightBranch, genieGlobalCVWeightBranch, 1.0);
        std::vector<TH1D> ppfxUniverses = createMultiSimUniverses(nominalHist, rawDataFrame, variableName, ppfxMultisimBranch, ppfxCVWeightBranch, ppfxGlobalCVWeightBranch, 1.0);
        std::vector<TH1D> reintUniverses = createMultiSimUniverses(nominalHist, rawDataFrame, variableName, reintMultisimBranch, reintCVWeightBranch, reintGlobalCVWeightBranch, 1.0);

        // Compute covariance matrices for each category
        TMatrixD genieCov = covarianceMatrixFromMultisims(genieUniverses, nominalHist);
        TMatrixD ppfxCov = covarianceMatrixFromMultisims(ppfxUniverses, nominalHist);
        TMatrixD reintCov = covarianceMatrixFromMultisims(reintUniverses, nominalHist);

        // Combine covariance matrices
        TMatrixD totalCov = genieCov + ppfxCov + reintCov;

        // Plot covariance matrices for debugging (optional)
        PlotMatrix(genieCov, "genie_cov");
        PlotMatrix(ppfxCov, "ppfx_cov");
        PlotMatrix(reintCov, "reint_cov");
        PlotMatrix(totalCov, "total_cov");

        // Plot fractional covariance matrix for debugging (optional)
        PlotFractionalCovarianceMatrix(genieCov, nominalHist, "genie_frac_cov");
        PlotFractionalCovarianceMatrix(ppfxCov, nominalHist, "ppfx_frac_cov");
        PlotFractionalCovarianceMatrix(reintCov, nominalHist, "reint_frac_cov");
        PlotFractionalCovarianceMatrix(totalCov, nominalHist, "total_frac_cov");

        // Prepare actual uncertainty histogram to return (diagonal elements of total covariance)
        TH1D totalVarHist(Form("%s_total_uncert", nominalHist.GetName()), Form("%s with Total Systematic Uncertainty;Variable;Events", nominalHist.GetTitle()), nominalHist.GetNbinsX(), nominalHist.GetXaxis()->GetXmin(), nominalHist.GetXaxis()->GetXmax());

        for (int i = 0; i < totalCov.GetNrows(); ++i)
        {
            double var_i = totalCov(i, i);
            //double uncert_i = std::sqrt(var_i);
            totalVarHist.SetBinContent(i + 1, var_i);
        }

        return totalVarHist;
    }

    TH1D SystematicsUtil::RunAllDetVarSystematics(
        const TH1D& nominalHist,
        const ROOT::RDF::RNode& rawDataFrame,
        const std::vector<ROOT::RDF::RNode>& detVarNodes,
        const std::vector<std::string>& detVarNames,
        const std::string& variableName,
        const std::vector<double>& detVarGlobalWeights,
        double nomHistScaleFactor,
        const std::string& weightColumn
    )
    {
        TMatrixD detVarCov = BuildDetVarCovariance(
            nominalHist,
            rawDataFrame,
            detVarNodes,
            detVarNames,
            variableName,
            detVarGlobalWeights,
            nomHistScaleFactor,
            weightColumn);

        // Plot covariance matrix for debugging (optional)
        PlotMatrix(detVarCov, "detvar_cov");
        PlotFractionalCovarianceMatrix(detVarCov, nominalHist, "detvar_frac_cov");

        // Prepare actual uncertainty histogram to return (diagonal elements of covariance)
        TH1D detVarUncertHist(Form("%s_detvar_uncert", nominalHist.GetName()), Form("%s with Detector Variation Uncertainty;Variable;Events", nominalHist.GetTitle()), nominalHist.GetNbinsX(), nominalHist.GetXaxis()->GetXmin(), nominalHist.GetXaxis()->GetXmax());

        for (int i = 0; i < detVarCov.GetNrows(); ++i)
        {
            double var_i = detVarCov(i, i);
            //double uncert_i = std::sqrt(var_i);
            detVarUncertHist.SetBinContent(i + 1, var_i);
        }
        return detVarUncertHist;
    }

} // namespace Analysis

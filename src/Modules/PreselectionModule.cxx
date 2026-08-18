#include "Modules/PreselectionModule.hxx"
#include "Utils/Plotter.hxx"
#include "Utils/SystematicsUtil.hxx"
#include "Utils/TimingUtils.hxx"

#include <TEnv.h>
#include <TDecompSVD.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TString.h>
#include <TH1D.h>
#include <TVectorD.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

using namespace Analysis;

namespace {

struct PlotRuntime {
    PlotConfig        config;
    std::vector<TH1D> histograms;
    std::vector<TMatrixD> systCovarianceMatrices;
    TMatrixD              detVarCovariance;
    bool                  hasDetVarCovariance = false;
};

std::string DefaultPlotWeightColumn(SampleType sampleType)
{
    return (IsOverlaySample(sampleType) || IsDirtSample(sampleType) || IsDetectorVariationInputSample(sampleType))
        ? "weight_cv"
        : "";
}

TMatrixD MakeZeroCovariance(const int nBins)
{
    TMatrixD covariance(nBins, nBins);
    covariance.Zero();
    return covariance;
}

TMatrixD MakeNormalisationCovariance(const TH1D& nominalHist)
{
    const int nBins = nominalHist.GetNbinsX();
    TMatrixD covariance(nBins, nBins);
    for (int i = 0; i < nBins; ++i) {
        const double contentI = nominalHist.GetBinContent(i + 1);
        for (int j = 0; j < nBins; ++j) {
            covariance(i, j) = contentI * nominalHist.GetBinContent(j + 1);
        }
    }
    return covariance;
}

struct ChiSquaredResult {
    double chiSquared = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
    int usableBins = 0;
    int droppedModes = 0;
    std::string error;
};

ChiSquaredResult CalculateCorrelatedChiSquared(
    const TVectorD& data,
    const TVectorD& prediction,
    const TMatrixD& covariance)
{
    ChiSquaredResult result;
    const int nBins = data.GetNrows();
    if (prediction.GetNrows() != nBins
        || covariance.GetNrows() != nBins
        || covariance.GetNcols() != nBins) {
        result.error = "dimension mismatch";
        return result;
    }

    std::vector<int> candidateBins;
    candidateBins.reserve(nBins);
    for (int i = 0; i < nBins; ++i) {
        if (std::isfinite(data(i))
            && std::isfinite(prediction(i))
            && std::isfinite(covariance(i, i))
            && covariance(i, i) > 0.0) {
            candidateBins.push_back(i);
        }
    }

    std::vector<int> usableBins;
    usableBins.reserve(candidateBins.size());
    for (const int i : candidateBins) {
        bool finiteRow = true;
        for (const int j : candidateBins) {
            if (!std::isfinite(covariance(i, j))) {
                finiteRow = false;
                break;
            }
        }
        if (finiteRow) usableBins.push_back(i);
    }

    result.usableBins = static_cast<int>(usableBins.size());
    if (usableBins.empty()) {
        result.error = "no bins with finite positive variance";
        return result;
    }

    const int dimension = result.usableBins;
    TMatrixD reducedCovariance(dimension, dimension);
    TVectorD residual(dimension);
    for (int i = 0; i < dimension; ++i) {
        residual(i) = data(usableBins[i]) - prediction(usableBins[i]);
        for (int j = 0; j < dimension; ++j) {
            // Numerical noise can make independently accumulated covariance
            // matrices slightly asymmetric. SVD should see a symmetric matrix.
            reducedCovariance(i, j) = 0.5 * (
                covariance(usableBins[i], usableBins[j])
                + covariance(usableBins[j], usableBins[i]));
        }
    }

    TDecompSVD decomposition(reducedCovariance);
    if (!decomposition.Decompose()) {
        result.error = "SVD decomposition failed";
        return result;
    }

    const TVectorD singularValues = decomposition.GetSig();
    const TMatrixD matrixU = decomposition.GetU();
    const TMatrixD matrixV = decomposition.GetV();
    double maximumSingularValue = 0.0;
    for (int i = 0; i < singularValues.GetNrows(); ++i) {
        maximumSingularValue = std::max(maximumSingularValue, singularValues(i));
    }
    if (!(maximumSingularValue > 0.0) || !std::isfinite(maximumSingularValue)) {
        result.error = "covariance has no positive singular values";
        return result;
    }

    constexpr double relativeTolerance = 1.0e-12;
    const double threshold = relativeTolerance * maximumSingularValue;
    TMatrixD pseudoInverse(dimension, dimension);
    pseudoInverse.Zero();
    for (int mode = 0; mode < singularValues.GetNrows(); ++mode) {
        const double singularValue = singularValues(mode);
        if (!std::isfinite(singularValue) || singularValue <= threshold) continue;
        ++result.ndf;
        const double inverseSingularValue = 1.0 / singularValue;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                pseudoInverse(i, j) += matrixV(i, mode)
                    * inverseSingularValue
                    * matrixU(j, mode);
            }
        }
    }
    result.droppedModes = dimension - result.ndf;
    if (result.ndf == 0) {
        result.error = "covariance has no retained SVD modes";
        return result;
    }

    const TVectorD weightedResidual = pseudoInverse * residual;
    result.chiSquared = residual * weightedResidual;
    if (!std::isfinite(result.chiSquared)) {
        result.error = "non-finite chi-squared";
    }
    return result;
}

} // namespace

//------------------------------------------------------------------------------
PreselectionModule::PreselectionModule(const TEnv& cfg)
    : Module(cfg)
    , fTreeName     (cfg.GetValue("Preselection.TreeName","nuselection/NeutrinoSelectionFilter"  ))
    , fRunLabel     (cfg.GetValue("Global.RunLabel","run_x") )
    , fMakePlots   (cfg.GetValue("Preselection.MakePlots", false))
{
    
    //----------------------------------------------------------------------
    // Parse the comma‑separated list of cuts, this is currently messy but was needed to ensure that the proper format of 
    // cutstring was passed to RDF::Filter
    //----------------------------------------------------------------------

    const std::string cutsString = cfg.GetValue("Preselection.Cuts", "");
    std::stringstream ssCuts{cutsString};
    std::string cutToken;

    while (std::getline(ssCuts, cutToken, ',')) {
        // trim leading/trailing whitespace (but keep quotes if present)
        cutToken.erase(0, cutToken.find_first_not_of(" \t\n\r"));
        cutToken.erase(cutToken.find_last_not_of(" \t\n\r") + 1);

        // Remove optional surrounding quotes so the expression is valid for RDataFrame::Filter
        if (cutToken.size() >= 2 &&
           ((cutToken.front() == '"'  && cutToken.back() == '"') ||
            (cutToken.front() == '\'' && cutToken.back() == '\'')))
        {
            cutToken = cutToken.substr(1, cutToken.size() - 2);
        }

        if (!cutToken.empty())
            cuts.push_back(cutToken);
    }
    
    if (cuts.empty()) {
        throw std::runtime_error("[Preselection] No cuts specified!");
    }

    std::stringstream ssKeep{cfg.GetValue("Preselection.Keep", "")};
    std::string keepItem;
    while (ssKeep >> keepItem) {
        if (keepItem.back()==',') keepItem.pop_back();
        fVarsToKeep.push_back(keepItem); 
    }

    std::stringstream ssInput{cfg.GetValue("Preselection.InputFiles", "")};
    std::string inputItem;
    while (ssInput >> inputItem) {
        if (inputItem.back()==',') inputItem.pop_back();
        fInputFiles.push_back(inputItem);
    }

    std::stringstream ssOut{cfg.GetValue("Preselection.Outputs", "")};
    std::string outItem;
    while (ssOut >> outItem) {
        if (outItem.back()==',') outItem.pop_back();
        fOutFiles.push_back(outItem);
    }

    const std::string labelsString = cfg.GetValue("Preselection.SampleLabels", "");
    std::stringstream ssLabels{labelsString};
    std::string labelToken;
    while (std::getline(ssLabels, labelToken, ',')) {
        // trim leading/trailing whitespace
        labelToken.erase(0, labelToken.find_first_not_of(" \t\n\r"));
        labelToken.erase(labelToken.find_last_not_of(" \t\n\r") + 1);

        // Remove optional surrounding quotes so labels like "Run 3 data" are kept intact
        if (labelToken.size() >= 2 &&
           ((labelToken.front() == '"'  && labelToken.back() == '"') ||
            (labelToken.front() == '\'' && labelToken.back() == '\'')))
        {
            labelToken = labelToken.substr(1, labelToken.size() - 2);
        }

        if (!labelToken.empty()) {
            fSampleLabels.push_back(labelToken);
        }
    }

    std::stringstream ssWeights{cfg.GetValue("Preselection.SampleWeights", "")};
    double weight;
    while (ssWeights >> weight) {
        fSampleWeights.push_back(weight);
    }

    const std::string sampleTypesKey = "Preselection.SampleTypes";
    fSampleTypes = ParseSampleTypes(RequireConfigValue(cfg, sampleTypesKey), sampleTypesKey);

    if (fVarsToKeep.empty()) {
        throw std::runtime_error("[Preselection] No variables to keep specified!");
    }

    const std::size_t nSamples = fInputFiles.size();
    if (nSamples == 0) {
        throw std::runtime_error("[Preselection] No input files specified!");
    }
    if (fOutFiles.size() != nSamples) {
        throw std::runtime_error("[Preselection] Outputs count (" + std::to_string(fOutFiles.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }
    if (fSampleLabels.size() != nSamples) {
        throw std::runtime_error("[Preselection] SampleLabels count (" + std::to_string(fSampleLabels.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }
    if (fSampleTypes.size() != nSamples) {
        throw std::runtime_error("[Preselection] SampleTypes count (" + std::to_string(fSampleTypes.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }
    if (fSampleWeights.size() != nSamples) {
        throw std::runtime_error("[Preselection] SampleWeights count (" + std::to_string(fSampleWeights.size())
                                 + ") does not match InputFiles count (" + std::to_string(nSamples) + ").");
    }

    const auto nDetVarCVSamples = std::count_if(
        fSampleTypes.begin(),
        fSampleTypes.end(),
        [](SampleType type) { return IsDetectorVariationCVSample(type); });
    const auto nDetVarSamples = std::count_if(
        fSampleTypes.begin(),
        fSampleTypes.end(),
        [](SampleType type) { return IsDetectorVariationSample(type); });

    if (nDetVarSamples > 0 && nDetVarCVSamples == 0) {
        throw std::runtime_error("[Preselection] SampleTypes contains detvar samples but no detvarcv sample.");
    }
    if (nDetVarCVSamples > 1) {
        throw std::runtime_error("[Preselection] SampleTypes must contain at most one detvarcv sample.");
    }
    if (nDetVarCVSamples == 1 && nDetVarSamples == 0) {
        std::cout << "[Preselection] Warning: detvarcv sample configured without detvar samples; "
                  << "detector variation systematics will be skipped.\n";
    }

    if (fMakePlots) {
        fPlotConfigs = ParsePlotConfigs(cfg, "Preselection");
        if (fPlotConfigs.empty()) {
            throw std::runtime_error("[Preselection] MakePlots is enabled but Preselection.Plots is empty.");
        }
    }
}

std::vector<std::unique_ptr<ROOT::RDataFrame>>
PreselectionModule::BuildDataFrames(const std::vector<std::string>& files,
                                     const std::string& treeName) const
{
    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfVec;
    for (const auto& fname : files) {
        auto file = TFile::Open(fname.c_str());
        if (!file || file->IsZombie()) {
            throw std::runtime_error("[Preselection] Cannot open file: " + fname);
        }
        auto tree = file->Get<TTree>(treeName.c_str());
        if (!tree) {
            throw std::runtime_error("[Preselection] Cannot find tree: " + treeName);
        }
        auto RDF = std::make_unique<ROOT::RDataFrame>(*tree);
        dfVec.push_back(std::move(RDF));
    }
    if (dfVec.empty())
        throw std::runtime_error("[Preselection] No data frames created!");

    return dfVec;
}

Long64_t PreselectionModule::EntryCount() const
{
    
    if (dfVec.size() == 0) {
        throw std::runtime_error("[Preselection] DataFrames not initialised!");
    }

    Long64_t totalEntries = 0;
    for (const auto& df : dfVec) {
        totalEntries += df->Count().GetValue();
    }
    return totalEntries;
}

void PreselectionModule::Initialise()
{
    dfVec = BuildDataFrames(fInputFiles, fTreeName);

     // One RNode per sample, initially pointing at the un-filtered DataFrame
     // This again feels messy because we initialise a set of pointers to rdataframes, but then RDF::Filter returns RNodes
     // and we have to keep track of those instead.
     
    std::vector<ROOT::RDF::RNode> nodes;
    nodes.reserve(dfVec.size());
    for (auto &dfPtr : dfVec) nodes.emplace_back(*dfPtr);

    // Apply every cut in sequence
    for (const auto &cut : cuts) {
        std::cout << "\n[Preselection] Cut: " << cut << '\n';

        for (std::size_t i = 0; i < nodes.size(); ++i) {
            auto before = nodes[i].Count().GetValue();
            std::cout << "    " << fSampleLabels[i] << " before: " << before << '\n';

            nodes[i] = nodes[i].Filter(cut);      
            auto after  = nodes[i].Count().GetValue();
            std::cout << "    " << fSampleLabels[i] << " after : " << after  << '\n';
        }
    }

    ROOT::RDF::RSnapshotOptions opt;
    opt.fMode = "RECREATE";
    opt.fCompressionAlgorithm = ROOT::kZLIB;
    opt.fCompressionLevel     = 4;

    // Test of whether any further timing alignment is needed. Use the timing utility function from TimingUtils to find the
    // mean timing offsets. Should be zero if no additional shift happens.
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        //Create vector of times from the dataframe
        std::vector<double> times = nodes[i].Take<double>("interaction_time_merged").GetValue();

        auto [A, mu, sigma, C, muError] = TimingUtils::WrappedGaussianFit(times, /*period=*/18.831, /*K=*/3);
        std::cout << "[Preselection] Timing fit results for sample " << fSampleLabels[i] << ":\n";
        std::cout << "  A     = " << A << "\n";
        std::cout << "  mu    = " << mu << "\n";
        std::cout << "  sigma = " << sigma << "\n";
        std::cout << "  C     = " << C << "\n";
        std::cout << "  muError = " << muError << "\n";

        //From this, there is a differenc of 0.5915 ns. For diagnostics, add another branch to the
        //RNode, adding this as a correction to the overlay sample only
        if (IsOverlaySample(fSampleTypes[i]) || IsDetectorVariationInputSample(fSampleTypes[i])) {
            nodes[i] = nodes[i].Redefine("interaction_time_merged",
                [mu](double t) {
                    //double corrected_time = t + 0.5915; // Apply the timing correction
                    double corrected_time = t - mu; // Apply the timing correction from the fit, which is more accurate than the fixed value
                    double remerged_time = std::fmod(corrected_time, 18.831); // Wrap around using the spill period
                    if (remerged_time < 0) remerged_time += 18.831; // Ensure non-negative
                    return remerged_time;
                },
                {"interaction_time_merged"}
            );
            std::cout << "[Preselection] Applied timing correction of " << mu << " ns to sample " << fSampleLabels[i] << ".\n";
        }
    }

    SystematicsConfig systConfig;
    systConfig.genieMultisimBranch = "weightsGenie";
    systConfig.genieCVWeightBranch = "weight_cv_untuned";
    systConfig.genieGlobalCVWeightBranch = "weight_cv";
    systConfig.ppfxMultisimBranch = "weightsPPFX";
    systConfig.ppfxCVWeightBranch = "weight_cv_noppfx";
    systConfig.ppfxGlobalCVWeightBranch = "weight_cv";
    systConfig.reintMultisimBranch = "weightsReint";
    systConfig.reintCVWeightBranch = "weight_cv";
    systConfig.reintGlobalCVWeightBranch = "weight_cv";

    std::vector<PlotRuntime> plots;
    plots.reserve(fPlotConfigs.size());
    for (const auto& plotConfig : fPlotConfigs) {
        plots.push_back({plotConfig, {}, {}});
    }

    std::vector<std::string> plotSampleLabels;
    std::vector<SampleType> plotSampleTypes;
    std::vector<double> plotSampleWeights;

    std::size_t detVarCVIndex = nodes.size();
    std::vector<std::size_t> detVarIndices;
    for (std::size_t i = 0; i < fSampleTypes.size(); ++i) {
        if (IsDetectorVariationCVSample(fSampleTypes[i])) {
            detVarCVIndex = i;
        } else if (IsDetectorVariationSample(fSampleTypes[i])) {
            detVarIndices.push_back(i);
        }
    }

    auto fillHistogramsForSample = [&](ROOT::RDF::RNode& node,
                                       const std::string& sampleLabel,
                                       SampleType sampleType,
                                       double sampleWeight) {
        for (auto& plot : plots) {
            const bool useFirstElement = plot.config.valueMode == "first_element";
            const std::string weightCol = plot.config.weightColumn.empty()
                ? DefaultPlotWeightColumn(sampleType)
                : plot.config.weightColumn;

            if (!weightCol.empty()) {
                std::cout << "    Applying " << weightCol
                          << " to plot " << plot.config.name
                          << " for sample: " << sampleLabel << "\n";
            }

            plot.histograms.push_back(
                Plotter::CreateTH1DFromRNode(
                    node,
                    plot.config.histNamePrefix + sampleLabel,
                    plot.config.column,
                    plot.config.xTitle,
                    plot.config.yTitle,
                    plot.config.nBins,
                    plot.config.xMin,
                    plot.config.xMax,
                    useFirstElement,
                    false,
                    weightCol));

            TH1D& nominalHist = plot.histograms.back();
            if (!plot.config.enableSystematics) {
                plot.systCovarianceMatrices.push_back(
                    MakeZeroCovariance(nominalHist.GetNbinsX()));
            } else if (IsOverlaySample(sampleType)) {
                std::cout << "    Computing systematic covariance matrix for sample: " << sampleLabel
                          << " and plot: " << plot.config.name << "\n";
                SystematicsUtil sysUtil;
                TMatrixD overlayCovariance = sysUtil.RunAllMultisimCovariance(
                    nominalHist,
                    node,
                    plot.config.column,
                    systConfig);
                plot.systCovarianceMatrices.push_back(std::move(overlayCovariance));
            } else if (IsDirtSample(sampleType)) {
                // A normalization uncertainty is one nuisance shared by all
                // bins, so the 100% dirt uncertainty is fully correlated.
                plot.systCovarianceMatrices.push_back(
                    MakeNormalisationCovariance(nominalHist));
            } else {
                plot.systCovarianceMatrices.push_back(
                    MakeZeroCovariance(nominalHist.GetNbinsX()));
            }
        }

        plotSampleLabels.push_back(sampleLabel);
        plotSampleTypes.push_back(sampleType);
        plotSampleWeights.push_back(sampleWeight);
    };

    auto fillDetectorVariationSystematics = [&]() {
        if (detVarCVIndex == nodes.size() || detVarIndices.empty()) {
            return;
        }

        for (auto& plot : plots) {
            if (!plot.config.enableSystematics) {
                continue;
            }

            const bool useFirstElement = plot.config.valueMode == "first_element";
            const std::string detVarColumn = useFirstElement
                ? plot.config.column + "_detvarFirstElement"
                : plot.config.column;
            const std::string weightCol = plot.config.weightColumn.empty()
                ? DefaultPlotWeightColumn(fSampleTypes[detVarCVIndex])
                : plot.config.weightColumn;

            ROOT::RDF::RNode detVarCVNode = nodes[detVarCVIndex];
            std::vector<ROOT::RDF::RNode> detVarNodes;
            std::vector<double> globalDetVarWeights; 
            std::vector<std::string> detVarNames;
            detVarNodes.reserve(detVarIndices.size());
            detVarNames.reserve(detVarIndices.size());
            globalDetVarWeights.reserve(detVarIndices.size());

            if (useFirstElement) {
                detVarCVNode = detVarCVNode.Define(
                    detVarColumn.c_str(),
                    [](const ROOT::VecOps::RVec<float>& vec) {
                        return vec.empty() ? -9999.0f : vec[0];
                    },
                    {plot.config.column.c_str()});
            }

            for (const std::size_t idx : detVarIndices) {
                ROOT::RDF::RNode detVarNode = nodes[idx];
                if (useFirstElement) {
                    detVarNode = detVarNode.Define(
                        detVarColumn.c_str(),
                        [](const ROOT::VecOps::RVec<float>& vec) {
                            return vec.empty() ? -9999.0f : vec[0];
                        },
                        {plot.config.column.c_str()});
                }
                detVarNodes.push_back(detVarNode);
                detVarNames.push_back(fSampleLabels[idx]);
                globalDetVarWeights.push_back(fSampleWeights[idx]);
            }

            double nomHistScaleFactor = fSampleWeights[detVarCVIndex];

            if (!weightCol.empty()) {
                std::cout << "    Applying " << weightCol
                          << " to detector variation systematics for plot: "
                          << plot.config.name << "\n";
            }

            TH1D detVarCVNominalHist = Plotter::CreateTH1DFromRNode(
                detVarCVNode,
                plot.config.histNamePrefix + fSampleLabels[detVarCVIndex] + "_detvarcv",
                detVarColumn,
                plot.config.xTitle,
                plot.config.yTitle,
                plot.config.nBins,
                plot.config.xMin,
                plot.config.xMax,
                false,
                false,
                weightCol);

            std::cout << "    Computing detector variation systematics for plot: "
                      << plot.config.name << " using CV sample: "
                      << fSampleLabels[detVarCVIndex] << "\n";

            SystematicsUtil sysUtil;
            const TMatrixD detVarCovariance = sysUtil.RunAllDetVarCovariance(
                detVarCVNominalHist,
                nodes[detVarCVIndex],
                detVarNodes,
                detVarNames,
                detVarColumn,
                globalDetVarWeights,
                nomHistScaleFactor,
                weightCol);
            plot.detVarCovariance.ResizeTo(detVarCovariance);
            plot.detVarCovariance = detVarCovariance;
            plot.hasDetVarCovariance = true;
        }
    };

    auto plotAllHistograms = [&]() {
        for (auto& plot : plots) {
            TH1D totalVarianceHist(
                plot.config.outputName.c_str(),
                (plot.config.outputName + " Total Variance").c_str(),
                plot.config.nBins,
                plot.config.xMin,
                plot.config.xMax);

            TMatrixD totalSystematicCovariance(plot.config.nBins, plot.config.nBins);
            totalSystematicCovariance.Zero();

            if (plot.systCovarianceMatrices.size() != plotSampleWeights.size()) {
                std::cerr << "[Preselection] ERROR: systCovarianceMatrices size (" << plot.systCovarianceMatrices.size()
                          << ") != plotted sample weights size (" << plotSampleWeights.size()
                          << ") for plot " << plot.config.outputName << "\n";
            }

            const size_t n = std::min(plot.systCovarianceMatrices.size(), plotSampleWeights.size());
            for (size_t i = 0; i < n; ++i) {
                const double weight = plotSampleWeights[i];
                TMatrixD covariance = plot.systCovarianceMatrices[i];
                covariance *= weight * weight;
                totalSystematicCovariance += covariance;
            }

            if (plot.hasDetVarCovariance) {
                totalSystematicCovariance += plot.detVarCovariance;
            }

            for (int bin = 0; bin < plot.config.nBins; ++bin) {
                totalVarianceHist.SetBinContent(
                    bin + 1,
                    totalSystematicCovariance(bin, bin));
            }

            std::cout << "Total variance for plot " << plot.config.outputName << ":\n";
            for (int bin = 1; bin <= totalVarianceHist.GetNbinsX(); ++bin) {
                std::cout << "  Bin " << bin << ": " << totalVarianceHist.GetBinContent(bin) << "\n";
            }

            TVectorD data(plot.config.nBins);
            TVectorD background(plot.config.nBins);
            data.Zero();
            background.Zero();
            TMatrixD chiSquaredCovariance = totalSystematicCovariance;
            bool hasData = false;
            bool hasBackground = false;

            for (std::size_t sample = 0; sample < plot.histograms.size(); ++sample) {
                const TH1D& histogram = plot.histograms[sample];
                const double weight = plotSampleWeights[sample];
                const double weightSquared = weight * weight;

                if (IsDataSample(plotSampleTypes[sample])) {
                    // Match FullDataMCSignalPlot's ratio behavior: use the
                    // first configured data histogram.
                    if (hasData) continue;
                    hasData = true;
                    for (int bin = 0; bin < plot.config.nBins; ++bin) {
                        data(bin) = weight * histogram.GetBinContent(bin + 1);
                        const double error = histogram.GetBinError(bin + 1);
                        chiSquaredCovariance(bin, bin) += weightSquared * error * error;
                    }
                } else if (!IsSignalSample(plotSampleTypes[sample])) {
                    hasBackground = true;
                    for (int bin = 0; bin < plot.config.nBins; ++bin) {
                        background(bin) += weight * histogram.GetBinContent(bin + 1);
                        const double error = histogram.GetBinError(bin + 1);
                        chiSquaredCovariance(bin, bin) += weightSquared * error * error;
                    }
                }
            }

            if (!hasData || !hasBackground) {
                std::cout << "[Preselection] Plot " << plot.config.name
                          << ": chi2 unavailable ("
                          << (!hasData ? "no data histogram" : "no background histogram")
                          << ")\n";
            } else {
                const ChiSquaredResult chiSquared = CalculateCorrelatedChiSquared(
                    data,
                    background,
                    chiSquaredCovariance);
                if (!chiSquared.error.empty()) {
                    std::cout << "[Preselection] Plot " << plot.config.name
                              << ": chi2 unavailable (" << chiSquared.error << ")\n";
                } else {
                    std::cout << "[Preselection] Plot " << plot.config.name
                              << ": chi2 = " << chiSquared.chiSquared
                              << ", ndf = " << chiSquared.ndf << "\n";
                    if (chiSquared.usableBins != plot.config.nBins) {
                        std::cout << "[Preselection] Plot " << plot.config.name
                                  << ": excluded "
                                  << (plot.config.nBins - chiSquared.usableBins)
                                  << " bin(s) without finite positive variance.\n";
                    }
                    if (chiSquared.droppedModes > 0) {
                        std::cout << "[Preselection] Plot " << plot.config.name
                                  << ": SVD pseudoinverse dropped "
                                  << chiSquared.droppedModes
                                  << " singular covariance mode(s) at relative tolerance 1e-12.\n";
                    }
                }
            }

            Plotter::FullDataMCSignalPlot(
                plot.histograms,
                plotSampleLabels,
                plotSampleTypes,
                plot.config.outputName,
                plot.config.logY,
                plotSampleWeights,
                0.5,
                1.5,
                &totalVarianceHist,
                plot.config.microBooNELabel,
                plot.config.showLowerCut,
                plot.config.showUpperCut);
        }
    };

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        std::cout << "\n[Preselection] Writing output for sample: " << fSampleLabels[i] << '\n';
        std::cout << "    to file: " << fOutFiles[i] << '\n';
        nodes[i].Snapshot(fTreeName, fOutFiles[i], fVarsToKeep, opt);

        //Test code for checking event matching in systematics utility. Provide overlay sample as both cv and the detvar,
        //should produce perfect matching histogram.

        if (IsOverlaySample(fSampleTypes[i])){
            SystematicsUtil sysUtil;
            ROOT::RDF::RNode cvNode = nodes[i];
            std::vector<ROOT::RDF::RNode> detVarNodes = {nodes[i]};
            std::vector<std::string> detVarNames = {"overlay"};
            std::string branchName = "NeutrinoEnergy2";
            std::string weightColumn = "weight_cv";
            TH1D cvNominalHist = Plotter::CreateTH1DFromRNode(
                cvNode,
                "cv_nominal",
                branchName,
                branchName,
                "Events",
                50,
                0.0,
                500.0,
                false,
                false,
                weightColumn);
            std::vector<TH1D> matchedHists = sysUtil.createMatchedDetVarHists(
                cvNode,
                detVarNodes,
                detVarNames,
                branchName,
                weightColumn,
                cvNominalHist);
            // Compare the matched histogram to the nominal histogram, they should be identical
            if (matchedHists.size() != 1) {
                std::cerr << "[Preselection] ERROR: Expected 1 matched histogram, got " << matchedHists.size() << "\n";
            } else {
                TH1D& matchedHist = matchedHists[0];
                bool identical = true;
                for (int bin = 1; bin <= cvNominalHist.GetNbinsX(); ++bin) {
                    double nomContent = cvNominalHist.GetBinContent(bin);
                    double matchedContent = matchedHist.GetBinContent(bin);
                    if (std::abs(nomContent - matchedContent) > 1e-6) {
                        identical = false;
                        std::cerr << "[Preselection] ERROR: Bin " << bin << " content mismatch: nominal = " << nomContent << ", matched = " << matchedContent << "\n";
                    }
                }
                if (!identical) {
                    std::cerr << "[Preselection] ERROR: Matched histogram does not match nominal histogram\n";
                }
                if (identical) {
                    std::cout << "[Preselection] SUCCESS: Matched histogram matches nominal histogram\n";
                }
            }
        }

        // End of test code for histogram matching in systematics utility

        // Test code to make a plot of the detvar cv and variations
        if (IsDetectorVariationCVSample(fSampleTypes[i])) {
            std::cout << "\n[Preselection] Creating test plot of detector variation CV and variations for sample: " << fSampleLabels[i] << '\n';
            TH1D detVarNominalTestHist = Plotter::CreateTH1DFromRNode(
                nodes[i],
                "detvar_test",
                "NeutrinoEnergy2",
                "NeutrinoEnergy2",
                "Events",
                30,
                0.0,
                500.0,
                false,
                false,
                "weight_cv");
            // Create matched histograms for all detvar variations using the systematics utility function
            SystematicsUtil sysUtil;
            std::vector<ROOT::RDF::RNode> detVarNodes;
            std::vector<std::string> detVarNames;
            for (std::size_t idx : detVarIndices) {
                detVarNodes.push_back(nodes[idx]);
                detVarNames.push_back(fSampleLabels[idx]);
            }
            std::vector<TH1D> detVarHists = sysUtil.createMatchedDetVarHists(
                nodes[i],
                detVarNodes,
                detVarNames,
                "NeutrinoEnergy2",
                "weight_cv",
                detVarNominalTestHist);

            //Output the number of entries in each histogram for diagnostics
            std::cout << "Entries in CV nominal histogram: " << detVarNominalTestHist.GetEntries()
                        << ", Entries in detvar histograms: ";
            for (const auto& hist : detVarHists) {
                std::cout << hist.GetEntries() << " ";
            }
            std::cout << "\n";
            // Plot the nominal histogram and all detvar variations on the same plot for comparison
            std::vector<TH1D> allHists = {detVarNominalTestHist};
            allHists.insert(allHists.end(), detVarHists.begin(), detVarHists.end());
            std::vector<std::string> allLabels = {fSampleLabels[i]};
            allLabels.insert(allLabels.end(), detVarNames.begin(), detVarNames.end());
            std::vector<SampleType> allTypes = {fSampleTypes[i]};
            allTypes.insert(allTypes.end(), detVarIndices.size(), SampleType::DetectorVariation);
            std::vector<double> allWeights(allLabels.size(), 1.0);
            TCanvas c("c", "c", 800, 600);
            for (size_t j = 0; j < allHists.size(); ++j) {
                allHists[j].SetLineColor(j + 1);
                allHists[j].SetMarkerColor(j + 1);
                //Standard root plot
                if (j == 0)
                    allHists[j].Draw("hist");
                else
                    allHists[j].Draw("same E1");
            }
            c.BuildLegend();
            c.SaveAs("detvar_comparison.png");
        }

        if (!fMakePlots) continue;

        if (!IsPlottableSample(fSampleTypes[i])) {
            std::cout << "[Preselection] Skipping normal plot histograms for hidden sample: "
                      << fSampleLabels[i] << " (" << SampleTypeName(fSampleTypes[i]) << ")\n";
            continue;
        }

        fillHistogramsForSample(nodes[i], fSampleLabels[i], fSampleTypes[i], fSampleWeights[i]);
    }

    if (!fMakePlots) return;
    fillDetectorVariationSystematics();
    plotAllHistograms();
}

void PreselectionModule::Finalise()
{
    // Nothing to do here
}

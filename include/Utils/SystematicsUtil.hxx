// SystematicsUtility.hxx
#pragma once

#include <map>
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <optional>

#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TDecompChol.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"

namespace Analysis {
    struct SystematicsConfig
    {
        std::string genieMultisimBranch;
        std::string genieCVWeightBranch;
        std::string genieGlobalCVWeightBranch;

        std::string ppfxMultisimBranch;
        std::string ppfxCVWeightBranch;
        std::string ppfxGlobalCVWeightBranch;

        std::string reintMultisimBranch;
        std::string reintCVWeightBranch;
        std::string reintGlobalCVWeightBranch;
    };
    class SystematicsUtil {
    public:
    explicit SystematicsUtil(std::string name = "SystematicsUtility");
    std::vector<TH1D> createMultiSimUniverses(
        const TH1D& nominalHist,
        const ROOT::RDF::RNode& dataFrame,
        const std::string& variableName,
        const std::string& multisimWeightColumn,
        const std::string& CVWeightColumn,
        const std::string& GlobalCVWeightColumn,
        const double& globalScaleFactor
    );

    TMatrixD covarianceMatrixFromMultisims(
        const std::vector<TH1D>& universeHists,
        const TH1D& nominalHist
    );

    TMatrixD combineCovarianceMatrices(
        const std::vector<TMatrixD>& matrices
    );

    std::pair<TMatrixD, TVectorD> EigenDecomposition(const TMatrixD& cov);

    void PlotMatrix(const TMatrixD& matrix, const std::string& name);

    void PlotFractionalCovarianceMatrix(
        const TMatrixD& cov,
        const TH1D& nominalHist,
        const std::string& outName = "fractional_covariance") const;

    std::vector<TVectorD> CreateNuisanceParams(const TMatrixD& cov);

    void PlotNuisanceParams(
        const std::vector<TVectorD>& nuisanceParams,
        const TH1D& nominalHist,
        const std::string& outPrefix = "nuisance",
        const int maxOverlay = 20) const;

    TH1D RunAllMultisimSystematics(
        const TH1D& nominalHist,
        const ROOT::RDF::RNode& rawDataFrame,
        const std::string& variableName,
        const SystematicsConfig& systConfig
    );

    private:
    };
}

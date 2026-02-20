#include "Utils/TimingUtils.hxx"

#include <TGraph.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <vector>
#include <tuple>
#include <set>
#include <unordered_map>


using namespace Analysis;

std::vector<double> TimingUtils::AddTimingOffset(
    const std::vector<double>& original_times,
    double offset)
{
    std::vector<double> corrected_times;
    corrected_times.reserve(original_times.size());

    for (const auto& t : original_times) {
        corrected_times.push_back(t + offset);
    }

    return corrected_times;
}

std::tuple<double,double,double,double>
TimingUtils::GaussianFit(const std::vector<double>& t)
{
    // Simple gaus, not really used
    std::vector<double> y(t.size(), 1.0);

    TGraph g(t.size(), t.data(), y.data());

    TF1 f("f", "[0]*exp(-0.5*((t-[1])/[2])^2) + [3]",
          *std::min_element(t.begin(),t.end()),
          *std::max_element(t.begin(),t.end()));

    f.SetParameters(1.0, 0.0, 1.0, 0.0);

    g.Fit(&f, "Q");

    return {
        f.GetParameter(0),
        f.GetParameter(1),
        f.GetParameter(2),
        f.GetParameter(3)
    };
}

TF1 TimingUtils::MakeWrappedGaussian(
    double period,
    int    K,
    double xmin,
    double xmax)
{
    auto wrapped_gauss = [period, K](double *x, double *p) {
        const double A     = p[0];
        const double mu    = p[1];
        const double sigma = p[2];
        const double C     = p[3];

        double sum = 0.0;
        for (int k = -K; k <= K; ++k) {
            double dx = x[0] - mu + k * period;
            sum += std::exp(-0.5 * dx * dx / (sigma * sigma));
        }
        return A * sum + C;
    };

    TF1 f("wrapped_gauss", wrapped_gauss, xmin, xmax, 4);
    f.SetParNames("A", "mu", "sigma", "C");

    return f;
}

std::tuple<double, double, double, double, double>
TimingUtils::WrappedGaussianFit(
    const std::vector<double>& t,
    double period,
    int    K)
{
    // K is the number of additional Gaussian "wraps" to include on each side
    std::vector<double> phase;
    phase.reserve(t.size());
    for (double ti : t)
        phase.push_back(std::fmod(ti, period));

    double xmin = 0.0;
    double xmax = period;

    TH1D h("h", "Spill timing;Offset [ns];Counts", 30, xmin, xmax);
    for (double x : phase) h.Fill(x);

    // Build wrapped Gaussian
    TF1 f = MakeWrappedGaussian(period, K, xmin, xmax);

    // Initial guesses
    f.SetParameters(
        h.GetMaximum(),        // A
        h.GetBinCenter(h.GetMaximumBin()), // mu
        0.1 * period,              // sigma
        h.GetMinimum()         // C
    );

    // Parameter bounds to help fit
    f.SetParLimits(1, 0.0, period);        // mu ∈ [0, period)
    f.SetParLimits(2, 0.001*period, 0.5*period);

    h.Fit(&f, "Q");

    //Plot fit for debugging
    TCanvas c("c", "c", 800, 600);
    h.Draw();
    f.Draw("same");
    c.SaveAs("wrapped_gaussian_fit.png");

    std::cout << "[TimingUtils] Wrapped Gaussian fit results:\n";
    std::cout << "  A     = " << f.GetParameter(0) << "\n";
    std::cout << "  mu    = " << f.GetParameter(1) << "\n";
    std::cout << "  sigma = " << f.GetParameter(2) << "\n";
    std::cout << "  C     = " << f.GetParameter(3) << "\n";
    std::cout << " Uncertainty on mu = " << f.GetParError(1) << "\n";

    return {
        f.GetParameter(0),
        f.GetParameter(1),
        f.GetParameter(2),
        f.GetParameter(3),
        f.GetParError(1)
    };
}

std::unordered_map<int, std::pair<double, double>>
TimingUtils::CreateRunOffsetMap(
    const std::vector<int>& run_numbers,
    const std::vector<double>& times,
    double period,
    int    K,
    int    runWindowSize)
{
    // Map of run number to (timing offset, uncertainty)
    std::cout << "[TimingUtils] Creating run offset map...\n";
    std::unordered_map<int, std::pair<double, double>> runTimeOffsetMap;

    int minRun = *std::min_element(run_numbers.begin(), run_numbers.end());
    int maxRun = *std::max_element(run_numbers.begin(), run_numbers.end());

    // Find the runs present in the data
    std::set<int> uniqueRuns(run_numbers.begin(), run_numbers.end());

    // Iterate over runs in windows of runWindowSize
    std::cout << "[TimingUtils] Creating run offset map from run " << minRun << " to " << maxRun << " in windows of " << runWindowSize << " runs.\n";
    for (int runStart = minRun; runStart <= maxRun; runStart += runWindowSize) {
        int runEnd = runStart + runWindowSize - 1;

        // Collect times for runs in this window
        std::vector<double> windowTimes;
        windowTimes.reserve(run_numbers.size());
        for (size_t i = 0; i < run_numbers.size(); ++i) {
            if (run_numbers[i] >= runStart && run_numbers[i] <= runEnd) {
                windowTimes.push_back(times[i]);
            }
        }

        // If we have collected any times, fit a Gaussian
        if (!windowTimes.empty()) {
            auto [A, mu, sigma, C, muError] = TimingUtils::WrappedGaussianFit(windowTimes, period, K);
            for (int run = runStart; run <= runEnd; ++run) {
                if (uniqueRuns.find(run) != uniqueRuns.end()) {
                    runTimeOffsetMap.emplace(run, std::make_pair(mu, sigma));
                }
            }
        }
    }

    return runTimeOffsetMap;
}

#ifndef TIMING_UTILS_HXX
#define TIMING_UTILS_HXX

#include <TGraph.h>
#include <TF1.h>
#include <TH1.h>
#include <vector>
#include <tuple>
#include <set>
#include <unordered_map>

/*--------------------------------------------------------------------------*
 *  Helper functions for ns timing corrections
 *--------------------------------------------------------------------------*/

 namespace Analysis {

class TimingUtils {
public:
    // ------------------------------------------------------------------
    //  Timing correction helpers
    // ------------------------------------------------------------------   
    static std::tuple<double, double, double, double> GaussianFit(const std::vector<double>& t);
    static std::vector<double> AddTimingOffset(const std::vector<double>& original_times, double offset);
    static std::tuple<double, double, double, double, double> WrappedGaussianFit(
        const std::vector<double>& t,
        double period,
        int    K);
    static std::unordered_map<int, std::pair<double, double>> CreateRunOffsetMap(
        const std::vector<int>& run_numbers,
        const std::vector<double>& times,
        double period,
        int    K,
        int    runWindowSize);

private:
    static TF1 MakeWrappedGaussian(
        double period,
        int    K,
        double xmin,
        double xmax);
};

} // namespace Analysis

#endif // TIMING_UTILS_HXX
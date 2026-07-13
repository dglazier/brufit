#pragma once

#include <TTree.h>
#include <RooArgSet.h>
#include <ROOT/RVec.hxx>
#include <vector>
#include <utility>

namespace HS {
namespace FIT {

    class BruMappedRhat {
    public:
        BruMappedRhat() = default;
        ~BruMappedRhat() = default;

        // --- PUBLIC CONSTANTS ---
        static constexpr Double_t kConvergenceFailure = 9999.0; 
        static constexpr Int_t    kMinRequiredEntries = 20;

        // Returns std::pair<R-hat, ESS>
        std::pair<Double_t, Double_t> CalculateDiagnostics(TTree* tree, const RooArgSet& activePars, Int_t maxPoints = 2500);

    private:
        std::vector<int> PM_NearestNeighborTour(const std::vector<ROOT::RVec<double>>& points);
        ROOT::RVec<double> PM_NearestNeighborProximityMap(const std::vector<ROOT::RVec<double>>& points);
        ROOT::RVec<double> Vehtari_RankNormalize(const ROOT::RVec<double>& x);
        
        std::pair<Double_t, Double_t> Vehtari_Diagnostics(const ROOT::RVec<double>& z);
    };

} // namespace FIT
} // namespace HS

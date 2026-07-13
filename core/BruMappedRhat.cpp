#include "BruMappedRhat.h"
#include <ROOT/RDataFrame.hxx>
#include <TMath.h>
#include <RooRealVar.h>
#include <cmath>
#include <limits>
#include <iostream>

namespace HS {
namespace FIT {

    std::vector<int> BruMappedRhat::PM_NearestNeighborTour(const std::vector<ROOT::RVec<double>>& points) {
        int N = points.size();
        std::vector<int> order;
        order.reserve(N);
        std::vector<bool> visited(N, false);

        int current = 0;
        order.push_back(current);
        visited[current] = true;

        for (int step = 1; step < N; ++step) {
            double best = std::numeric_limits<double>::max();
            int bestIdx = -1;
            for (int j = 0; j < N; ++j) {
                if (visited[j]) continue;
                auto diff = points[current] - points[j];
                double d2 = ROOT::VecOps::Sum(diff * diff); 
                if (d2 < best) { best = d2; bestIdx = j; }
            }
            order.push_back(bestIdx);
            visited[bestIdx] = true;
            current = bestIdx;
        }
        return order;
    }

    ROOT::RVec<double> BruMappedRhat::PM_NearestNeighborProximityMap(const std::vector<ROOT::RVec<double>>& points) {
        int N = points.size();
        if (N < 2) return ROOT::RVec<double>(N, 0.0);

        std::vector<int> tour = PM_NearestNeighborTour(points);
        ROOT::RVec<double> edge(N, 0.0);
        
        for (int k = 0; k < N; ++k) {
            auto diff = points[tour[k]] - points[tour[(k + 1) % N]];
            edge[k] = std::sqrt(ROOT::VecOps::Sum(diff * diff)); 
        }

        std::vector<double> prefix(2 * N + 1, 0.0);
        for (int k = 0; k < 2 * N; ++k) prefix[k + 1] = prefix[k] + edge[k % N];

        std::vector<int> posInTour(N, -1);
        for (int m = 0; m < N; ++m) posInTour[tour[m]] = m;

        double bestD = std::numeric_limits<double>::max();
        ROOT::RVec<double> bestMapped(N, 0.0);

        for (int m = 0; m < N; ++m) {
            ROOT::RVec<double> mapped(N, 0.0);
            for (int i = 0; i < N; ++i) {
                int steps = (posInTour[i] - m + N) % N;
                mapped[i] = prefix[m + steps] - prefix[m];
            }
            double Dm = 0.0;
            for (int i = 1; i < N; ++i) Dm += std::fabs(mapped[i] - mapped[i - 1]);
            if (Dm < bestD) { bestD = Dm; bestMapped = mapped; }
        }
        return bestMapped;
    }

    ROOT::RVec<double> BruMappedRhat::Vehtari_RankNormalize(const ROOT::RVec<double>& x) {
        int n = x.size();
        auto indices = ROOT::VecOps::Argsort(x); 
        ROOT::RVec<double> z(n);
        for(int i = 0; i < n; ++i) {
            double rank = i + 1.0;
            z[indices[i]] = TMath::NormQuantile((rank - 0.375) / (n + 0.25));
        }
        return z;
    }

    std::pair<Double_t, Double_t> BruMappedRhat::Vehtari_Diagnostics(const ROOT::RVec<double>& z) {
        int n = z.size();
        int half = n / 2;
        
        if (half < 2) return {kConvergenceFailure, 0.0}; 

        ROOT::RVec<double> z1(z.begin(), z.begin() + half);
        ROOT::RVec<double> z2(z.begin() + half, z.end());

        double m1 = ROOT::VecOps::Mean(z1);
        double m2 = ROOT::VecOps::Mean(z2);
        
        double v1 = ROOT::VecOps::Sum((z1 - m1) * (z1 - m1)) / (half - 1.0);
        double v2 = ROOT::VecOps::Sum((z2 - m2) * (z2 - m2)) / (half - 1.0);

        double W = 0.5 * (v1 + v2);
        if (W <= 0.0) return {kConvergenceFailure, 0.0}; 

        double grand_mean = (m1 + m2) / 2.0;
        double B = half * ((m1 - grand_mean)*(m1 - grand_mean) + (m2 - grand_mean)*(m2 - grand_mean));
        double var_plus = ((half - 1.0) / half) * W + B / half;
        double rhat = std::sqrt(var_plus / W);

        // --- NEW: Calculate Effective Sample Size (ESS) via Variogram ---
        double sum_rho = 0.0;
        for (int t = 1; t < half; ++t) {
            double vt = 0.0;
            for(int i = t; i < half; ++i) {
                vt += (z1[i] - z1[i-t])*(z1[i] - z1[i-t]);
                vt += (z2[i] - z2[i-t])*(z2[i] - z2[i-t]);
            }
            vt /= (2.0 * (half - t));
            
            double rho_t = 1.0 - (vt / (2.0 * var_plus));
            if (rho_t < 0.0) break; // Geyer's initial positive sequence estimator truncation
            
            sum_rho += rho_t;
        }
        
        double ess = (2.0 * half) / (1.0 + 2.0 * sum_rho);
        return {rhat, ess};
    }

    std::pair<Double_t, Double_t> BruMappedRhat::CalculateDiagnostics(TTree* tree, const RooArgSet& activePars, Int_t maxPoints) {
        if (!tree) return {kConvergenceFailure, 0.0};
        
        Long64_t nEntries = tree->GetEntries();
        if (nEntries < kMinRequiredEntries) return {kConvergenceFailure, 0.0};

        Int_t nPars = activePars.getSize();
        if (nPars == 0) return {kConvergenceFailure, 0.0};
        
        std::vector<Double_t> vals(nPars, 0.0);
        
        tree->ResetBranchAddresses();
        int idx = 0;
        for (auto* arg : activePars) {
            tree->SetBranchAddress(arg->GetName(), &vals[idx]);
            idx++;
        }

        int stride = (nEntries > maxPoints) ? std::ceil((double)nEntries / maxPoints) : 1;
        std::vector<std::vector<double>> rawPoints;
        rawPoints.reserve(nEntries / stride + 1);

        for (Long64_t i = 0; i < nEntries; i += stride) {
            tree->GetEntry(i);
            rawPoints.emplace_back(vals.begin(), vals.end()); 
        }
        tree->ResetBranchAddresses();

        int numPoints = rawPoints.size();

        // Z-SCORE NORMALIZATION
        if (numPoints > 0 && nPars > 0) {
            std::vector<double> means(nPars, 0.0);
            std::vector<double> stdDevs(nPars, 0.0);

            for (int p = 0; p < nPars; ++p) {
                for (int i = 0; i < numPoints; ++i) means[p] += rawPoints[i][p];
                means[p] /= numPoints;
            }
            for (int p = 0; p < nPars; ++p) {
                for (int i = 0; i < numPoints; ++i) stdDevs[p] += (rawPoints[i][p] - means[p]) * (rawPoints[i][p] - means[p]);
                stdDevs[p] = std::sqrt(stdDevs[p] / numPoints);
                if (stdDevs[p] < 1e-12) stdDevs[p] = 1.0; 
            }
            for (int i = 0; i < numPoints; ++i) {
                for (int p = 0; p < nPars; ++p) rawPoints[i][p] = (rawPoints[i][p] - means[p]) / stdDevs[p];
            }
        }
        
        std::vector<ROOT::RVec<double>> points;
        points.reserve(numPoints);
        for (int i = 0; i < numPoints; ++i) points.emplace_back(rawPoints[i].data(), nPars);

        ROOT::RVec<double> mapped = PM_NearestNeighborProximityMap(points);
        ROOT::RVec<double> z = Vehtari_RankNormalize(mapped);
        
        return Vehtari_Diagnostics(z);
    }

} // namespace FIT
} // namespace HS

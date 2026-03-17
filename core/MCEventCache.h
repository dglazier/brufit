/**
 * @file MCEventCache.h
 * @brief Standalone data cache for Monte Carlo integration.
 */

#pragma once

#include <TTree.h>
#include <TString.h>
#include "Weights.h"
#include <vector>

namespace bru {

    class MCEventCache {
    public:
        MCEventCache() = default;
        ~MCEventCache() = default;

        // Disallow copying to guarantee memory is never accidentally duplicated
        MCEventCache(const MCEventCache&) = default;
        MCEventCache& operator=(const MCEventCache&) = default;

        /** * @brief Loads and flattens a TTree into contiguous memory vectors.
         * @param fCut Passed by reference so it can be amended if branches are missing.
         */
        bool LoadTree(TTree* tree, TString& fCut,
                      const std::vector<TString>& varNames,
                      const std::vector<TString>& catNames,
                      const TString& fTruthPrefix,
                      const HS::FIT::WeightsConfig& wgtsConf,
                      TTree* mcGenTree = nullptr);

        // Extracted Contiguous Data Arrays
        std::vector<Float_t> _vecReal;
        std::vector<Float_t> _vecRealGen;
        std::vector<Float_t> _vecRealMCGen;

        std::vector<Int_t> _vecCat;
        std::vector<Int_t> _vecCatGen;
        std::vector<Int_t> _vecCatMCGen;

        std::vector<Float_t> _EvWeights;
        std::vector<Long64_t> _TreeEntryNumber;

        // Metadata
        Long64_t _NTreeEntries = 0;
        Long64_t _NMCGenTreeEntries = 0;
        size_t _Nvars = 0;
        size_t _Ncats = 0;
        bool _HasMCGenTree = false;
        bool _UseEvWeights = false;

        inline double GetWeight(Long64_t entry) const {
            return _UseEvWeights ? _EvWeights[entry] : 1.0;
        }
    };

} // namespace bru

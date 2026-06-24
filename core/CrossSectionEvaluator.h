////////////////////////////////////////////////////////////////
///
/// Class:       CrossSectionEvaluator
/// Description: Engine for evaluating yields, acceptances, and 
///              generic cross sections using user-injected 
///              lambda functions over RooDataSets.
///

#pragma once

#include <utility>
#include <vector>
#include <map>
#include <functional>
#include <memory>
#include <TTree.h>
#include <TString.h>

#include "FitManager.h"
#include "Setup.h"

namespace HS {
namespace FIT {

    struct CSData {
        Double_t yield = 0.;
        Double_t yield_err = 0.;
        Double_t acceptance = 0.;
        Double_t acceptance_err = 0.;
        std::map<TString, BinVolume> kinematics;
    };

    enum class AcceptanceMode {
        kFullPDF,         
        kPhaseSpaceScale  
    };

    class CrossSectionEvaluator : public FitManager {

    public:
        using CSFormula_t = std::function<std::pair<Double_t, Double_t>(const CSData&)>;

        CrossSectionEvaluator() = default;
        CrossSectionEvaluator(const CrossSectionEvaluator&) = default;
        CrossSectionEvaluator(const FitManager& fm, TString outDir="", TString resultFile="")
            : FitManager(fm), fResultDir(std::move(outDir)), fResultFileName(std::move(resultFile)) {};
        CrossSectionEvaluator(CrossSectionEvaluator&&) = delete;
        ~CrossSectionEvaluator() override = default;
        CrossSectionEvaluator& operator=(const CrossSectionEvaluator& other) = default;
        CrossSectionEvaluator& operator=(CrossSectionEvaluator&& other) = delete;

        Bool_t Run() override;
        void SaveResults() override;
        void InitTree();

        void SetCrossSectionFormula(CSFormula_t formula) { fCSFormula = std::move(formula); }
        void SetAcceptanceMode(AcceptanceMode mode) { fAccMode = mode; }
        
        // Scale factor for generated events (e.g. 10.0 if generated is 10x smaller than luminosity)
        void SetGeneratedScale(Double_t scale) { fGenScale = scale; }

        void SetNThreads(Int_t nThreads) { fNThreads = nThreads; }

        void SampleAcceptance(Bool_t b = kTRUE, Int_t nMinuitSamples = 100) { 
            fSampleAcceptance = b; 
            fNSamplesMinuit = nMinuitSamples;
        }

        void SetResultDir(TString name) { fResultDir = std::move(name); }
        void SetResultFileName(TString name) { fResultFileName = std::move(name); }

    private:

        Bool_t LoadDefaultBins();
        Bool_t LoadFitResult(Int_t globalBinIndex);
        void CalcYield(Int_t globalBinIndex, CSData& binData);
        void CalcAcceptance(Int_t globalBinIndex, CSData& binData);

        TString fResultDir;
        TString fResultFileName;
        CSFormula_t fCSFormula;
        
        AcceptanceMode fAccMode = AcceptanceMode::kFullPDF;
        Double_t fGenScale = 1.0;
        Int_t fNThreads = 4;
        
        Bool_t fSampleAcceptance = kFALSE;
        Int_t fNSamplesMinuit = 100;

        TTree* fOutTree = nullptr;
        Double_t fOutCrossSection = 0.;
        Double_t fOutCrossSectionErr = 0.;
        Double_t fOutYield = 0.;
        Double_t fOutYieldErr = 0.;
        Double_t fOutAcceptance = 0.;
        Double_t fOutAcceptanceErr = 0.;
        Int_t fOutGlobalBin = 0;

        std::vector<Double_t> fBinCenters;
        std::vector<Double_t> fBinWidths;

        ClassDefOverride(HS::FIT::CrossSectionEvaluator, 1);
    }; 

} // namespace FIT
} // namespace HS

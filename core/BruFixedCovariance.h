#pragma once

#include "BruMcmc.h" // Ensures BruMcmcCovariance is available
#include <TString.h>
#include <TMatrixDSym.h>
#include <TList.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <vector>
#include <string>

namespace HS {
namespace FIT {

// ========================================================================
    // Helper Class: Validates, aligns, and loads external Covariance Matrices
    // ========================================================================
    class BruCovarianceReader {
    public:
      
        BruCovarianceReader() = default;
        ~BruCovarianceReader() = default;

        // Loads the matrix and its associated parameter names from the file
        Bool_t Load(const TString& filePath, const TString& matrixName, const TString& listName);

        // Validates dimensions and physically reorders the matrix to match the target RooFit order
        Bool_t AlignAndValidate(const RooArgSet& currentPars);

        // Returns a const reference to the fully aligned matrix object
        const TMatrixDSym& GetMatrix() const { return fMatrix; }

        // Generic, non-RooFit dependent method to save a matrix and its parameters
        static Bool_t SaveThis(const TString& filePath, const TString& matrixName, 
                               const TMatrixDSym& mat, const std::vector<std::string>& params);

    private:
        TMatrixDSym fMatrix;
        std::vector<std::string> fSavedNames;
    };  
    // ========================================================================
    // MCMC Class: Phase 1 (Seq Burn-in) -> Phase 3 (Fixed Covariance)
    // ========================================================================
    class BruMcmcFixedCovariance : public BruMcmcCovariance {
    public:
      BruMcmcFixedCovariance(const TString& covFilePath, const TString& matrixName, 
			     std::vector<Int_t> Niters = {1000, 10000}, Int_t Nburn = 10, Float_t norm = 0.01, 
			     float target = 0.234, float accmin = 0.15, float accmax = 0.35);
      BruMcmcFixedCovariance() : BruMcmcCovariance() {};
      
        ~BruMcmcFixedCovariance() override = default;

        void Run(Setup& setup, RooAbsData& fitdata) override;
        
        void SetCovarianceList(const TString& listName) { fListName = listName; }
        void SetShrinkage(Double_t shrinkage) { _shrinkage = shrinkage; }

    private:
        BruCovarianceReader _covReader;

        TString fCovFilePath;
        TString fMatrixName;
        TString fListName;       

        Double_t _shrinkage = 0.15; 
        
        ClassDefOverride(HS::FIT::BruMcmcFixedCovariance, 1);
    };

} // namespace FIT
} // namespace HS

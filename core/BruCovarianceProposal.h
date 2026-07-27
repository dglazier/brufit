////////////////////////////////////////////////////////////////
///
///Class:           BruCovarianceProposal    
///Description:  Fast Multivariate Gaussian proposal using Cholesky Decomposition
///           
#pragma once

#include "BruSequentialProposal.h"
#include <TMath.h>
#include <RooArgSet.h>
#include <RooMsgService.h>
#include <RooRealVar.h>
#include <TMatrixDSym.h>
#include <TMatrixD.h>
#include <vector>

namespace HS{
  namespace FIT{

    class BruCovarianceProposal : public BruSequentialProposal{

    public:
      BruCovarianceProposal() : BruSequentialProposal{} {}
      
      BruCovarianceProposal(float scale,float target=0.234,float accmin=0.15,float accmax=0.35);
      
      void Propose(RooArgSet& xPrime, RooArgSet& x) override;
      bool IsSymmetric(RooArgSet& x1, RooArgSet& x2) override ;
      double GetProposalDensity(RooArgSet& x1, RooArgSet& x2) override;
 
      virtual void Reset() {
        _baseMatrix.Clear(); _baseMatrix.ResizeTo(0,0);
        _covMatrix.Clear(); _covMatrix.ResizeTo(0,0);
        _lMatrix.Clear(); _lMatrix.ResizeTo(0,0);
        _xVec.clear();
        _varCache.clear();
        _primeCache.clear(); // Clear the fast pointer cache
      }

      void SetCovariance(const TMatrixDSym& mat, const RooArgSet& vars);
      void UpdateCholesky();
      
      TMatrixDSym GetCovariance() { return _covMatrix; }
      TMatrixDSym GetBaseMatrix() { return _baseMatrix; }

      // Cleanly updates the scale and rebuilds the matrix without losing base covariance
      void ApplyNewScale(Float_t newScale) {
          SetScale(newScale);
          UpdateCholesky();
      }

      Bool_t CheckStepSize(Float_t acceptance) override {
        if(_tuneCovStep == kFALSE) return kTRUE;
        auto doExit = BruSequentialProposal::CheckStepSize(acceptance);
        UpdateCholesky();
        return doExit;
      }
      
      void TuneCovarianceStep(Bool_t tune) { _tuneCovStep = tune; }

      // Set the yields so they are masked from global step scaling
      void SetYields(const RooArgList& yields) { 
          _yieldPars.removeAll(); 
          _yieldPars.add(yields); 
      }

    private:
      TMatrixDSym _baseMatrix; 
      TMatrixDSym _covMatrix;  
      TMatrixD _lMatrix;       

      RooArgList _xVec;
      RooArgList _yieldPars; // Track yield parameters
      std::vector<RooRealVar*> _varCache;   // Fast cache for base variables
      std::vector<RooRealVar*> _primeCache; // Fast cache for xPrime targets
       
      Bool_t _tuneCovStep = kFALSE;

      ClassDefOverride(BruCovarianceProposal,1) 
    };

  }
}

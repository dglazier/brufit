#include "BruCovarianceProposal.h"
#include <RooStats/RooStatsUtils.h>
#include <RooRandom.h>
#include <TVectorD.h>
#include <TDecompChol.h>
#include <TError.h>
#include <cmath>

namespace HS{
  namespace FIT{

    BruCovarianceProposal::BruCovarianceProposal(float scale,float target,float accmin,float accmax) :
      BruSequentialProposal(scale, target, accmin, accmax)
    {
    }

    void BruCovarianceProposal::SetCovariance(const TMatrixDSym& mat, const RooArgSet& vars){
      Reset();

      if (vars.getSize() != mat.GetNcols()) {
        std::cerr << "BruCovarianceProposal::SetCovariance: Variables not same as matrix!" << std::endl;
        exit(0);
      }
      
      _baseMatrix.ResizeTo(mat.GetNrows(), mat.GetNcols());
      _baseMatrix = mat;

      if(_xVec.getSize() == 0) {
        for (auto *r : static_range_cast<RooRealVar *>(vars)){
          _xVec.add(*r);
          _varCache.push_back(dynamic_cast<RooRealVar*>(r));
        }
      }
      
      UpdateCholesky();
    }

    void BruCovarianceProposal::UpdateCholesky() {
      if (_baseMatrix.GetNrows() == 0) return;

      int d = _baseMatrix.GetNrows();
      double gelmanScale = (2.382 * 2.382) / (double)d; 
 
      _covMatrix.ResizeTo(_baseMatrix.GetNrows(), _baseMatrix.GetNcols());
      
      // 1. Apply the global MCMC tuning scale to the entire matrix
      _covMatrix = _baseMatrix * static_cast<Double_t>(StepSizeFactor() * gelmanScale);
      
      // 2. THE YIELD MASK: Force yields to ignore the global scale and off-diagonals
      // A pure 1D Gelman jump is strictly (2.382^2) / 1.0 = 5.673924
      double yieldGelmanScale = 5.673924;
      
      for (size_t i = 0; i < _varCache.size(); ++i) {
          if (_yieldPars.contains(*_varCache[i])) {
              // Restore the raw diagonal variance and apply the 1D scale
              _covMatrix(i, i) = _baseMatrix(i, i) * yieldGelmanScale;
              
              // Zero out any residual off-diagonal correlations just to be absolutely safe
              for (size_t j = 0; j < _varCache.size(); ++j) {
                  if (i != j) {
                      _covMatrix(i, j) = 0.0;
                      _covMatrix(j, i) = 0.0;
                  }
              }
          }
      }
      
      // 3. Adaptive Jitter for numerical stability during Cholesky
      double maxDiag = 1e-12;
      for (int i = 0; i < _covMatrix.GetNrows(); ++i) {
          if (_covMatrix(i, i) > maxDiag) maxDiag = _covMatrix(i, i);
      }
      double jitter = maxDiag * 1e-8; 
      
      bool decomposed = false;
      
      int oldLevel = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kFatal;
      
      for (int attempts = 0; attempts < 10; attempts++) {
          TDecompChol chol(_covMatrix);
          if (chol.Decompose()) {
              _lMatrix.ResizeTo(_covMatrix.GetNrows(), _covMatrix.GetNcols());
              _lMatrix = chol.GetU(); 
              _lMatrix.T(); 
              decomposed = true;
              break; 
          } else {
              for (int i = 0; i < _covMatrix.GetNrows(); ++i) {
                  _covMatrix(i, i) += jitter;
              }
              jitter *= 10.0; 
          }
      }
      
      gErrorIgnoreLevel = oldLevel;

      if (!decomposed) {
          std::cerr << "BruCovarianceProposal: FATAL - Matrix is fundamentally singular even with maximum jitter!" << std::endl;
      }
    }

    void BruCovarianceProposal::Propose(RooArgSet& xPrime, RooArgSet& x )
    {
      RooStats::SetParameters(&x, &xPrime); // Start at current location
      int n = _varCache.size();

      // --- OPTIMIZATION: Smart Pointer Cache ---
      // We do a single, ultra-fast validation check on the first parameter.
      // If the memory address doesn't match, we know MakeChain() launched 
      // a new sequence with fresh variables, and we rebuild the cache safely.
      bool rebuildCache = _primeCache.empty();
      if (!rebuildCache && n > 0) {
          if (_primeCache[0] != xPrime.find(_varCache[0]->GetName())) {
              rebuildCache = true;
          }
      }

      if (rebuildCache) {
          _primeCache.clear();
          _primeCache.reserve(n);
          for (int i = 0; i < n; ++i) {
              _primeCache.push_back(dynamic_cast<RooRealVar*>(xPrime.find(_varCache[i]->GetName())));
          }
      }

      // 1. Generate independent standard normal steps
      TVectorD z(n);
      for (int i = 0; i < n; ++i) {
          z[i] = RooRandom::gaussian();
      }

      // 2. Instantly correlate them via Cholesky decomposition (L * z)
      TVectorD step = _lMatrix * z;

      // 3. Apply steps to parameters safely
      for (int i = 0; i < n; ++i) {
          RooRealVar* primeVar = _primeCache[i]; // Completely bypasses string lookup!
          if (!primeVar) continue;

          double val = primeVar->getVal();
          double max = primeVar->getMax();
          double min = primeVar->getMin();
          double len = max - min;
          double s = step[i];

          if (fCyclicPars.contains(*primeVar)) {
              // --- CYCLIC WRAP AROUND ---
              double center = min + (len / 2.0);
              primeVar->setVal(center + std::remainder(val + s - center, len));
          } else {
              // --- HARD BOUNDARY REFLECTION ---
              while ((val + s > max) || (val + s < min)) {
                  if (val + s > max) s = 2 * (max - val) - s;
                  if (val + s < min) s = 2 * (min - val) - s;
              }
              primeVar->setVal(val + s);
          }
      }
    }
    
    bool BruCovarianceProposal::IsSymmetric(RooArgSet& , RooArgSet& ) {
      return true;
    }
 
    double BruCovarianceProposal::GetProposalDensity(RooArgSet&, RooArgSet&) {
      return 1.0; 
    }
  }
}

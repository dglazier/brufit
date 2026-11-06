#include "BruCovarianceProposal.h"
#include <RooStats/RooStatsUtils.h>
#include <RooRandom.h>
#include <TVectorD.h>
#include <TDecompChol.h>
#include <TError.h>
#include <cmath>
#include <algorithm>

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

      int nTotal = _baseMatrix.GetNrows();
      int d = (_gibbsBlockSize > 0 && _gibbsBlockSize < nTotal) 
              ? _gibbsBlockSize : nTotal;
      double gelmanScale = (2.382 * 2.382) / (double)d; 
      // Yields are updated on every step even in Gibbs mode, so keep yield scale locked to nTotal
      double gelmanScaleYield = (2.382 * 2.382) / (double)nTotal;
 
      _covMatrix.ResizeTo(nTotal, nTotal);
      
      // 1. Apply the global MCMC tuning scale to the physics matrix
      double scale = static_cast<Double_t>(StepSizeFactor() * gelmanScale);
      _covMatrix = _baseMatrix * scale;
      
      // 2. THE YIELD MASK: Lock yields to gelmanScaleYield (ignoring StepSizeFactor) and zero off-diagonals
      for (size_t i = 0; i < _varCache.size(); ++i) {
          if (_yieldPars.contains(*_varCache[i])) {
              for (size_t j = 0; j < _varCache.size(); ++j) {
                  if (i != j) {
                      _covMatrix(i, j) = 0.0;
                      _covMatrix(j, i) = 0.0;
                  } else {
                      _covMatrix(i, i) = _baseMatrix(i, i) * gelmanScaleYield;
                  }
              }
          }
      }
      
      // 3. Adaptive Jitter for numerical stability
      double relJitter = 1e-6; 
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
                  _covMatrix(i, i) += std::abs(_covMatrix(i, i)) * relJitter + 1e-12;
              }
              relJitter *= 10.0; 
          }
      }
      
      gErrorIgnoreLevel = oldLevel;

      if (!decomposed) {
          std::cerr << "BruCovarianceProposal: FATAL - Matrix is fundamentally singular even with maximum jitter!" << std::endl;
      }
    }

    void BruCovarianceProposal::Propose(RooArgSet& xPrime, RooArgSet& x )
    {
      RooStats::SetParameters(&x, &xPrime); 
      int n = _varCache.size();

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

      // 2. Correlate them via Cholesky decomposition (L * z)
      TVectorD step = (_lMatrix * z);
      
      // 3. Determine active Gibbs block if _gibbsBlockSize > 0
      int startIdx = 0;
      int endIdx = n;
      if (_gibbsBlockSize > 0 && _gibbsBlockSize < n) {
          int numBlocks = (n + _gibbsBlockSize - 1) / _gibbsBlockSize;
          startIdx = RooRandom::integer(numBlocks) * _gibbsBlockSize;
          endIdx = std::min(startIdx + _gibbsBlockSize, n);
      }

      // 4. Apply steps only to the active block (plus yields)
      for (int i = 0; i < n; ++i) {
          bool isYieldPar = _yieldPars.contains(*_varCache[i]);
          bool inActiveBlock = (i >= startIdx && i < endIdx) || isYieldPar;
          if (!inActiveBlock) continue;

          RooRealVar* primeVar = _primeCache[i];
          if (!primeVar) continue;

          double val = primeVar->getVal();
          double max = primeVar->getMax();
          double min = primeVar->getMin();
          double len = max - min;
          double s = step[i];

          if (fCyclicPars.contains(*primeVar)) {
              double center = min + (len / 2.0);
              primeVar->setVal(center + std::remainder(val + s - center, len));
          } else {
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

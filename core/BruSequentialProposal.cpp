#include "BruSequentialProposal.h"
#include <RooStats/RooStatsUtils.h>
#include <RooRandom.h>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace HS{
  namespace FIT{

    BruSequentialProposal::BruSequentialProposal(float scale,float target,float accmin,float accmax) :
      ProposalFunction{},
      fScale{scale},
      fTargetAcc{target},
      fMinAcc{accmin},
      fMaxAcc{accmax}
    {
    }

    void BruSequentialProposal::SetParStepScales(const RooArgList& pars, const std::vector<Double_t>& sigmas) {
      _parStepNames.clear();
      _parStepSigmas.clear();
      _cachedStepWidths.clear();

      int n = std::min(static_cast<int>(pars.getSize()), static_cast<int>(sigmas.size()));
      for (int i = 0; i < n; ++i) {
        _parStepNames.emplace_back(pars.at(i)->GetName());
        _parStepSigmas.push_back(sigmas[i]);
      }
    }
 
    void BruSequentialProposal::Propose(RooArgSet& xPrime, RooArgSet& x )
    {
      RooStats::SetParameters(&x, &xPrime);
      int n = xPrime.getSize();
      
      // 1. FAST C++ POINTER & STEP-WIDTH CACHE
      if (_varCache.empty() || _varCache.size() != (size_t)n) {
          _varCache.clear();
          _cachedStepWidths.clear();
          for (auto *var : static_range_cast<RooRealVar *>(xPrime)) {
              _varCache.push_back(var);
              double len = var->getMax() - var->getMin();
              double width = len;
              if (!_parStepNames.empty()) {
                  std::string vname = var->GetName();
                  for (size_t k = 0; k < _parStepNames.size(); ++k) {
                      if (_parStepNames[k] == vname && _parStepSigmas[k] > 0.0) {
                          width = std::clamp(_parStepSigmas[k], 1e-4 * len, 0.5 * len);
                          break;
                      }
                  }
              }
              _cachedStepWidths.push_back(width);
          }
      }

      int j = int( floor(RooRandom::uniform()*n) );
      
      // 2. THE VALLEY TRAP MIXTURE
      bool force1D = (!_isNotSequential) || (_inValley && RooRandom::uniform() < 0.5);

      // 3. DETERMINE THE GIBBS BLOCK
      int startIdx = 0;
      int endIdx = n;
      
      if (!force1D && _gibbsBlockSize > 0) {
          int numBlocks = (n + _gibbsBlockSize - 1) / _gibbsBlockSize;
          startIdx = RooRandom::integer(numBlocks) * _gibbsBlockSize;
          endIdx = std::min(startIdx + _gibbsBlockSize, n);
      }

      // 4. APPLY THE PROPOSAL
      for (int i = 0; i < n; ++i) {
        
        bool updateThis = false;
        if (force1D) {
            updateThis = (i == j); // Strict 1D Step
        } else if (_gibbsBlockSize > 0) {
            updateThis = (i >= startIdx && i < endIdx); // Block-Wise Gibbs Step
        } else {
            updateThis = true; // Global ND Step
        }

        if (updateThis) {
          RooRealVar* var = _varCache[i];
          double val = var->getVal(), max = var->getMax(), min = var->getMin(), len = max - min;
          double baseWidth = (i < (int)_cachedStepWidths.size()) ? _cachedStepWidths[i] : len;
          bool isCyclic = fCyclicPars.contains(*var);
          
          auto step = RooRandom::gaussian() * baseWidth * fScale;
          
          if (isCyclic) {
              // --- CYCLIC WRAP AROUND ---
              double center = min + (len / 2.0);
              var->setVal(center + std::remainder(val + step - center, len));
          } else {
              // --- LEGACY STICKY BOUNDARY (RE-DRAW) ---
              int redraws = 0;
              while (((val + step > max) || (val + step < min)) && redraws < 100) {
                  step = RooRandom::gaussian() * baseWidth * fScale;
                  redraws++;
              }
              if (val + step > max) step = max - val;
              if (val + step < min) step = min - val;
              var->setVal(val + step);
          }
        }
      }
    }

    Bool_t BruSequentialProposal::CheckStepSize(Float_t acceptance){
      if(acceptance < fMinAcc || acceptance > fMaxAcc){
        Double_t oldScale = fScale; 
        
        // Damped scaling prevents 40D-77D acceptance ping-ponging during Phase 2
        if (acceptance <= 0.005f) {
            fScale *= 0.25f;
        } else {
            Double_t ratio = std::sqrt(static_cast<Double_t>(acceptance) / static_cast<Double_t>(fTargetAcc));
            ratio = std::clamp(ratio, 0.4, 2.0);
            fScale *= static_cast<Float_t>(ratio);
        }
        
        if (fScale > 3.0f) {
            fScale = 3.0f; 
        }

        if(fScale < fMinScale){
          fScale = fMinScale; 
          fNminScale++;
          
          if (_isNotSequential) {
              _inValley = kTRUE;
              std::cout << "  BruSequentialProposal: Valley Detected! Mixing in 1D steps to recover." << std::endl;
          }
        } else {
            _inValley = kFALSE;
        }
        
        std::cout<<"BruSequentialProposal::CheckStepSize Changed to "<<fScale<<" (min "<<fMinScale<<") from "<< oldScale <<std::endl;

        if(fNminScale > 10){
          fNminScale=0;
          std::cerr<<"BruSequentialProposal::CheckStepSize cannot get within allowed acceptance limits. Must be stuck. Will exit."<<std::endl;
          return kFALSE;
        }
        return kTRUE;
      }
      
      _inValley = kFALSE; 
      fNminScale = 0; 
      return kTRUE;
    }
    
    bool BruSequentialProposal::IsSymmetric(RooArgSet& , RooArgSet& ) {
      return true;
    }
 
    double BruSequentialProposal::GetProposalDensity(RooArgSet& , RooArgSet& )
    {
      return 1.0; 
    }

  }
}

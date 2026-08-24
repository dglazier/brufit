#include "BruSequentialProposal.h"
#include <RooStats/RooStatsUtils.h>
#include <RooRandom.h>
#include <cmath>
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
 
    void BruSequentialProposal::Propose(RooArgSet& xPrime, RooArgSet& x )
    {
      RooStats::SetParameters(&x, &xPrime);
      int n = xPrime.getSize();
      
      // 1. FAST C++ POINTER CACHE
      if (_varCache.empty() || _varCache.size() != (size_t)n) {
          _varCache.clear();
          for (auto *var : static_range_cast<RooRealVar *>(xPrime)) {
              _varCache.push_back(var);
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
          bool isCyclic = fCyclicPars.contains(*var);
          
          auto step = RooRandom::gaussian() * len * fScale;
          
          if (isCyclic) {
              // --- CYCLIC WRAP AROUND ---
              double center = min + (len / 2.0);
              var->setVal(center + std::remainder(val + step - center, len));
          } else {
              // --- TRUE BILLIARD REFLECTION (CPU FIX) ---
              // Replaces the infinite re-roll 'while' loop with a deterministic mathematical bounce
              double newVal = val + step;
              while (newVal > max || newVal < min) {
                  if (newVal > max) newVal = 2.0 * max - newVal; 
                  if (newVal < min) newVal = 2.0 * min - newVal; 
              }
              var->setVal(newVal);
          }
        }
      }
    }

    Bool_t BruSequentialProposal::CheckStepSize(Float_t acceptance){
      if(acceptance < fMinAcc || acceptance > fMaxAcc){
        Double_t oldScale = fScale; 
        
        Double_t acc = acceptance > 0 ? acceptance : 0.01;
        fScale *= (acc)/(fTargetAcc);
        
        // --- MAXIMUM SCALE CLAMP (ASYMPTOTE FIX) ---
        // If the proposal spans the entire boundary volume, shrinking it further is impossible.
        // Cap the multiplier at 3.0 to prevent numbers exploding to 300+.
        if (fScale > 3.0) {
            fScale = 3.0; 
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

#include "BruSequentialProposal.h"
#include <RooStats/RooStatsUtils.h>
#include <RooRandom.h>
#include <cmath>

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

      for (int i = 0; i < n; ++i) {
        if (i == j || !force1D) {
          RooRealVar* var = _varCache[i];
          double val = var->getVal(), max = var->getMax(), min = var->getMin(), len = max - min;
          bool isCyclic = fCyclicPars.contains(*var);
          
          auto step = RooRandom::gaussian() * len * fScale;
          
          if (isCyclic) {
              // --- CYCLIC WRAP AROUND ---
              // Centers the math around the parameter's physical center (e.g. 0 for -pi to pi)
              double center = min + (len / 2.0);
              var->setVal(center + std::remainder(val + step - center, len));
          } else {
              // --- STANDARD BOUNCE ---
              while ((val + step > max) || (val + step < min)) {
                  step = RooRandom::gaussian() * len * fScale;
              }
              var->setVal(val + step);
          }
        }
      }
    }

    Bool_t BruSequentialProposal::CheckStepSize(Float_t acceptance){
      if(acceptance < fMinAcc || acceptance > fMaxAcc){
        Double_t oldScale = fScale; 
        
        Double_t acc = acceptance > 0 ? acceptance : 0.01;
        fScale *= (acc)/(fTargetAcc);
        
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
}// #include "BruSequentialProposal.h"
// #include <RooStats/RooStatsUtils.h>
// #include <RooRandom.h>

// namespace HS{
//   namespace FIT{

//     BruSequentialProposal::BruSequentialProposal(float scale,float target,float accmin,float accmax) :
//       ProposalFunction{},
//       fScale{scale},
//       fTargetAcc{target},
//       fMinAcc{accmin},
//       fMaxAcc{accmax}
//     {
//     }
 
//     void BruSequentialProposal::Propose(RooArgSet& xPrime, RooArgSet& x )
//     {
//       RooStats::SetParameters(&x, &xPrime);
//       int n = xPrime.getSize();
      
//       // 1. FAST C++ POINTER CACHE
//       // Build the cache only on the first pass (or if the parameter set changes size)
//       if (_varCache.empty() || _varCache.size() != (size_t)n) {
//           _varCache.clear();
//           for (auto *var : static_range_cast<RooRealVar *>(xPrime)) {
//               _varCache.push_back(var);
//           }
//       }

//       int j = int( floor(RooRandom::uniform()*n) );
      
//       // 2. THE VALLEY TRAP MIXTURE
//       // If we are doing N-D steps (_isNotSequential == true) BUT we are stuck in a valley,
//       // we flip a coin. 50% chance we temporarily force a 1D step to crawl along the valley floor.
//       bool force1D = (!_isNotSequential) || (_inValley && RooRandom::uniform() < 0.5);

//       // Fast iteration using the pre-built C++ vector
//       for (int i = 0; i < n; ++i) {
//         if (i == j || !force1D) {
//           RooRealVar* var = _varCache[i];
//           double val = var->getVal(), max = var->getMax(), min = var->getMin(), len = max - min;
          
//           auto step = RooRandom::gaussian() * len * fScale;
//           while ((val + step > max) || (val + step < min)) {
//               step = RooRandom::gaussian() * len * fScale;
//           }

//           var->setVal(val + step);
//         }
//       }
//     }

//     Bool_t BruSequentialProposal::CheckStepSize(Float_t acceptance){
//       if(acceptance < fMinAcc || acceptance > fMaxAcc){
//         Double_t oldScale = fScale; // <--- ADD THIS: Save the real previous scale
        
//         Double_t acc = acceptance > 0 ? acceptance : 0.01;
//         fScale *= (acc)/(fTargetAcc);
        
//         if(fScale < fMinScale){
//           fScale = fMinScale; 
//           fNminScale++;
          
//           // ALGORITHM UPGRADE: If we hit the floor, we are stuck in a valley!
//           if (_isNotSequential) {
//               _inValley = kTRUE;
//               std::cout << "  BruSequentialProposal: Valley Detected! Mixing in 1D steps to recover." << std::endl;
//           }
//         } else {
//             // If scale bounces back, we have escaped the valley
//             _inValley = kFALSE;
//         }
        
// 	std::cout<<"BruSequentialProposal::CheckStepSize Changed to "<<fScale<<" (min 1E-6) from "<< oldScale <<std::endl;

//         if(fNminScale > 10){
//           fNminScale=0;
//           std::cerr<<"BruSequentialProposal::CheckStepSize cannot get within allowed acceptance limits. Must be stuck. Will exit."<<std::endl;
//           return kFALSE;
//         }
//         return kTRUE;
//       }
      
//       // If acceptance is healthy, ensure valley mode is off
//       _inValley = kFALSE; 
//       fNminScale = 0; // Reset counter on healthy step
//       return kTRUE;
//     }
    
//     bool BruSequentialProposal::IsSymmetric(RooArgSet& , RooArgSet& ) {
//       return true;
//     }
 
//     double BruSequentialProposal::GetProposalDensity(RooArgSet& , RooArgSet& )
//     {
//       return 1.0; 
//     }

//   }
// }

#include "BruMetropolisHastings.h"
#include <RooStats/RooStatsUtils.h>
#include <RooStats/MarkovChain.h>
#include <RooStats/PdfProposal.h>
#include <RooStats/SequentialProposal.h>
#include <RooStats/ProposalHelper.h>
#include <TRandom.h>
#include <algorithm>
#include "BruSequentialProposal.h"

namespace HS{
  namespace FIT{
    
    RooStats::MarkovChain* BruMetropolisHastings::ConstructChain()
    {
      if (fParameters.getSize() == 0 || !fPropFunc || !fFunction) {
        coutE(Eval) << "Critical members unintialized: parameters, proposal function, or (log) likelihood function" << std::endl;
        return nullptr;
      }
  
      if (fChainParams.getSize() == 0) fChainParams.add(fParameters);
      
      RooArgSet x;
      RooArgSet xPrime;
      x.addClone(fParameters);
      xPrime.addClone(fParameters);
   
      auto* chain = new RooStats::MarkovChain();
      chain->SetParameters(fChainParams);

      // --- 1. FAST MEMORY MAPPING ---
      std::vector<RooRealVar*> curVars, propVars, masterVars;
      for (auto *var : static_range_cast<RooRealVar*>(x)) curVars.push_back(var);
      for (auto *var : static_range_cast<RooRealVar*>(xPrime)) propVars.push_back(var);
      for (auto *var : static_range_cast<RooRealVar*>(fParameters)) masterVars.push_back(var);
      size_t nVars = curVars.size();

      // --- 2. CHAIN BUFFER ---
      std::vector<std::vector<Double_t>> chainBufferVals;
      std::vector<Double_t> chainBufferNLL;
      std::vector<Double_t> chainBufferWeights;
      chainBufferVals.reserve(fNumIters);
      chainBufferNLL.reserve(fNumIters);
      chainBufferWeights.reserve(fNumIters);

      Int_t weight = 0;
      Double_t xL = 0.0, xPrimeL = 0.0, a = 0.0;

      RooFit::MsgLevel oldMsgLevel = RooMsgService::instance().globalKillBelow();
      RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);
      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
      RooAbsReal::clearEvalErrorLog();

      bool hadEvalError = true;
      ooccoutP((TObject *)nullptr, Generation) << "Metropolis-Hastings progress: ";

      int icount = 0;
      int totcount = 0;
      int snapcount = 0;
      int havePrinted = 0;

      // Allow down to 0.1% acceptance (1000 proposals per desired accepted step)
      // Use 64-bit integer so large Phase 2 chains never overflow 32-bit Int_t
      const Long64_t maxProposals = std::max<Long64_t>(
          static_cast<Long64_t>(fBailoutThreshold) * 2LL,
          static_cast<Long64_t>(fNumIters) * 1000LL);
      const Int_t maxConsecutiveRejects = 50000;
      Long64_t totalProposals = 0;
      Int_t consecutiveRejects = 0;
      fEvalErrorCount = 0;

      auto bseqprop = dynamic_cast<BruSequentialProposal*>(fPropFunc);
      if(bseqprop == nullptr){
        std::cerr<<"ERROR BruMetropolisHastings need a BruSequentialProposal "<<std::endl;
        delete chain;
        return nullptr;
      }
      bseqprop->ResetCounter();

      // --- INITIAL EVALUATION ---
      for (size_t i = 0; i < nVars; ++i) curVars[i]->setVal(masterVars[i]->getVal());
      xL = fFunction->getVal();

      std::cout<<"INITIAL STATUS "<<std::endl;
      for (size_t i = 0; i < nVars; ++i) {
        std::cout<< masterVars[i]->GetName()<<" "<<masterVars[i]->getVal()<<std::endl;
      }
      std::cout<< xL<<std::endl;

      // --- THE BARE-METAL MCMC LOOP ---
      while (icount < fNumIters) {

        // --- STOCHASTIC HOT-SWAP ---
        if (fBatchedNLLs.size() > 0 && fSwapFreq > 0 && totalProposals > 0 && (totalProposals % fSwapFreq == 0)) {
          fCurrentBatch = (fCurrentBatch + 1) % fBatchedNLLs.size();
          fFunction = fBatchedNLLs[fCurrentBatch]; 
          std::cout << "Shaking the landscape! Swapped to Batch " << fCurrentBatch << std::endl;
          
          for (size_t i = 0; i < nVars; ++i) {
            masterVars[i]->setVal(curVars[i]->getVal());
          }
          
          xL = fFunction->getVal(); 
        }

        totcount++; 
        totalProposals++;  
        hadEvalError = false;

        // ==========================================================
        // --- FAST-FAIL BAILOUT & HARD PROPOSAL CEILING ---
        // ==========================================================
        if (totalProposals % 1000 == 0) {
            double currentAcc = (double)icount / totalProposals;
            
            if (totalProposals >= fBailoutThreshold && currentAcc < fBailoutAcceptance) {
                std::cout << "\nBruMetropolisHastings: WARNING - Fast Bailout Triggered! "
                          << "Acceptance is critically low (" << currentAcc * 100 << "%) after " 
                          << totalProposals << " proposals. Aborting chain to allow immediate retuning." << std::endl;
                for (size_t i = 0; i < nVars; ++i) {
                  std::cout<< masterVars[i]->GetName()<<" "<<masterVars[i]->getVal()<<std::endl;
                }
                std::cout<< xL<<" " << xPrimeL<<std::endl;
                break;
            }

            if (totalProposals >= maxProposals) {
                std::cout << "\nBruMetropolisHastings: WARNING - Reached hard proposal ceiling ("
                          << totalProposals << " proposals, " << icount << "/" << fNumIters
                          << " accepted, Acc = " << currentAcc * 100.0 << "%). Terminating chain safely." << std::endl;
                break;
            }
        }

        // ==========================================================
        // Acceptance Tuning Logic
        if (totcount % 1000 == 0) {
          fAcceptance = ((Double_t)snapcount) / totcount;
          std::cout << "  BruMetropolisHastings accepted " << snapcount << " out of " << totcount << " for acceptance " << fAcceptance << std::endl;
          
          RooMsgService::instance().setGlobalKillBelow(oldMsgLevel);
          auto ok = bseqprop->CheckStepSize(fAcceptance);
          if (ok == kFALSE) {
            delete chain;
            return nullptr; 
          }
          snapcount = 0;
          totcount = 0;
          continue;
        }

        // Progress Printing
        if (icount % 100 == 0 && havePrinted == 0) {
          ooccoutP((TObject*)nullptr, Generation) << " " << icount << "/" << fNumIters;
          havePrinted = 1;
        }
        if (icount % 100 == 1) havePrinted = 0;

        // Propose new step
        fPropFunc->Propose(xPrime, x);
        
        // Fast push proposed values to the master NLL parameters
        for (size_t i = 0; i < nVars; ++i) masterVars[i]->setVal(propVars[i]->getVal());
        xPrimeL = fFunction->getVal();
        
        if (wasEvalErrors()) {
          xPrimeL = RooNumber::infinity();
          hadEvalError = true;
          fEvalErrorCount++;
        }

        a = xPrimeL - xL;
       
        if (!hadEvalError && !fPropFunc->IsSymmetric(xPrime, x)) {
          Double_t xPrimePD = fPropFunc->GetProposalDensity(xPrime, x);
          Double_t xPD      = fPropFunc->GetProposalDensity(x, xPrime);
          a += TMath::Log(xPrimePD) - TMath::Log(xPD);
        }
       
        if (!hadEvalError && ShouldTakeStep(a)) {
          consecutiveRejects = 0;
          // ACCEPTED: Record CURRENT state (x) with its accumulated weight
          if (weight != 0.0) {
            std::vector<Double_t> stepVals(nVars);
            for (size_t i = 0; i < nVars; ++i) stepVals[i] = curVars[i]->getVal();
              
            chainBufferVals.push_back(std::move(stepVals));
            chainBufferNLL.push_back(CalcNLL(xL));
            chainBufferWeights.push_back((Double_t)weight);
              
            icount++;
            snapcount++;
          }
          
          // Jump to new state (xPrime)
          weight = 1;
          for (size_t i = 0; i < nVars; ++i) curVars[i]->setVal(propVars[i]->getVal());
          xL = xPrimeL;
        } else {
          // REJECTED: Stay at current state, increment weight
          weight++;
          consecutiveRejects++;
          if (consecutiveRejects >= maxConsecutiveRejects) {
            std::cout << "\nBruMetropolisHastings: WARNING - Chain deadlocked for " << consecutiveRejects
                      << " consecutive proposals (" << icount << "/" << fNumIters
                      << " accepted). Terminating chain safely." << std::endl;
            break;
          }
        }
      }

      // Record the final point
      if (weight != 0.0) {
        std::vector<Double_t> stepVals(nVars);
        for (size_t i = 0; i < nVars; ++i) stepVals[i] = curVars[i]->getVal();
        chainBufferVals.push_back(std::move(stepVals));
        chainBufferNLL.push_back(CalcNLL(xL));
        chainBufferWeights.push_back((Double_t)weight);
      }
      
      ooccoutP((TObject *)nullptr, Generation) << std::endl;

      // --- 3. FLUSH BUFFER TO CHAIN ---
      for(size_t i = 0; i < chainBufferVals.size(); ++i) {
        for (size_t v = 0; v < nVars; ++v) masterVars[v]->setVal(chainBufferVals[i][v]);
        chain->Add(fParameters, chainBufferNLL[i], chainBufferWeights[i]);
      }

      RooMsgService::instance().setGlobalKillBelow(oldMsgLevel);

      Int_t numAccepted = chain->Size();
      fAcceptance = (totalProposals > 0) ? (((Double_t)icount) / totalProposals) : 0.0;
      
      coutI(Eval) << "Proposal acceptance rate: " << fAcceptance * 100.0 << "%" << std::endl;
      coutI(Eval) << "Number of steps in chain: " << numAccepted << std::endl;
 
      return chain;
    }

    Bool_t BruMetropolisHastings::wasEvalErrors(){
      if(RooAbsReal::numEvalErrors()){
        RooAbsReal::clearEvalErrorLog();
        return kTRUE;
      }
      return kFALSE;
    }

    Bool_t BruMetropolisHastings::CheckForBurnIn(RooStats::MarkovChain* chain){
      return true;
    } 
  }//namespace FIT
}//namespace HS

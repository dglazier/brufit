#include "AmpMcmc.h"
#include "AmpConfigure.h"
#include "BruMappedRhat.h"
#include <TMath.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include "AmpMcmc.h"
#include "AmpConfigure.h"
#include "BruMappedRhat.h"
#include <TMath.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooStats/RooStatsUtils.h>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace HS{
  namespace FIT{
    
    AmpMcmc::AmpMcmc(AmpConfigure* configure, std::vector<Double_t> Nuser, UInt_t nrefits,Int_t Nburn, Float_t norm,float target,float accmin,float accmax, Bool_t nozeroinit) 
      : BruMcmcCovariance(Nuser,Nburn,norm,target,accmin,accmax), fNFits{nrefits}, fNoZeroInitialVal{nozeroinit}, _ampHelper{configure} {
      SetNameTitle("HSAmpMcmc","Mcmc multi fit for amplitudes");
    }

    AmpMcmc::AmpMcmc(AmpConfigure* configure, std::vector<Int_t> Niters, UInt_t nrefits,Int_t Nburn, Float_t norm,float target,float accmin,float accmax, Bool_t nozeroinit) 
      : BruMcmcCovariance(Niters,Nburn,norm,target,accmin,accmax), fNFits{nrefits}, fNoZeroInitialVal{nozeroinit}, _ampHelper{configure} {
      SetNameTitle("HSAmpMcmc","Mcmc multi fit for amplitudes");
    }

    void AmpMcmc::Run(Setup &setup, RooAbsData &fitdata){
      TString baseName = GetName();
      std::cout << "AmpMcmc::Run " << baseName << " | Fits to attempt: " << fNFits << std::endl;
      fSetup = &setup;
      fData = &fitdata;

      _ampHelper.ConfigAmps(fSetup);
 
      std::vector<ChainSummary> summaries;
      summaries.reserve(fNFits);

      UInt_t nrefit = 0;
      UInt_t failureCount = 0;
      const UInt_t maxFailures = std::max<UInt_t>(3, fNFits * 2);

      while (nrefit < fNFits) {
        nrefit++;
        SetName(Form("%s_%u_", baseName.Data(), nrefit));
        std::cout << "\n=======================================================" << std::endl;
        std::cout << " AmpMcmc::Run - Starting Refit " << nrefit << " of " << fNFits << " (" << GetName() << ")" << std::endl;
        std::cout << "=======================================================" << std::endl;
        
        AmpMcmc::RandomiseParameters();
        
        // Reset the base class phase tracker before launching the multi-phase run
        ResetNIter();
        
        BruMcmcCovariance::Run(setup, fitdata);

        if (Success() == kFALSE || !fChainData || fChainData->numEntries() == 0) {
          failureCount++;
          std::cout << "AmpMcmc::Run Refit " << nrefit << " Failed (" << failureCount << "/" << maxFailures << ")." << std::endl;
          if (failureCount < maxFailures) {
            std::cout << " -> Retrying refit " << nrefit << "..." << std::endl;
            nrefit--;
            continue;
          } else {
            std::cerr << "AmpMcmc::Run WARNING: Maximum refit failures reached. Moving on." << std::endl;
            continue;
          }
        }

        // Clone fChainData BEFORE SaveInfo() deletes fTreeMCMC and closes fOutFile
        std::unique_ptr<RooDataSet> chainDataCopy(
            dynamic_cast<RooDataSet*>(fChainData->Clone(Form("mcmcChain_refit_%u", nrefit))));

        // SaveInfo() calls Result(), which computes the weighted posterior means and sigmas
        // and writes Results...HSAmpMcmc_i_.root to disk.
        SaveInfo();

        ChainSummary summary;
        summary.chainIdx    = nrefit;
        summary.chainName   = GetName();
        summary.meanNLL     = GetMeanNLL();
        summary.sigmaNLL    = GetSigmaNLL();
        summary.minNLL      = GetMinNLL();
        summary.maxRhat     = GetMaxRhat();
        summary.minESS      = GetMinESS();
        summary.meanESS     = GetMeanESS();
        summary.acceptance  = fChainAcceptance;
        summary.wasRescued  = WasRescued();
        summary.chainData   = std::move(chainDataCopy);
        summary.posteriorCov.ResizeTo(GetLastPosteriorCov().GetNrows(), GetLastPosteriorCov().GetNcols());
        summary.posteriorCov = GetLastPosteriorCov();
        summary.rhatVals    = GetLastRhatVals();
        summary.essVals     = GetLastEssVals();

        RooArgSet allParsAndYields(fSetup->ParsAndYields());
        summary.paramsSnapshot.reset( dynamic_cast<RooArgSet*>(allParsAndYields.snapshot()) );

        std::cout << "--> Refit " << nrefit << " Summary:"
                  << " | Mean NLL = " << Form("%.2f +/- %.2f (Min: %.2f)", summary.meanNLL, summary.sigmaNLL, summary.minNLL)
                  << " | Worst R-hat = " << Form("%.4f", summary.maxRhat)
                  << " | Worst ESS = " << Form("%.1f", summary.minESS)
                  << " | Mean ESS = " << Form("%.1f", summary.meanESS) << std::endl;

        summaries.push_back(std::move(summary));
      }

      if (summaries.empty()) {
        std::cerr << "AmpMcmc::Run FATAL: No valid chains completed!" << std::endl;
        SetName(baseName);
        return;
      }

      // =========================================================================
      // MULTI-CHAIN SELECTION:
      // 1. Find the lowest-likelihood basin (using Mean NLL to avoid single-sample
      //    fluctuations, preferring converged chains if available).
      // 2. Identify all chains statistically consistent with the lowest Mean NLL
      //    within natural MCMC thermal fluctuations.
      // 3. Among consistent chains, pick the chain with the best R-hat and ESS.
      // =========================================================================
      const Double_t nPars = static_cast<Double_t>(fSetup->NonConstParsAndYields().size());
      const Double_t thermalSigmaFloor = std::sqrt( std::max(1.0, nPars / 2.0) );

      int refIdx = -1;
      for (size_t i = 0; i < summaries.size(); ++i) {
        bool isConverged = (summaries[i].maxRhat != BruMappedRhat::kConvergenceFailure &&
                            summaries[i].maxRhat <= _rhatTarget &&
                            !summaries[i].wasRescued);
        if (isConverged) {
          if (refIdx < 0 || summaries[i].meanNLL < summaries[refIdx].meanNLL) {
            refIdx = static_cast<int>(i);
          }
        }
      }
      if (refIdx < 0) {
        for (size_t i = 0; i < summaries.size(); ++i) {
          if (summaries[i].maxRhat != BruMappedRhat::kConvergenceFailure) {
            if (refIdx < 0 || summaries[i].meanNLL < summaries[refIdx].meanNLL) {
              refIdx = static_cast<int>(i);
            }
          }
        }
      }
      if (refIdx < 0) refIdx = 0;

      const Double_t refMeanNLL  = summaries[refIdx].meanNLL;
      const Double_t refSigmaNLL = std::max(summaries[refIdx].sigmaNLL, thermalSigmaFloor);

      int bestIdx = refIdx;
      bool bestIsConverged = false;
      Double_t bestDiagScore = 1e18;

      Double_t targetESS = (_fNumItersVec.size() >= 3 ? _fNumItersVec[2] : fNumIters) * _essFraction;
      if (targetESS < 50.0) targetESS = 50.0;

      std::vector<bool> isConsistent(summaries.size(), false);
      std::vector<Double_t> diagScores(summaries.size(), 1e18);

      for (size_t i = 0; i < summaries.size(); ++i) {
        const auto& s = summaries[i];
        Double_t chainSigma = std::max(s.sigmaNLL, thermalSigmaFloor);
        Double_t nllTolerance = std::max(2.0 * std::sqrt(refSigmaNLL * refSigmaNLL + chainSigma * chainSigma),
                                         std::sqrt(2.0 * std::max(1.0, nPars)));

        isConsistent[i] = (s.meanNLL <= refMeanNLL + nllTolerance) &&
                          (s.maxRhat != BruMappedRhat::kConvergenceFailure);

        Double_t rhatTerm    = std::max(0.0, s.maxRhat - 1.0) / std::max(1e-3, _rhatTarget - 1.0);
        Double_t minEssTerm  = targetESS / std::max(1.0, s.minESS);
        Double_t meanEssTerm = 0.25 * (targetESS / std::max(1.0, s.meanESS));
        Double_t rescuePen   = s.wasRescued ? 10.0 : 0.0;
        diagScores[i]        = rhatTerm + minEssTerm + meanEssTerm + rescuePen;

        if (!isConsistent[i]) continue;

        bool sConverged = (s.maxRhat <= _rhatTarget && !s.wasRescued);
        if (sConverged && !bestIsConverged) {
          bestIsConverged = true;
          bestDiagScore   = diagScores[i];
          bestIdx         = static_cast<int>(i);
        } else if (sConverged == bestIsConverged) {
          if (diagScores[i] < bestDiagScore) {
            bestDiagScore = diagScores[i];
            bestIdx       = static_cast<int>(i);
          }
        }
      }

      std::cout << "\n=========================================================================================================" << std::endl;
      std::cout << " AmpMcmc Multi-Chain Selection Audit (Lowest Mean NLL Ref = " << Form("%.2f +/- %.2f", refMeanNLL, refSigmaNLL) << ")" << std::endl;
      std::cout << "---------------------------------------------------------------------------------------------------------" << std::endl;
      for (size_t i = 0; i < summaries.size(); ++i) {
        const auto& s = summaries[i];
        std::cout << Form(" Chain %u (%s) | MeanNLL: %9.2f +/- %5.2f (Min: %9.2f) | R-hat: %7.4f | MinESS: %6.1f | MeanESS: %6.1f | NLL-OK: %s %s",
                          s.chainIdx, s.chainName.Data(), s.meanNLL, s.sigmaNLL, s.minNLL,
                          s.maxRhat, s.minESS, s.meanESS,
                          (isConsistent[i] ? "YES" : " NO"),
                          (static_cast<int>(i) == bestIdx ? "<-- [SELECTED BEST]" : "")) << std::endl;
      }
      std::cout << "=========================================================================================================\n" << std::endl;

      auto& winner = summaries[bestIdx];
      SetName(baseName);
      SetTag("");
      SetupBasicUsage();

      if (winner.paramsSnapshot) {
        for (auto* arg : fSetup->ParsAndYields()) {
          auto* masterVar = dynamic_cast<RooRealVar*>(arg);
          auto* snapVar   = dynamic_cast<RooRealVar*>(winner.paramsSnapshot->find(arg->GetName()));
          if (masterVar && snapVar) {
            masterVar->setVal(snapVar->getVal());
            masterVar->setError(snapVar->getError());
          }
        }
      }

      fChainAcceptance = winner.acceptance;
      fMeanNLL         = winner.meanNLL;
      fSigmaNLL        = winner.sigmaNLL;
      fMinNLL          = winner.minNLL;
      fMaxRhat         = winner.maxRhat;
      fMinESS          = winner.minESS;
      fMeanESS         = winner.meanESS;
      fWasRescued      = winner.wasRescued;

      fChainData = std::move(winner.chainData);
      if (fTreeMCMC) { delete fTreeMCMC; fTreeMCMC = nullptr; }

      if (fChainData) {
        auto saveDir = gDirectory;
        if (fOutFile) fOutFile->cd();
        fTreeMCMC = RooStats::GetAsTTree("MCMCTree", "MCMCTree", *fChainData);
        WriteDiagnosticsAndCovariance(winner.posteriorCov, winner.rhatVals, winner.essVals);
        if (saveDir) saveDir->cd();
      }
    }

    void AmpMcmc::RandomiseParameters(){
      std::cout<<"AmpMcmc::RandomiseParameters() "<<fSetup<<std::endl;
      Double_t intensity0=0;
 
      auto Nrand=10000;
      for(Int_t irand=0;irand<Nrand;++irand){

        _ampHelper.RandomiseFitParameters();
        
        if(fNoZeroInitialVal==kFALSE)
          break;
        
        intensity0 = fSetup->Model()->getVal();
        std::cout<<"AmpMcmc::RandomiseParameters() rand "<<irand<<" "<<intensity0<<std::endl;
        if (intensity0!=0) break;
        if(irand==(Nrand-1)){
          std::cerr<<"AmpMcmc::Run FATAL : error tried to randomise parameters 9999 times without a non-zero intensity and you specified not have a non zero inital intensity. You may try removing the true from the Minuit constructor if you do not need this. Note on irand = "<<irand<<std::endl;
          exit(0);
        }
      }
    }
  }
}

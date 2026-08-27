#include "AmpMcmcFixedCovariance.h"
#include "AmpConfigure.h"
#include <TMath.h>
#include <RooArgSet.h>
#include <iostream>

namespace HS{
  namespace FIT{
    
    // Mode 1: Pass static file config to the fixed covariance base class
    AmpMcmcFixedCovariance::AmpMcmcFixedCovariance(AmpConfigure* configure, const TString& covFilePath, const TString& matrixName, std::vector<Int_t> Niters, UInt_t nrefits, Int_t Nburn, Float_t norm, float target, float accmin, float accmax, Bool_t nozeroinit) 
      : BruMcmcFixedCovariance(covFilePath, matrixName, Niters, Nburn, norm, target, accmin, accmax), fNFits{nrefits}, fNoZeroInitialVal{nozeroinit}, _ampHelper{configure} {
      SetNameTitle("HSAmpMcmcFixedCovariance","Mcmc multi fit for amplitudes using Fixed Covariance");
    }

    // Mode 2 & 3: Pass dynamic config to the fixed covariance base class
    AmpMcmcFixedCovariance::AmpMcmcFixedCovariance(AmpConfigure* configure, CovLoadMode mode, const TString& pathStr1, const TString& pathStr2, const TString& matrixName, std::vector<Int_t> Niters, UInt_t nrefits, Int_t Nburn, Float_t norm, float target, float accmin, float accmax, Bool_t nozeroinit) 
      : BruMcmcFixedCovariance(mode, pathStr1, pathStr2, matrixName, Niters, Nburn, norm, target, accmin, accmax), fNFits{nrefits}, fNoZeroInitialVal{nozeroinit}, _ampHelper{configure} {
      SetNameTitle("HSAmpMcmcFixedCovariance","Mcmc multi fit for amplitudes using Fixed Covariance");
    }

    void AmpMcmcFixedCovariance::Run(Setup &setup, RooAbsData &fitdata){
      std::cout << "AmpMcmcFixedCovariance::Run " << GetName() << " | Fits to attempt: " << fNFits << std::endl;
      fSetup = &setup;
      fData = &fitdata;

      _ampHelper.ConfigAmps(fSetup);
 
      UInt_t nrefit = 0;
      Double_t bestNLL = 1e12; // Track the lowest NLL
      std::unique_ptr<RooArgSet> bestParams; // Snapshot tracker for best variables

      while(nrefit++ < fNFits){
        SetName(Form("HSAmpMcmcFixedCovariance_%d_", nrefit));
        std::cout << "\n=======================================================" << std::endl;
        std::cout << " AmpMcmcFixedCovariance::Run - Starting Refit " << nrefit << " of " << fNFits << " (" << GetName() << ")" << std::endl;
        std::cout << "=======================================================" << std::endl;
        
        RandomiseParameters();
        
        // CRITICAL FIX: Reset the base class phase tracker before launching the burn-in
        ResetNIter();
        
        // CRITICAL FIX: Call the Fixed Covariance base class Run method so matrices are loaded!
        BruMcmcFixedCovariance::Run(setup, fitdata);

        if(Success() == kFALSE){
          // failed so try again
          std::cout << "AmpMcmcFixedCovariance::Run Refit " << nrefit << " Failed. Retrying..." << std::endl;
          nrefit--;
          continue;
        }
        
        // --- Track the Best Minimum ---
        Double_t currentNLL = NLL();
        std::cout << "--> Refit " << nrefit << " completed with NLL = " << currentNLL << std::endl;

        if (currentNLL < bestNLL) {
            bestNLL = currentNLL;
            std::cout << "    (*** New Best Minimum Found! ***)" << std::endl;
            
            // Snapshot the winning parameters so we can restore them later
            if (bestParams) {
                bestParams->assignValueOnly(fSetup->Parameters());
            } else {
                bestParams.reset((RooArgSet*)fSetup->Parameters().snapshot());
            }
        }

        // Save intermediate chain trees
        if(nrefit < fNFits) SaveInfo();
      }

      // --- Restore the Best Minimum ---
      if (bestParams) {
          std::cout << "\n=======================================================" << std::endl;
          std::cout << " AmpMcmcFixedCovariance: All fits complete. Restoring Best Minimum (NLL = " << bestNLL << ")" << std::endl;
          std::cout << "=======================================================\n" << std::endl;
          
          fSetup->Parameters().assignValueOnly(*bestParams);
      }
    }

    void AmpMcmcFixedCovariance::RandomiseParameters(){
      std::cout << "AmpMcmcFixedCovariance::RandomiseParameters() " << fSetup << std::endl;
      Double_t intensity0=0;
 
      //Try 10000 times to get a physical starting value
      auto Nrand=10000;
      for(Int_t irand=0;irand<Nrand;++irand){

        _ampHelper.RandomiseFitParameters();
        
        if(fNoZeroInitialVal==kFALSE) //if don't care about 0 iniatal intensity break
          break;
        
        // if we do care keep trying until we get non-zero intensity
        intensity0 = fSetup->Model()->getVal();
        std::cout << "AmpMcmcFixedCovariance::RandomiseParameters() rand " << irand << " " << intensity0 << std::endl;
        
        if (intensity0!=0) break;
        if(irand==(Nrand-1)){
          std::cerr << "AmpMcmcFixedCovariance::Run FATAL : error tried to randomise parameters 9999 times without a non-zero intensity and you specified not have a non zero inital intensity. You may try removing the true from the Minuit constructor if you do not need this. Note on irand = " << irand << std::endl;
          exit(0);
        }
      }
    }
  } //namespace FIT
} //namespace HS

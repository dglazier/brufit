#include "AmpMcmc.h"
#include "AmpConfigure.h"
#include <TMath.h>
#include <RooArgSet.h>

namespace HS{
  namespace FIT{
    
    // Pass the vector of iterations to the base class, and hardcode norm to 0.01 since the tuner overwrites it anyway
    AmpMcmc::AmpMcmc(AmpConfigure* configure, std::vector<Int_t> Niters, UInt_t nrefits,Int_t Nburn, Float_t norm,float target,float accmin,float accmax, Bool_t nozeroinit) 
      : BruMcmcCovariance(Niters,Nburn,norm,target,accmin,accmax), fNFits{nrefits}, fNoZeroInitialVal{nozeroinit}, _ampHelper{configure} {
      SetNameTitle("HSAmpMcmc","Mcmc multi fit for amplitudes");
    }

    void AmpMcmc::Run(Setup &setup, RooAbsData &fitdata){
      cout << "AmpMcmc::Run " << GetName() << " | Fits to attempt: " << fNFits << endl;
      fSetup = &setup;
      fData = &fitdata;

      _ampHelper.ConfigAmps(fSetup);
 
      UInt_t nrefit = 0;
      Double_t bestNLL = 1e12; // Track the lowest NLL
      std::unique_ptr<RooArgSet> bestParams; // Snapshot tracker for best variables

      while(nrefit++ < fNFits){
        SetName(Form("HSAmpMcmc_%d_", nrefit));
        cout << "\n=======================================================" << endl;
        cout << " AmpMcmc::Run - Starting Refit " << nrefit << " of " << fNFits << " (" << GetName() << ")" << endl;
        cout << "=======================================================" << endl;
        
        AmpMcmc::RandomiseParameters();
        
        // CRITICAL FIX: Reset the base class phase tracker before launching the 3-phase burn in!
        ResetNIter();
        
        BruMcmcCovariance::Run(setup, fitdata);

        if(Success() == kFALSE){
          // failed so try again
          cout << "AmpMcmc::Run Refit " << nrefit << " Failed. Retrying..." << endl;
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
                // FIX: Use assignValueOnly to safely copy values from the list to the set
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
          std::cout << " AmpMcmc: All fits complete. Restoring Best Minimum (NLL = " << bestNLL << ")" << std::endl;
          std::cout << "=======================================================\n" << std::endl;
          
          // FIX: Safely copy the winning values back into the master physics parameters
          fSetup->Parameters().assignValueOnly(*bestParams);
      }
    }

    void AmpMcmc::RandomiseParameters(){
      cout<<"AmpMcmc::RandomiseParameters() "<<fSetup<<endl;//exit(0);
      Double_t intensity0=0;
 
      //Try 10000 times to get a physical starting value
      auto Nrand=10000;
      for(Int_t irand=0;irand<Nrand;++irand){

        _ampHelper.RandomiseFitParameters();
        
        if(fNoZeroInitialVal==kFALSE) //if don't care about 0 iniatal intensity break
          break;
        
        // if we do care keep trying until we get non-zero intensity
        intensity0 = fSetup->Model()->getVal();
        std::cout<<"AmpMcmc::RandomiseParameters() rand "<<irand<<" "<<intensity0<<endl;
        if (intensity0!=0) break;
        if(irand==(Nrand-1)){
          std::cerr<<"AmpMcmc::Run FATAL : error tried to randomise parameters 9999 times without a non-zero intensity and you specified not have a non zero inital intensity. You may try removing the true from the Minuit constructor if you do not need this. Note on irand = "<<irand<<std::endl;
          exit(0);
        }
      }
    }
  }
}

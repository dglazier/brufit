////////////////////////////////////////////////////////////////
///
///Class:               AmpMcmc
///Description:
///           

#include "BruMcmc.h"
#include "PhotoTwoSpin0Amps.h"
#include "AmpHelpers.h"


#pragma once
namespace HS{
  namespace FIT{


    class AmpMcmc  : public BruMcmcCovariance {
      
    public:

      AmpMcmc(AmpConfigure* configure,Int_t Niter=100,Int_t Nburn=10, Float_t norm=0.1,UInt_t nrefits=0,Bool_t nozeroinit=kFALSE);
      // AmpMcmc(const AmpMcmc&)=default;
      //AmpMcmc(AmpMcmc&&)=default;
      ~AmpMcmc() override =default;
      //AmpMcmc& operator=(const AmpMcmc& other)=default;
      // AmpMcmc& operator=(AmpMcmc&& other) = default;  

      void RandomiseParameters();
      void Run(Setup &setup,RooAbsData &fitdata) override;

      //void ConfigAmps(const PhotoTwoSpin0Amps& config);
      //void ConfigAmps(AmpConfigure* config);
      void  CopyToMomentPars();
      void  CopyToAmpPars();
      
    private:
      UInt_t fIFit=0;
      UInt_t fNFits=1;
      Bool_t fNoZeroInitialVal=kFALSE;
      Bool_t _IsAmplitudes=kTRUE;

      // Setup _ampSetup;
      //PhotoTwoSpin0Amps* _ampConfig={nullptr};
      AmpHelpers _ampHelper;
      
      ClassDefOverride(HS::FIT::AmpMcmc,1);
      
    };

  }//namespaces
}

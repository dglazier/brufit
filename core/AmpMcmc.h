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

      // Updated: Replaced single Int_t with vector, removed Float_t norm
      AmpMcmc(AmpConfigure* configure, std::vector<Int_t> Niters, UInt_t nrefits,Int_t Nburn=10, Float_t norm=0.01,float target=0.234,float accmin=0.16,float accmax=0.3, Bool_t nozeroinit=false);
      ~AmpMcmc() override = default;

      void RandomiseParameters();
      void Run(Setup &setup, RooAbsData &fitdata) override;

      void  CopyToMomentPars();
      void  CopyToAmpPars();
      
    private:
      UInt_t fIFit=0;
      UInt_t fNFits=1;
      Bool_t fNoZeroInitialVal=kFALSE;
      Bool_t _IsAmplitudes=kTRUE;

      AmpHelpers _ampHelper;
      
      ClassDefOverride(HS::FIT::AmpMcmc,1);
      
    };

  }//namespaces
}

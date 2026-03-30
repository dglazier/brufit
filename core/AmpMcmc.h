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
      AmpMcmc(AmpConfigure* configure, std::vector<Int_t> Niters, Int_t Nburn=10, UInt_t nrefits=0, Bool_t nozeroinit=kFALSE);
      
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

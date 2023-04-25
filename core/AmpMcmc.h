////////////////////////////////////////////////////////////////
///
///Class:               AmpMcmc
///Description:
///           

#include "RooMcmc.h"


#pragma once
namespace HS{
  namespace FIT{


    class AmpMcmc  : public RooMcmcSeqHelper {
      
    public:

      AmpMcmc(Int_t Niter=100,Int_t Nburn=10, Float_t norm=0.1,UInt_t nrefits=0,Bool_t nozeroinit=kFALSE);
      // AmpMcmc(const AmpMcmc&)=default;
      //AmpMcmc(AmpMcmc&&)=default;
      ~AmpMcmc() override =default;
      //AmpMcmc& operator=(const AmpMcmc& other)=default;
      // AmpMcmc& operator=(AmpMcmc&& other) = default;  

      void RandomiseParameters();
      void Run(Setup &setup,RooAbsData &fitdata) override;

      UInt_t fIFit=0;
      UInt_t fNFits=1;
      Bool_t fNoZeroInitialVal=kFALSE;

      ClassDefOverride(HS::FIT::AmpMcmc,1);
      
    };

  }//namespaces
}

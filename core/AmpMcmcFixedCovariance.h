////////////////////////////////////////////////////////////////
///
/// Class:       AmpMcmcFixedCovariance
/// Description: Multi-chain MCMC fitter for amplitudes utilizing
///              injected Fixed Covariance matrices.
///           

#pragma once

#include "BruFixedCovariance.h"
#include "PhotoTwoSpin0Amps.h"
#include "AmpHelpers.h"

namespace HS{
  namespace FIT{

    class AmpMcmcFixedCovariance  : public BruMcmcFixedCovariance {
      
    public:

      // Mode 1 Constructor: Static file for all bins
      AmpMcmcFixedCovariance(AmpConfigure* configure, const TString& covFilePath, const TString& matrixName, std::vector<Int_t> Niters, UInt_t nrefits, Int_t Nburn=10, Float_t norm=0.01, float target=0.234, float accmin=0.16, float accmax=0.3, Bool_t nozeroinit=false);
      
      // Mode 2 & 3 Constructor: Dynamic loading per bin
      AmpMcmcFixedCovariance(AmpConfigure* configure, CovLoadMode mode, const TString& pathStr1, const TString& pathStr2, const TString& matrixName, std::vector<Int_t> Niters, UInt_t nrefits, Int_t Nburn=10, Float_t norm=0.01, float target=0.234, float accmin=0.16, float accmax=0.3, Bool_t nozeroinit=false);
      
      ~AmpMcmcFixedCovariance() override = default;

      void RandomiseParameters();
      void Run(Setup &setup, RooAbsData &fitdata) override;
      
    private:
      UInt_t fIFit=0;
      UInt_t fNFits=1;
      Bool_t fNoZeroInitialVal=kFALSE;
      Bool_t _IsAmplitudes=kTRUE;

      AmpHelpers _ampHelper;
      
      ClassDefOverride(HS::FIT::AmpMcmcFixedCovariance, 1);
      
    };

  }//namespace FIT
}//namespace HS

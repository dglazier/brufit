////////////////////////////////////////////////////////////////
///
///Class:               AmpMcmc
///Description:
///           

#pragma once

#include "BruMcmc.h"
#include "PhotoTwoSpin0Amps.h"
#include "AmpHelpers.h"
#include <TMatrixDSym.h>
#include <memory>
#include <vector>

namespace HS{
  namespace FIT{

    class AmpMcmc  : public BruMcmcCovariance {
      
    public:

      AmpMcmc(AmpConfigure* configure, std::vector<Double_t> Nuser, UInt_t nrefits,Int_t Nburn=10, Float_t norm=0.01,float target=0.234,float accmin=0.16,float accmax=0.3, Bool_t nozeroinit=false);
      AmpMcmc(AmpConfigure* configure, std::vector<Int_t> Niters, UInt_t nrefits,Int_t Nburn=10, Float_t norm=0.01,float target=0.234,float accmin=0.16,float accmax=0.3, Bool_t nozeroinit=false);
      ~AmpMcmc() override = default;

      void RandomiseParameters();
      void Run(Setup &setup, RooAbsData &fitdata) override;

      void CopyToMomentPars();
      void CopyToAmpPars();
      
    private:
      struct ChainSummary {
        UInt_t chainIdx = 0;
        TString chainName;
        Double_t meanNLL = 1e12;
        Double_t sigmaNLL = 0.0;
        Double_t minNLL = 1e12;
        Double_t maxRhat = 9999.0;
        Double_t minESS = 0.0;
        Double_t meanESS = 0.0;
        Double_t acceptance = 0.0;
        Bool_t wasRescued = kFALSE;
        std::unique_ptr<RooArgSet> paramsSnapshot;
        std::unique_ptr<RooDataSet> chainData;
        TMatrixDSym posteriorCov;
        std::vector<Double_t> rhatVals;
        std::vector<Double_t> essVals;
      };

      UInt_t fIFit=0;
      UInt_t fNFits=1;
      Bool_t fNoZeroInitialVal=kFALSE;
      Bool_t _IsAmplitudes=kTRUE;

      AmpHelpers _ampHelper;
      
      ClassDefOverride(HS::FIT::AmpMcmc,1);
      
    };

  }//namespaces
}

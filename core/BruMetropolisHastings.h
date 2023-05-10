////////////////////////////////////////////////////////////////
///
///Class:               BruMetropolisHastings
///Description:
///           

#pragma once

#include <RooStats/MetropolisHastings.h>
#include <vector>

namespace HS{
  namespace FIT{

    
    class BruMetropolisHastings : public RooStats::MetropolisHastings {
      
    public:
      BruMetropolisHastings(){
	SetType(kLog);
	SetSign(kNegative);
      }
      BruMetropolisHastings(const BruMetropolisHastings&)=default;
      BruMetropolisHastings(BruMetropolisHastings&&)=default;
      ~BruMetropolisHastings() override =default;
      BruMetropolisHastings& operator=(const BruMetropolisHastings& other)=default;
      BruMetropolisHastings& operator=(BruMetropolisHastings&& other) = default;

      RooStats::MarkovChain* ConstructChain() override;

      Double_t GetAcceptance()const {return fAcceptance;}

      void Help(){fTryHelp=kTRUE;}
      Bool_t CheckForBurnIn(RooStats::MarkovChain* chain);

      void SetAcceptanceRange(Double_t min,Double_t max){
	fMinAcc=min;
	fMaxAcc=max;
      }

      void SetBalance(RooAbsPdf* bal){fBalancePDF=bal;}

      void SetNorm(Float_t norm){fNorm=norm;}
      void SetTargetAccept(Float_t target){fTargetAcc=target;};
    protected:
      Bool_t wasEvalErrors();
    private:
      RooAbsPdf* fBalancePDF=nullptr;
      
      // Bool_t fRandomiseStart=kTRUE;
      Bool_t fTryHelp=kFALSE;
      Int_t fNWorse=0;
      Int_t fNBetter=0;
      Int_t fNWorseThanSave=0;
      Int_t fNBetterThanSave=0;
      Double_t fSaveNLL=0;
      Double_t fAcceptance=0;
      Double_t fLastEntries=0;
      Double_t fMinAcc=0.15;
      Double_t fMaxAcc=0.3;
      Double_t fOldBalance=1.;
      Float_t fNorm=0;
      Double_t fTargetAcc=0.234;
 
      std::vector<Double_t> fMeans;
      std::vector<Double_t> fSigmas;
      ClassDefOverride(HS::FIT::BruMetropolisHastings,1);
     };
    
  }//namespace FIT
}//namespace HS


////////////////////////////////////////////////////////////////
///
///Class:               HSMetropolisHastings
///Description:
///           

#pragma once

#include <RooStats/MetropolisHastings.h>
#include <vector>

namespace HS{
  namespace FIT{

    
    class HSMetropolisHastings : public RooStats::MetropolisHastings {
      
    public:
      HSMetropolisHastings()=default;
      HSMetropolisHastings(const HSMetropolisHastings&)=default;
      HSMetropolisHastings(HSMetropolisHastings&&)=default;
      ~HSMetropolisHastings() override =default;
      HSMetropolisHastings& operator=(const HSMetropolisHastings& other)=default;
      HSMetropolisHastings& operator=(HSMetropolisHastings&& other) = default;

      void SetKeepStart(){fRandomiseStart=kFALSE;}

      RooStats::MarkovChain* ConstructChain() override;

      Double_t GetAcceptance()const {return fAcceptance;}

      void Help(){fTryHelp=kTRUE;}
      Bool_t CheckForBurnIn(RooStats::MarkovChain* chain);
      
    protected:
      Bool_t wasEvalErrors();
    private:
      Bool_t fRandomiseStart=kTRUE;
      Bool_t fTryHelp=kFALSE;
      Int_t fNWorse=0;
      Int_t fNBetter=0;
      Int_t fNWorseThanSave=0;
      Int_t fNBetterThanSave=0;
      Double_t fSaveNLL=0;
      Double_t fAcceptance=0;
      Double_t fLastEntries=0;
      std::vector<Double_t> fMeans;
      std::vector<Double_t> fSigmas;
      ClassDefOverride(HS::FIT::HSMetropolisHastings,1);
     };
    
  }//namespace FIT
}//namespace HS


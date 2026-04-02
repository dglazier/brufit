////////////////////////////////////////////////////////////////
///
///Class:           BruSequentialProposal    
///Description:  Just sequential proposal with adaptove step size to
///              keep in some requested range
///           

#pragma once

#include <RooStats/ProposalFunction.h>
#include <RooRealVar.h>
#include <RooArgList.h>
#include <vector>

namespace HS{
  namespace FIT{

    class BruSequentialProposal : public RooStats::ProposalFunction{

    public:
      BruSequentialProposal() : RooStats::ProposalFunction() {}
      
      BruSequentialProposal(float scale,float target=0.234,float accmin=0.15,float accmax=0.35) ;
      
      /// Populate xPrime with a new proposed point
      void Propose(RooArgSet& xPrime, RooArgSet& x) override;
 
      /// Determine whether or not the proposal density is symmetric for
      /// points x1 and x2 - that is, whether the probability of reaching x2
      /// from x1 is equal to the probability of reaching x1 from x2
      bool IsSymmetric(RooArgSet& x1, RooArgSet& x2) override ;
 
      /// Return the probability of proposing the point x1 given the starting
      /// point x2
      double GetProposalDensity(RooArgSet& x1, RooArgSet& x2) override;
 
      ~BruSequentialProposal() override =default;
 
      void SetAcceptanceRange(Double_t min,Double_t max){
        fMinAcc=min;
        fMaxAcc=max;
      }
      void SetTargetAccept(Float_t target){fTargetAcc=target;};
      void SetMinScale(Float_t ms){fMinScale=ms;};
      void SetScale(Float_t scale) { fScale = scale; }
      virtual Bool_t CheckStepSize(Float_t acceptance);

      Float_t StepSizeFactor() const {return fScale;}
      void SetIsSequential(Bool_t isit){_isNotSequential= (!isit);}

      virtual void ResetCounter(){
        fNminScale = 0;
        _inValley = kFALSE;
        _varCache.clear(); 
      }    

      // --- CYCLIC PARAMETERS ---
      void SetCyclicParameters(const RooArgList& cyclics) { 
          fCyclicPars.removeAll(); 
          fCyclicPars.add(cyclics);
	  std::cout<< "BruSequentialProposal  SetCyclicParameters :"<<std::endl;cyclics.Print();
      }
    protected:
     RooArgList fCyclicPars; // List of parameters that should wrap around boundaries

    private:
 
      Float_t fScale=1.;
      Float_t fMinScale=1E-6;
      Float_t fMinAcc=0.15;
      Float_t fMaxAcc=0.3;
      Float_t fTargetAcc=0.234;
      UShort_t fNminScale = 0;
      Bool_t _isNotSequential=kFALSE;
      
      Bool_t _inValley = kFALSE; 
      std::vector<RooRealVar*> _varCache; 
      
 
      ClassDefOverride(BruSequentialProposal,1) 

    };

  }
}

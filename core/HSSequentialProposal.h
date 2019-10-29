#pragma once

#include <RooStats/SequentialProposal.h>
#include <RooArgSet.h>

namespace HS{
  namespace FIT{
    
    
    class HSSequentialProposal	 : public RooStats::SequentialProposal {      
   
    public:	

    HSSequentialProposal(double divisor):fDivisor(1./divisor){}

      void   SetDivisor(double divisor){fDivisor = 1./divisor;}
      /// Populate xPrime with a new proposed point
      void Propose(RooArgSet& xPrime, RooArgSet& x) override;
      
      /// Determine whether or not the proposal density is symmetric for
      /// points x1 and x2 - that is, whether the probabilty of reaching x2
      /// from x1 is equal to the probability of reaching x1 from x2
      Bool_t IsSymmetric(RooArgSet& x1, RooArgSet& x2) override;
      
      /// Return the probability of proposing the point x1 given the starting
      /// point x2
      Double_t GetProposalDensity(RooArgSet& x1, RooArgSet& x2) override;
      
      
    private:
      
      Double_t fDivisor=1;
      
      ClassDefOverride(HS::FIT::HSSequentialProposal,1); // A concrete implementation of ProposalFunction, that uniformly samples the parameter space.
    };
  }
}

#include "BruSequentialProposal.h"

namespace HS{
  namespace FIT{

    ////////////////////////////////////////////////////////////////////////////////
 
    BruSequentialProposal::BruSequentialProposal(double divisor) :
      ProposalFunction(),
      fDivisor(1./divisor)
    {
    }
 
    ////////////////////////////////////////////////////////////////////////////////
    /// Populate xPrime with a new proposed point
 
    void BruSequentialProposal::Propose(RooArgSet& xPrime, RooArgSet& x )
    {
      RooStats::SetParameters(&x, &xPrime);
      int n = xPrime.getSize();
      int j = int( floor(RooRandom::uniform()*n) );
      int i = 0;
      for (auto *var : static_range_cast<RooRealVar *>(xPrime)) {
	if (i == j) {
	  double val = var->getVal(), max = var->getMax(), min = var->getMin(), len = max - min;
	  val += RooRandom::gaussian() * len * fDivisor;
	  while (val > max) val -= len;
	  while (val < min) val += len;
	  var->setVal(val);
	  //std::cout << "Proposing a step along " << var->GetName() << std::endl;
	}
	++i;
      }
    }
 
    ////////////////////////////////////////////////////////////////////////////////
    /// Return the probability of proposing the point x1 given the starting
    /// point x2
 
    bool BruSequentialProposal::IsSymmetric(RooArgSet& , RooArgSet& ) {
      return true;
    }
 
    double BruSequentialProposal::GetProposalDensity(RooArgSet& ,
						     RooArgSet& )
    {
      return 1.0; // should not be needed
    }
 

  }
}

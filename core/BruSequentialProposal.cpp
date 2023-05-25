#include "BruSequentialProposal.h"
#include <RooStats/RooStatsUtils.h>
#include <RooRandom.h>

namespace HS{
  namespace FIT{

    ////////////////////////////////////////////////////////////////////////////////
 
    BruSequentialProposal::BruSequentialProposal(float scale,float target,float accmin,float accmax) :
      ProposalFunction{},
      fScale{scale},
      fTargetAcc{target},
      fMinAcc{accmin},
      fMaxAcc{accmax}
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
      // std::cout<<"BruSequentialProposal::Propose  static_range_cast< does not work until 6.28 c++14"<<std::endl;
      //exit(0);
       for (auto *var : static_range_cast<RooRealVar *>(xPrime)) {
      	if (i == j || _isNotSequential) {
      	  // std::cout<<"Propose "<<i<<" "<<j<<std::endl;
      	  double val = var->getVal(), max = var->getMax(), min = var->getMin(), len = max - min;
      	  auto step = RooRandom::gaussian() * len * fScale;//scale=#number of sigma
      	  while ((val+step > max) || (val+step <min) ) step = RooRandom::gaussian() * len * fScale;

      	  //val += RooRandom::gaussian() * len * fScale;
      	  //	  while (val > max)  val -= len;//should try val=max-val
      	  //while (val < min) val += len;

      	  var->setVal(val+step);
      	  //std::cout << "Proposing a step along " << var->GetName() << std::endl;
      	}
      	++i;
      }
    }
    Bool_t  BruSequentialProposal::CheckStepSize(Float_t acceptance){
      if(acceptance<fMinAcc||acceptance>fMaxAcc){
	//in case no event accepted start with correcting for 1%
	Double_t acc = acceptance >0 ? acceptance:0.01;
	fScale *= (acc)/(fTargetAcc);
	if(fScale<1E-6) fScale = 1E-6; //lower limit to scale
 	std::cout<<"BruSequentialProposal::CheckStepSize Changed to "<<fScale<<"(min 1E-6)) from "<<(fTargetAcc)/acc*fScale<<std::endl;
	return kTRUE;
      }
      return kFALSE;
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

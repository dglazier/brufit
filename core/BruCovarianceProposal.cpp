#include "BruCovarianceProposal.h"
#include <RooStats/RooStatsUtils.h>
#include <RooRandom.h>

namespace HS{
  namespace FIT{

    ///////////////////////////////////////////////////////////////////////////////
  
     BruCovarianceProposal::BruCovarianceProposal(float scale,float target,float accmin,float accmax) :
      BruSequentialProposal(scale, target, accmin, accmax),
      fPdf{nullptr},
      fCacheSize{1},
      fCachePosition{0},
      fCache{nullptr}
    {

    };
    ////////////////////////////////////////////////////////////////////////////////
    /// Populate xPrime with a new proposed point
 
    void BruCovarianceProposal::Propose(RooArgSet& xPrime, RooArgSet& x )
    {
      /*
      if (fLastX.empty()) {
	// fLastX not yet initialized
	fLastX.addClone(x);
	// generate initial cache
	RooStats::SetParameters(&x, &fMaster);
	if (fMap.size() > 0) {
	  for (fIt = fMap.begin(); fIt != fMap.end(); fIt++)
            fIt->first->setVal(fIt->second->getVal(&x));
	}
	fCache.reset(fPdf->generate(xPrime, fCacheSize));
      }
      */
      // bool moved = false;
      if (fMap.size() > 0) {
	//Just assume moved
	RooStats::SetParameters(&x, &fMaster);
	
	for (fIt = fMap.begin(); fIt != fMap.end(); fIt++)
	  fIt->first->setVal(fIt->second->getVal(&x));

	  /*
	moved = !Equals(fLastX, x);

	// if we've moved, set the values of the variables in the PDF to the
	// corresponding values of the variables in x, according to the
	// mappings (i.e. let the variables in x set the given values for the
	// PDF that will generate xPrime)
	if (moved) {
	  // update the pdf parameters
	  RooStats::SetParameters(&x, &fMaster);

	  for (fIt = fMap.begin(); fIt != fMap.end(); fIt++)
            fIt->first->setVal(fIt->second->getVal(&x));

	  // save the new x in fLastX
	  RooStats::SetParameters(&x, &fLastX);
	  
	}
	  */
      }

    /*
      // generate new cache if necessary
      if (moved || fCachePosition >= fCacheSize) {
	delete fCache;
	fCache = fPdf->generate(xPrime, fCacheSize);
	fCachePosition = 0;
      }

      const RooArgSet* proposal = fCache->get(fCachePosition);
      fCachePosition++;
    */


    fCache.reset(fPdf->generate(xPrime, 1));
    const RooArgSet* proposal = fCache->get(0);
    RooStats::SetParameters(proposal, &xPrime);
     
    }
    
 
    bool BruCovarianceProposal::IsSymmetric(RooArgSet& , RooArgSet& ) {
      return true;
    }
 
    ////////////////////////////////////////////////////////////////////////////////
    /// Return the probability of proposing the point x1 given the starting
    /// point x2
    double BruCovarianceProposal::GetProposalDensity(RooArgSet& ,
						     RooArgSet& )
    {
      return 1.0; // should not be needed
    }
 
    ////////////////////////////////////////////////////////////////////////////////
    /// specify a mapping between a parameter of the proposal function and
    /// a parameter of interest.  this mapping is used to set the value of
    /// proposalParam equal to the value of update to determine the
    /// proposal function.
    /// proposalParam is a parameter of the proposal function that must
    /// be set to the value of update (from the current point) in order to
    /// propose a new point.

    void BruCovarianceProposal::AddMapping(RooRealVar& proposalParam, RooAbsReal& update)
    {
      fMaster.add(*update.getParameters(static_cast<RooAbsData const*>(nullptr)));
      if (update.getParameters(static_cast<RooAbsData const*>(nullptr))->empty())
	fMaster.add(update);
      fMap.insert(std::pair<RooRealVar*, RooAbsReal*>(&proposalParam, &update));
    }

    void BruCovarianceProposal::SetCovariance(const TMatrixDSym& mat,const RooArgSet& vars){
      //Based on void ProposalHelper::CreatePdf()
      if (vars.getSize()!=mat.GetNcols()) {
	std::cerr << "BruCovarianceProposal::SetCovariance: " <<
	  "Variables to create proposal function for are not same as matrix" << std::endl;
	exit(0);
      }
      if(_covMatrix.GetNrows()==0)_covMatrix.ResizeTo(mat.GetNrows(),mat.GetNcols());
      _covMatrix=mat*static_cast<Double_t>(StepSizeFactor());
      _covMatrix.Print();
         
      if(_xVec.getSize()==0){
	//_vars = RooArgSet(vars);
	 	for (auto *r : static_range_cast<RooRealVar *> (vars)){
	  //make an offset variable mu for each var
	  _xVec.add(*r);
	  
	  TString cloneName = TString::Format("%s%s", "mu__", r->GetName());
	  auto clone = static_cast<RooRealVar*>(r->clone(cloneName.Data()));
	  _muVec.add(*clone);
	  AddMapping(*clone, *r);
	}
      }

        
      fPdf.reset(new RooMultiVarGaussian("mvg", "MVG Proposal", _xVec, _muVec, _covMatrix));
      }
  }
}

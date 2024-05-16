////////////////////////////////////////////////////////////////
///
///Class:           BruCovarianceProposal    
///Description:  Just sequential proposal with adaptove step size to
///              keep in some requested range
///           

#pragma once

#include "BruSequentialProposal.h"

#include <TMath.h>
#include <RooArgSet.h>
#include <RooMsgService.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooMultiVarGaussian.h>
#include <vector>
#include <map>

namespace HS{
  namespace FIT{

    class BruCovarianceProposal : public  BruSequentialProposal{

    public:
      BruCovarianceProposal() : BruSequentialProposal{} {}
      
      BruCovarianceProposal(float scale,float target=0.234,float accmin=0.15,float accmax=0.35);
      
      /// Populate xPrime with a new proposed point
      void Propose(RooArgSet& xPrime, RooArgSet& x) override;
 
      /// Determine whether or not the proposal density is symmetric for
      /// points x1 and x2 - that is, whether the probability of reaching x2
      /// from x1 is equal to the probability of reaching x1 from x2
      bool IsSymmetric(RooArgSet& x1, RooArgSet& x2) override ;
 
      /// Return the probability of proposing the point x1 given the starting
      /// point x2
      double GetProposalDensity(RooArgSet& x1, RooArgSet& x2) override;
 
      //  ~BruCovarianceProposal() override =default;

      /// Set the PDF to be the proposal density function
      virtual void AddMapping(RooRealVar& proposalParam, RooAbsReal& update);
      virtual void Reset()
      {
	fCache.reset();
	fCachePosition = 0;
	fLastX.removeAll();

	fMaster.clear();
	fMap.clear();

	_covMatrix.Clear();
	_covMatrix.ResizeTo(0,0);
	_covMatrix = TMatrixDSym();
	_xVec.clear();
	_muVec.clear();
       
 
      }

      virtual void printMappings()
      {
         std::map<RooRealVar*, RooAbsReal*>::iterator it;
         for (it = fMap.begin(); it != fMap.end(); it++)
            std::cout << it->first->GetName() << " => " << it->second->GetName() << std::endl;
      }

      /// Set how many points to generate each time we propose from a new point
      /// Default (and minimum) is 1
      virtual void SetCacheSize(Int_t size)
      {
         if (size > 0)
            fCacheSize = size;
         else
            coutE(Eval) << "Warning: Requested non-positive cache size: " <<
               size << ". Cache size unchanged." << std::endl;
      }

      /// set whether we own the PDF that serves as the proposal density function
      /// By default, when constructed, PdfProposal does NOT own the PDF.
  
      void SetCovariance(const TMatrixDSym& mat,const RooArgSet& vars);

      Bool_t CheckStepSize(Float_t acceptance) override{
	//changing step size with covariance gives a "blotchy" posterior
	if(_tuneCovStep==kFALSE){
	  //SetCovariance(_covMatrix*static_cast<Double_t>(TMath::Sqrt(StepSizeFactor())),_xVec);
	  return kTRUE;
	}
	auto doExit = BruSequentialProposal::CheckStepSize(acceptance);
	SetCovariance(_covMatrix*static_cast<Double_t>(TMath::Sqrt(StepSizeFactor())),_xVec);
	std::cout<< "CheckStepSize "<<_tuneCovStep<<" "<<StepSizeFactor()<<std::endl;
	//if tuning covariance matrix start again with new step size
	return kFALSE;

    }
      void TuneCovarianceStep(Bool_t tune){_tuneCovStep=tune;}
 
    private:
 
      std::unique_ptr<RooMultiVarGaussian> fPdf={nullptr}; //! the proposal density function
      std::map<RooRealVar*, RooAbsReal*> fMap; //! map of values in pdf to update
      std::map<RooRealVar*, RooAbsReal*>::iterator fIt; //! pdf iterator
      RooArgSet fLastX; //! the last point we were at
      Int_t fCacheSize; /// how many points to generate each time
      Int_t fCachePosition; /// our position in the cached proposal data set
      std::unique_ptr<RooDataSet> fCache={nullptr}; //! the cached proposal data set
      RooArgSet fMaster; //! pointers to master variables needed for updates
 
      TMatrixDSym _covMatrix;//!

      RooArgList _xVec;
      RooArgList _muVec;
       
      Bool_t _tuneCovStep=kFALSE;

      ClassDefOverride(BruCovarianceProposal,1) // A concrete implementation of ProposalFunction, that uniformly samples the parameter space.
 


    };

  }
}

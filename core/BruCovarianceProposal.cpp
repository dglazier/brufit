#include "BruCovarianceProposal.h"
#include <RooStats/RooStatsUtils.h>
#include <RooRandom.h>
#include <TVectorD.h>
#include <TDecompChol.h>
#include <cmath>

namespace HS{
  namespace FIT{

    BruCovarianceProposal::BruCovarianceProposal(float scale,float target,float accmin,float accmax) :
      BruSequentialProposal(scale, target, accmin, accmax)
    {
    }

    void BruCovarianceProposal::SetCovariance(const TMatrixDSym& mat, const RooArgSet& vars){
      Reset();

      if (vars.getSize() != mat.GetNcols()) {
        std::cerr << "BruCovarianceProposal::SetCovariance: Variables not same as matrix!" << std::endl;
        exit(0);
      }
      
      _baseMatrix.ResizeTo(mat.GetNrows(), mat.GetNcols());
      _baseMatrix = mat;

      if(_xVec.getSize() == 0) {
        for (auto *r : static_range_cast<RooRealVar *>(vars)){
          _xVec.add(*r);
          _varCache.push_back(dynamic_cast<RooRealVar*>(r));
        }
      }
      
      UpdateCholesky();
    }

void BruCovarianceProposal::UpdateCholesky() {
      if (_baseMatrix.GetNrows() == 0) return;

      int d = _baseMatrix.GetNrows();
      // The universally optimal Gelman MCMC scaling factor for d dimensions
      double gelmanScale = (2.38 * 2.38) / (double)d; 
 
      _covMatrix.ResizeTo(_baseMatrix.GetNrows(), _baseMatrix.GetNcols());
      // Apply the user's StepSizeFactor on top of the native Gelman scaling
      _covMatrix = _baseMatrix * static_cast<Double_t>(StepSizeFactor() * gelmanScale);
      
      double jitter = 1e-12; 
      bool decomposed = false;
      
      for (int attempts = 0; attempts < 10; attempts++) {
          TDecompChol chol(_covMatrix);
          if (chol.Decompose()) {
              _lMatrix.ResizeTo(_covMatrix.GetNrows(), _covMatrix.GetNcols());
              _lMatrix = chol.GetU(); 
              _lMatrix.T(); 
              decomposed = true;
              break; 
          } else {
              for (int i = 0; i < _covMatrix.GetNrows(); ++i) {
                  _covMatrix(i, i) += jitter;
              }
              jitter *= 10.0; 
          }
      }
        if (!decomposed) {
          std::cerr << "BruCovarianceProposal: FATAL - Matrix is fundamentally singular even with maximum jitter!" << std::endl;
      }
   }
    
    void BruCovarianceProposal::Propose(RooArgSet& xPrime, RooArgSet& x )
    {
      RooStats::SetParameters(&x, &xPrime); // Start at current location
      int n = _varCache.size();

      // 1. Generate independent standard normal steps
      TVectorD z(n);
      for (int i = 0; i < n; ++i) {
          z[i] = RooRandom::gaussian();
      }

      // 2. Instantly correlate them via Cholesky decomposition (L * z)
      TVectorD step = _lMatrix * z;

      // 3. Apply steps to parameters safely
      for (int i = 0; i < n; ++i) {
          RooRealVar* baseVar = _varCache[i];
          RooRealVar* primeVar = dynamic_cast<RooRealVar*>(xPrime.find(baseVar->GetName()));
          if (!primeVar) continue;

          double val = primeVar->getVal();
          double max = primeVar->getMax();
          double min = primeVar->getMin();
          double len = max - min;
          double s = step[i];

          if (fCyclicPars.contains(*primeVar)) {
              // --- CYCLIC WRAP AROUND ---
              double center = min + (len / 2.0);
              primeVar->setVal(center + std::remainder(val + s - center, len));
          } else {
              // --- HARD BOUNDARY REFLECTION ---
              // Bouncing preserves the local detailed balance of the matrix!
              while ((val + s > max) || (val + s < min)) {
                  if (val + s > max) s = 2 * (max - val) - s;
                  if (val + s < min) s = 2 * (min - val) - s;
              }
              primeVar->setVal(val + s);
          }
      }
    }
    
    bool BruCovarianceProposal::IsSymmetric(RooArgSet& , RooArgSet& ) {
      return true;
    }
 
    double BruCovarianceProposal::GetProposalDensity(RooArgSet&, RooArgSet&) {
      return 1.0; 
    }
  }
}// #include "BruCovarianceProposal.h"
// #include <RooStats/RooStatsUtils.h>
// #include <RooRandom.h>

// namespace HS{
//   namespace FIT{

//     ///////////////////////////////////////////////////////////////////////////////
  
//      BruCovarianceProposal::BruCovarianceProposal(float scale,float target,float accmin,float accmax) :
//       BruSequentialProposal(scale, target, accmin, accmax),
//       fPdf{nullptr},
//       fCacheSize{1},
//       fCachePosition{0},
//       fCache{nullptr}
//     {

//     };
//     ////////////////////////////////////////////////////////////////////////////////
//     /// Populate xPrime with a new proposed point
 
//     void BruCovarianceProposal::Propose(RooArgSet& xPrime, RooArgSet& x )
//     {
//       /*
//       if (fLastX.empty()) {
// 	// fLastX not yet initialized
// 	fLastX.addClone(x);
// 	// generate initial cache
// 	RooStats::SetParameters(&x, &fMaster);
// 	if (fMap.size() > 0) {
// 	  for (fIt = fMap.begin(); fIt != fMap.end(); fIt++)
//             fIt->first->setVal(fIt->second->getVal(&x));
// 	}
// 	fCache.reset(fPdf->generate(xPrime, fCacheSize));
//       }
//       */
//       // bool moved = false;
//       if (fMap.size() > 0) {
// 	//Just assume moved
// 	RooStats::SetParameters(&x, &fMaster);
	
// 	for (fIt = fMap.begin(); fIt != fMap.end(); fIt++)
// 	  fIt->first->setVal(fIt->second->getVal(&x));

// 	  /*
// 	moved = !Equals(fLastX, x);

// 	// if we've moved, set the values of the variables in the PDF to the
// 	// corresponding values of the variables in x, according to the
// 	// mappings (i.e. let the variables in x set the given values for the
// 	// PDF that will generate xPrime)
// 	if (moved) {
// 	  // update the pdf parameters
// 	  RooStats::SetParameters(&x, &fMaster);

// 	  for (fIt = fMap.begin(); fIt != fMap.end(); fIt++)
//             fIt->first->setVal(fIt->second->getVal(&x));

// 	  // save the new x in fLastX
// 	  RooStats::SetParameters(&x, &fLastX);
	  
// 	}
// 	  */
//       }

//     /*
//       // generate new cache if necessary
//       if (moved || fCachePosition >= fCacheSize) {
// 	delete fCache;
// 	fCache = fPdf->generate(xPrime, fCacheSize);
// 	fCachePosition = 0;
//       }

//       const RooArgSet* proposal = fCache->get(fCachePosition);
//       fCachePosition++;
//     */


//     fCache.reset(fPdf->generate(xPrime, 1));
//     const RooArgSet* proposal = fCache->get(0);
//     RooStats::SetParameters(proposal, &xPrime);
     
//     }
    
 
//     bool BruCovarianceProposal::IsSymmetric(RooArgSet& , RooArgSet& ) {
//       return true;
//     }
 
//     ////////////////////////////////////////////////////////////////////////////////
//     /// Return the probability of proposing the point x1 given the starting
//     /// point x2
//     double BruCovarianceProposal::GetProposalDensity(RooArgSet& ,
// 						     RooArgSet& )
//     {
//       return 1.0; // should not be needed
//     }
 
//     ////////////////////////////////////////////////////////////////////////////////
//     /// specify a mapping between a parameter of the proposal function and
//     /// a parameter of interest.  this mapping is used to set the value of
//     /// proposalParam equal to the value of update to determine the
//     /// proposal function.
//     /// proposalParam is a parameter of the proposal function that must
//     /// be set to the value of update (from the current point) in order to
//     /// propose a new point.

//     void BruCovarianceProposal::AddMapping(RooRealVar& proposalParam, RooAbsReal& update)
//     {
//       fMaster.add(*update.getParameters(static_cast<RooAbsData const*>(nullptr)));
//       if (update.getParameters(static_cast<RooAbsData const*>(nullptr))->empty())
// 	fMaster.add(update);
//       fMap.insert(std::pair<RooRealVar*, RooAbsReal*>(&proposalParam, &update));
//     }

//     void BruCovarianceProposal::SetCovariance(const TMatrixDSym& mat,const RooArgSet& vars){
//       Reset();

//       //Based on void ProposalHelper::CreatePdf()
//       if (vars.getSize()!=mat.GetNcols()) {
// 	std::cerr << "BruCovarianceProposal::SetCovariance: " <<
// 	  "Variables to create proposal function for are not same as matrix" << std::endl;
// 	exit(0);
//       }
//       if(_covMatrix.GetNrows()==0)_covMatrix.ResizeTo(mat.GetNrows(),mat.GetNcols());
//       _covMatrix=mat*static_cast<Double_t>(StepSizeFactor());
         
//       if(_xVec.getSize()==0){
// 	//	std::cout<<" BruCovarianceProposal::SetCovariance  static_range_cast does not work until 6.28 c++14"<<std::endl;exit(0);
// 	for (auto *r : static_range_cast<RooRealVar *> (vars)){
// 	  //make an offset variable mu for each var
// 	  _xVec.add(*r);
	  
// 	  TString cloneName = TString::Format("%s%s", "mu__", r->GetName());
// 	  auto clone = static_cast<RooRealVar*>(r->clone(cloneName.Data()));
// 	  _muVec.addOwned(*clone);
// 	  AddMapping(*clone, *r);
// 	}
	
//       }

        
//       fPdf.reset(new RooMultiVarGaussian("mvg", "MVG Proposal", _xVec, _muVec, _covMatrix));
//       }
//   }
// }

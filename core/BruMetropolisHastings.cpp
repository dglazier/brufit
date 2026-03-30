#include "BruMetropolisHastings.h"
#include <RooStats/RooStatsUtils.h>
#include <RooStats/MarkovChain.h>
#include <RooStats/PdfProposal.h>
#include <RooStats/SequentialProposal.h>
#include <RooStats/ProposalHelper.h>
#include <TRandom.h>
#include "BruSequentialProposal.h"

namespace HS{
  namespace FIT{
    
   
    RooStats::MarkovChain* BruMetropolisHastings::ConstructChain()
    {
      if (fParameters.getSize() == 0 || !fPropFunc || !fFunction) {
        coutE(Eval) << "Critical members unintialized: parameters, proposal function, or (log) likelihood function" << std::endl;
        return nullptr;
      }
  
      if (fChainParams.getSize() == 0) fChainParams.add(fParameters);
      
      RooArgSet x;
      RooArgSet xPrime;
      x.addClone(fParameters);
      xPrime.addClone(fParameters);
   
      auto* chain = new RooStats::MarkovChain();
      chain->SetParameters(fChainParams);

      // --- 1. FAST MEMORY MAPPING ---
      // Map the RooRealVars into fast C++ vectors to bypass string-lookups during the loop
      std::vector<RooRealVar*> curVars, propVars, masterVars;
      for (auto *var : static_range_cast<RooRealVar*>(x)) curVars.push_back(var);
      for (auto *var : static_range_cast<RooRealVar*>(xPrime)) propVars.push_back(var);
      for (auto *var : static_range_cast<RooRealVar*>(fParameters)) masterVars.push_back(var);
      size_t nVars = curVars.size();

      // --- 2. CHAIN BUFFER ---
      // Buffer the accepted steps in RAM to avoid RooDataSet overhead during the loop
      std::vector<std::vector<Double_t>> chainBufferVals;
      std::vector<Double_t> chainBufferNLL;
      std::vector<Double_t> chainBufferWeights;
      chainBufferVals.reserve(fNumIters);
      chainBufferNLL.reserve(fNumIters);
      chainBufferWeights.reserve(fNumIters);

      Int_t weight = 0;
      Double_t xL = 0.0, xPrimeL = 0.0, a = 0.0;

      RooFit::MsgLevel oldMsgLevel = RooMsgService::instance().globalKillBelow();
      RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);
      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
      RooAbsReal::clearEvalErrorLog();

      bool hadEvalError = true;
      ooccoutP((TObject *)nullptr, Generation) << "Metropolis-Hastings progress: ";

      int icount = 0;
      int totcount = 0;
      int snapcount = 0;
      int havePrinted = 0;
      Long64_t totalProposals = 0;
      
      auto bseqprop = dynamic_cast<BruSequentialProposal*>(fPropFunc);
      if(bseqprop == nullptr){
        std::cerr<<"ERROR BruMetropolisHastings need a BruSequentialProposal "<<std::endl;
        delete chain;
        return nullptr;
      }
      bseqprop->ResetCounter();

      // --- INITIAL EVALUATION ---
      for (size_t i = 0; i < nVars; ++i) curVars[i]->setVal(masterVars[i]->getVal());
      xL = fFunction->getVal();

      // --- THE BARE-METAL MCMC LOOP ---
      while (icount < fNumIters) {

	// --- STOCHASTIC HOT-SWAP ---
	if (fBatchedNLLs.size() > 0 && fSwapFreq > 0 && totalProposals > 0 && (totalProposals % fSwapFreq == 0)) {
          fCurrentBatch = (fCurrentBatch + 1) % fBatchedNLLs.size();
          fFunction = fBatchedNLLs[fCurrentBatch]; 
          std::cout << "Shaking the landscape! Swapped to Batch " << fCurrentBatch << std::endl;
          
          // =======================================================
          // THE FIX: Reset the Baseline NLL safely
          // =======================================================
          // 1. Force the master RooRealVars back to the current accepted 
          //    state (curVars). If the last step was rejected, masterVars 
          //    are currently stuck at the bad proposed values!
          for (size_t i = 0; i < nVars; ++i) {
	    masterVars[i]->setVal(curVars[i]->getVal());
          }
          
          // 2. Re-evaluate the baseline NLL on the new batch
          xL = fFunction->getVal(); 
	}

        totcount++; 
        totalProposals++;  
        hadEvalError = false;

	// ==========================================================
	// --- FAST-FAIL BAILOUT ---
	// If we are grinding infinitely without accepts, abort the 
	// chain early to allow the TuneCovarianceStep to rescue the scale.
	// ==========================================================

        if (totalProposals % 1000 == 0) {
            double currentAcc = (double)icount / totalProposals;
            
            if (totalProposals >= 3000 && currentAcc < 0.01) {
                std::cout << "\nBruMetropolisHastings: WARNING - Fast Bailout Triggered! "
                          << "Acceptance is critically low (" << currentAcc * 100 << "%) after " 
                          << totalProposals << " proposals. Aborting chain to allow immediate retuning." << std::endl;
                break;
            }
        }

	// ==========================================================
        // Acceptance Tuning Logic
        if (totcount % 1000 == 0) {
          fAcceptance = ((Double_t)snapcount) / totcount;
          std::cout << "  BruMetropolisHastings accepted " << snapcount << " out of " << totcount << " for acceptance " << fAcceptance << std::endl;
          
          RooMsgService::instance().setGlobalKillBelow(oldMsgLevel);
          auto ok = bseqprop->CheckStepSize(fAcceptance);
          if (ok == kFALSE) {
            delete chain;
            return nullptr; 
          }
          snapcount = 0;
          totcount = 0;
          continue;
        }

        // Progress Printing
        if (icount % 100 == 0 && havePrinted == 0) {
          ooccoutP((TObject*)nullptr, Generation) << " " << icount << "/" << fNumIters;
          havePrinted = 1;
        }
        if (icount % 100 == 1) havePrinted = 0;

        // Propose new step
        fPropFunc->Propose(xPrime, x);
        
        // Fast push proposed values to the master NLL parameters
        for (size_t i = 0; i < nVars; ++i) masterVars[i]->setVal(propVars[i]->getVal());
        xPrimeL = fFunction->getVal();
        
        if (wasEvalErrors()) {
          xPrimeL = RooNumber::infinity();
          hadEvalError = true;
        }

        a = xPrimeL - xL;
       
        if (!hadEvalError && !fPropFunc->IsSymmetric(xPrime, x)) {
          Double_t xPrimePD = fPropFunc->GetProposalDensity(xPrime, x);
          Double_t xPD      = fPropFunc->GetProposalDensity(x, xPrime);
          a += TMath::Log(xPrimePD) - TMath::Log(xPD);
        }
       
        if (!hadEvalError && ShouldTakeStep(a)) {
          // ACCEPTED: Record CURRENT state (x) with its accumulated weight
          if (weight != 0.0) {
	    std::vector<Double_t> stepVals(nVars);
	    for (size_t i = 0; i < nVars; ++i) stepVals[i] = curVars[i]->getVal();
              
	    chainBufferVals.push_back(std::move(stepVals));
	    chainBufferNLL.push_back(CalcNLL(xL));
	    chainBufferWeights.push_back((Double_t)weight);
              
	    icount++;
	    snapcount++;
          }
          
          // Jump to new state (xPrime)
          weight = 1;
          for (size_t i = 0; i < nVars; ++i) curVars[i]->setVal(propVars[i]->getVal());
          xL = xPrimeL;
        } else {
          // REJECTED: Stay at current state, increment weight
          weight++;
        }
      }

      // Record the final point
      if (weight != 0.0) {
	std::vector<Double_t> stepVals(nVars);
	for (size_t i = 0; i < nVars; ++i) stepVals[i] = curVars[i]->getVal();
	chainBufferVals.push_back(std::move(stepVals));
	chainBufferNLL.push_back(CalcNLL(xL));
	chainBufferWeights.push_back((Double_t)weight);
      }
      
      ooccoutP((TObject *)nullptr, Generation) << std::endl;

      // --- 3. FLUSH BUFFER TO CHAIN ---
      // Now we pay the RooDataSet memory penalty exactly once at the very end
      for(size_t i = 0; i < chainBufferVals.size(); ++i) {
	for (size_t v = 0; v < nVars; ++v) masterVars[v]->setVal(chainBufferVals[i][v]);
	chain->Add(fParameters, chainBufferNLL[i], chainBufferWeights[i]);
      }

      RooMsgService::instance().setGlobalKillBelow(oldMsgLevel);

      Int_t numAccepted = chain->Size();
      fAcceptance = ((Double_t)icount) / totalProposals;
      
      coutI(Eval) << "Proposal acceptance rate: " << ((float)icount)/totalProposals * 100 << "%" << std::endl;
      coutI(Eval) << "Number of steps in chain: " << numAccepted << std::endl;
 
      return chain;
    }

    Bool_t BruMetropolisHastings::wasEvalErrors(){
      if(RooAbsReal::numEvalErrors()){
	RooAbsReal::clearEvalErrorLog();
	return kTRUE;
      }
      return kFALSE;
	    
    }

    Bool_t BruMetropolisHastings::CheckForBurnIn(RooStats::MarkovChain* chain){
      return true;
      //      Check RMS and mean for last 100 events
      auto Nentries=chain->Size();
      auto vars = chain->Get();
      RooDataSet current("current","current",*vars);
     
      Int_t Ntests=vars->getSize()*50;
      
      if(Ntests>Nentries){
	//	std::cout<<" ntests "<<Ntests<<" "<<Nentries<<std::endl; return kTRUE;
      }//not enough events Ntests=Nentries/2;
      if(Ntests>Nentries-fLastEntries){std::cout<<" not enough entries "<<Ntests<<" > "<<Nentries-fLastEntries<<" "<<Nentries<<std::endl; return kTRUE;}
 
      for(Int_t i=Nentries-1;i>Nentries-100;--i){
	if(gRandom->Uniform()<1./chain->Weight(i))current.add(*chain->Get(i));
      }

      if(fMeans.empty()==true){
	fMeans.resize(vars->getSize());
	fSigmas.resize(vars->getSize());
      }
      //loop over variables and get current mean and rms
      Int_t iv=0;
      Int_t Npass=0;
      for(auto& var: *vars){
	RooRealVar* rvar=static_cast<RooRealVar*>(var);
	Double_t mean=current.mean(*rvar);
	Double_t sigma=current.sigma(*rvar);
	//std::cout<<"BruMetropolisHastings CheckForBurnIn "<<var->GetName()<<" mean "<<mean<<" "<<sigma<<" out of "<<current.numEntries()<<" diff "<< (mean-fMeans[iv])/TMath::Sqrt(fSigmas[iv]*fSigmas[iv]+sigma*sigma)<<std::endl;

	auto biggestSigma=fSigmas[iv]>sigma? fSigmas[iv]:sigma;
	if(biggestSigma==0)biggestSigma=1;
	//std::cout<<"sigma "<<biggestSigma<<" "<<fSigmas[iv]<<" "<<fMeans[iv]<<" "<<TMath::Abs((mean-fMeans[iv])/biggestSigma)<<std::endl;
	
	if(TMath::Abs((mean-fMeans[iv])/biggestSigma)<3)++Npass;
	

	if(var==(*vars)[vars->getSize()-1]){
	  std::cout<<"NLL diff "<<mean<<" "<<fMeans[iv]<<" diff "<<mean-fMeans[iv]<<" "<<TMath::Abs((mean-fMeans[iv])/biggestSigma)<<" "<<std::endl;

	  if(mean-fMeans[iv]>0){
	    ++fNWorse;
	    if(fSaveNLL==0)
	      fSaveNLL=mean;
	  }
	  
	  else if(fNWorse>0) ++fNBetter;

	  
	  if(mean-fSaveNLL>0)
	    ++fNWorseThanSave;
	  if(mean-fSaveNLL<0)
	    ++fNBetterThanSave;
	  
	}
	
	fMeans[iv]=mean;
	fSigmas[iv]=sigma;
	++iv;
      }
      std::cout<<Npass<<" out of "<<vars->getSize()<<" up to "<<Nentries<<" and Worse "<<fNWorse<<" and Better "<<fNBetter<<" and Worse than save "<<fNWorseThanSave<<" and Better "<<fNBetterThanSave<<" which is " <<fSaveNLL<<std::endl;
      fLastEntries=Nentries;
      return kTRUE;
    } 
  }//namespace FIT

}//namespace HS

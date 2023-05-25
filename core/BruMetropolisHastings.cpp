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
	coutE(Eval) << "Critical members unintialized: parameters, proposal " <<
	  " function, or (log) likelihood function" << std::endl;
	return nullptr;
      }
  
      if (fChainParams.getSize() == 0) fChainParams.add(fParameters);
      
      // RooRealVar NSinceIntChanged("NSinceIntChanged","NSinceIntChanged",0,0,1E12);
      //fChainParams.add(NSinceIntChanged);
      
      RooArgSet x;
      RooArgSet xPrime;
      x.addClone(fParameters);
      xPrime.addClone(fParameters);
   
      //      xChain.add(x);
      //      xChain.add(NSinceIntChanged);
      
      
      auto* chain = new RooStats::MarkovChain();
      // only the POI will be added to the chain
      chain->SetParameters(fChainParams);

      Int_t weight = 0;
      Double_t xL = 0.0, xPrimeL = 0.0, a = 0.0;

      // ibucur: i think the user should have the possibility to display all the message
      //    levels should they want to; maybe a setPrintLevel would be appropriate
      //    (maybe for the other classes that use this approach as well)?
      RooFit::MsgLevel oldMsgLevel = RooMsgService::instance().globalKillBelow();
      RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);

      // We will need to check if log-likelihood evaluation left an error status.
      // Now using faster eval error logging with CountErrors.
      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
      //N.B: need to clear the count in case of previous errors !
      // the clear needs also to be done after calling setEvalErrorLoggingMode
      RooAbsReal::clearEvalErrorLog();

      bool hadEvalError = true;

    
      ooccoutP((TObject *)nullptr, Generation) << "Metropolis-Hastings progress: ";

      // do main loop
      //for (i = 0; i < fNumIters; i++) {
      int icount=0;
      int totcount=0;
      int snapcount=0;
      int NsinceIntegralChanged=0;
      int havePrinted=0;
      //x.Print("v");
      
      //Check if we have PDF proposal
      
      RooStats::PdfProposal* pdfProposal=dynamic_cast<RooStats::PdfProposal*>(fPropFunc);
      std::unique_ptr<RooStats::SequentialProposal> sp{};
      //std::cout<<"DEBUG RooStats::MarkovChain* BruMetropolisHastings::ConstructChain"<<std::endl; fParameters.Print();

      Long64_t totalProposals=0;
      while (icount <fNumIters) {
	totcount++; //can be reset 
	totalProposals++;  //will not be reset
	
	//std::cout<<"DEBUG RooStats::MarkovChain another one "<<std::endl;
	//std::cout<<"        METHAST "<<totcount<<std::endl;
	// reset error handling flag
	hadEvalError = false;

	//recalc XL with sampled integral
	/*	auto xLtest = fFunction->getVal();//DEBUGGING
	Double_t newBalance = fBalancePDF==nullptr?1: fBalancePDF->getVal();
	if(TMath::Abs(newBalance-fOldBalance)>1E-8){
	  xL=xLtest;
	  NSinceIntChanged.setVal(0);
	  //std::cout<<"Changed integral for CHAIN "<<newBalance <<" "<<fOldBalance<<std::endl;
	}
	else{
	  NSinceIntChanged.setVal(NSinceIntChanged.getVal()+1);
	  //std::cout<<"Increment for CHAIN"<<newBalance <<" "<<fOldBalance<<std::endl;

	}
	*/
	
	// print a dot every 1% of the chain construction
	if (totcount%1000 == 0){

	  
	  fAcceptance = ((Double_t)snapcount)/totcount;
	  std::cout<<"  BruMetropolisHastings accepted "<<snapcount<<" out of "<<totcount<<" for acceptance "<<fAcceptance<<std::endl;
	  // CheckForBurnIn(chain);
	  //check acceptance, if too low exit
	  //if(fTryHelp==kTRUE){
	  //if(fAcceptance<fMinAcc||fAcceptance>fMaxAcc){
	  RooMsgService::instance().setGlobalKillBelow(oldMsgLevel);
	  //	  std::cout<<"WARNING BruMetropolisHastings acceptance not optimal exiting..." <<" current Likelihood "<<CalcNLL(xL)<<" "<<fPropFunc->ClassName()<<" "<<dynamic_cast<RooStats::SequentialProposal*>(fPropFunc)<<std::endl;

	  // x.Print("v");
	  auto bseqprop=dynamic_cast<BruSequentialProposal*>(fPropFunc);
	  if(bseqprop!=nullptr){
	    //	    sp.reset(new RooStats::SequentialProposal(fNorm));
	    //std::cout<<"new seq "<<sp.get()<<std::endl;
	    //SetProposalFunction(*sp.get());
	    //std::cout<<"new seq "<<sp.get()<<std::endl;
	    bseqprop->CheckStepSize(fAcceptance);
	    snapcount =0;
	    totcount=0;
	    continue;
	  }
	  else{
	    std::cerr<<"ERROR BruMetropolisHastings need a BruSequentialProposal "<<std::endl;
	    // std::cout<<"x "<<std::endl;
	    // x.Print("v");
	    // std::cout<<"xprime "<<std::endl;
	    // xPrime.Print("v");
	    //	exit(0);
	    delete chain;
	    return nullptr;
	  }
	
	  

	}
	if (icount%100 == 0&&havePrinted==0){
	  ooccoutP((TObject*)nullptr, Generation) << " "<<icount<<"/"<<fNumIters;
	  havePrinted=1;
	}
	if (icount%100 == 1) havePrinted=0;

	//PdfProposal uses a cache, but it samples from
	//the wrong Pdf parameters so best to reset every proposal
	if(pdfProposal){
	   pdfProposal->Reset();
	}

	if(xL==0)    {
	  RooStats::SetParameters(&x, &fParameters);
	  xL = fFunction->getVal();
	}
	fPropFunc->Propose(xPrime, x);
	RooStats::SetParameters(&xPrime, &fParameters);
	xPrimeL = fFunction->getVal();
	
	//	std::cout<<"********************************************L"<<xPrimeL<<" "<<xL<<std::endl;

	// check if log-likelihood for xprime had an error status
	if (wasEvalErrors()) {
	  xPrimeL = RooNumber::infinity();
	  hadEvalError = true;
	  
	  // --totcount;//dont't count for acceptance purposes
	  // std::cout<<"DEBUG  BruMetropolisHastings was errors "<<std::endl;
	  // for(Int_t ix =0;ix<x.getSize();ix++){
	  //   if(x.getRealValue(x[ix]->GetName())!=xPrime.getRealValue(x[ix]->GetName())) std::cout<<"\t"<<x[ix]->GetName()<<"\t"<<x.getRealValue(x[ix]->GetName())<<" "<<xPrime.getRealValue(x[ix]->GetName());
	  // }
	  // std::cout<<std::endl;
	}
	else{
	  // std::cout<<"**************DEBUG  BruMetropolisHastings no errors "<<std::endl;
	  // for(Int_t ix =0;ix<x.getSize();ix++){
	  //   if(x.getRealValue(x[ix]->GetName())!=xPrime.getRealValue(x[ix]->GetName())) std::cout<<"\t"<<x[ix]->GetName()<<"\t"<<x.getRealValue(x[ix]->GetName())<<" "<<xPrime.getRealValue(x[ix]->GetName());

	  //  }
	}

	// why evaluate the last point again, can't we cache it?
	// kbelasco: commenting out lines below to add/test caching support
	//RooStats::SetParameters(&x, &fParameters);
	//xL = fFunction->getVal();


       a = xPrimeL - xL;//**DO this one

       //std::cout<<"DEBUG RooStats::MarkovChain a "<<a<<" "<<xL<<" "<<xPrimeL<<" "<<icount<<" "<<snapcount<<" weight "<<weight<<" should "<<ShouldTakeStep(a)<<" error "<<hadEvalError<<std::endl;
       
       if (!hadEvalError && !fPropFunc->IsSymmetric(xPrime, x)) {
	  //Sequential proposal has density 1 so this does nothing
	  Double_t xPrimePD = fPropFunc->GetProposalDensity(xPrime, x);
	  Double_t xPD      = fPropFunc->GetProposalDensity(x, xPrime);

	  a += TMath::Log(xPrimePD) - TMath::Log(xPD);

	  if(TMath::IsNaN(a))std::cerr<<"BruMetropolisHastings::MakeChain() whn making symmtric s is nan"<<std::endl; exit(0);
	  
	  //std::cout<<"DEBUG RooStats::MarkovChain symmetric  "<<a<<ShouldTakeStep(a)<<" "<<TMath::Log(xPrimePD) - TMath::Log(xPD)<<" x'pd "<<xPrimePD<<" xpd "<<xPD<<std::endl;
	  //if(pdfProposal){
	  // std::cout<<"DEBUG RooStats::MarkovChain pdf "<<pdfProposal->GetPdf()->getVal(x)<<" "<<pdfProposal->GetPdf()->getVal(xPrime)<<std::endl;
	  //}
	

	}
       
	if (!hadEvalError && ShouldTakeStep(a)) {
	  // go to the proposed point xPrime
	  // add the current point with the current weight
	  // ? dglazier should this not have xPrimeL the current val
	  if (weight != 0.0){
	    //  std::cout<<"DEBUG RooStats::MarkovChain Add to chain "<<a<<std::endl;
	    //xPrime.Print("v");
	    chain->Add(xPrime, CalcNLL(xPrimeL),(Double_t)weight);
	    //	    chain->Add(xChain, CalcNLL(xPrimeL),(Double_t)weight);
	    //      fOldBalance=newBalance;
	    icount++;
	    snapcount++;
	  }
	  
	  // reset the weight and go to xPrime
	  weight = 1;
	  RooStats::SetParameters(&xPrime, &x);
	  xL = xPrimeL;
	  

	} else {
	  // stay at the current point
	  weight++;
	}
      }

      // make sure to add the last point
      if (weight != 0.0)
	chain->Add(x, CalcNLL(xPrimeL), (Double_t)weight);
      //chain->Add(xChain, CalcNLL(xPrimeL), (Double_t)weight);
      ooccoutP((TObject *)nullptr, Generation) << std::endl;

      RooMsgService::instance().setGlobalKillBelow(oldMsgLevel);

      Int_t numAccepted = chain->Size();
      coutI(Eval) << "Proposal acceptance rate: " <<
	((float)icount)/totalProposals * 100 << "%" << std::endl;
      coutI(Eval) << "Number of steps in chain: " << numAccepted << std::endl;

      //clear burnin checker
      // fMeans.clear();
      //fSigmas.clear();
      //fNWorse=0;
 
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

	  /*  if(fNBetterThanSave>fNWorseThanSave+5)
	    {
	      fNWorse=0;
	      fNBetter=0;
	      fNWorseThanSave=0;
	      fNBetterThanSave=0;
	      fSaveNLL=0;
	    }
	  */
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

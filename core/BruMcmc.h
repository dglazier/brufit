////////////////////////////////////////////////////////////////
///
///Class:               BruMcmc
///Description:
///           


#pragma once


#include "Minimiser.h"
#include "BruSequentialProposal.h"
#include "BruCovarianceProposal.h"
#include <RooAbsData.h>
#include <TTree.h>
#include <RooStats/ModelConfig.h>
#include <RooStats/MarkovChain.h>
#include <RooStats/ProposalFunction.h>
#include <RooFitResult.h>

namespace HS{
  namespace FIT{

    class BruMcmc  : public Minimiser {
      
    public:

      BruMcmc(Int_t Niter=100,Int_t Nburn=10, Float_t norm=0.1): fNumIters(Niter),fNumBurnInSteps(Nburn),fNorm(norm){
	SetNameTitle("HSBruMcmc","BruMcmc minimiser");
      }
      // BruMcmc(const BruMcmc&)=default;
      //BruMcmc(BruMcmc&&)=default;
      ~BruMcmc() override;
      //BruMcmc& operator=(const BruMcmc& other)=default;
      //BruMcmc& operator=(BruMcmc&& other) = default;  

      void Run(Setup &setup,RooAbsData &fitdata) override;

      void FitTo();
      
      file_uptr SaveInfo() override;
      void AddFormulaToMCMCTree();
      
      Bool_t MakeChain();
      TMatrixDSym MakeMinuitCovarianceMatrix();
      TMatrixDSym MakeMcmcCovarianceMatrix(TTree* tree,size_t burnin);
      TMatrixDSym MakeMcmcNonYieldCovarianceMatrix(TTree* tree,size_t burnin);
      TMatrixDSym MakeMcmcPrincipalCovarianceMatrix(TTree* tree,size_t burnin);
      TTree* GetTree(){return fTreeMCMC;}
      Double_t SumWeights();
      Double_t SumWeights2();
      
      void NoWeightCorrection(){fCorrectForWeights=kFALSE;}
      
      void Result();
      Double_t NLL(){return fChain->NLL();}



      //MCMCCalculator

      void SetModel(RooStats::ModelConfig *model);

      /// Set the DataSet if not already there
      void SetData(RooAbsData& data) { fData = &data; }

      /// Set the Pdf if not already there
      void SetPdf(RooAbsPdf& pdf) { fPdf = &pdf; }

      /// Set the Prior Pdf if not already there
      void SetPriorPdf(RooAbsPdf& pdf) { fPriorPdf = &pdf; }

      /// specify the parameters of interest in the interval
      void SetParameters(const RooArgSet& set) { fPOI.removeAll(); fPOI.add(set); }

      /// specify the parameters to store in the Markov chain
      /// By default all the parameters are stored
      void SetChainParameters(const RooArgSet & set) { fChainParams.removeAll(); fChainParams.add(set); }

      /// specify the nuisance parameters (eg. the rest of the parameters)
      void SetNuisanceParameters(const RooArgSet& set) {fNuisParams.removeAll(); fNuisParams.add(set);}

      /// set the conditional observables which will be used when creating the NLL
      /// so the pdf's will not be normalized on the conditional observables when computing the NLL
      void SetConditionalObservables(const RooArgSet& set) {fConditionalObs.removeAll(); fConditionalObs.add(set);}

      /// set the global observables which will be used when creating the NLL
      /// so the constraint pdf's will be normalized correctly on the global observables when computing the NLL
      void SetGlobalObservables(const RooArgSet& set) {fGlobalObs.removeAll(); fGlobalObs.add(set);}

      /// set the proposal function for suggesting new points for the MCMC
      void SetProposalFunction(RooStats::ProposalFunction& proposalFunction)
      { fPropFunc = &proposalFunction; }

      /// set the number of iterations to run the metropolis algorithm
      void SetNumIters(Int_t numIters)
      { fNumIters = numIters; }

      /// set the number of steps in the chain to discard as burn-in,
      /// starting from the first
      void SetNumBurnInSteps(Int_t numBurnInSteps)
      { fNumBurnInSteps = numBurnInSteps; }

      void SetupBasicUsage();
      void SetKeepStart(Bool_t keep=kTRUE){fKeepStart=keep;}

      virtual Int_t GetNumBurnInSteps()const {return fNumBurnInSteps;}

      void SetDesiredAcceptance(Double_t min,Double_t max,Double_t target=0){
	fMinAcc=min;
	fMaxAcc=max;
	if(target)
	  fTargetAcc=target;
	else
	  fTargetAcc = (max-min)/2;
      }
      void SetUncorrelateYields(Int_t un){fUncorrelateYields=un;}
      void SetParVals(RooArgSet* toThesePars);

      Bool_t Success(){return fChain != nullptr;}
    protected :
      void AddEntryBranch();
      void CleanMakeChain();
      
      RooStats::MarkovChain* fChain =nullptr; //!
      RooDataSet* fChainData=nullptr;//!
      TTree* fTreeMCMC=nullptr;//!
      Bool_t fCorrectForWeights=kTRUE;
      RooArgSet* fParams=nullptr;//!
      std::shared_ptr<TFile> fTempFile;//!
      file_uptr fOutFile;//!
      
      Bool_t fKeepStart=kFALSE; //randomise starting values
      Bool_t fMCMCHelp=kFALSE;//automate acceptance etc.
      //MCMCCalculator
      RooStats::ModelConfig *fModelConfig=nullptr;
      
      RooArgSet   fPOI;        //! parameters of interest for interval
      RooArgSet   fNuisParams; //! nuisance parameters for interval (not really used)
      RooArgSet   fChainParams; //! parameters to store in the chain (if not specified they are all of them )
      RooArgSet   fConditionalObs; //! conditional observables
      RooArgSet   fGlobalObs;     //! global observables
      RooStats::ProposalFunction* fPropFunc{}; //! Proposal function for MCMC integration
      RooAbsPdf * fPdf=nullptr;        //! pointer to common PDF (owned by the workspace)
      RooAbsPdf * fPriorPdf=nullptr;   //! pointer to prior  PDF (owned by the workspace)
      //RooAbsData * fData=nullptr;     //! pointer to the data (owned by the workspace)
      Int_t fNumIters; // number of iterations to run metropolis algorithm
      Int_t fNumBurnInSteps; // number of iterations to discard as burn-in, starting from the first

      Int_t fNumBins{}; // set the number of bins to create for each
      Int_t fWarmup{}; //ignore these events
      Float_t fNorm=1;
      Int_t fNumBurnInStepsCov; //Number of steps to remove from chain to make covariance matrix for proposal function

      vector<Double_t> _formVals;//(formulas.getSize(),0);
      vector<TBranch*> _formBranches;//(formulas.getSize(),nullptr);

      Double_t fChainAcceptance=0;//!
      Double_t fMinAcc=0.15;
      Double_t fMaxAcc=0.3;
      Double_t fTargetAcc=0.234;
      Int_t  fUncorrelateYields=0;
      
      ClassDefOverride(HS::FIT::BruMcmc,1);
      
     };

    class BruMcmcSeq  : public BruMcmc {
      
    public:

      BruMcmcSeq(Int_t Niter=100,Int_t Nburn=10, Float_t norm=0.1):
	BruMcmc(Niter,Nburn,norm){
	SetNameTitle("BruMcmcSeq","BruMcmcSeq minimiser");
      }
      //  BruMcmcSeq(const BruMcmcSeq&)=default;
      // BruMcmcSeq(BruMcmcSeq&&)=default;
      // ~BruMcmcSeq() override =default;
      // BruMcmcSeq& operator=(const BruMcmcSeq& other)=default;
      //BruMcmcSeq& operator=(BruMcmcSeq&& other) = default;  

      void Run(Setup &setup,RooAbsData &fitdata) override;

      ClassDefOverride(HS::FIT::BruMcmcSeq,1);
   };

    class BruMcmcSeqHelper  : public BruMcmc {
      
    public:
      
      BruMcmcSeqHelper(Int_t Niter=100,Int_t Nburn=10, Float_t norm=0.1,float target=0.234,float accmin=0.15,float accmax=0.35):BruMcmc(Niter,Nburn,norm),
	_proposal{norm,target,accmin,accmax}{
	SetNameTitle("HSBruMcmcSeqHelper","BruMcmcSeqHelper minimiser");
      }
      //~BruMcmcSeqHelper() override =default;
  
      void Run(Setup &setup,RooAbsData &fitdata) override;

 
    private:
      BruSequentialProposal _proposal;
      
      ClassDefOverride(HS::FIT::BruMcmcSeqHelper,1);
   };

   class BruMcmcCovariance  : public BruMcmc {
      
    public:
      
      BruMcmcCovariance(Int_t Niter=100,Int_t Nburn=10, Float_t norm=0.01,float target=0.234,float accmin=0.15,float accmax=0.35):
	BruMcmc(Niter,Nburn,norm),
	_propSeq{norm,target,accmin,accmax},
	//	_propSeq{norm,0.1,0.15,0.5}, //try for largish inial steps
	_propCov{norm,target,accmin,accmax}//1 for norm if covariance is correct
     {
       SetNameTitle("BruMcmcCovariance","BruMcmcCovariance minimiser");
     }
  
      void Run(Setup &setup,RooAbsData &fitdata) override;


     void TurnOffSequential(){_doSeq=kFALSE;}
     void TurnOffNDStep(){_doND=kFALSE;}
     void TurnOffCovariance(){_doCov=kFALSE;}
     void TuneCovarianceStep(){_tuneCovStep=kTRUE;}
     
    private:

     BruSequentialProposal _propSeq;
     BruCovarianceProposal _propCov;


     Bool_t _doSeq=kTRUE;
     Bool_t _doND=kTRUE;
     Bool_t _doCov=kTRUE;
     Bool_t _tuneCovStep=kFALSE;
     
      ClassDefOverride(HS::FIT::BruMcmcCovariance,1);
   };

  }//namespace FIT
}//namespace HS


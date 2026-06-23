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
#include <RooStats/MarkovChain.h>
#include <RooStats/ProposalFunction.h>
#include <RooFitResult.h>
#include <RooProduct.h>
#include <memory> 

namespace HS{
  namespace FIT{

    // --- BATCHED MCMC SUPPORT ---
      struct NLLCache {
          std::unique_ptr<RooAbsPdf> localProdPdf;
          std::unique_ptr<RooAbsReal> baseNll;
          std::unique_ptr<RooRealVar> alphaVar;
          std::unique_ptr<RooProduct> correctedNll;
          RooAbsReal* finalNll = nullptr;
      };

    
    class BruMcmc  : public Minimiser {
      
    public:

      BruMcmc(Int_t Niter=100,Int_t Nburn=10, Float_t norm=0.1): fNumIters(Niter),fNumBurnInSteps(Nburn),fNorm(norm){
        SetNameTitle("HSBruMcmc","BruMcmc minimiser");
      }
      ~BruMcmc() override = default;

      void Run(Setup &setup,RooAbsData &fitdata) override;

      void FitTo();
      
      file_uptr SaveInfo() override;
      void AddFormulaToMCMCTree();
      
      Bool_t MakeChain();
      TMatrixDSym MakeMinuitCovarianceMatrix();
      TMatrixDSym MakeMcmcCovarianceMatrix(TTree* tree,size_t burnin, Bool_t decoupleYields = kFALSE);
      TTree* GetTree(){return fTreeMCMC;}
      Double_t SumWeights();
      Double_t SumWeights2();
      
      void NoWeightCorrection(){fCorrectForWeights=kFALSE;}
      
      void Result();
      Double_t NLL(){return fChain->NLL();}

      //MCMCCalculator
      void InitModel(); 
      
      void SetData(RooAbsData& data) { fData = &data; }
      void SetPdf(RooAbsPdf& pdf) { fPdf = &pdf; }
      void SetPriorPdf(RooAbsPdf& pdf) { fPriorPdf = &pdf; }
      void SetParameters(const RooArgSet& set) { fPOI.removeAll(); fPOI.add(set); }
      void SetChainParameters(const RooArgSet & set) { fChainParams.removeAll(); fChainParams.add(set); }
      void SetNuisanceParameters(const RooArgSet& set) {fNuisParams.removeAll(); fNuisParams.add(set);}
      void SetConditionalObservables(const RooArgSet& set) {fConditionalObs.removeAll(); fConditionalObs.add(set);}
      void SetGlobalObservables(const RooArgSet& set) {fGlobalObs.removeAll(); fGlobalObs.add(set);}
      void SetProposalFunction(RooStats::ProposalFunction& proposalFunction) { fPropFunc = &proposalFunction; }
      void SetNumIters(Int_t numIters) { fNumIters = numIters; }
      void SetNumBurnInSteps(Int_t numBurnInSteps) { fNumBurnInSteps = numBurnInSteps; }

      void SetupBasicUsage();
      void SetKeepStart(Bool_t keep=kTRUE){fKeepStart=keep;}
      void SetTuneCovariance(Bool_t keep=kTRUE){fTuneCovStep=keep;}

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
 
      void SetTag(const TString& tag){fFileTag=tag;}
      const TString& GetTag()const {return fFileTag;}
      void SaveStepInfo();

      void SetCyclicParameters(const RooArgList& cyclics) { fCyclicPars.removeAll(); fCyclicPars.add(cyclics); }


      
   protected :
     // --- NEW MODULAR METHODS ---
      bool RunHastings(RooAbsReal* nll);
      void SaveChainToTree();
      // ---------------------------

      void AddEntryBranch();
      void CleanMakeChain(){};
      
      // --- UPGRADED TO SMART POINTERS ---
      std::unique_ptr<RooStats::MarkovChain> fChain; //!
      std::unique_ptr<RooDataSet> fChainData;        //!
      std::unique_ptr<RooArgSet> fParams;            //!
      // ----------------------------------

  
      void BuildBatchedNLLs(int numBatches);
      void ClearBatches();
      void SetStochasticSwapping(int freq) { fBatchSwapFreq = freq; }

      std::vector<std::unique_ptr<RooAbsData>> fBatchedData;
      std::vector<std::unique_ptr<NLLCache>> fBatchedNLLCache;
      std::vector<RooAbsReal*> fBatchedNLLPointers;
      int fBatchSwapFreq = 0;
      // ----------------------------

      // Modify the signature to take a data pointer!
      RooAbsReal* BuildNLL(RooAbsData* data, 
                           std::unique_ptr<RooAbsPdf>& localProdPdf,
                           std::unique_ptr<RooAbsReal>& baseNll,
                           std::unique_ptr<RooRealVar>& alphaVar,
                           std::unique_ptr<RooProduct>& correctedNll);

      
      TTree* fTreeMCMC=nullptr;//! ROOT manages this, DO NOT use unique_ptr
      Bool_t fCorrectForWeights=kTRUE;
      
      std::shared_ptr<TFile> fTempFile;//!
      file_uptr fOutFile;//!
      TString fFileTag;
      
      Bool_t fKeepStart=kFALSE; //randomise starting values
      Bool_t fMCMCHelp=kFALSE;//automate acceptance etc.
      Bool_t fTuneCovStep=kTRUE;
      
      RooArgSet   fPOI;        //! parameters of interest for interval
      RooArgSet   fNuisParams; //! nuisance parameters for interval (not really used)
      RooArgSet   fChainParams; //! parameters to store in the chain (if not specified they are all of them )
      RooArgSet   fConditionalObs; //! conditional observables
      RooArgSet   fGlobalObs;     //! global observables
      RooStats::ProposalFunction* fPropFunc{}; //! Proposal function for MCMC integration
      RooAbsPdf * fPdf=nullptr;        //! pointer to common PDF (owned by the workspace)
      RooAbsPdf * fPriorPdf=nullptr;   //! pointer to prior  PDF (owned by the workspace)
      Int_t fNumIters; // number of iterations to run metropolis algorithm
      Int_t fNumBurnInSteps; // number of iterations to discard as burn-in, starting from the first

      Int_t fNumBins{}; // set the number of bins to create for each
      Int_t fWarmup{}; //ignore these events
      Float_t fNorm=1;
      Int_t fNumBurnInStepsCov; //Number of steps to remove from chain to make covariance matrix for proposal function

      std::vector<Double_t> _formVals;
      std::vector<TBranch*> _formBranches;

      Double_t fChainAcceptance=0;//!
      Double_t fMinAcc=0.15;
      Double_t fMaxAcc=0.3;
      Double_t fTargetAcc=0.234;
      Int_t  fUncorrelateYields=0;

      RooArgList fCyclicPars;
      
      ClassDefOverride(HS::FIT::BruMcmc,1);
      
     };

    class BruMcmcSeq  : public BruMcmc {
      
    public:

      BruMcmcSeq(Int_t Niter=100,Int_t Nburn=10, Float_t norm=0.1):
        BruMcmc(Niter,Nburn,norm){
        SetNameTitle("BruMcmcSeq","BruMcmcSeq minimiser");
      }

      void Run(Setup &setup,RooAbsData &fitdata) override;

      ClassDefOverride(HS::FIT::BruMcmcSeq,1);
   };

    class BruMcmcSeqHelper  : public BruMcmc {
      
    public:
      
      BruMcmcSeqHelper(Int_t Niter=100,Int_t Nburn=10, Float_t norm=0.1,float target=0.234,float accmin=0.15,float accmax=0.35):BruMcmc(Niter,Nburn,norm),
        _proposal{norm,target,accmin,accmax}{
        SetNameTitle("HSBruMcmcSeqHelper","BruMcmcSeqHelper minimiser");
      }
  
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
        _propCov{norm,target,accmin,accmax}
     {
       SetNameTitle("BruMcmcCovariance","BruMcmcCovariance minimiser");
     }
     BruMcmcCovariance(std::vector<Int_t> Niters,Int_t Nburn=10, Float_t norm=0.01,float target=0.234,float accmin=0.16,float accmax=0.3):
        BruMcmc(Niters[0],Nburn,norm),
        _propSeq{norm,target,accmin,accmax},
        _propCov{norm,target,accmin,accmax},
        _fNumItersVec{Niters}
     {
       SetNameTitle("BruMcmcCovariance","BruMcmcCovariance minimiser");
     }
  
      void Run(Setup &setup,RooAbsData &fitdata) override;

     void TurnOffSequential(){_doSeq=kFALSE;}
     void TurnOffNDStep(){_doND=kFALSE;}
     void TurnOffCovariance(){_doCov=kFALSE;}
     void TuneCovarianceStep(){_tuneCovStep=kTRUE;}

     void ChangeNIter(){
       if(_iNIter==_fNumItersVec.size()) return;
        SetNumIters(_fNumItersVec[_iNIter]);
       std::cout<<" change nuiters "<<_fNumItersVec[_iNIter]<<" "<<_iNIter<<" "<<_fNumItersVec.size()<<std::endl;
       _iNIter++;
    }
     void ResetNIter() { _iNIter = 0; }
     
    private:

     BruSequentialProposal _propSeq;
     BruCovarianceProposal _propCov;

     std::vector<Int_t>_fNumItersVec;
     UInt_t _iNIter=0;
     
     Bool_t _doSeq=kTRUE;
     Bool_t _doND=kTRUE;
     Bool_t _doCov=kTRUE;
     Bool_t _tuneCovStep=kFALSE;
     
      ClassDefOverride(HS::FIT::BruMcmcCovariance,1);
   };

  }//namespace FIT
}//namespace HS

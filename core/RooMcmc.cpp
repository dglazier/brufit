#include "RooMcmc.h"
#include "HSMetropolisHastings.h"
#include <TIterator.h>
#include <TLeaf.h>
#include <TTreeIndex.h>
#include <RooStats/UniformProposal.h>
#include <RooStats/SequentialProposal.h>
#include <RooStats/ProposalHelper.h>
#include <TRobustEstimator.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>

namespace HS{
  namespace FIT{

    using namespace RooFit;
    using namespace RooStats;

    
    RooMcmc::~RooMcmc(){
      if(!_formBranches.empty()){
	for(auto br:_formBranches)
	  delete br;
      }
    }
    
    void RooMcmc::Run(Setup &setup,RooAbsData &fitdata){
      //initialise MCMCCalculator
      fSetup=&setup;
      fData=&fitdata;
      cout<<"RooMcmc::Run"<<endl;
      SetModel(fSetup->GetModelConfig());
      SetupBasicUsage();

      MakeChain();
      
    }
  
    ///////////////////////////////////////////
    void RooMcmc::MakeChain()
    {
      cout<<"HSMCMC::MakeChain() "<<fData<<" "<<fPdf<<" "<<fPOI.getSize()<<endl;
      if (!fData || !fPdf   ) return;
      if (fPOI.getSize() == 0) return;
     std::cout<<"proceed"<<endl;
 
   
      // if a proposal function has not been specified create a default one
      bool useDefaultPropFunc = (fPropFunc == nullptr);
      bool usePriorPdf = (fPriorPdf != nullptr);
      if (useDefaultPropFunc) fPropFunc = new UniformProposal();

      // if prior is given create product
      RooAbsPdf * prodPdf = fPdf;
      if (usePriorPdf) {
	TString prodName = TString("product_") + TString(fPdf->GetName()) + TString("_") + TString(fPriorPdf->GetName() );
	prodPdf = new RooProdPdf(prodName,prodName,RooArgList(*fPdf,*fPriorPdf) );
      }
 
      RooArgSet* constrainedParams = prodPdf->getParameters(*fData);
      RooAbsReal* nll = prodPdf->createNLL(*fData, Constrain(*constrainedParams),ConditionalObservables(fConditionalObs));
      //RooAbsReal* nll = prodPdf->createNLL(*fData); 
      delete constrainedParams;

      nll->constOptimizeTestStatistic(RooAbsArg::Activate,false) ;
      
      // add in sumw/sumw2 term
      if(fData->isNonPoissonWeighted()&&fCorrectForWeights){
      	Double_t SumW=SumWeights();
      	Double_t SumW2=SumWeights2();
	TString NllName=nll->GetName();
	NllName.ReplaceAll("-","m");
	NllName.ReplaceAll("+","p");
	nll->SetName(NllName);
      	RooFormulaVar *alphanll=new RooFormulaVar("alphanll",Form("%lf*%s",SumW/SumW2,nll->GetName()),RooArgSet(*nll));
      	nll=alphanll;
      }
    
      fParams = nll->getParameters(*fData);
      RemoveConstantParameters(fParams);
      std::cout<<"metropolis"<<endl;
      HSMetropolisHastings mh;
      if(fKeepStart) mh.SetKeepStart();
      mh.SetFunction(*nll);
      mh.SetType(MetropolisHastings::kLog);
      mh.SetSign(MetropolisHastings::kNegative);
      mh.SetParameters(*fParams);
      if (fChainParams.getSize() > 0) mh.SetChainParameters(fChainParams);
      mh.SetProposalFunction(*fPropFunc);
      mh.SetNumIters(fNumIters);

      if(fChain){ delete fChain; fChain=nullptr;}
      
      fChain= mh.ConstructChain(); //mh is still owner and will delete
      
      if(fChainData){ delete fChainData; fChainData=nullptr;}
      fChainData=fChain->GetAsDataSet(EventRange(0, fChain->Size()));

      if(fChainData){
	if(fTreeMCMC){ delete fTreeMCMC; fTreeMCMC=nullptr;}
	
 	fTreeMCMC=RooStats::GetAsTTree("MCMCTree","MCMCTree",*fChainData);
	delete fChainData; fChainData=nullptr;
      }  
      if(fChain->Size()>fNumBurnInSteps)
	fChainData=fChain->GetAsDataSet(EventRange(fNumBurnInSteps, fChain->Size()));

 
      nll->constOptimizeTestStatistic(RooAbsArg::DeActivate,false) ;
      if (useDefaultPropFunc) delete fPropFunc;
      if (usePriorPdf) delete prodPdf;
      delete nll;
      
      return;
    }
    ////////////////////////////////////////////////////////

    TMatrixDSym RooMcmc::MakeMcmcCovarianceMatrix(TTree* tree,size_t burnin){
          
      auto pars = fSetup->NonConstParsAndYields();
      Int_t Npars = pars.size();
      Int_t Nentries = tree->GetEntries()-burnin;
      Int_t param_index=0;     
      vector<Double_t> params(Npars);
      Double_t data[Npars];
      //Int_t NburnC = fNumBurnInStepsCov;
  
     
      //Loop over parameters of the model and set values from the tree
      //Needed for RobustEstimator
      //Only needed once
      int pindex=0;
      for(RooAbsArg* ipar : pars)
	{
	  if(ipar->isConstant()) continue;
	  if(tree->SetBranchAddress(ipar->GetName(), &params[pindex])==0){
	    pindex++;
	  }
	}
      Npars=pindex; //should be = number of branches in tree, protects for constant pars
      cout<<"Robust "<<Nentries<<" "<<Npars<<" "<<tree->GetEntries()<<" "<<burnin<<endl;
      //Create instance of TRobustEstimator
      TRobustEstimator r(Nentries,Npars);

      //Loop over entries of the tree to 'AddRow' of data to RobustEstimator
      // for (int ientry = 0; ientry<Nentries; ientry++)
      for (int ientry = burnin; ientry<Nentries+burnin; ientry++)
	{//Loop over entries of the tree
	  tree->GetEntry(ientry);
	 
	  for (int param_index = 0; param_index<Npars; param_index++)
	    { //Loop over parameters of the model
	      //And set 'data' element
	      data[param_index]=params[param_index];
	    }
	 
	  r.AddRow(data);//Appends data to RE
	}
     
      r.Evaluate(); //Necessary to calculate RE properly
      const TMatrixDSym* covMatSym=nullptr;
      covMatSym = r.GetCovariance();
      covMatSym->Print();
      //covMatSym is the symmetric covariance matrix to be used in the proposal function
     
      TMatrixDSym covMatSymNorm=*covMatSym;
     
      tree->ResetBranchAddresses();
      return covMatSymNorm;
    }

    /////////////////////////////////////////////////////////

    void RooMcmc::Result(){
      //Add entry branch to mcmc tree for easy cutting on BurnIn
      //fMCMCtree contains all events
      Long64_t entry=0;
      auto entryBranch=fTreeMCMC->Branch("entry",&entry,"entry/L");
      for(entry=0;entry<fTreeMCMC->GetEntries();entry++)
	entryBranch->Fill();
      //Add any formulas
      //Need to get a copy of variables first or setting
      //the means as parameter values does not seem to work...
      RooArgList saveFloatFinalList(*fChainData->get()) ;

      AddFormulaToMCMCTree();
  
      //set paramters to mean values of post burn in distributions
      //     RooArgList saveFloatFinalList(*fChainData->get()) ;
      for(Int_t i=0;i<fParams->getSize();i++){

	auto* var=dynamic_cast<RooRealVar*>(saveFloatFinalList.at(i));
      	cout<<var->GetName()<<" "<<var->getVal()<<" "<<fChainData->mean(*var)<<" +- "<<fChainData->sigma(*var)<<endl;
	auto var2=dynamic_cast<RooRealVar*>(fParams->find(var->GetName()));
	var2->setVal(fChainData->mean(*var));
	var2->setError(fChainData->sigma(*var));
      }
       
      //  fChainData->covarianceMatrix()->Print(); //crashin 6.20
 
 
      //look for the best likelihood

      //It is not recommended to use the best likelihood
      //code is lef here as an example
      // fTreeMCMC->BuildIndex(TString("1E6*nll_MarkovChain_local_"));
      // TTreeIndex *tindex = (TTreeIndex*)fTreeMCMC->GetTreeIndex();

      // auto index=tindex->GetIndex();
      // auto values=tindex->GetIndexValues();
      //  for(Int_t i=0;i<fParams->getSize();i++){
      //   RooArgList saveFloatMaxLikeList(*fChainData->get(index[0])) ;

      // 	RooRealVar* var=dynamic_cast<RooRealVar*>(saveFloatMaxLikeList.at(i));
      // 	Double_t val=var->getVal();
      // 	auto var2=dynamic_cast<RooRealVar*>(fParams->find(var->GetName()));
      // 	var2->setVal(val);
      // }
      
     
    }
    void RooMcmc::AddFormulaToMCMCTree(){
      //avoid future memory leaks/crashes...
      fTreeMCMC->ResetBranchAddresses();


      auto formulas=fSetup->ParameterFormulas(); //formulas that just depend on parameters, not variables/observables
      if(!formulas.getSize()) return;

      _formVals.reserve(formulas.getSize());
      _formBranches.reserve(formulas.getSize());

      TIter iter=formulas.createIterator();
      Int_t iform=0;

      //getLeaves before extra branches
      auto parLeaves=fTreeMCMC->GetListOfLeaves();
      
      while(auto* formu=dynamic_cast<RooFormulaVar*>(iter())){
	TString formuName=formu->GetName();
	_formVals[iform]=0;
	_formBranches[iform]=nullptr;
	_formBranches[iform]=fTreeMCMC->Branch(formuName,&_formVals[iform],formuName+"/D");
	iform++;
      }

      Long64_t Nmcmc=fTreeMCMC->GetEntries();
      Int_t Nleaf=parLeaves->GetEntries();
 
      for(Int_t entry=0;entry<Nmcmc;entry++){
	
	fTreeMCMC->GetEntry(entry);
 
	//Set value of parameters to value in tree for this event
	for(Int_t ibr=0;ibr<Nleaf;ibr++){
	  auto *leaf=dynamic_cast<TLeaf*>(parLeaves->At(ibr));	
	  auto* brVar=dynamic_cast<RooRealVar*>(fParams->find(leaf->GetName()));
	  if(brVar!=nullptr) brVar->setVal(leaf->GetValue());
	    
	}
	//now calculate value of formula for these parameters
	iter.Reset();
	iform=0;
	while(auto* formu=dynamic_cast<RooFormulaVar*>(iter())){
	  
	  _formVals[iform]=formu->getValV();
	  _formBranches[iform]->Fill();
	  iform++;

	}
      }  
    }
    ///////////////////////////////////////////////
    Double_t  RooMcmc::SumWeights(){
      // Otherwise sum the weights in the event
      Double_t sumw(0), carry(0);
      Int_t i ;
      for (i=0 ; i<fData->numEntries() ; i++) {
	fData->get(i) ;
 
	Double_t y = fData->weight() - carry;
	Double_t t = sumw + y;
	carry = (t - sumw) - y;
	sumw = t;
      }
      return sumw;
    }
    ///////////////////////////////////////////////////////
    Double_t  RooMcmc::SumWeights2(){
      // Otherwise sum the weights in the event
      Double_t sumw(0), carry(0);
      Int_t i ;
      for (i=0 ; i<fData->numEntries() ; i++) {
	fData->get(i) ;
 
	Double_t y = fData->weight()*fData->weight() - carry;
	Double_t t = sumw + y;
	carry = (t - sumw) - y;
	sumw = t;
      }
      return sumw;
    }
    //FRom MCMCCalculator
    void RooMcmc::SetModel( ModelConfig*  model) {
      cout<<"RooMcmc::SetModel"<<endl;
      // set the model
      fModelConfig=model;
      //fPdf = fModelConfig->GetPdf();
      fPdf = fSetup->Model();
      fPriorPdf = fModelConfig->GetPriorPdf();
      fPOI.removeAll();
      fNuisParams.removeAll();
      fConditionalObs.removeAll();
      fGlobalObs.removeAll();
      if (fModelConfig->GetParametersOfInterest())
	fPOI.add(*fModelConfig->GetParametersOfInterest());
      if (fModelConfig->GetNuisanceParameters())
	fNuisParams.add(*fModelConfig->GetNuisanceParameters());
      if (fModelConfig->GetConditionalObservables())
	fConditionalObs.add( *(fModelConfig->GetConditionalObservables() ) );
      if (fModelConfig->GetGlobalObservables())
	fGlobalObs.add( *(fModelConfig->GetGlobalObservables() ) );
      
    }
    /////////////////////////////////////////////////////
    void RooMcmc::SetupBasicUsage()
    {
      fPropFunc = nullptr;
      // fNumIters = 100;
      //fNumBurnInSteps = 10;
      //fWarmup=fNumBurnInSteps;
      
     }
    ///////////////////////////////////////////////////////////////
    file_uptr RooMcmc::SaveInfo(){
      
      TString fileName=fSetup->GetOutDir()+fSetup->GetName()+"/Results"+fSetup->GetTitle()+GetName()+".root";
      file_uptr file(TFile::Open(fileName,"recreate"));
      Result();
      
      fTreeMCMC->Write();
      //save paramters and chi2s in  dataset (for easy merging)
      //RooArgSet saveArgs(*fParams);
      RooArgSet saveArgs(fSetup->Parameters());
      saveArgs.add(fSetup->Yields());
      
      RooRealVar Nllval("NLL","NLL",NLL());
      saveArgs.add(Nllval);
     
      RooDataSet saveDS(FinalParName(),TString(GetName())+"Results",saveArgs);
      saveDS.add(saveArgs);
      saveDS.Write();
      TTree* treeDS=RooStats::GetAsTTree(ResultTreeName(),ResultTreeName(),saveDS);
      treeDS->Write();
      delete treeDS;treeDS=nullptr;

      return std::move(file);
    }
     //////////////////////////////////////////////////////////////

    
   void RooMcmcSeq::Run(Setup &setup,RooAbsData &fitdata){
     fSetup=&setup;
    fData=&fitdata;
    //initialise MCMCCalculator
    SetData(fitdata);
    SetModel(setup.GetModelConfig());
    SetupBasicUsage();
     
    RooStats::SequentialProposal sp(fNorm);
    SetProposalFunction(sp);
    fKeepStart=kTRUE; //start values from previous
    MakeChain();
    
   }

    //////////////////////////////////////////////////

    void RooMcmcSeqCov::Run(Setup &setup, RooAbsData &fitdata){
      
      fSetup=&setup;
      fData=&fitdata;
      //initialise MCMCCalculator
      SetData(fitdata);
      SetModel(setup.GetModelConfig());
      SetupBasicUsage();

      /*
      // Want to get the covariance matrix from another run of the mcmc (RooMcmcSeq) and use this to generate a proposal function in RooMcmcSeqCov.
      See below for how to run with covMatrix from the same run (RooMcmcSeqThenCov)
      */

//Check if ResultsHSRooMcmcSeq.root file exists
      if(gSystem->AccessPathName(fSetup->GetOutDir() + "/ResultsHSRooMcmcSeq.root"))
	{std::cout<<"\n\n ResultsHSRooMcmcSeq.root file not found. \n\n Are you sure you ran RooMcmcSeq? \n"<<std::endl;
	  exit(-1);}

      auto saveDir = gDirectory;
      TFile *file = TFile::Open(fSetup->GetOutDir() + "/ResultsHSRooMcmcSeq.root");
      saveDir->cd();
      auto tree=dynamic_cast<TTree*>(file->Get("MCMCTree"));
 
      
      TMatrixDSym covMat =  MakeMcmcCovarianceMatrix(tree,fNumBurnInStepsForCov);
      delete file; //Must be at end of func

      //scale covariance matrix by Norm
      
      auto divideNorm = 1./fNorm;
      covMat*= divideNorm;
 
      ProposalHelper ph;
      ph.SetVariables(fSetup->ParsAndYields());
      ph.SetUpdateProposalParameters(true); // auto-create mean vars and add mappings
      ph.SetCacheSize(1);
      ph.SetCovMatrix(covMat);
      ProposalFunction* pf = ph.GetProposalFunction();
      SetProposalFunction(*pf);
     	
      fKeepStart=kTRUE; //start values from previous
      MakeChain();


    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 void RooMcmcSeqThenCov::Run(Setup &setup, RooAbsData &fitdata){
      
      fSetup=&setup;
      fData=&fitdata;
      //initialise MCMCCalculator
      SetData(fitdata);
      SetModel(setup.GetModelConfig());
      SetupBasicUsage();

      /*
      //A class that 
1.Runs through a number of burn in events
2.Runs a sequential proposal mcmc
3.Finds the covariance matrix from the seq prop run
4.Uses the cov mat to generate new prop func and run
      */

      RooStats::SequentialProposal sp(fNorm);
      SetProposalFunction(sp);
      fKeepStart=kTRUE; //start values from previous
      MakeChain();

      auto saveN = fNumIters;
      fNumIters = fNumItersThenCov;
      auto saveNorm=fNorm;
      fNorm=fNormThenCov;
      auto saveBurn=fNumBurnInSteps;
      fNumBurnInSteps=fNumBurnInStepsThenCov;

      //      TMatrixDSym covMat =  MakeMcmcCovarianceMatrix(fTreeMCMC,fNumBurnInSteps);
      TMatrixDSym covMat =  MakeMcmcCovarianceMatrix(fTreeMCMC,saveBurn);
      auto divideNorm = 1./fNorm;
      covMat*= divideNorm;
 
        
      ProposalHelper ph;
      ph.SetVariables(fSetup->NonConstParsAndYields());
      ph.SetUpdateProposalParameters(true); // auto-create mean vars and add mappings


      //nite we reset the cache every proposal as was
      //causing the chain to get stuck
      ph.SetCacheSize(1);
      ph.SetCovMatrix(covMat);
      ProposalFunction* pf = ph.GetProposalFunction();
      SetProposalFunction(*pf);
     	
      fKeepStart=kTRUE; //start values from previous
      MakeChain();
      
      fNumIters = saveN; //switch back to seq Niters for next bin
      fNorm=saveNorm;
      fNumBurnInSteps=saveBurn;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void RooMcmcUniform2Seq::Run(Setup &setup,RooAbsData &fitdata){
     fSetup=&setup;
     fData=&fitdata;
     //initialise MCMCCalculator
     SetData(fitdata);
     SetModel(fSetup->GetModelConfig());
     SetupBasicUsage();
      
     auto NumIters0=fNumIters;
     fNumIters=2*fPOI.getSize();

     MakeChain();//default uniform for burn-in

     fNumIters=NumIters0;
     fWarmup=5;

     RooStats::SequentialProposal sp(fNorm);
     SetProposalFunction(sp);
     fKeepStart=kTRUE; //start values from previous
     MakeChain();
     
 
   }

    void RooMcmcGaus::Run(Setup &setup,RooAbsData &fitdata){
     fSetup=&setup;
     fData=&fitdata;

     //initialise MCMCCalculator
     SetData(fitdata);
     SetModel(setup.GetModelConfig());
     SetupBasicUsage();
     
     ProposalHelper ph;
     ph.SetVariables(fPOI);

     auto npars=fPOI.getSize();
    
     ph.SetWidthRangeDivisor(fNorm);
    
     ph.SetUpdateProposalParameters(true);
     ph.SetCacheSize(1);
     fPropFunc = (ph.GetProposalFunction());
     
     
     MakeChain();
     
      
   }

 
    
  }
}


    

#include "BruMcmc.h"
#include "BruMetropolisHastings.h"

#include <TROOT.h>
#include <TIterator.h>
#include <TDirectory.h>
#include <TLeaf.h>
#include <TTreeIndex.h>
#include <RooStats/UniformProposal.h>
#include <RooStats/SequentialProposal.h>
#include <RooStats/ProposalHelper.h>
#include <TRobustEstimator.h>
#include <TPrincipal.h>
#include <TMatrixD.h>
#include <TH1D.h>
#include <TMatrixDSym.h>
#include <RooGlobalFunc.h>
namespace HS{
  namespace FIT{

    using namespace RooFit;
    using namespace RooStats;

    
    BruMcmc::~BruMcmc(){
      if(!_formBranches.empty()){
	for(auto br:_formBranches)
	  delete br;
      }
      if(fChain!=nullptr){delete fChain;fChain =nullptr;}
      if(fChainData!=nullptr) {delete fChainData;fChainData=nullptr;}
    }
    
    void BruMcmc::Run(Setup &setup,RooAbsData &fitdata){
      //initialise MCMCCalculator
      fSetup=&setup;
      fData=&fitdata;
      cout<<"BruMcmc::Run"<<endl;
      SetModel(fSetup->GetModelConfig());
      SetupBasicUsage();

      MakeChain();
      
    }
  
    ///////////////////////////////////////////
    Bool_t BruMcmc::MakeChain()
    {
      cout<<"HSMCMC::MakeChain() "<<fData<<" "<<fPdf<<" "<<fPOI.getSize()<<endl;
      if (!fData || !fPdf   ) return kFALSE;
      if (fPOI.getSize() == 0) return kFALSE;
     std::cout<<"proceed"<<endl;
 
     ProcInfo_t info;
     gSystem->GetProcInfo(&info);
     cout<<"1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
     cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
     cout<<"PRocInfo "<<0.000001*info.fMemResident<<" "<<0.000001*info.fMemVirtual<<endl;
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

      //   RooAbsReal* nll = prodPdf->createNLL(*fData, Constrain(*constrainedParams),ConditionalObservables(fConditionalObs));
      auto foptions = fSetup->FitOptions();

      //remove not applicable options      
      TObject* opt=nullptr;
      if((opt=foptions.find("Save"))!=nullptr){
	foptions.Remove(opt);
      }
       if((opt=foptions.find("SumW2Error"))!=nullptr){
	foptions.Remove(opt);
      }

       
      auto cmd = ConditionalObservables(fConditionalObs);
      foptions.Add(dynamic_cast<RooCmdArg*>(&cmd));
      cmd=Constrain(*constrainedParams);
      foptions.Add(dynamic_cast<RooCmdArg*>(&cmd));
 
      RooAbsReal* nll = prodPdf->createNLL(*fData,foptions);
      //RooAbsReal* nll = prodPdf->createNLL(*fData); 
      delete constrainedParams;

      
      nll->constOptimizeTestStatistic(RooAbsArg::Activate,false) ;
      
      // add in sumw/sumw2 term
      RooAbsReal* delnll=nullptr;
      if(fData->isNonPoissonWeighted()&&fCorrectForWeights){
      	Double_t SumW=SumWeights();
      	Double_t SumW2=SumWeights2();
	TString NllName=nll->GetName();
	NllName.ReplaceAll("-","m");
	NllName.ReplaceAll("+","p");
	nll->SetName(NllName);
      	RooFormulaVar *alphanll=new RooFormulaVar("alphanll",Form("%lf*%s",SumW/SumW2,nll->GetName()),RooArgSet(*nll));
	delnll=nll;//keep a pointer for deleting
      	nll=alphanll;
      }
 
      fParams = nll->getParameters(*fData);
      RemoveConstantParameters(fParams);

      //if using smapling integral need to add params ?
      //probably do not need this now, as just ampling in EVentsPDF class
      //auto* sampPDF= dynamic_cast<RooAbsPdf*>(fSetup->Constraints().at(0));    
      //if(sampPDF)fParams->add(*(sampPDF->getParameters(*fData)));
      //RemoveConstantParameters(fParams);
    gSystem->GetProcInfo(&info);
     cout<<"2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
     cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
     cout<<"PRocInfo "<<0.000001*info.fMemResident<<" "<<0.000001*info.fMemVirtual<<endl;
       fParams->Print("v");
      BruMetropolisHastings mh;
      // if(fKeepStart) mh.SetKeepStart();
      // if(fMCMCHelp){
      // 	mh.Help();
      // 	mh.SetAcceptanceRange(fMinAcc,fMaxAcc);
      // 	mh.SetTargetAccept(fTargetAcc);
      // 	mh.SetNorm(fNorm);
      // }
      
      mh.SetFunction(*nll);
      // mh.SetType(MetropolisHastings::kLog);
      // mh.SetSign(MetropolisHastings::kNegative);
      mh.SetParameters(*fParams);
      if (fChainParams.getSize() > 0) mh.SetChainParameters(fChainParams);
      mh.SetProposalFunction(*fPropFunc);
      mh.SetNumIters(fNumIters);


      //cout<<"number of constraints "<<sampPDF->GetName()<<endl;
      //try balancing integral sampling
       // mh.SetBalance(sampPDF);
      
      if(fChain){ delete fChain; fChain=nullptr;}
      gSystem->GetProcInfo(&info);
     cout<<"3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
     cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
     cout<<"PRocInfo "<<0.000001*info.fMemResident<<" "<<0.000001*info.fMemVirtual<<endl;
     
      fChain= mh.ConstructChain(); //mh is still owner and will delete
      cout<<"DEBUG "<<" Got chain "<<fChain<<endl;
      if(fChain==nullptr){
	if (useDefaultPropFunc) delete fPropFunc;
	if (usePriorPdf) delete prodPdf;
	if(fTreeMCMC!=nullptr){ delete fTreeMCMC; fTreeMCMC=nullptr;}

	fChainAcceptance=mh.GetAcceptance();
	
	delete nll;
	if(delnll) delete delnll;

	return kFALSE; //unsuccessful
      }
     gSystem->GetProcInfo(&info);
     cout<<"4~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
     cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
     cout<<"PRocInfo "<<0.000001*info.fMemResident<<" "<<0.000001*info.fMemVirtual<<endl;
       cout<<"DEBUG "<<" Got chain data 1 "<<fChainData<<endl;
     
      if(fChainData!=nullptr){ delete fChainData; fChainData=nullptr;}
      cout<<"DEBUG "<<" Got chain data 2 "<<fChainData<<" "<<fChain->Size()<<endl;
     
      fChainData=fChain->GetAsDataSet(EventRange(0, fChain->Size()));

      cout<<"DEBUG "<<" Got chain data 3 "<<fChainData<<" "<<fTreeMCMC<<endl;
     if(fChainData!=nullptr){
	if(fTreeMCMC!=nullptr){ delete fTreeMCMC; fTreeMCMC=nullptr;}
	cout<<"Get tree from chains "<<endl;//(*gDirectory).GetName()<<endl;
	auto saveDir=gDirectory;
	// if((*gDirectory).IsWritable()==false){
	//   TString saveName=fSetup->GetOutDir()+fSetup->GetName()+"/MCMCTemp.root";
	//   fTempFile=std::make_shared<TFile>(saveName,"recreate");
	// }
	//	fChainData->Print("v");
	fOutFile->cd();
	fTreeMCMC=RooStats::GetAsTTree("MCMCTree","MCMCTree",*fChainData);
	saveDir->cd();
 	delete fChainData; fChainData=nullptr;
      }
     cout<<"DEBUG "<<" Got chain size  "<< fChain->Size() <<" burnin "<< fNumBurnInSteps<<endl;

     if(fChain->Size()>fNumBurnInSteps){
       fChainData=fChain->GetAsDataSet(EventRange(fNumBurnInSteps, fChain->Size()));
     }
     
 
      nll->constOptimizeTestStatistic(RooAbsArg::DeActivate,false) ;

      CleanMakeChain();
      if (useDefaultPropFunc) delete fPropFunc;
      if (usePriorPdf) delete prodPdf;
      delete nll;
      if(delnll) delete delnll;
      
     gSystem->GetProcInfo(&info);
     cout<<"5~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
     cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
     cout<<"PRocInfo "<<0.000001*info.fMemResident<<" "<<0.000001*info.fMemVirtual<<endl;
     
  
      return kTRUE;
    }
    void BruMcmc::CleanMakeChain(){
    }
    ////////////////////////////////////////////////////////
    
    TMatrixDSym BruMcmc::MakeMcmcCovarianceMatrix(TTree* tree,size_t burnin){
          
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
      covMatSym = r.GetCovariance(); //actually returns pointer to reference so do not delete
      covMatSym->Print();
      //covMatSym is the symmetric covariance matrix to be used in the proposal function
     
      TMatrixDSym covMatSymNorm=*covMatSym;
     
      
      tree->ResetBranchAddresses();

      // TString saveName=fSetup->GetOutDir()+fSetup->GetName()+"/MCMCSeq.root";
      
      // TFile* saveSeq=new TFile(saveName,"recreate");

      // cout<<"Set tree file "<<tree->GetDirectory()<<endl;
      // // auto saveDir=tree->GetDirectory();
      // //tree->SetDirectory(saveSeq);
      // AddEntryBranch();
      // //tree->Write();
      // saveSeq->WriteObject(tree,tree->GetName());
      // // tree->SetDirectory(saveDir);
      // delete saveSeq;
      //      tree->SetDirectory(nullptr);
      
      return covMatSymNorm;
    }
    ////////////////////////////////////////////////////////
    
    TMatrixDSym BruMcmc::MakeMcmcNonYieldCovarianceMatrix(TTree* tree,size_t burnin){
          
      //auto pars = fSetup->NonConstParsAndYields();
      auto pars = fSetup->Parameters();
      Int_t Npars = pars.size();
      auto yields = fSetup->Yields();
      Int_t Nyields = yields.size();


      Int_t Nentries = tree->GetEntries()-burnin;
      vector<Double_t> params(Npars);
      vector<Double_t> byields(Nyields);
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
      int yindex=0;
      std::vector<TH1D> yield_hists;
      for(RooAbsArg* iyield : yields)
	{
	  if(iyield->isConstant()) continue;
	  
	  auto realvar=dynamic_cast<RooRealVar*>(iyield);

	  if(tree->SetBranchAddress(iyield->GetName(), &byields[yindex])==0){
	    yield_hists.push_back(std::move(TH1D(iyield->GetName(),iyield->GetName(),1000,realvar->getMax(),realvar->getMin())));
	    yindex++;
	  }
	}
      Nyields=yindex; //should be = number of branches in tree, protects for constant pars
      cout<<" MakeMcmcNonYieldCovarianceMatrix "<<Nentries<<" "<<Npars<<" "<<tree->GetEntries()<<" "<<burnin<<endl;
      //Create instance of TRobustEstimator
      TRobustEstimator r(Nentries,Npars+Nyields);
      vector<Double_t> data(Npars+Nyields);
     

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
	 
	  r.AddRow(data.data());//Appends data to RE
	  
	  for (int yield_index = 0; yield_index<Nyields; yield_index++)
	    { //Loop over yields of the model
	      //And fill histogram to get rms
	      yield_hists[yield_index].Fill(byields[yield_index]);
	      data[Npars+yield_index]=0;
	    }


	}

      
     
      r.Evaluate(); //Necessary to calculate RE properly
      const TMatrixDSym* covMatSym=nullptr;
      covMatSym = r.GetCovariance(); //actually returns pointer to reference so do not delete
       //covMatSym is the symmetric covariance matrix to be used in the proposal function
     
      TMatrixDSym covMatSymNorm=*covMatSym;

      //add final rows with yeilds uncorrelated
      for(int iy=0;iy<Nyields;++iy){
	//zero all but last elements of row
	for(Int_t id=0;id<data.size();++id)
	  data[id]=0;
	//add new row for his yield

	cout<<"MakeMcmcNonYieldCovarianceMatrix add yeilds "<<yield_hists[iy].GetName()<<" covariance "<<yield_hists[iy].GetRMS()*yield_hists[iy].GetRMS()<<endl;
	//data.push_back(yield_hists[iy].GetRMS()*yield_hists[iy].GetRMS());
	//covMatSymNorm.AddRow(data.data());
	//	covMatSymNorm.ResizeTo(covMatSymNorm.GetNrows()+1,covMatSymNorm.GetNcols()+1);
	//	covMatSymNorm.Use(covMatSymNorm.GetNrows()-(Nyields-iy),data.data());
	auto row=covMatSymNorm.GetNrows()-(Nyields-iy);
	cout<<"change "<<row<<" row from "<<covMatSymNorm(row,row)<<" to ";
	covMatSymNorm(row,row)=yield_hists[iy].GetRMS()*yield_hists[iy].GetRMS();
	cout<<covMatSymNorm(row,row)<<endl;
	
      }
      covMatSymNorm.Print();
 
      tree->ResetBranchAddresses();

      TString saveName=fSetup->GetOutDir()+fSetup->GetName()+"/MCMCSeq.root";
      
      TFile* saveSeq=new TFile(saveName,"recreate");
      AddEntryBranch();
      tree->Write();
      delete saveSeq;
      
      return covMatSymNorm;
    }
    TMatrixDSym BruMcmc::MakeMcmcPrincipalCovarianceMatrix(TTree* tree,size_t burnin){
          
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
      //TRobustEstimator r(Nentries,Npars);
      TPrincipal principal(Npars,"ND");
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
	  principal.AddRow(data);
	  //	  r.AddRow(data);//Appends data to RE
	}
      // Do the actual analysis
      //principal.MakePrincipals();
      //principal.Print();
 
 
      // r.Evaluate(); //Necessary to calculate RE properly
      const TMatrixD* covMat=nullptr;
      covMat = principal.GetCovarianceMatrix(); //actually returns pointer to reference so do not delete
      //covMat->Print();
      //covMatSym is the symmetric covariance matrix to be used in the proposal function
      // const TMatrixD* m = p.GetCovarianceMatrix();
     TMatrixD mt = *covMat; mt.T();
     TMatrixDDiag d(mt); d = 0;
     TMatrixD tempMat = *covMat+mt;
     //tempMat.Print();
     TMatrixDSym  covMatSym(0,Npars-1);
     covMatSym.SetMatrixArray(tempMat.GetMatrixArray());
     covMatSym.Print();
     tree->ResetBranchAddresses();

      TString saveName=fSetup->GetOutDir()+fSetup->GetName()+"/MCMCSeq.root";
       
      TFile* saveSeq=new TFile(saveName,"recreate");
      AddEntryBranch();
      tree->Write();
      delete saveSeq;
      
      return covMatSym;
    }

    /////////////////////////////////////////////////////////
    void BruMcmc::AddEntryBranch(){
    std::cout<<"BruMcmc::AddEntryBranch()"<<std::endl;
       //Add entry branch to mcmc tree for easy cutting on BurnIn
      //fMCMCtree contains all events
      Long64_t entry=0;
      auto entryBranch=fTreeMCMC->Branch("entry",&entry,"entry/L");
      for(entry=0;entry<fTreeMCMC->GetEntries();entry++)
	entryBranch->Fill();
      
      std::cout<<"BruMcmc::AddEntryBranch() Done"<<std::endl;
    }
    void BruMcmc::Result(){
      AddEntryBranch();
      //Add any formulas
      //Need to get a copy of variables first or setting
      //the means as parameter values does not seem to work...
      RooArgList saveFloatFinalList(*fChainData->get()) ;

      AddFormulaToMCMCTree();
  
      //set paramters to mean values of post burn in distributions
      //     RooArgList saveFloatFinalList(*fChainData->get()) ;
      for(Int_t i=0;i<fParams->getSize();i++){

	auto* var=dynamic_cast<RooRealVar*>(saveFloatFinalList.at(i));
      	cout<<var->GetName()<<" "<<fChainData->mean(*var)<<" +- "<<fChainData->sigma(*var)<<endl;
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
    void BruMcmc::AddFormulaToMCMCTree(){
      //avoid future memory leaks/crashes...
      fTreeMCMC->ResetBranchAddresses();

      std::cout<<"BruMcmc::AddFormulaToMCMCTree()"<<std::endl;
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
    std::cout<<"BruMcmc::AddFormulaToMCMCTree() done"<<std::endl;
     }
    ///////////////////////////////////////////////
    Double_t  BruMcmc::SumWeights(){
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
    Double_t  BruMcmc::SumWeights2(){
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
    void BruMcmc::SetModel( ModelConfig*  model) {
      cout<<"BruMcmc::SetModel"<<endl;
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
    void BruMcmc::SetupBasicUsage()
    {
      fPropFunc = nullptr;
      // fNumIters = 100;
      //fNumBurnInSteps = 10;
      //fWarmup=fNumBurnInSteps;
     
      TString fileName=fSetup->GetOutDir()+fSetup->GetName()+"/Results"+fSetup->GetTitle()+GetName()+".root";
      //TString fileName=fSetup->GetOutDir()+fSetup->GetName()+"/"+FileName();

      fOutFile.reset(TFile::Open(fileName,"recreate"));
       
     }
    ///////////////////////////////////////////////////////////////
    file_uptr BruMcmc::SaveInfo(){
      auto saveDir= gDirectory;
      fOutFile->cd();
      fTreeMCMC->SetDirectory(fOutFile.get());
      Result();
      fTreeMCMC->Write();
    
      delete fTreeMCMC; fTreeMCMC=nullptr;//or else crashes in destructor
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

      std::cout<<"BruMcmc::SaveInfo() Done to "<<fOutFile->GetName()<<std::endl;
      saveDir->cd();
      return std::move(fOutFile);
    }
     //////////////////////////////////////////////////////////////

    void BruMcmc::SetParVals(RooArgSet* toThesePars){
      for( auto &pory: fSetup->ParsAndYields()){
	if( dynamic_cast<RooRealVar*>(pory)){
	  dynamic_cast<RooRealVar*>(pory)->setVal(dynamic_cast<RooRealVar*>(toThesePars->find(pory->GetName()))->getVal());
	}
      }
    }
 
   void BruMcmcSeq::Run(Setup &setup,RooAbsData &fitdata){
     fSetup=&setup;
    fData=&fitdata;
    //initialise MCMCCalculator
    SetData(fitdata);
    SetModel(setup.GetModelConfig());
    SetupBasicUsage();
     
    RooStats::SequentialProposal sp(fNorm);
    SetProposalFunction(sp);
    MakeChain();
    
   }

   void BruMcmcSeqHelper::Run(Setup &setup,RooAbsData &fitdata){

     fSetup=&setup;
     fData=&fitdata;
     //initialise MCMCCalculator
     SetData(fitdata);
     SetModel(setup.GetModelConfig());
     SetupBasicUsage();

     SetProposalFunction(_proposal);
     MakeChain();

   }

  void BruMcmcCovariance::Run(Setup &setup,RooAbsData &fitdata){

    //auto foptions = fSetup->FitOptions();
    /* //Prefit with reduced stats, needs work....remove for now
    std::cout<<"BruMcmcCovariance::Run check for prefit "<<dynamic_cast<RooDataSet*>(&fitdata)<<std::endl;
    RooDataSet tiny("tiny", "tiny", *fitdata.get(),
		    fitdata.isWeighted() ? RooFit::WeightVar(dynamic_cast<RooDataSet*>(&fitdata)->weightVar()->GetName())  : RooCmdArg());

    std::cout<<"BruMcmcCovariance::Run check for prefit "<<std::endl;
    auto tinyFraction = 0.1;
    if(0){
      
      auto step =(Int_t) 1/tinyFraction; //for 10%
      for (int i=0; i<fitdata.numEntries(); i+=step)
	{
	  const RooArgSet *event = fitdata.get(i);
	  tiny.add(*event, fitdata.weight());
	}
      fData=&tiny;

      std::cout<<"BruMcmcCovariance::Run made tiny datset for prefit!!"<<std::endl;
      fData->Print("v");
    }
    else{
      fData=&fitdata;
    }
    */
    
    fData=&fitdata;
    fSetup=&setup;
    //initialise MCMCCalculator
    
    SetModel(setup.GetModelConfig());
    SetupBasicUsage();
    
    //find a region of hgh likelihood
    if(_doSeq==kTRUE){
      SetProposalFunction(_propSeq);
       MakeChain();
     }
     
     //now move in all parameters simultaneosuly to
     //give chain for covaiance matrix
     if(_doND==kTRUE){
       _propSeq.SetIsSequential(kFALSE);
       MakeChain();
     }

     //Now find accurate covariance matrix for final sampling
     if(fTreeMCMC!=nullptr){

       // if(0){
       // 	 auto& yields = fSetup->Yields();
       // 	 for(auto yld:static_range_cast<RooRealVar *>(yields)){
       // 	   yld->setVal(yld->getVal()/tinyFraction);
       // 	 }
       // 	 fData=&fitdata; //make sure have full dataset
       // }

       if(_doCov==kTRUE){
	 std::cout<<" BruMcmcCovariance::Run "<<_doCov<<std::endl;
	 std::unique_ptr<TMatrixDSym> covMat; 
	 covMat.reset(new TMatrixDSym(MakeMcmcCovarianceMatrix(fTreeMCMC,fNumBurnInSteps)));
	 _propCov.SetCovariance(*covMat.get(),fSetup->NonConstParsAndYields());
	 SetProposalFunction(_propCov);
	 MakeChain();
       }
       
     }
     
  }
  

   
  }
}


    

#include "ToyManager.h"
#include "RooHSEventsPDF.h"
#include "RooHSEventsHistPDF.h"
#include "RooComponentsPDF.h"
#include <TH1.h>
#include <TDirectory.h>
#include <RooRandom.h>
#include <RooStats/RooStatsUtils.h>

#include <memory>

#include <utility>

namespace HS{
  namespace FIT{

    Bool_t ToyManager::Run(){


      LoadResult();
      
      if(GetCurrToy()==0){
	InitSummary(); 
      }
      
      CreateCurrSetup();

      //Look for Special case of RooHSEventsPDFs
      FillEventsPDFs();
     
      //create extended max likelihood pdf
      fCurrSetup->TotalPDF();
      
      Generate();

      return kTRUE;
    }
    
    void  ToyManager::Generate(){
      cout<<"ToyManager::Generate()"<<endl;
      Long64_t geni=0;

      auto *model=fCurrSetup->Model();

      auto fitvars=fCurrSetup->FitVarsAndCats();
  
     RooRandom::randomGenerator()->SetSeed(0);//random seed
     //RooRandom::randomGenerator()->SetSeed(111);//random seed

      // Long64_t nexp=RooRandom::randomGenerator()->Poisson(model->expectedEvents(fitpars));
 
 
      fitvars.Print("v");
      auto iter=fitvars.iterator();
      const RooAbsArg *tmp=nullptr;
      while ((tmp = dynamic_cast<RooAbsArg*>(iter.Next()))){
	auto arg=fitvars.find(tmp->GetName());
	//cout<<"ToyManager::Generate() "<<tmp->GetName()<<" "<<model->isDirectGenSafe(*arg)<<endl;
      }
      while(fToyi<fNToys){//Note we do not parallelise toy generation, just run sequentially here
	cout<<"ToyManager::Generate() "<<fToyi<<" of "<<fNToys<<endl;
	//use number of events set,
	//or number =yields (e.g. from previous) fits if not.
	Long64_t NtoGen = fNEvents==-1 ? fCurrSetup->SumOfYields() : fNEvents;
	
	Long64_t nexp=RooRandom::randomGenerator()->Poisson(NtoGen);
	
	model->Print();
	fGenData=model->generate(fitvars,nexp);
	fGenData->SetName("ToyData");
	SaveResults();
	fToyi++;
      }
      fToyi=0;
      
    }
    void ToyManager::SaveResults(){
      if(!fGenData) return;
       TString fileName=Form("%s%s/%s.root",fCurrSetup->GetOutDir().Data(),GetCurrName().Data(),GetCurrTitle().Data());
      
      fToyFileNames.push_back(fileName);
      
      auto outfile=std::unique_ptr<TFile>(new TFile{fileName,"recreate"});
      //convert dataset to a tree for saving
      TTree* tree=RooStats::GetAsTTree("ToyData","ToyData",*fGenData);

      auto numEntries=fGenData->numEntries();
      auto cats=fCurrSetup->Cats();

      //add categories to tree! not done by GetAsTTree
      TIter iter=cats.createIterator();
      vector<TBranch*> branches(cats.getSize());
      vector<Int_t> branchVal(cats.getSize());
      Int_t ib=0;
      while(auto* arg=dynamic_cast<RooCategory*>(iter())){	
	TString catName=arg->GetName();
	branches[ib]=tree->Branch(catName,&branchVal[ib],catName+"/I");
	ib++;
      }
      //add ID branch
      auto id = fCurrSetup->WS().var(fCurrSetup->GetIDBranchName());
      TBranch *idbranch=nullptr;
      if(id!=nullptr)
	idbranch = tree->Branch(fCurrSetup->GetIDBranchName(),&fIDval,fCurrSetup->GetIDBranchName()+"/D");
      
      //Now loop over dataset
      for(Int_t entry=0;entry<numEntries;entry++){
	auto vars=fGenData->get(entry);
	ib=0;
	//extra for categories
	for(auto& branch: branches){
	  branchVal[ib++]=vars->getCatIndex(branch->GetName());
	  branch->Fill();
	}
	//extra for new ID branch
	if(id!=nullptr){
	  idbranch->Fill();
	  ++fIDval;
	}
      }
      tree->Write();
       // for(auto& branch: branches){
      // 	delete branch;
      // }

      //Use generated entry list to save full tree variables
      //Not just those used in generation
      auto PDFs = SetUp().PDFs();
      for(auto* pdf:PDFs){
	if(auto evPdf=dynamic_cast<RooHSEventsPDF*> (pdf) ) {
	  TFile entryFile(TString("entryFile_")+evPdf->GetName()+".root");
	  if(entryFile.IsOpen()==kFALSE) continue; //no events tree or entry list, will just have used Accept or Reject
	  auto entryList=dynamic_cast<TEntryList*>(entryFile.Get("GenEvents"));

	  auto idata=GetDataBin(GetFiti()); //bin for this generation
	  //get the tree used in RooHEEventsPDF::SetTree
	  auto filetree=FiledTree::Read(Bins().TreeName(evPdf->GetName()),
					Bins().FileNames(evPdf->GetName())[idata]);
	  
	  auto pdftree=filetree->Tree();
	  pdftree->SetEntryList(entryList);
	  outfile->cd();//new tree in outfile
	  auto generatedTree=pdftree->CopyTree("","");
	  generatedTree->Write();
	  delete generatedTree;
	}
      }
      
      delete tree;
      delete fGenData; fGenData=nullptr;
    }
    ///////////////////////////////////////////////////////////////
    void ToyManager::InitSummary(){
       
      std::unique_ptr<TFile> resFile{TFile::Open(SetUp().GetOutDir()+GetCurrName()+"/ToySummary.root","recreate")};
       
      
      auto initpars=*(dynamic_cast<RooArgSet*>(SetUp().ParsAndYields().selectByAttrib("Constant",kFALSE)));
      cout<<"ToyManager::InitSummary()"<<endl;
      initpars.Print("v");
      initpars.setName(InitialParsName());
      //initpars.Write();//commented out for amplitudes test as crashing?!>
      //resFile->WriteObject(&initpars,InitialParsName());//crahss when RooFormulaVar is used....
      //cant save pars to file directly if got formulavar connected 
      //dont know why
      RooDataSet d(InitialParsName(),InitialParsName(),initpars);
      d.fill();
      d.Print("v");
      d.Write();
     
    }
    ////////////////////////////////////////////////////////////////
    void ToyManager::PreRun(){
  
      //if using bins make sure directories are made
      Bins().GetBins().SetOutDir(SetUp().GetOutDir());
      Bins().GetBins().MakeDirectories();
    }
    ////////////////////////////////////////////////////////////////
    std::unique_ptr<FitManager> ToyManager::Fitter(){
      std::cout<<"ToyManager::Fitter() "<<std::endl;
      std::unique_ptr<FitManager> fit{new FitManager(*this)};      
      std::cout<<"ToyManager::Fitter() "<<fit.get()<<std::endl;
      //      fit->LoadData("ToyData",fToyFileNames);
      //if we have a eventpdf generation, use the original tree
      //filtered by the entry list. If we just use ToyData then
      //there can be some events outwith Dataset limits due to
      //resolution effects and generating with truth
      //Note the entry list probably does not work in the case of
      //multiple PDFS, might need to improve this
      if(SetUp().PDFs().getSize()==1){
	std::cout<<"ToyManager::Fitter() "<<std::endl;
 	auto evHPdf=dynamic_cast<RooHSEventsHistPDF*> (&(SetUp().PDFs()[0]));
	auto evPdf=dynamic_cast<RooHSEventsPDF*> (&(SetUp().PDFs()[0]));
	if(evHPdf!=nullptr){ //check if want to generate from histogram template
	  if(evHPdf->UsingHistGenerator()==kTRUE){ //use ToyData not eventree
	    evPdf=nullptr;
	  }
	}
	
	if(evPdf!=nullptr)
	  fit->LoadData(Bins().TreeName(evPdf->GetName()),fToyFileNames);
	else
	  fit->LoadData("ToyData",fToyFileNames);
      }
      else //multiple PDFs must use ToyData
	fit->LoadData("ToyData",fToyFileNames);

      std::cout<<"ToyManager::Fitter() "<<std::endl;
  
      fit->Data().Toys(fNToys);
      return fit;
      //      return std::move(fit);
    }
    //  std::shared_ptr<FitManager> ToyManager::Fitter(){
    //   std::shared_ptr<FitManager> fit{new FitManager(*this)};      
    //   fit->LoadData("ToyData",fToyFileNames);
    //   fit->Data().Toys(fNToys);
      
    //   return std::move(fit);
    // }
  
    ///////////////////////////////////////////////////////////////
    std::shared_ptr<ToyManager> ToyManager::GetFromFit(Int_t N,const TString& filename,const TString& resultFile){
       std::unique_ptr<TFile> fitFile{TFile::Open(filename)};
       // auto fit=dynamic_cast<FitManager*>( fitFile->Get("HSFit") );
       std::unique_ptr<FitManager> fit{dynamic_cast<FitManager*>( fitFile->Get("HSFit"))};
       return GetFromFit(N,*fit.get(),(resultFile));
       //       return GetFromFit(*fit,result);
     }
    ///////////////////////////////////////////////////////////////
    std::shared_ptr<ToyManager> ToyManager::GetFromFit(Int_t N,const std::shared_ptr<FitManager>& fit,const TString& resultFile){
      return GetFromFit(N,*fit.get(),(resultFile));
    }
    ////////////////////////////////////////////////////////////////
    std::shared_ptr<ToyManager> ToyManager::GetFromFit(Int_t N,FitManager& fit,const TString& resultFile){

      
      std::shared_ptr<ToyManager> toy{new ToyManager(N,fit,fit.SetUp().GetOutDir(),(resultFile))};

      return std::move(toy);
    }

    void ToyManager::UseMyToyData(FitManager& fitter){
      //Get data file names from summary file
      auto summaryFile = std::unique_ptr<TFile> (TFile::Open(SetUp().GetOutDir()+"/ToySummary.root") );
      std::vector<TString> *fnames=nullptr;
      summaryFile->GetObject("ToyFiles", fnames);
      fToyFileNames =*fnames;
      fitter.LoadData("ToyData",fToyFileNames);
      fNToys= fToyFileNames.size();

      //set number of toys
      fitter.Data().Toys(fNToys);

      //cp my summary file to fitter directory
      gSystem->Exec(Form("cp %s %s",(SetUp().GetOutDir()+"/ToySummary.root").Data(),(fitter.SetUp().GetOutDir()+"ToySummary.root").Data()));

      //Take the fitter output directory
      SetUp().SetOutDir(fitter.SetUp().GetOutDir());

    }
    
    void ToyManager::LoadResult(){
      if(fResultFileName==TString())
	return;
      cout<<"ToyManager::LoadResult()  "<<fResultFileName<<endl;
      TString resultFile=fResultOutDir+Bins().BinName(GetDataBin(GetFiti()))+"/"+fResultFileName;
      std::unique_ptr<TFile> fitFile{TFile::Open(resultFile)};
      std::unique_ptr<RooDataSet> result{dynamic_cast<RooDataSet*>( fitFile->Get(Minimiser::FinalParName()))};
      //Set the values of the paramteres to those in the given result
      if(result.get()){
	auto newPars = SetUp().ParsAndYields();
	auto* resAll = result->get(); //get all result info
	auto* resPars=resAll->selectCommon(newPars); //just select pars and yieds
	newPars.assignFast(*resPars); //set values to results
       	cout<<"ToyManager::LoadResult setting values from fit results "<<resultFile<<" : "<<endl;
	newPars.Print("v");
      }
    }
    //////////////////////////////////////////////////////////////////////
    void ToyManager::Summarise(){
      if(!Bins().GetSize())
	Summarise(0);
  
      for(UInt_t i=0;i<Bins().GetSize();i++)
	Summarise(i);
      
    }
    //////////////////////////////////////////////////////////////////////
    void ToyManager::Summarise(Int_t ibin){
      cout<<"Summarise "<<Minimiser::ResultTreeName() <<" "<<SetUp().GetOutDir()+Bins().BinName(ibin)<<endl;
      TChain resChain(Minimiser::ResultTreeName());
      resChain.Add(SetUp().GetOutDir()+Bins().BinName(ibin)+"/Results*.root");
      std::unique_ptr<TFile> resFile{TFile::Open(SetUp().GetOutDir()+Bins().BinName(ibin)+"/ToySummary.root","update")};

      resFile->WriteObject(&fToyFileNames,"ToyFiles");
 
      auto tree=resChain.CloneTree();
       
      //Loop over all parameters
      std::unique_ptr<RooDataSet> parsData{dynamic_cast<RooDataSet*>( resFile->Get(InitialParsName()))};
      cout<<" ToyManager::Summarise() Initial Parameters"<<endl;
      if(!parsData.get()) return;
      auto pars=parsData->get();
      // pars->Print("v");
      // auto pars = SetUp().ParsAndYields();   
      TIter iter=pars->createIterator();
      while(auto* arg=dynamic_cast<RooRealVar*>(iter())){	
	TString parName=arg->GetName();
	Double_t val=arg->getValV();
	//in casre a -ve sign on the name replace with minus for drawing
	if(parName.Contains("-")){
	  auto oldName=parName;
	  parName.ReplaceAll("-","minus");
	  tree->GetBranch(oldName)->SetName(parName);
	  tree->GetBranch(oldName+"_err")->SetName(parName+"_err");
	  
	}
	
	tree->Draw(parName,"","goff");
	//parameter
	auto *hpar = dynamic_cast<TH1F*>((gDirectory->FindObject("htemp"))->Clone(parName));
	Double_t mean=hpar->GetMean();
	Double_t rms=hpar->GetRMS();

	//parameter error
	tree->Draw(parName+"_err","","goff");
	//parameter
	TH1F *hparErr = dynamic_cast<TH1F*>((gDirectory->FindObject("htemp"))->Clone(parName+"_err"));
	Double_t meanErr=hparErr->GetMean();
	Double_t rmsErr=hparErr->GetRMS();

	//pull distribution
	tree->Draw(Form("(%s-%lf)/%s_err",parName.Data(),mean,parName.Data()),"","goff");
	TH1F *hpull = dynamic_cast<TH1F*>((gDirectory->FindObject("htemp"))->Clone(parName+"_pull"));
	Double_t meanPull=hpull->GetMean();
	Double_t rmsPull=hpull->GetRMS();


	//bias
	tree->Draw(Form("(%s-%lf)",parName.Data(),val),"","goff");
	TH1F *hbias = dynamic_cast<TH1F*>((gDirectory->FindObject("htemp"))->Clone(parName+"_bias"));
	Double_t meanBias=hbias->GetMean();
	
	//bias pull
	tree->Draw(Form("(%s-%lf)/%s_err",parName.Data(),val,parName.Data()),"","goff");
	TH1F *hbiasPull = dynamic_cast<TH1F*>((gDirectory->FindObject("htemp"))->Clone(parName+"_biasPull"));
	Double_t meanBPull=hbiasPull->GetMean();
	Double_t rmsBPull=hbiasPull->GetRMS();

	

	hpar->SetNameTitle(parName,parName);
	hpar->Write(nullptr,TFile::kOverwrite);
	hparErr->SetNameTitle(parName+"_err",parName+"_err");
	hparErr->Write(nullptr,TFile::kOverwrite);
	hpull->SetNameTitle(parName+"_pull",parName+"_pull");
	hpull->Write(nullptr,TFile::kOverwrite);
	hbias->SetNameTitle(parName+"_bias",parName+"_bias");
	hpull->Write(nullptr,TFile::kOverwrite);
	hbiasPull->SetNameTitle(parName+"_biasPull",parName+"_biasPull");
	hbiasPull->Write(nullptr,TFile::kOverwrite);
	
	cout<<endl<<parName <<" "<<mean<<" +- "<<meanErr<<" sigma "<<rms<<" meanPull "<<meanPull<<" sigmaPull "<<rmsPull<<"\n      bias "<< meanBias << " bias Pull "<<meanBPull<< " sigma "<<rmsBPull<<endl<<endl;
	
      }

    }

   
   
  }
}	

#include "FitManager.h"

#include <memory>
#include "RooHSEventsPDF.h"
#include "RooHSEventsHistPDF.h"
#include "RooComponentsPDF.h"
#include "TSystem.h"


namespace HS{
  namespace FIT{

 
    FitManager::FitManager(const FitManager& other):TNamed(other.fName,other.fName){
      fSetup=other.fSetup;
      fBinner=other.fBinner;
      //LoadData(other.GetDataTreeName(),other.GetDataFileNames());
  
    }

    FitManager&  FitManager::operator=(const FitManager& other){
      cout<<"=============FitManager"<<endl;
      fSetup=other.fSetup;
      fBinner=other.fBinner;
      //LoadData(other.fData.ParentTreeName(),other.fData.FileNames());
      return *this;
    }
    
    void FitManager::Run(){
      
      CreateCurrSetup();
     
      //get dataset fFiti
      fCurrDataSet=std::move(Data().Get(fFiti));

      //Look for Special case of RooHSEventsPDFs
      FillEventsPDFs();
   
      //Add fit constraints
      fCurrSetup->AddFitOption(RooFit::ExternalConstraints
			       (fCurrSetup->Constraints()));
      //initialise yields
      if(fCurrSetup->Yields().getSize()==1){//special case only 1 yield)
	Double_t yld=fCurrDataSet->sumEntries();
	SetAllValLimits(fCurrSetup->Yields(),
			yld,0,1.2*yld);
      }
      else
	SetAllValLimits(fCurrSetup->Yields(),
			fCurrDataSet->sumEntries()/2,0,fCurrDataSet->sumEntries()*1.2);
      //create extended max likelihood pdf
      //fCurrSetup->Parameters().Print("v");
      fCurrSetup->TotalPDF();
      FitTo();
    }
    void FitManager::CreateCurrSetup(){
      fCurrSetup = std::make_unique<Setup>(fSetup); //Copy setup from template
      fCurrSetup->SetName(GetCurrName());
      fCurrSetup->SetTitle(GetCurrTitle());
      //make sure we take current setup values
      //If not it will use the string from Factory() etc,
      fCurrSetup->ParsAndYields().assignFast(fSetup.ParsAndYields());

      //Look to see if taking previous fit results as initial pars
      if(fUsePrevResult){
	LoadPrevResult(fPrevResultDir,fPrevResultMini);
      }
    }
    /////////////////////////////////////////////////////////////
    void FitManager::RunAll(){

      PreRun();

      UInt_t Nf=GetN();
      for(UInt_t i=0;i<Nf;i++){
	 RunOne(i);
      }
    }

    ////////////////////////////////////////////////////////////
    void FitManager::FitTo(){
      if(!fMinimiser.get()) SetMinimiser(new HS::FIT::Minuit2());
      fMinimiser->Run(*fCurrSetup,*fCurrDataSet);
      
      ///////////////////////////
      //Plot best fit and return
      PlotDataModel();
	  CalcAcceptanceCorrection();

    }
    void FitManager::RunOne(Int_t ifit){
      fFiti=ifit;
      if(fRedirect) RedirectOutput(fSetup.GetOutDir()+Form("logRooFit%d.txt",fFiti));
      Run();
      if(fRedirect) RedirectOutput();
      SaveResults();
      
      Reset();
    }
    
    void FitManager::FillEventsPDFs(){
      UInt_t idata=GetDataBin(fFiti);

      auto pdfs=fCurrSetup->PDFs();

      auto savedir=gDirectory;
      
      for(Int_t ip=0;ip<pdfs.getSize();ip++){
	auto pdf=dynamic_cast<RooHSEventsPDF*>( &pdfs[ip]);

	if(pdf){
	  //  pdf->SetConstInt();
	  if(fBinner.FileNames(pdf->GetName()).size()==0)
	    continue;
	  auto filetree=FiledTree::
	    Read(fBinner.TreeName(pdf->GetName()),
		 fBinner.FileNames(pdf->GetName())[idata]);
	  auto tree=filetree->Tree();
	
	  auto mcgenfiletree= (fBinner.FileNames(pdf->GetName()+TString("__MCGen")).empty() ? nullptr : FiledTree::Read(fBinner.TreeName(pdf->GetName()+TString("__MCGen")),fBinner.FileNames(pdf->GetName()+TString("__MCGen"))[idata]));
	  auto mcgentree=(mcgenfiletree ? mcgenfiletree->Tree() : nullptr);
	  
	  savedir->cd();
 	  if(!tree.get()){
	    cout<<"WARNING FitManager::FillEventsPDFs :"<<
	      "    No tree data found for EventPDF "<<pdf->GetName()<<endl;
	    continue;
	  }
	  //if too few events remove this PDF
	  if(!tree->GetEntries()||!pdf->IsValid()){
	    fCurrSetup->Yields().remove(fCurrSetup->Yields()[ip]);
	    pdfs.remove(*pdf);
	    ip--;
	  }
	  else{ //use it and give it the simulated tree
	    pdf->SetInWeights(fCurrSetup->GetPDFInWeights(pdf->GetName()));
		pdf->SetEvTree(tree.get(),fCurrSetup->Cut(),mcgentree.get());

	    //See if data to load for proto data
	    if(!fCurrDataSet.get())
	      fCurrDataSet=std::move(Data().Get(idata));

	    if(fCurrDataSet.get())pdf->AddProtoData(fCurrDataSet.get());
	    RooHSEventsHistPDF* histspdf=nullptr;
	    if((histspdf=dynamic_cast<RooHSEventsHistPDF*>(pdf))){
	      histspdf->CreateHistPdf();
	      fCurrSetup->AddGausConstraint(histspdf->AlphaConstraint());
	      fCurrSetup->AddGausConstraint(histspdf->OffConstraint());
	      fCurrSetup->AddGausConstraint(histspdf->ScaleConstraint());
	    }
	  }
	  //keep the simulated tree alive until Reset()
	  fFiledTrees.push_back(std::move(filetree));	
	  if(mcgenfiletree)
		fFiledTrees.push_back(std::move(mcgenfiletree));	  
	}
      }
      savedir->cd();
    }
    void FitManager::SaveSetup(){
      auto file=TFile::Open(fSetup.GetOutDir()+"HSSetup.root","recreate");
      fSetup.Write("HSSetup");
      delete file;
    }
    void FitManager::LoadSetup(const TString& dir){
      auto file=TFile::Open(dir+"/HSSetup.root");
      fSetup=*(dynamic_cast<HS::FIT::Setup*>(file->Get("HSSetup")));
      delete file;
    }
    //Read in paarameters from previous fit
    void FitManager::InitPrevResult(const TString& resultDir,const TString& resultMinimiser){
      fUsePrevResult=kTRUE;
      
      if(resultDir==TString()) fPrevResultDir=fSetup.GetOutDir(); //use current
      else fPrevResultDir=resultDir;

      if(resultMinimiser==TString())fPrevResultMini=fMinimiser->GetName();
      else fPrevResultMini=resultMinimiser;
      
    }
    
    void FitManager::LoadPrevResult(const TString& resultDir,const TString& resultMinimiser){
      TString resultFile=resultDir+"/"+fCurrSetup->GetName()+"/Results"+fCurrSetup->GetTitle()+resultMinimiser+".root";

     
      std::unique_ptr<TFile> fitFile{TFile::Open(resultFile)};
      std::unique_ptr<RooDataSet> result{dynamic_cast<RooDataSet*>( fitFile->Get(Minimiser::FinalParName()))};
      //fitFile.reset();
      //      auto result=dynamic_cast<RooDataSet*>( fitFile->Get(Minimiser::FinalParName())->Clone());//**
      //Set the values of the paramteres to those in the given result
      if(result.get()){
      //if(result){
	auto newPars = fCurrSetup->ParsAndYields();
	auto* resAll = result->get(); //get all result info
	auto* resPars=resAll->selectCommon(newPars); //just select pars and yieds
       	newPars.assignFast(*resPars); //set values to results
	cout<<"FitManager::LoadResult setting values from fit results "<<resultFile<<" : "<<endl;
	newPars.Print("v");
	//	delete result;result=nullptr;
      }
    }
    void FitManager::WriteThis(){
      auto file=TFile::Open(fSetup.GetOutDir()+"HSFit.root","recreate");
      if(!fMinimiser.get()) SetMinimiser(new HS::FIT::Minuit2());
      file->WriteObject(this,"HSFit");
      file->WriteObject(fMinimiser.get(),fMinimiserType);

      if(fCompiledMacros.size()){
	auto* macList=new TList();
	//	macList->SetName("HS_COMPILEDMACROS");
	macList->SetOwner();
	for(auto& macro : fCompiledMacros)
	  macList->Add(new TObjString(macro));
	file->WriteObject(macList,"HS_COMPILEDMACROS");
      }
      delete file;
    }
    void FitManager::RedirectOutput(const TString& log){
      Info("FitManager::RedirectOutput",Form("text ouput will be sent to file %s",log.Data()));
      if(log==TString(""))
	gSystem->RedirectOutput(nullptr,"w");
      else
	gSystem->RedirectOutput(log.Data(),"w");
      
    }
    
	void FitManager::CalcAcceptanceCorrection(){
		UInt_t idata=GetDataBin(fFiti);
		auto pdfs=fCurrSetup->PDFs();
		
		fAcceptanceTree = new TTree("acceptanceCorrection","acceptanceCorrection");
		
		cout<< "FitManager::CalcAcceptanceCorrection() Found " << pdfs.getSize() << "pdfs." << endl;
		for(Int_t ip=0;ip<pdfs.getSize();ip++){
			auto pdf=dynamic_cast<RooHSEventsPDF*>( &pdfs[ip]);
			if(pdf){
				if(fBinner.FileNames(pdf->GetName()).size()==0)
					continue;
// 				
				Double_t integralAccepted=pdf->unnormalisedIntegral(1,"");
				Double_t integralGenerated=pdf->unnormalisedIntegral(2,"");
				Double_t sumofweightsData=fData.Get(idata)->sumEntries();
				Double_t acceptanceCorrectedYield=0;
				if(integralGenerated)
					acceptanceCorrectedYield=sumofweightsData/(integralAccepted/integralGenerated);
				
				if(integralGenerated)
					cout << "FitManager::CalcAcceptanceCorrection() accepted=" << integralAccepted << " generated=" << integralGenerated << " ratio=" << integralAccepted/integralGenerated << endl;
				else
					cout << "FitManager::CalcAcceptanceCorrection() accepted=" << integralAccepted << " generated=" << integralGenerated << " Can't calculate acceptance!!!" << endl;
				cout << "FitManager::CalcAcceptanceCorrection() data yield=" << sumofweightsData << endl;
				cout << "FitManager::CalcAcceptanceCorrection() acceptance corrected data yield=" << acceptanceCorrectedYield << endl;
				
				fAcceptanceTree = new TTree("acceptanceCorrection","acceptanceCorrection");
				fAcceptanceTree->Branch(pdf->GetName()+TString("_acc"),&integralAccepted);
				fAcceptanceTree->Branch(pdf->GetName()+TString("_gen"),&integralGenerated);
				fAcceptanceTree->Branch(pdf->GetName()+TString("_yld"),&sumofweightsData);
				fAcceptanceTree->Branch(pdf->GetName()+TString("_yldcorr"),&acceptanceCorrectedYield);
				fAcceptanceTree->Fill();
			}
		}
	}

    void FitManager::SaveResults(){
     
      auto outFile=fMinimiser->SaveInfo();
      if(fPlots.size())fPlots.back()->Write(); //just save the last one
      if(fAcceptanceTree) fAcceptanceTree->Write(); //save acceptance TTree, if generated MC not passed to FitManager this is just a nullptr
      //outfile is unique_ptr so will be deleted and saved here
    }

  }//namespace FIT
}//namespace HS

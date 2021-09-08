#include "sPlot.h"

#include <memory>
#include "TDirectory.h"
#include "TBenchmark.h"

namespace HS{
  namespace FIT{


  
    Bool_t sPlot::Run(){
      cout<<"HS::FIT::sPlot::Do prelimanry fits "<<gDirectory->GetName()<<endl;

      if(FitManager::Run()==kFALSE) return kFALSE;

      //Zero yield check
      //returns true if single yield left and weights generated
      //else if zero yields it will remove them for sPlot
      if(ZeroYieldCheck()==kTRUE) return kTRUE;
   
 
      //Note sPlot is much (10X) faster with tree store
      //Normal fit is 2X faster with vector...
      //RooAbsData::setDefaultStorageType(RooAbsData::Tree);
      //auto* dataset =dynamic_cast<RooDataSet*>( fCurrDataSet->emptyClone());
      //dataset->append(*fCurrDataSet.get());
      //RooAbsData::setDefaultStorageType(RooAbsData::Vector);

      //Note at tested again with 6.24 and Vector store is now faster...
       RooDataSet* dataset =fCurrDataSet.get();

      
      auto *model=fCurrSetup->Model();
      
      ////////////////////////////////////////////////////////
       //sPlot
       cout<<"HS::FIT::sPlot::Run create sWeights "<<endl;
       fCurrSetup->Parameters().setAttribAll("Constant");
       fSPlot.reset(new RooStats::SPlot{"splot_fit", "sPlot Fit",
	     *dataset,model,fCurrSetup->Yields()});
      
       fCurrSetup->Parameters().setAttribAll("Constant",kFALSE);

       CreateWeights();
       
       // PlotDataModel();

       //delete dataset;

       return kTRUE;
   }
    
    void sPlot::CreateWeights(){
      //If single yield it gets new weight =1; other species weight 0
      //if(fSingleYield!=TString("")) ExportWeights();
      //Check that the fit was succesfull
      Double_t TotalYield=0;
      auto yields=fCurrSetup->Yields();
      for(Int_t iy=0;iy<yields.getSize();iy++)
    	TotalYield+=(dynamic_cast<RooRealVar*>(&yields[iy]))->getVal();

      if(TotalYield>0){ //got some weights
    	fWeights=std::make_shared<Weights>("HSsWeights");//initialise weights
    	fWeights->SetIDName(fCurrSetup->GetIDBranchName());
    	fWeights->SetTitle(fCurrSetup->GetName());
    	fWeights->SetFile(fCurrSetup->GetOutDir()+TString("Weights")+fCurrSetup->GetName()+fCurrSetup->GetTitle()+".root");
    	ExportWeights();
 	//fWeights->PrintWeight();
	fWeights->SortWeights();

	}
    
      else Warning("sPlot::sPlot()"," total weights 0, fit did not converge. Make sure the non-sweight fit to fix parameters was succesful. No weights will be assigned for these events");
      
    }

    void sPlot::ExportWeights(){
      cout<<"sPlot::ExportWeights "<<endl;
      const TString idname=fCurrSetup->GetIDBranchName();
      // cout<<"sPlot::ExportWeights()  "<<idname<<endl;
      const RooArgSet* vars=fCurrDataSet->get(0);
      Bool_t gotID=kFALSE;
      if(vars->find(idname))
	gotID=kTRUE;
      
      auto yields=fCurrSetup->Yields();
      
      Int_t NSpecies=yields.getSize();
      TVectorD eventW(NSpecies); //initialise weights vector
      for(Int_t iw=0;iw<NSpecies;iw++){//set name for each species, 
	fWeights->SetSpecies(TString(yields.at(iw)->GetName()).Remove(0,4));
      }
      //include species where the yeild was found to be zero, i.e. with 0 weight
      for(const auto& zeroSpecies:fZeroYields){
	fWeights->SetSpecies(TString(zeroSpecies).Remove(0,4));
	eventW.ResizeTo(eventW.GetNrows()+1);
      }
      
        
      //Loop over all events and asign weights
      for(Long64_t ev=0;ev<fCurrDataSet->numEntries();ev++){//loop over events
	//Include special case of single species
	if(NSpecies==1)	{
	  fCurrDataSet->get(ev); //move to this event
	  eventW[0]=fCurrDataSet->weight(); //get input weight, which as no other species must also be the output weight
	}
	else{//normal fill with sWeights
	
	  for(Int_t iw=0;iw<NSpecies;iw++){//loop over species
	    eventW[iw]=fSPlot->GetSWeight(ev,yields.at(iw)->GetName());//get weight for this species
	  }
	}
	
	//Include special case of weights when zero yield
	UInt_t iz=NSpecies; //now count on from NSpecies
	for(const auto& zeroSpecies:fZeroYields){
	  //give 0 weight to species with 0 yield
	  eventW[iz++]=0;
	}

       

	if(gotID){//use ID from initial tree
	  auto* vars=fCurrDataSet->get(ev);
	  fWeights->FillWeights((Long64_t)vars->getRealValue(idname),eventW);
	} //ID not defined just use entry number in dataset
	else fWeights->FillWeights(ev,eventW);
      }
    }
    
    weights_uptr sPlot::MergeWeights(){
      //in addition combine the weights into 1 and load them
      weights_uptr wts(new Weights("HSsWeights"));
      //Note the output file cannot contain the word Weights (because of Merge), hence Tweights!
      wts->Merge(SetUp().GetOutDir()+"/Weights",
		 SetUp().GetOutDir()+"/"+SetUp().GetName()+"Tweights.root",
		 "HSsWeights");
      //wts->Save();

      //reset to save and reopen
      wts.reset(new Weights{});
      wts->LoadSaved(SetUp().GetOutDir()+TString("Tweights.root"),"HSsWeights");
      
      return std::move(wts);
    }

    void sPlot::WeightedTree(){
      if(!fWeights.get()){
	if(Bins().GetSize()>0)
	  fWeights = MergeWeights();
	else{
	  fWeights.reset(new Weights());
	  fWeights->LoadSaved(SetUp().GetOutDir()+TString("Weights")+SetUp().GetName()+".root","HSsWeights");
	}
      }
      TDirectory* saveDir=gDirectory;
      
      //Open tree file
      auto ftree=FiledTree::Read(Data().ParentTreeName(),
      				 Data().ParentName());
      //create a copy in a new file to append the weights to
      //Keep it in file as large trees can use too much memory
      fWeightedFiledTree=(FiledTree::CloneFull(ftree->Tree(),SetUp().GetOutDir()+"DataWeightedTree.root"));
  
      //delete original
      ftree.reset();
      //Add weights to tree
      fWeights->AddToTree(fWeightedFiledTree->Tree().get());
      fWeightedFiledTree->Tree()->SetBranchStatus("*",true);

      saveDir->cd();
    }

    void sPlot::DrawWeighted(const TString& var,const TString& wname,TString cut, const TString& opt){
      if(!fWeightedFiledTree.get())
	WeightedTree();

      if(!fWeightedFiledTree.get())
	return;

 
      if(cut==TString()) cut="1";
     
      fWeightedFiledTree->Tree()->Draw(var,wname+"*("+cut+")",opt);

    }

    Bool_t sPlot::ZeroYieldCheck(){
   auto& yields=fCurrSetup->Yields(); //get reference to yields
      auto& pdfs=fCurrSetup->PDFs(); //get reference to yields
      bool removedPdf=false;    
      for(Int_t iy=0;iy<yields.getSize();iy++){
	auto checkYield=dynamic_cast<RooRealVar*>(&yields[iy]);
	Double_t  thisYield=checkYield->getVal();
	if(thisYield<1E-2){
	  //Need to remove this pdf all weights weights will not be written for this species
	  Warning("sPlot::sPlot()",Form("Found zero yield for %s, will remove from sPlot, weights for this species will be set to 0",checkYield->GetName()),"");
	  pdfs.remove(pdfs[iy]);
	  //yields->remove(*(fCurrSetup->WS().var(yield0->GetName())));
	  yields.remove(yields[iy]);

	  removedPdf=true;
	  fZeroYields.push_back(checkYield->GetName());
	}

      }
      //remake sum of species without zero yied pdf
      if(removedPdf)fCurrSetup->TotalPDF();

      if(yields.getSize()==1){//Only 1 species all weights==inWeights
	fSingleYield=yields[0].GetName();
	CreateWeights();
	return kTRUE;
      }
      return kFALSE;
    }

    
  }//namespace FIT 
}//namespace HS

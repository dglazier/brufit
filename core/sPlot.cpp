/**
 * @file sPlot.cpp
 * @brief Implementation of the sPlot class.
 */

#include "sPlot.h"

#include <memory>
#include <TDirectory.h>
#include <TBenchmark.h>
#include <iostream>

namespace HS {
namespace FIT {
  // Copy Constructor
  sPlot::sPlot(const sPlot& other) : FitManager(other) {
    fSingleYield = other.fSingleYield;
    fZeroYields = other.fZeroYields;
    // Notice we intentionally DO NOT copy fSPlot, fWeights, or fWeightedFiledTree
  }

    // Assignment Operator
    sPlot& sPlot::operator=(const sPlot& other) {
        if (this != &other) {
            FitManager::operator=(other);
            fSingleYield = other.fSingleYield;
            fZeroYields = other.fZeroYields;
        }
        return *this;
    }
    // ========================================================================
    // Core Execution
    // ========================================================================

    Bool_t sPlot::Run() {
        std::cout << "HS::FIT::sPlot::Run() Do preliminary fits in " << gDirectory->GetName() << std::endl;

        // Perform the standard unbinned maximum likelihood fit first via the base class
        if (!FitManager::Run()) return kFALSE;

        // Clean up zero yields to prevent matrix inversion errors in sPlot
        fZeroYields.clear();
        fSingleYield = "";
        if (ZeroYieldCheck()) return kTRUE; // If only one species is left, weights are trivial (1.0 or 0.0)
   
        // --- SAFE DOWNCAST ---
        // FitManager holds fCurrDataSet as a generic RooAbsData (could be unbinned or binned).
        // RooStats::SPlot strictly requires an unbinned RooDataSet. We cast and verify here.
        auto* dataset = dynamic_cast<RooDataSet*>(fCurrDataSet.get());
        if (!dataset) {
            std::cerr << "ERROR (sPlot::Run): sPlot requires an unbinned RooDataSet. "
                      << "Binned fits (RooDataHist) are not supported by RooStats::SPlot." << std::endl;
            return kFALSE;
        }

        auto* model = fCurrSetup->Model();
      
        // ========================================================
        // sPlot Calculation
        // ========================================================
        std::cout << "HS::FIT::sPlot::Run() Create sWeights " << std::endl;
        
        // Fix shape parameters to calculate weights based strictly on the fitted yields
        fCurrSetup->Parameters().setAttribAll("Constant");
       
        // Instantiate the SPlot object (RAII handled by unique_ptr)
        fSPlot = std::make_unique<RooStats::SPlot>("splot_fit", "sPlot Fit", *dataset, model, fCurrSetup->Yields());
      
        // Release parameters back to floating status for subsequent operations
        fCurrSetup->Parameters().setAttribAll("Constant", kFALSE);

        CreateWeights();
       
        return kTRUE;
    }
    
    // ========================================================================
    // Weight Processing
    // ========================================================================

    void sPlot::CreateWeights() {
        Double_t TotalYield = 0;
        auto& yields = fCurrSetup->Yields();
        
        // Accumulate total yield across all active species
        for (auto* arg : yields) {
            if (auto* y = dynamic_cast<RooRealVar*>(arg)) {
                TotalYield += y->getVal();
            }
        }

        // Only proceed if the fit actually found signal/background
        if (TotalYield > 0) { 
            fWeights = std::make_shared<Weights>("HSsWeights");
            fWeights->SetIDName(fCurrSetup->GetIDBranchName());
            fWeights->SetTitle(fCurrSetup->GetName());
            fWeights->SetFile(fCurrSetup->GetOutDir() + "Weights" + fCurrSetup->GetName() + fCurrSetup->GetTitle() + ".root");
            
            ExportWeights();
            fWeights->SortWeights(); // Organize weights by event ID for efficient merging
        } else {
            Warning("sPlot::CreateWeights()", 
                    "Total weights 0, fit did not converge. Make sure the non-sweight fit to fix parameters was successful. No weights will be assigned for these events");
        }
    }

    void sPlot::ExportWeights() {
        std::cout << "sPlot::ExportWeights" << std::endl;
        const TString idname = fCurrSetup->GetIDBranchName();
        
        // Check if the dataset contains a unique event ID variable
        auto* vars = fCurrDataSet->get(0);
        Bool_t gotID = (vars && vars->find(idname));
      
        auto& yields = fCurrSetup->Yields();
        Int_t NSpecies = yields.getSize();
        TVectorD eventW(NSpecies); 
        
        // Register the names of the active species in the Weights object
        for (Int_t iw = 0; iw < NSpecies; iw++) {
            fWeights->SetSpecies(TString(yields[iw].GetName()).Remove(0, 4)); // Remove "Yld_" prefix
        }
        
        // Register species that collapsed to zero (they will receive a forced 0.0 weight)
        for (const auto& zeroSpecies : fZeroYields) {
            fWeights->SetSpecies(TString(zeroSpecies).Remove(0, 4));
            eventW.ResizeTo(eventW.GetNrows() + 1);
        }
              
        // Iterate through the unbinned dataset to extract weights event-by-event
        for (Long64_t ev = 0; ev < fCurrDataSet->numEntries(); ev++) {
            
            if (NSpecies == 1) { 
                // Trivial case: Only 1 valid species means it gets 100% of the input weight
                fCurrDataSet->get(ev); 
                eventW[0] = fCurrDataSet->weight(); 
            } else { 
                // Normal case: Extract calculated sWeights from the SPlot object
                for (Int_t iw = 0; iw < NSpecies; iw++) {
                    eventW[iw] = fSPlot->GetSWeight(ev, yields[iw].GetName());
                }
            }
            
            // Assign a hard 0.0 weight to species that collapsed during the fit
            UInt_t iz = NSpecies; 
            for (size_t z = 0; z < fZeroYields.size(); z++) {
                eventW[iz++] = 0;
            }

            // Save the weight vector mapped to the correct event ID
            if (gotID) { 
                auto* eventVars = fCurrDataSet->get(ev);
                fWeights->FillWeights((Long64_t)eventVars->getRealValue(idname), eventW);
            } else { 
                // Fallback: Use the sequential entry index if no ID variable exists
                fWeights->FillWeights(ev, eventW);
            }
        }
        
        // Final sanity check to ensure no events were skipped
        if (fWeights->Size() != fCurrDataSet->numEntries()) {
            std::cout << "FATAL (sPlot::ExportWeights): Mismatch between number of weights (" 
                      << fWeights->Size() << ") and number of data events (" 
                      << fCurrDataSet->sumEntries() << ")." << std::endl;
            yields.Print("v");
            exit(0); // Better to hard-exit than proceed with corrupted weights
        }
    }
    
    weights_uptr sPlot::MergeWeights() {
        std::cout << "sPlot::MergeWeights()" << std::endl;
        auto wts = std::make_unique<Weights>("HSsWeights");
        
        // Merges chunked weight files into one consolidated file
        wts->Merge(SetUp().GetOutDir() + "/Weights",
                   SetUp().GetOutDir() + "/" + SetUp().GetName() + "Tweights.root",
                   "HSsWeights");

        // Reload the newly merged file to verify integrity and return it
        wts = std::make_unique<Weights>();
        wts->LoadSaved(SetUp().GetOutDir() + "Tweights.root", "HSsWeights");
        
        return wts;
    }

    // ========================================================================
    // TTree operations
    // ========================================================================

    void sPlot::WeightedTree() {
        // Retrieve or merge weights if not already loaded into memory
        if (!fWeights) {
            if (Bins().GetSize() > 0) {
                fWeights = MergeWeights();
            } else {
                fWeights = std::make_shared<Weights>();
                fWeights->LoadSaved(SetUp().GetOutDir() + "Weights" + SetUp().GetName() + ".root", "HSsWeights");
            }
        }
        
        TDirectory* saveDir = gDirectory;
      
        // Read the original unweighted parent tree
        auto ftree = FiledTree::Read(Data().ParentTreeName(), Data().ParentName());
        
        // Deep clone the tree to a new file so we can safely inject the new weight branches
        fWeightedFiledTree = FiledTree::CloneFull(ftree->Tree(), SetUp().GetOutDir() + "DataWeightedTree.root");
  
        // Attach the calculated weights to the cloned tree as new branches
        fWeights->AddToTree(fWeightedFiledTree->Tree().get());
        fWeightedFiledTree->Tree()->SetBranchStatus("*", true);

        saveDir->cd();
    }

    void sPlot::DrawWeighted(const TString& var, const TString& wname, TString cut, const TString& opt) {
        if (!fWeightedFiledTree) WeightedTree();
        if (!fWeightedFiledTree) return; // Fail gracefully if tree cloning failed

        if (cut == "") cut = "1";
     
        // Draw using standard ROOT syntax, scaling the variable by the designated weight species
        fWeightedFiledTree->Tree()->Draw(var, wname + "*(" + cut + ")", opt);
    }

    Bool_t sPlot::ZeroYieldCheck() {
        auto& yields = fCurrSetup->Yields(); 
        auto& pdfs = fCurrSetup->PDFs(); 
        bool removedPdf = false;    
        
        // --- SAFE REMOVAL LOOP ---
        // We iterate backwards through the ROOT collection. If we iterate forward (0 -> size) 
        // and remove index 0, what used to be index 1 shifts down to 0, and the next loop 
        // iteration (i=1) will skip it entirely. Backward iteration prevents this shifting bug.
        for (Int_t iy = yields.getSize() - 1; iy >= 0; iy--) {
            auto* checkYield = dynamic_cast<RooRealVar*>(&yields[iy]);
            if (!checkYield) continue;

            // If a yield collapsed to near-zero, remove the associated PDF completely.
            // sPlot relies on covariance matrix inversion, which fails (singular matrix) 
            // if a component has zero yield.
            if (checkYield->getVal() < 1E-2) {
	      Warning("sPlot::ZeroYieldCheck()", 
                        "%s", Form("Found zero yield for %s, removing from sPlot. Weights will be 0.", checkYield->GetName()));
	      
                pdfs.remove(pdfs[iy]);
                yields.remove(yields[iy]);

                removedPdf = true;
                fZeroYields.push_back(checkYield->GetName());
            }
        }
        
        // Rebuild the sum of species without the zero-yield PDF
        if (removedPdf) fCurrSetup->TotalPDF();

        // If all but one species collapsed, calculating weights is trivial
        if (yields.getSize() == 1) { 
            fSingleYield = yields[0].GetName();
            CreateWeights();
            return kTRUE;
        }
        return kFALSE;
    }

} // namespace FIT 
} // namespace HS

// #include "sPlot.h"

// #include <memory>
// #include "TDirectory.h"
// #include "TBenchmark.h"

// namespace HS{
//   namespace FIT{


  
//     Bool_t sPlot::Run(){
//       cout<<"HS::FIT::sPlot::Do prelimanry fits "<<gDirectory->GetName()<<endl;

//       if(FitManager::Run()==kFALSE) return kFALSE;

//       //Zero yield check
//       //returns true if single yield left and weights generated
//       //else if zero yields it will remove them for sPlot
//       fZeroYields.clear();
//       if(ZeroYieldCheck()==kTRUE) return kTRUE;
   
 
//       RooDataSet* dataset =fCurrDataSet.get();

      
//       auto *model=fCurrSetup->Model();
      
//       ////////////////////////////////////////////////////////
//        //sPlot
//        cout<<"HS::FIT::sPlot::Run create sWeights "<<endl;
//        fCurrSetup->Parameters().setAttribAll("Constant");
//        fSPlot.reset(new RooStats::SPlot{"splot_fit", "sPlot Fit",
// 	     *dataset,model,fCurrSetup->Yields()});
      
//        fCurrSetup->Parameters().setAttribAll("Constant",kFALSE);

//        CreateWeights();
       
//        // cout<<"HS::FIT::sPlot::Run Done "<<endl;
//        return kTRUE;
//    }
    
//     void sPlot::CreateWeights(){
//       //If single yield it gets new weight =1; other species weight 0
//       //if(fSingleYield!=TString("")) ExportWeights();
//       //Check that the fit was succesfull
//       Double_t TotalYield=0;
//       auto yields=fCurrSetup->Yields();
//       for(Int_t iy=0;iy<yields.getSize();iy++)
//     	TotalYield+=(dynamic_cast<RooRealVar*>(&yields[iy]))->getVal();

//       if(TotalYield>0){ //got some weights
//     	fWeights=std::make_shared<Weights>("HSsWeights");//initialise weights
//     	fWeights->SetIDName(fCurrSetup->GetIDBranchName());
//     	fWeights->SetTitle(fCurrSetup->GetName());
//     	fWeights->SetFile(fCurrSetup->GetOutDir()+TString("Weights")+fCurrSetup->GetName()+fCurrSetup->GetTitle()+".root");
//     	ExportWeights();
//  	//fWeights->PrintWeight();
// 	fWeights->SortWeights();

// 	}
    
//       else Warning("sPlot::sPlot()"," total weights 0, fit did not converge. Make sure the non-sweight fit to fix parameters was succesful. No weights will be assigned for these events");
      
//     }

//     void sPlot::ExportWeights(){
//       cout<<"sPlot::ExportWeights "<<endl;
//       const TString idname=fCurrSetup->GetIDBranchName();
//       // cout<<"sPlot::ExportWeights()  "<<idname<<endl;
//       const RooArgSet* vars=fCurrDataSet->get(0);
//       Bool_t gotID=kFALSE;
//       if(vars->find(idname))
// 	gotID=kTRUE;
      
//       auto yields=fCurrSetup->Yields();
      
//       Int_t NSpecies=yields.getSize();
//       TVectorD eventW(NSpecies); //initialise weights vector
//       for(Int_t iw=0;iw<NSpecies;iw++){//set name for each species, 
// 	fWeights->SetSpecies(TString(yields.at(iw)->GetName()).Remove(0,4));
//       }
//       //include species where the yeild was found to be zero, i.e. with 0 weight
//       for(const auto& zeroSpecies:fZeroYields){
// 	fWeights->SetSpecies(TString(zeroSpecies).Remove(0,4));
// 	eventW.ResizeTo(eventW.GetNrows()+1);
//       }
      
        
//       //Loop over all events and asign weights
//       for(Long64_t ev=0;ev<fCurrDataSet->numEntries();ev++){//loop over events
// 	//Include special case of single species
// 	if(NSpecies==1)	{
// 	  fCurrDataSet->get(ev); //move to this event
// 	  eventW[0]=fCurrDataSet->weight(); //get input weight, which as no other species must also be the output weight
// 	}
// 	else{//normal fill with sWeights
	
// 	  for(Int_t iw=0;iw<NSpecies;iw++){//loop over species
// 	    eventW[iw]=fSPlot->GetSWeight(ev,yields.at(iw)->GetName());//get weight for this species
// 	  }
// 	}
	
// 	//Include special case of weights when zero yield
// 	UInt_t iz=NSpecies; //now count on from NSpecies
// 	for(const auto& zeroSpecies:fZeroYields){
// 	  //give 0 weight to species with 0 yield
// 	  eventW[iz++]=0;
// 	}

       

// 	if(gotID){//use ID from initial tree
// 	  auto* vars=fCurrDataSet->get(ev);
// 	  fWeights->FillWeights((Long64_t)vars->getRealValue(idname),eventW);
// 	} //ID not defined just use entry number in dataset
// 	else fWeights->FillWeights(ev,eventW);
//       }
//       if(fWeights->Size()!=fCurrDataSet->numEntries()){
// 	cout<<"sPlot::ExportWeights Done but mismatch between number of weights and number of data events"<<fWeights->Size()<<" "<<fCurrDataSet->sumEntries()<<" ..... exiting"<<endl;
// 	yields.Print("v");
// 	exit(0);
//       }
//     }
    
//     weights_uptr sPlot::MergeWeights(){
//       std::cout<<"sPlot::MergeWeights() "<<std::endl;
//       //in addition combine the weights into 1 and load them
//       weights_uptr wts(new Weights("HSsWeights"));
//       //Note the output file cannot contain the word Weights (because of Merge), hence Tweights!
//       wts->Merge(SetUp().GetOutDir()+"/Weights",
// 		 SetUp().GetOutDir()+"/"+SetUp().GetName()+"Tweights.root",
// 		 "HSsWeights");
//       //wts->Save();

//       //reset to save and reopen
//       wts.reset(new Weights{});
//       wts->LoadSaved(SetUp().GetOutDir()+TString("Tweights.root"),"HSsWeights");
      
//       return std::move(wts);
//     }

//     void sPlot::WeightedTree(){
//        if(!fWeights.get()){
// 	if(Bins().GetSize()>0)
// 	  fWeights = MergeWeights();
// 	else{
// 	  fWeights.reset(new Weights());
// 	  fWeights->LoadSaved(SetUp().GetOutDir()+TString("Weights")+SetUp().GetName()+".root","HSsWeights");
// 	}
//       }
//       TDirectory* saveDir=gDirectory;
      
//       //Open tree file
//       auto ftree=FiledTree::Read(Data().ParentTreeName(),
//       				 Data().ParentName());
//       //create a copy in a new file to append the weights to
//       //Keep it in file as large trees can use too much memory
//       fWeightedFiledTree=(FiledTree::CloneFull(ftree->Tree(),SetUp().GetOutDir()+"DataWeightedTree.root"));
  
//       //delete original
//       ftree.reset();
//       //Add weights to tree
//       fWeights->AddToTree(fWeightedFiledTree->Tree().get());
//       fWeightedFiledTree->Tree()->SetBranchStatus("*",true);

//       saveDir->cd();
//     }

//     void sPlot::DrawWeighted(const TString& var,const TString& wname,TString cut, const TString& opt){
//       if(!fWeightedFiledTree.get())
// 	WeightedTree();

//       if(!fWeightedFiledTree.get())
// 	return;

 
//       if(cut==TString()) cut="1";
     
//       fWeightedFiledTree->Tree()->Draw(var,wname+"*("+cut+")",opt);

//     }

//     Bool_t sPlot::ZeroYieldCheck(){
      
//       auto& yields=fCurrSetup->Yields(); //get reference to yields
//       auto& pdfs=fCurrSetup->PDFs(); //get reference to yields
//       bool removedPdf=false;    
//       for(Int_t iy=0;iy<yields.getSize();iy++){
// 	auto checkYield=dynamic_cast<RooRealVar*>(&yields[iy]);
// 	Double_t  thisYield=checkYield->getVal();
// 	if(thisYield<1E-2){
// 	  //Need to remove this pdf all weights weights will not be written for this species
// 	  Warning("sPlot::sPlot()",Form("Found zero yield for %s, will remove from sPlot, weights for this species will be set to 0",checkYield->GetName()),"");
// 	  pdfs.remove(pdfs[iy]);
// 	  //yields->remove(*(fCurrSetup->WS().var(yield0->GetName())));
// 	  yields.remove(yields[iy]);

// 	  removedPdf=true;
// 	  fZeroYields.push_back(checkYield->GetName());
// 	}

//       }
//       //remake sum of species without zero yied pdf
//       if(removedPdf)fCurrSetup->TotalPDF();

//       if(yields.getSize()==1){//Only 1 species all weights==inWeights
// 	fSingleYield=yields[0].GetName();
// 	CreateWeights();
// 	return kTRUE;
//       }
//       return kFALSE;
//     }

    
//   }//namespace FIT 
// }//namespace HS

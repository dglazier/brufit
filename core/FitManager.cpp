/**
 * @file FitManager.cpp
 * @brief Implementation of the FitManager class.
 * Refactored for C++17, removing manual memory management 
 * and improving ROOT I/O safety.
 */

#include "FitManager.h"

// --- NEW BRU NAMESPACE INCLUDES ---
#include "BruEventsPDF.h"
#include "BruEventsHistPDF.h"
//#include "BruComponentsPDF.h"
//#include "RooHSEventsPDF.h"
//#include "RooHSEventsHistPDF.h"
//#include "RooComponentsPDF.h"
#include "AmpMinuit2.h" // Ensure the default minimizer is included

#include <TSystem.h>
#include <TFile.h>
#include <TList.h>
#include <TObjString.h>

#include <iostream>
#include <memory>

namespace HS {
namespace FIT {

    // ========================================================================
    // Constructors & Assignment Operators
    // ========================================================================

    FitManager::FitManager(const FitManager& other) : TNamed(other.fName, other.fTitle) {
        fSetup = other.fSetup;
        fBinner = other.fBinner;
        fData = other.fData;
        fPlotOptions = other.fPlotOptions;
        fUsePrevResult = other.fUsePrevResult;
        fPrevResultDir = other.fPrevResultDir;
        fPrevResultMini = other.fPrevResultMini;
        fYldMaxFactor = other.fYldMaxFactor;
        fuseBinnedFit = other.fuseBinnedFit;
        
        // --- FIXED: Copy the missing persistent configuration state ---
        fMinimiserType = other.fMinimiserType;
        fCompiledMacros = other.fCompiledMacros;
        fTruthPrefix = other.fTruthPrefix;
        fDoPlotting = other.fDoPlotting;
        fRedirect = other.fRedirect;

        // --- FIXED: Safely Deep-Copy the Minimiser ---
        // If the master has a minimizer configured, clone its configuration 
        // to the worker node without sharing the memory pointer.
        if (other.fMinimiser) {
            // Assumes Minimiser has a Clone() method (via TObject/TNamed)
            fMinimiser.reset(dynamic_cast<Minimiser*>(other.fMinimiser->Clone()));
        }
    }

    FitManager& FitManager::operator=(const FitManager& other) {
        if (this != &other) {
            TNamed::operator=(other); 
            fSetup = other.fSetup;
            fBinner = other.fBinner;
            fData = other.fData;
            fPlotOptions = other.fPlotOptions;
            fUsePrevResult = other.fUsePrevResult;
            fPrevResultDir = other.fPrevResultDir;
            fPrevResultMini = other.fPrevResultMini;
            fYldMaxFactor = other.fYldMaxFactor;
            fuseBinnedFit = other.fuseBinnedFit;

            // --- FIXED: Copy missing state ---
            fMinimiserType = other.fMinimiserType;
            fCompiledMacros = other.fCompiledMacros;
            fTruthPrefix = other.fTruthPrefix;
            fDoPlotting = other.fDoPlotting;
            fRedirect = other.fRedirect;

            // --- FIXED: Safely Deep-Copy the Minimiser ---
            if (other.fMinimiser) {
                fMinimiser.reset(dynamic_cast<Minimiser*>(other.fMinimiser->Clone()));
            } else {
                fMinimiser.reset();
            }
        }
        return *this;
    }
  
    // ========================================================================
    // Core Execution Loop
    // ========================================================================

    Bool_t FitManager::Run() {
        fSetup.RequiredFitOptions();
        std::cout << " FitManager::Run()  "<<fFiti<<std::endl;
        CreateCurrSetup();
    
        // Get dataset for the current bin (fFiti)
        fCurrDataSet = std::move(Data().Get(fFiti));
        std::cout << " FitManager::Run()  "<< fCurrDataSet.get() <<std::endl;
 
        std::cout << "fCurrDataSet->numEntries() = " << fCurrDataSet->numEntries() << std::endl;
        
        // Use higher threshold to skip bins that will fail to fit
        if (fCurrDataSet->numEntries() < 10) {  
            std::cout << "WARNING FitManager::Run: <10 entries in dataset for this bin, moving to next..." << std::endl;
            return kFALSE;
        }
        if (fCurrDataSet->sumEntries() <= 0) {
            std::cout << "WARNING FitManager::Run: weighted entries <= 0 (" 
                      << fCurrDataSet->sumEntries() << "), moving to next..." << std::endl;
            return kFALSE;
        }
 
        // Look for Special case of BruEventsPDFs and initialize them
        FillEventsPDFs();
 
        // Add external fit constraints
        fCurrSetup->AddFitOption(RooFit::ExternalConstraints(fCurrSetup->Constraints()));
        
        // Initialise species yields based on dataset entries
        if (fCurrSetup->Yields().getSize() == 1) { 
            // Special case: only 1 yield
            Double_t yld = fCurrDataSet->sumEntries();
            SetAllValLimits(fCurrSetup->Yields(), yld, 0, fYldMaxFactor * yld);
        } else {
            SetAllValLimits(fCurrSetup->Yields(), 
                            fCurrDataSet->sumEntries() / 2, 0, 
                            fCurrDataSet->sumEntries() * fYldMaxFactor);
        }
        
        // Create extended maximum likelihood PDF
        fCurrSetup->TotalPDF();
        
        // Execute the minimization
        FitTo(); 

        return kTRUE;
    }

    void FitManager::CreateCurrSetup() {
        // Safely create a transient copy of the setup for this specific bin
        fCurrSetup = std::make_unique<Setup>(fSetup);
        fCurrSetup->SetName(GetCurrName());
        fCurrSetup->SetTitle(GetCurrTitle());
        
        // Ensure we take current setup values (ranges etc.)
        auto& currpy = fCurrSetup->ParsAndYields();
        currpy.assign(fSetup.ParsAndYields());
        
        // Range-based for-loop over the parameters
        for (auto* par : currpy) {
            auto* orig = dynamic_cast<RooRealVar*>(fSetup.ParsAndYields().find(par->GetName()));
            if (orig != nullptr) {
                dynamic_cast<RooRealVar*>(par)->setMin(orig->getMin());
                dynamic_cast<RooRealVar*>(par)->setMax(orig->getMax());
            }
        }

        // Seed with previous results if requested
        if (fUsePrevResult) {
            LoadPrevResult(fPrevResultDir, fPrevResultMini);
        }
    }

    void FitManager::RunAll() {
        PreRun();
        UInt_t Nf = GetN();
        for (UInt_t i = 0; i < Nf; i++) {
            RunOne(i);
        }
    }

    void FitManager::RunOne(Int_t ifit) {
        fFiti = ifit;
        if (fRedirect) RedirectOutput(fSetup.GetOutDir() + Form("logRooFit%d.txt", fFiti));
        
        auto success = Run();
        
        if (fRedirect) RedirectOutput();
        if (success) SaveResults();
        
        Reset();
    }
    
    void FitManager::FitTo() {
        // Fallback to Minuit2 if no minimizer was specified
        if (!fMinimiser) SetMinimiser(std::make_unique<HS::FIT::Minuit2>());

        if (fuseBinnedFit == kFALSE) {
            fMinimiser->Run(*fCurrSetup, *fCurrDataSet);
        } else {
            // binnedClone() is a RooDataSet method. We must downcast first.
            if (auto* unbinnedData = dynamic_cast<RooDataSet*>(fCurrDataSet.get())) {
                std::unique_ptr<RooAbsData> binnedData(unbinnedData->binnedClone());
                fMinimiser->Run(*fCurrSetup, *binnedData);
            } else {
                // If the cast fails, it's likely already a RooDataHist (binned data)
                fMinimiser->Run(*fCurrSetup, *fCurrDataSet);
            }
        }
        
        // Plot best fit and return
        if (fDoPlotting) PlotDataModel();
    }

    // ========================================================================
    // Special PDF Handling & Plotting
    // ========================================================================

    void FitManager::FillEventsPDFs() {
        UInt_t idata = GetDataBin(fFiti);
        auto& pdfs = fCurrSetup->PDFs();
        auto savedir = gDirectory;
        
        for (Int_t ip = 0; ip < pdfs.getSize(); ip++) {
            
            // --- UPDATED: Use bru::BruEventsPDF ---
            auto pdf = dynamic_cast<bru::BruEventsPDF*>(&pdfs[ip]);
	  // auto pdf = dynamic_cast<RooHSEventsPDF*>(&pdfs[ip]);
            
            if (pdf != nullptr) {
                // Set truth prefix for MC
                pdf->SetTruthPrefix(fTruthPrefix);
                
                if (fBinner.FileNames(pdf->GetName()).empty()) continue;
                
                // Open tree files for getting events
                auto filetree = FiledTree::Read(
                    fBinner.TreeName(pdf->GetName()),
                    fBinner.FileNames(pdf->GetName())[idata]
                );
                auto tree = filetree->Tree();
    
                auto mcgenfiletree = (fBinner.FileNames(pdf->GetName() + TString("__MCGen")).empty() ? 
                    nullptr : 
                    FiledTree::Read(fBinner.TreeName(pdf->GetName() + TString("__MCGen")), fBinner.FileNames(pdf->GetName() + TString("__MCGen"))[idata])
                );
                auto mcgentree = (mcgenfiletree ? mcgenfiletree->Tree() : nullptr);
      
                savedir->cd();
                
                if (!tree.get()) {
                    std::cout << "WARNING FitManager::FillEventsPDFs: No tree data found for EventPDF " 
                              << pdf->GetName() << std::endl;
                    continue;
                }
                
                // If too few events, remove this PDF from the fit
                if (!tree->GetEntries() || !pdf->IsValid()) {
                    std::cout << "WARNING FitManager::FillEventsPDFs: Too few events for EventPDF " 
                              << pdf->GetName() << std::endl;
                    fCurrSetup->Yields().remove(fCurrSetup->Yields()[ip]);
                    pdfs.remove(*pdf);
                    ip--;
                } else { 
                    // Use it and pass the simulated tree
                    // (SetEvTree automatically extracts to MCEventCache under the hood)
                    pdf->SetInWeights(fCurrSetup->GetPDFInWeights(pdf->GetName()));
                    pdf->SetEvTree(tree.get(), fCurrSetup->Cut(), mcgentree.get());

                    // See if data to load for proto data
                    if (!fCurrDataSet) fCurrDataSet = std::move(Data().Get(idata));
                    
                    if (fCurrDataSet) {
                        if (auto* unbinnedData = dynamic_cast<RooDataSet*>(fCurrDataSet.get())) {
                            pdf->AddProtoData(unbinnedData);
                        } else {
                            std::cerr << "WARNING: AddProtoData requires an unbinned RooDataSet." << std::endl;
                        }
                    }
                    
                    // --- UPDATED: Use bru::BruEventsHistPDF ---
                    bru::BruEventsHistPDF* histspdf = nullptr;
                    if ((histspdf = dynamic_cast<bru::BruEventsHistPDF*>(pdf))) {
		      //RooHSEventsHistPDF* histspdf = nullptr;
		      //if ((histspdf = dynamic_cast<RooHSEventsHistPDF*>(pdf))) {
                        histspdf->CreateHistPdf();
                        fCurrSetup->AddGausConstraint(histspdf->AlphaConstraint());
                        fCurrSetup->AddGausConstraint(histspdf->OffConstraint());
                        fCurrSetup->AddGausConstraint(histspdf->ScaleConstraint());
                    }

                    pdf->MakeAssertPostiveData();
                    pdf->AssertPositivePDF(); // cache it
                }
                
                // Keep the simulated tree alive until Reset() is called
                fFiledTrees.push_back(std::move(filetree));    
                if (mcgenfiletree) {
                    fFiledTrees.push_back(std::move(mcgenfiletree));    
                }
            }
        }
        savedir->cd();
    }

    void FitManager::PlotDataModel() {
        // Cast the dataset. If it's a binned RooDataHist, this becomes nullptr,
        // so ensure your Plot classes check for null!
        auto* dataset = dynamic_cast<RooDataSet*>(fCurrDataSet.get());

        if (auto* mcmc = dynamic_cast<RooMcmc*>(fMinimiser.get())) {
            fPlots.push_back(std::make_unique<MCMCPlotResults>(
                fCurrSetup.get(), dataset, GetCurrName() + GetCurrTitle(), mcmc, fPlotOptions));
                
        } else if (auto* bruMcmc = dynamic_cast<BruMcmc*>(fMinimiser.get())) {
            fPlots.push_back(std::make_unique<MCMCPlotResults>(
                fCurrSetup.get(), dataset, GetCurrName() + GetCurrTitle(), bruMcmc, fPlotOptions));
                
        } else {
            fPlots.push_back(std::make_unique<PlotResults>(
                fCurrSetup.get(), dataset, GetCurrName() + GetCurrTitle(), fPlotOptions));
        }
    }
  
    // ========================================================================
    // File I/O and Persistence (RAII Protected)
    // ========================================================================

    void FitManager::SaveSetup() {
        // unique_ptr ensures the file is closed correctly even if exceptions are thrown
        std::unique_ptr<TFile> file(TFile::Open(fSetup.GetOutDir() + "HSSetup.root", "recreate"));
        if (file && !file->IsZombie()) {
            fSetup.Write("HSSetup");
        }
    }

    void FitManager::LoadSetup(const TString& dir) {
        std::unique_ptr<TFile> file(TFile::Open(dir + "/HSSetup.root"));
        if (file && !file->IsZombie()) {
            auto* setupPtr = dynamic_cast<HS::FIT::Setup*>(file->Get("HSSetup"));
            if (setupPtr) fSetup = *setupPtr;
        }
    }

    void FitManager::InitPrevResult(const TString& resultDir, const TString& resultMinimiser) {
        fUsePrevResult = kTRUE;
        fMinimiserType = resultMinimiser;

        if (resultDir == TString("")) fPrevResultDir = fSetup.GetOutDir(); 
        else fPrevResultDir = resultDir;

        if (resultMinimiser == TString("")) fPrevResultMini = fMinimiser->GetName();
        else fPrevResultMini = resultMinimiser;
    }
    
    void FitManager::LoadPrevResult(const TString& resultDir, const TString& resultMinimiser) {
        TString resultFile = resultDir + "/" + fCurrSetup->GetName() + "/Results" + resultMinimiser + ".root";
        std::cout << "FitManager::LoadPrevResult open file " << resultFile << std::endl;
        
        std::unique_ptr<TFile> fitFile(TFile::Open(resultFile));
        if (!fitFile || fitFile->IsZombie()) return;

        // Retrieve the result dataset safely
        std::unique_ptr<RooDataSet> result(dynamic_cast<RooDataSet*>(fitFile->Get(Minimiser::FinalParName())));
        
        if (result) {
            auto newPars = fCurrSetup->ParsAndYields();
            auto* resAll = result->get(); // get all result info
            auto* resPars = resAll->selectCommon(newPars); // just select pars and yields
            
            newPars.assign(*resPars); // set values to results
            
            std::cout << "FitManager::LoadResult setting values from fit results " << resultFile << " : " << std::endl;
            newPars.Print("v");
            
            delete resPars; // Clean up the caller-owned set returned by selectCommon
        }
    }

    void FitManager::WriteThis() {
        std::cout << "FitManager::WriteThis() to " << fSetup.GetOutDir() + "HSFit.root" << std::endl;
        
        std::unique_ptr<TFile> file(TFile::Open(fSetup.GetOutDir() + "HSFit.root", "recreate"));
        if (!file || file->IsZombie()) return;

        if (!fMinimiser) SetMinimiser(std::make_unique<HS::FIT::Minuit2>());
        
        file->WriteObject(this, "HSFit");
        file->WriteObject(fMinimiser.get(), fMinimiserType);
  
        if (!fCompiledMacros.empty()) {
            // Memory safe list construction for file writing
            auto macList = std::make_unique<TList>();
            macList->SetOwner(kTRUE);
            for (const auto& macro : fCompiledMacros) {
                macList->Add(new TObjString(macro));
            }
            file->WriteObject(macList.get(), "HS_COMPILEDMACROS");
        }
    }

    void FitManager::RedirectOutput(const TString& log) {
      //const char* mess = Form("text output will be sent to file %s", log.Data());
      // std::cout << "FitManager::RedirectOutput " << mess << std::endl;
        
        if (log == TString("")) {
            gSystem->RedirectOutput(nullptr, "w");
        } else {
            gSystem->RedirectOutput(log.Data(), "w");
        }
    }

    void FitManager::SaveResults() {
        auto saveDir = gDirectory;
        // fMinimiser->SaveInfo() implicitly creates/modifies files, ensure it manages its own handles
        auto outFile = fMinimiser->SaveInfo(); 
        
        if (!fPlots.empty() && fPlots.back()) {
            fPlots.back()->Write(); // Just save the last one
        }
        
        saveDir->cd();
    }

} // namespace FIT
} // namespace HS

// #include "FitManager.h"

// #include <memory>
// #include "RooHSEventsPDF.h"
// #include "RooHSEventsHistPDF.h"
// #include "RooComponentsPDF.h"
// #include "TSystem.h"


// namespace HS{
//   namespace FIT{

 
//     FitManager::FitManager(const FitManager& other):TNamed(other.fName,other.fName){
//       fSetup=other.fSetup;
//       fBinner=other.fBinner;
//       //LoadData(other.GetDataTreeName(),other.GetDataFileNames());
//       fPlotOptions=other.fPlotOptions;
//       fUsePrevResult=other.fUsePrevResult;
//       fPrevResultDir=other.fPrevResultDir;
//       fPrevResultMini=other.fPrevResultMini;
//       fYldMaxFactor=other.fYldMaxFactor;
//       fuseBinnedFit=other.fuseBinnedFit;
      
//       //fIsSamplingIntegrals=other.fIsSamplingIntegrals;
//     }

//     FitManager&  FitManager::operator=(const FitManager& other){
//       cout<<"=============FitManager"<<endl;
//       fSetup=other.fSetup;
//       fBinner=other.fBinner;
//       fPlotOptions=other.fPlotOptions;
//       fUsePrevResult=other.fUsePrevResult;
//       fPrevResultDir=other.fPrevResultDir;
//       fPrevResultMini=other.fPrevResultMini;
//       fYldMaxFactor=other.fYldMaxFactor;
//       fuseBinnedFit=other.fuseBinnedFit;
//       //fIsSamplingIntegrals=other.fIsSamplingIntegrals;
  
//       return *this;
//     }
    
//     Bool_t FitManager::Run(){
//       fSetup.RequiredFitOptions();
      
//       CreateCurrSetup();
     
//       //get dataset fFiti
//       fCurrDataSet=std::move(Data().Get(fFiti));

//       cout<<"fCurrDataSet->numEntries() = "<<fCurrDataSet->numEntries()<<endl;
//       // if(fCurrDataSet->numEntries()==0){
//       if(fCurrDataSet->numEntries()<10){  // use higher threshold to remove bins that will fail for sure
// 	cout<<"WARNING FitManager::Run no entries in dataset for this bin will move to next...."<<endl;
// 	return kFALSE;
//       }
//       if(fCurrDataSet->sumEntries()<=0){
// 	cout<<"WARNING FitManager::Run weighted entries <=0 actually, "<<fCurrDataSet->sumEntries()<<" in dataset for this bin will move to next...."<<endl;
// 	return kFALSE;
//       }
 
//       //Look for Special case of RooHSEventsPDFs
//       FillEventsPDFs();
 
//       //Add fit constraints
//       fCurrSetup->AddFitOption(RooFit::ExternalConstraints
//       			       (fCurrSetup->Constraints()));
      
//       //initialise yields
//       if(fCurrSetup->Yields().getSize()==1){//special case only 1 yield)
// 	Double_t yld=fCurrDataSet->sumEntries();
// 	SetAllValLimits(fCurrSetup->Yields(),
// 			yld,0,fYldMaxFactor*yld);
//       }
//       else{
// 	SetAllValLimits(fCurrSetup->Yields(),
// 			fCurrDataSet->sumEntries()/2,0,
// 			fCurrDataSet->sumEntries()*fYldMaxFactor);
//       }
      
//       //create extended max likelihood pdf
//       //std::cout<<"DEBUG FitManager::Run()"<<std::endl;fCurrSetup->Parameters().Print("v");
//       fCurrSetup->TotalPDF();
//       FitTo(); 

//       return kTRUE;
//     }
//     void FitManager::CreateCurrSetup(){
//       fCurrSetup = std::unique_ptr<Setup>(new Setup{fSetup}); //Copy setup from template
//       fCurrSetup->SetName(GetCurrName());
//       fCurrSetup->SetTitle(GetCurrTitle());
//       //make sure we take current setup values
//       //If not it will use the string from Factory() etc,
//       auto& currpy = fCurrSetup->ParsAndYields();
//       currpy.assign(fSetup.ParsAndYields());
//       for(auto& par:currpy){//assignFast doesnt do ranges...
// 	auto* orig=dynamic_cast<RooRealVar*>(fSetup.ParsAndYields().find(par->GetName()));
// 	if(orig!=nullptr){
// 	  dynamic_cast<RooRealVar*>(par)->setMin(orig->getMin());
// 	  dynamic_cast<RooRealVar*>(par)->setMax(orig->getMax());
// 	}
//       }

//       //Look to see if taking previous fit results as initial pars
//       if(fUsePrevResult){
// 	LoadPrevResult(fPrevResultDir,fPrevResultMini);
//       }
//     }

//     /////////////////////////////////////////////////////////////
//     void FitManager::RunAll(){

//       PreRun();

//       UInt_t Nf=GetN();
//       for(UInt_t i=0;i<Nf;i++){
// 	 RunOne(i);
//       }
//     }

//     ////////////////////////////////////////////////////////////
//     void FitManager::FitTo(){
//       //      std::cout<<"DEBUG FitManager::FitTo()"<<std::endl;
//       if(!fMinimiser.get()) SetMinimiser(new HS::FIT::Minuit2());

//       if(fuseBinnedFit==kFALSE){
// 	fMinimiser->Run(*fCurrSetup,*fCurrDataSet);
//       }
//       else{
// 	auto binnedData = fCurrDataSet->binnedClone();
// 	fMinimiser->Run(*fCurrSetup,*binnedData);
// 	delete binnedData;
	
//       }
//       ///////////////////////////
//       //Plot best fit and return
//       if(fDoPlotting) PlotDataModel();

//     }
//     void FitManager::RunOne(Int_t ifit){
//       fFiti=ifit;
//       if(fRedirect) RedirectOutput(fSetup.GetOutDir()+Form("logRooFit%d.txt",fFiti));
//       auto success=Run();
//       if(fRedirect) RedirectOutput();

//       if(success)SaveResults();
      
//       Reset();
//     }
    
//     void FitManager::FillEventsPDFs(){
//       UInt_t idata=GetDataBin(fFiti);
//       auto& pdfs=fCurrSetup->PDFs();
 
//       auto savedir=gDirectory;
      
//       for(Int_t ip=0;ip<pdfs.getSize();ip++){
// 	auto pdf=dynamic_cast<RooHSEventsPDF*>( &pdfs[ip]);
// 	if(pdf!=nullptr){
// 	  //SetTruthprefix
// 	  pdf->SetTruthPrefix(fTruthPrefix);
// 	  if(fBinner.FileNames(pdf->GetName()).size()==0)
// 	    continue;
// 	  //Open tree files for getting events
// 	  auto filetree=FiledTree::
// 	    Read(fBinner.TreeName(pdf->GetName()),
// 		 fBinner.FileNames(pdf->GetName())[idata]);
// 	  auto tree=filetree->Tree();
	
// 	  auto mcgenfiletree= (fBinner.FileNames(pdf->GetName()+TString("__MCGen")).empty() ? nullptr : FiledTree::Read(fBinner.TreeName(pdf->GetName()+TString("__MCGen")),fBinner.FileNames(pdf->GetName()+TString("__MCGen"))[idata]));
// 	  auto mcgentree=(mcgenfiletree ? mcgenfiletree->Tree() : nullptr);
	  
// 	  savedir->cd();
//  	  if(!tree.get()){
// 	    cout<<"WARNING FitManager::FillEventsPDFs :"<<
// 	      "    No tree data found for EventPDF "<<pdf->GetName()<<endl;
// 	    continue;
// 	  }
// 	  //if too few events remove this PDF
// 	  if(!tree->GetEntries()||!pdf->IsValid()){
// 	    cout<<"WARNING FitManager::FillEventsPDFs :"<<
// 	      "    too few events for for EventPDF "<<pdf->GetName()<<endl;
// 	    fCurrSetup->Yields().remove(fCurrSetup->Yields()[ip]);
// 	    pdfs.remove(*pdf);
// 	    ip--;
// 	  }
// 	  else{ //use it and give it the simulated tree
	    
// 	    pdf->SetInWeights(fCurrSetup->GetPDFInWeights(pdf->GetName()));
// 	    pdf->SetEvTree(tree.get(),fCurrSetup->Cut(),mcgentree.get());

// 	    //See if data to load for proto data
// 	    if(!fCurrDataSet.get())
// 	      fCurrDataSet=std::move(Data().Get(idata));

// 	    if(fCurrDataSet.get())pdf->AddProtoData(fCurrDataSet.get());
// 	    RooHSEventsHistPDF* histspdf=nullptr;
// 	    if((histspdf=dynamic_cast<RooHSEventsHistPDF*>(pdf))){
// 	      histspdf->CreateHistPdf();
// 	      fCurrSetup->AddGausConstraint(histspdf->AlphaConstraint());
// 	      fCurrSetup->AddGausConstraint(histspdf->OffConstraint());
// 	      fCurrSetup->AddGausConstraint(histspdf->ScaleConstraint());
// 	    }
// 	    //cout<<"FitManager IsSamplingIntegrals "<<fIsSamplingIntegrals<<" "<<endl;
// 	    //if sampling PDF create constraint for fit
// 	    //if(fIsSamplingIntegrals==kTRUE){
// 	    // pdf->SetIsSamplingIntegral();
// 	      /////fCurrSetup->AddGausConstraint(pdf->GetIntegralPDF()->getPDF());
// 	    //}

// 	    pdf->MakeAssertPostiveData();
// 	    pdf->AssertPositivePDF();//cache it

// 	  }
// 	  //keep the simulated tree alive until Reset()
// 	  fFiledTrees.push_back(std::move(filetree));	
// 	  if(mcgenfiletree)
// 	    fFiledTrees.push_back(std::move(mcgenfiletree));	  
// 	}
//       }
//       savedir->cd();
//     }
//     void FitManager::SaveSetup(){
//       auto file=TFile::Open(fSetup.GetOutDir()+"HSSetup.root","recreate");
//       fSetup.Write("HSSetup");
//       delete file;
//     }
//     void FitManager::LoadSetup(const TString& dir){
//       auto file=TFile::Open(dir+"/HSSetup.root");
//       fSetup=*(dynamic_cast<HS::FIT::Setup*>(file->Get("HSSetup")));
//       delete file;
//     }
//     //Read in paarameters from previous fit
//     void FitManager::InitPrevResult(const TString& resultDir,const TString& resultMinimiser){
//       fUsePrevResult=kTRUE;
//       fMinimiserType=resultMinimiser;

//       if(resultDir==TString()) fPrevResultDir=fSetup.GetOutDir(); //use current
//       else fPrevResultDir=resultDir;

//       if(resultMinimiser==TString())fPrevResultMini=fMinimiser->GetName();
//       else fPrevResultMini=resultMinimiser;
      
//     }
    
//     void FitManager::LoadPrevResult(const TString& resultDir,const TString& resultMinimiser){
//       //TString resultFile=resultDir+"/"+fCurrSetup->GetName()+"/Results"+fCurrSetup->GetTitle()+resultMinimiser+".root";
//       TString resultFile=resultDir+"/"+fCurrSetup->GetName()+"/Results"+resultMinimiser+".root";

//       cout<<" FitManager::LoadPrevResult open file "<<resultFile<<endl;
//       std::unique_ptr<TFile> fitFile{TFile::Open(resultFile)};
//       std::unique_ptr<RooDataSet> result{dynamic_cast<RooDataSet*>( fitFile->Get(Minimiser::FinalParName()))};
//       //fitFile.reset();
//       //      auto result=dynamic_cast<RooDataSet*>( fitFile->Get(Minimiser::FinalParName())->Clone());//**
//       //Set the values of the paramteres to those in the given result
//       if(result.get()){
//       //if(result){
// 	auto newPars = fCurrSetup->ParsAndYields();
// 	auto* resAll = result->get(); //get all result info
// 	auto* resPars=resAll->selectCommon(newPars); //just select pars and yieds
//        	newPars.assign(*resPars); //set values to results
// 	cout<<"FitManager::LoadResult setting values from fit results "<<resultFile<<" : "<<endl;
// 	newPars.Print("v");
// 	//	delete result;result=nullptr;
//       }
//     }
//     void FitManager::WriteThis(){
//       cout<<"FitManager::WriteThis() to "<<fSetup.GetOutDir()+"HSFit.root"<<endl;
//       auto file=TFile::Open(fSetup.GetOutDir()+"HSFit.root","recreate");
//       if(!fMinimiser.get()) SetMinimiser(new HS::FIT::Minuit2());
//       file->WriteObject(this,"HSFit");
//       file->WriteObject(fMinimiser.get(),fMinimiserType);
  
//       if(fCompiledMacros.size()){
// 	auto* macList=new TList();
// 	//	macList->SetName("HS_COMPILEDMACROS");
// 	macList->SetOwner();
// 	for(auto& macro : fCompiledMacros)
// 	  macList->Add(new TObjString(macro));
// 	file->WriteObject(macList,"HS_COMPILEDMACROS");
//       }
//       delete file;
//     }
//     void FitManager::RedirectOutput(const TString& log){
//       const char* mess=Form("text ouput will be sent to file %s",log.Data());
//       cout<<"FitManager::RedirectOutput "<<mess<<endl;
//       if(log==TString(""))
// 	gSystem->RedirectOutput(nullptr,"w");
//       else
// 	gSystem->RedirectOutput(log.Data(),"w");
      
//     }

//     void FitManager::SaveResults(){
//       auto saveDir = gDirectory;
//       auto outFile=fMinimiser->SaveInfo();
//       if(fPlots.size())fPlots.back()->Write(); //just save the last one
//       saveDir->cd();
//     }


//   }//namespace FIT
// }//namespace HS

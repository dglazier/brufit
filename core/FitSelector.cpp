#define FitSelector_cxx
// The class definition in FitSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("FitSelector.C")
// root> T->Process("FitSelector.C","some options")
// root> T->Process("FitSelector.C+")
//


#include "FitSelector.h"
#include "FitManager.h"
#include <TClassTable.h>
#include <TH2.h>
#include <TStyle.h>
#include <chrono>
#include <thread>
#include <iostream>

// #include "Data.h"
// #include "Binner.h"
// #include "Setup.h"
// #include "FiledTree.h"
namespace HS{
  namespace FIT{
    void FitSelector::Begin(TTree * /*tree*/)
    {
      // The Begin() function is called at the start of the query.
      // When running with PROOF Begin() is only called on the client.
      // The tree argument is deprecated (on PROOF 0 is passed).

      TString option = GetOption();

       
      if(!fInput) fInput=new TList();

      auto outdir=dynamic_cast<TNamed*>(fInput->FindObject("HSOUTDIR"));
      if(outdir) {
        // reuse already existing entry for output directory
        outdir->SetTitle(fFitManager->SetUp().GetOutDir().Data());
      } else {
        // create new entry for output directory
        outdir=new TNamed("HSOUTDIR",fFitManager->SetUp().GetOutDir().Data());
        fInput->Add(outdir);
      }

      fFitManager->PreRun();
    }

    void FitSelector::SlaveBegin(TTree * /*tree*/)
    {
      //PROOF does not load the namespaces properly when loading classes
      //from shard library
      TClassTable::AddAlternate("HS::FIT::RooHSEventsPDF","RooHSEventsPDF");
      TClassTable::AddAlternate("HS::FIT::RooHSEventsHistPDF","RooHSEventsHistPDF");
      TClassTable::AddAlternate("HS::FIT::BruEventsHistPeakPDF","BruEventsHistPeakPDF");
      TClassTable::AddAlternate("HS::FIT::RooComponentsPDF","RooComponentsPDF");
      TClassTable::AddAlternate("HS::FIT::RooHSSphHarmonic","RooHSSphHarmonic");
      TClassTable::AddAlternate("HS::FIT::RooHSSphHarmonicIm","RooHSSphHarmonicIm");
      TClassTable::AddAlternate("HS::FIT::RooHSSphHarmonicRe","RooHSSphHarmonicRe");
      TClassTable::AddAlternate("HS::FIT::RelBreitWigner","RelBreitWigner");
      TClassTable::AddAlternate("HS::FIT::RooHSComplexSumSqdTerm","RooHSComplexSumSqdTerm");

      TClassTable::AddAlternate("HS::FIT::RooHSDWigner","RooHSDWigner");
      TClassTable::AddAlternate("HS::FIT::RooHSDWignerIm","RooHSDWignerIm");
      TClassTable::AddAlternate("HS::FIT::RooHSDWignerRe","RooHSDWignerRe");
       TClassTable::AddAlternate("HS::FIT::RooHSDWignerProduct","RooHSDWignerProduct");
      TClassTable::AddAlternate("HS::FIT::RooHSDWignerProductIm","RooHSDWignerProductIm");
      TClassTable::AddAlternate("HS::FIT::RooHSDWignerProductRe","RooHSDWignerProductRe");
  
      // TClassTable::AddAlternate("HS::FIT::","");

      TString option = GetOption();
      fInput->Print();

      //Get the outpur directory where HSFit.root should reside
      auto outdir=dynamic_cast<TNamed*>(fInput->FindObject("HSOUTDIR"));
      TString outdirstr=TString(outdir->GetTitle());
         if(outdirstr!=TString(""))
	outdirstr.Append('/');

      //Get the fitmanager
      fFitfile=TFile::Open(outdirstr+"HSFit.root");
      
      fFitManager=dynamic_cast<FitManager*>( fFitfile->Get("HSFit")->Clone() );

      if(fFitManager->GetMinimiserType()!=TString())
	fFitManager->SetMinimiser(dynamic_cast<Minimiser*>( fFitfile->Get(fFitManager->GetMinimiserType())->Clone() ));
  
      delete fFitfile; //close file to stop memeory resident issue!
      
      fFitManager->Data().LoadSetup(&fFitManager->SetUp());
      cout<<"FitSelector::SlaveBegin( "<<fFitManager->SetUp().GetOutDir()<<endl;
    }

    Bool_t FitSelector::Process(Long64_t entry)
    {
      cout<<"FitSelector::Process Run entry "<<entry<<endl; 
      fFitManager->RunOne(entry);

  
      return kTRUE;
    }

    void FitSelector::SlaveTerminate()
    {
      // The SlaveTerminate() function is called after all entries or objects
      // have been processed. When running with PROOF SlaveTerminate() is called
      // on each slave server.
    }

    void FitSelector::Terminate()
    {
      // The Terminate() function is the last function to be called during
      // a query. It always runs on the client, it can be used to present
      // the results graphically or save the results to file.

      auto outdir=dynamic_cast<TNamed*>(fInput->FindObject("HSOUTDIR"));
      if (outdir)
        fInput->Remove(outdir);

      //HERE we should gather information for slaves if we want to
      //pass to another fit manager, i.e. produce toys in proof
      //then pass file names to fitter
    }
    void FitSelector::Init(TTree *tree)
    {
      // The Init() function is called when the selector needs to initialize
      // a new tree or chain. Typically here the reader is initialized.
      // It is normally not necessary to make changes to the generated
      // code, but the routine can be extended by the user if needed.
      // Init() will be called many times when running on PROOF
      // (once per file to be processed).

      // fReader.SetTree(tree);
    }

    Bool_t FitSelector::Notify()
    {
      // The Notify() function is called when a new file is opened. This
      // can be either for a new TTree in a TChain or when when a new TTree
      // is started when using PROOF. It is normally not necessary to make changes
      // to the generated code, but the routine can be extended by the
      // user if needed. The return value is currently not used.

      return kTRUE;
    }

  }
}

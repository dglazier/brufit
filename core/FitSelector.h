//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 25 13:34:24 2019 by ROOT version 6.14/04
// from TTree bins/A FiledTree
// found on file: BinIndices.root
//////////////////////////////////////////////////////////

#ifndef FitSelector_h
#define FitSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
/* #include <TTreeReader.h> */
/* #include <TTreeReaderValue.h> */
/* #include <TTreeReaderArray.h> */
#include <TString.h>
#include <vector>

// Headers needed by this particular selector
#include "FitManager.h"

namespace HS{
  namespace FIT{

    class FitManager;
    
    class FitSelector : public TSelector {
      public :
      /* TTreeReader     fReader;  //!the tree reader */
      /* TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain */

      /* // Readers to access the data (delete the ones you do not need). */
      /* TTreeReaderValue<Int_t> bindex = {fReader, "bindex"}; */


      FitSelector(TTree * /*tree*/ =nullptr) { }
      ~FitSelector() override = default;
      Int_t   Version() const override { return 2; }
      void    Begin(TTree *tree) override;
      void    SlaveBegin(TTree *tree) override;
      void    Init(TTree *tree) override;
      Bool_t  Notify() override;
      Bool_t  Process(Long64_t entry) override;
      //virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEnqtry(entry, getall) : 0; }
      void    SetOption(const char *option) override { fOption = option; }
      void    SetObject(TObject *obj) override { fObject = obj; }
      void    SetInputList(TList *input) override { fInput = input; }
      TList  *GetOutputList() const override { return fOutput; }
      void    SlaveTerminate() override;
      void    Terminate() override;


      void SetFitManager(FitManager *fm){fFitManager=fm;}
    private:
      std::vector<TString> fFitFileNames;

      FitManager *fFitManager{};
      TFile* fFitfile=nullptr;
      
      ClassDefOverride(HS::FIT::FitSelector,0);

    };
  }
}

#endif // #ifdef FitSelector_cxx

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 30 09:47:45 2023 by ROOT version 6.26/04
// from TTree ResultTree/ResultTree
// found on file: Malte1Fits4Waves/ResultsHSAmpMinuit2.root
//////////////////////////////////////////////////////////

#ifndef PlotAmbigs_h
#define PlotAmbigs_h

#include <TROOT.h>
#include <vector>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class PlotAmbigs : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Double_t> NLL = {fReader, "NLL"};
   TTreeReaderValue<Double_t> Yld_PWA = {fReader, "Yld_PWA"};
   TTreeReaderValue<Double_t> Yld_PWA_err = {fReader, "Yld_PWA_err"};
   TTreeReaderValue<Double_t> a_0_0 = {fReader, "a_0_0"};
   TTreeReaderValue<Double_t> a_0_0_err = {fReader, "a_0_0_err"};
   TTreeReaderValue<Double_t> a_2_m1 = {fReader, "a_2_-1"};
   TTreeReaderValue<Double_t> a_2_m1_err = {fReader, "a_2_-1_err"};
   TTreeReaderValue<Double_t> a_2_0 = {fReader, "a_2_0"};
   TTreeReaderValue<Double_t> a_2_0_err = {fReader, "a_2_0_err"};
  //TTreeReaderValue<Double_t> a_2_1 = {fReader, "a_2_1"};
  // TTreeReaderValue<Double_t> a_2_1_err = {fReader, "a_2_1_err"};
   TTreeReaderValue<Double_t> aphi_0_0 = {fReader, "aphi_0_0"};
   TTreeReaderValue<Double_t> aphi_2_m1 = {fReader, "aphi_2_-1"};
   TTreeReaderValue<Double_t> aphi_2_m1_err = {fReader, "aphi_2_-1_err"};
   TTreeReaderValue<Double_t> aphi_2_0 = {fReader, "aphi_2_0"};
   TTreeReaderValue<Double_t> aphi_2_0_err = {fReader, "aphi_2_0_err"};
   TTreeReaderValue<Double_t> aphi_2_1 = {fReader, "aphi_2_1"};
   TTreeReaderValue<Double_t> aphi_2_1_err = {fReader, "aphi_2_1_err"};
   TTreeReaderValue<Double_t> fitstatus = {fReader, "status"};
  
   PlotAmbigs(TTree * /*tree*/ =0) { }
   virtual ~PlotAmbigs() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

  TMultiGraph* _multiMags = {new TMultiGraph()};
  TMultiGraph* _multiPhases= {new TMultiGraph()};
  std::vector<TGraphErrors> _graphMags=std::vector<TGraphErrors>(4);
  std::vector<TGraphErrors> _graphPhases=std::vector<TGraphErrors>(4);
  Int_t _ipoint=0;
  
  ClassDef(PlotAmbigs,0);

};

#endif

#ifdef PlotAmbigs_cxx
void PlotAmbigs::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t PlotAmbigs::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef PlotAmbigs_cxx

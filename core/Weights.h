#pragma once

#include "FiledTree.h"
#include <TTree.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TList.h>
#include <TVectorD.h>
#include <ROOT/RDataFrame.hxx>
#include <map>
#include <utility>
#include <utility>
#include <vector>
#include <iostream>
 

namespace HS{
  namespace FIT{
    
    using std::vector;
    using std::map;
    using std::pair;
    using std::cout;
    using std::endl;
    
    using  StrIntMap_t = map<TString, Int_t >;
    using DF_uptr=std::unique_ptr<ROOT::RDataFrame>;

    
    class Weights : public TNamed{
    
    public:
      Weights() =default;
      Weights(const TString& name);
      ~Weights() override;
    
      TTree* GetIDTree(){return fIDTree;};
      void SetIDTree(TTree* tree){fIDTree=tree;}
      TTree* GetTree(){return fWTree;};
      void SetTree(TTree* tree){fWTree=tree;}
      void FillWeights(Long64_t ev,const TVectorD& wgt){ fID=ev; fWVals=wgt; fWTree->Fill();fIDTree->Fill();fN++;}
      void FillWeight(Long64_t ev,Double_t wgt){if(GetNSpecies()==1){ fID=ev; fWVals[0]=wgt; fWTree->Fill();fIDTree->Fill();fN++;}}//Special case of single species!!!!
    
      void GetEntry(Long64_t ent){fWTree->GetEntry(ent);fIDTree->GetEntry(ent);}; 
      Bool_t GetEntryFast(Long64_t id); //use id branch with sorted tree
      Bool_t GetEntrySlow(Long64_t id); //use id branch
      Bool_t GetEntryBinarySearch(Long64_t id); //use binary search to give faster entrys when used by unsorted trees
      Double_t GetWeight(const TString& spe){if(fSpecies.count(spe))return GetWeight(fSpecies[spe]);return 1;}
      Double_t GetWeight(Int_t ispe){
	if(fGotEntry)
	  return fWVals[ispe];
	else return 0; //entry for id not found
      }
      Long64_t GetID(){return fID;}
      Long64_t* GetIDi(){return fIDi;};
      void SetCurrEntry(Long64_t ent){fCurrEntry=ent;}
      Bool_t GotEntry(){return fGotEntry;}
      Bool_t IsSorted(){return fIsSorted;}
      Long64_t GetCurrEntry(){return fCurrEntry;}
      Long64_t Size(){if(!fWTree) return 0;return fWTree->GetEntries();}
      void Add(Weights* wm);
      // void Multiply(Weights* other,TString species);
      void SetSpecies(TString name);
      Int_t GetNSpecies(){return fSpecies.size();}
      StrIntMap_t GetSpecies(){return fSpecies;}
      StrIntMap_t* GetSpeciesp(){return &fSpecies;}
      TString GetSpeciesName(UInt_t isp);
      void SetSpecies(StrIntMap_t species){fSpecies=std::move(species);};
      Int_t GetSpeciesID(const TString& name){if(fSpecies.count(name))return fSpecies[name]; else {cout<<"Weights:: GetSpeciesID species not found in "<<GetName()<<" species= "<<name<<endl;return -1;}}
      TString GetIDName(){return fIDName;}
      void SetIDName(TString name){fIDName=std::move(name);}
    
      TList* GetWeightList(){return fWeightList;}
      void PrintWeight();
      Long64_t Merge(const TString& tempName,const TString& outName="",const TString& wmName="WeightMap");
      void SortWeights();
      void BuildIndex();
      void SetFile(const TString& filename);
      void Save();
      void LoadSaved(const TString& fname,const TString& wname);
      void LoadSavedDisc(const TString& fname,const TString& wname);
      void WeightBySelection(TTree* tree,const TCut& cut,Double_t wgt);
      void WeightBySelection(TTree* tree,const TCut& cut,const TString& wgt);

      void AddToTree(TTree* tree);
      void AddToTreeDisc(TTree* tree,const TString& fileName);
      void ImportanceSampling(TTree* MCTree, TTree* dataTree, TH1* weightHist, TString var);
      void Draw1DWithWeights(TTree* tree,TH1* his,TString var,TString species="");

      filed_uptr DFAddToTree(const TString& wname,const TString& outfname,const TString& tname,const TString& infname);
    private:
      TTree *fWTree=nullptr;  //! not saved tree of weights, branchname = species
      TTree *fIDTree=nullptr;  //! not saved tree of ids, branchname = species
      TList* fWeightList=nullptr; //list of weight bins which have been merged to make this
      TFile* fFile=nullptr;
      TFile* fBranchFile=nullptr;
      TVectorD fWVals;
      Long64_t fID{};
      Long64_t fCurrEntry{};
      Long64_t fN{};
      Long64_t *fIDi=nullptr;//!
      Long64_t *fIDv=nullptr;//!
      StrIntMap_t fSpecies;//names of species with index in map
      TString fIDName; //name of tree branch with event ID
      Bool_t fGotEntry{};
      Bool_t fIsSorted{};
    
      ClassDefOverride(HS::FIT::Weights, 2);  // Writeble Weight map  class
    };

    class WeightsConfig{
      
    public:
      WeightsConfig()=default;
      WeightsConfig(const TString& wopt){
	auto opts=wopt.Tokenize(",");
	if(opts->GetEntries()!=3)
	  cout<<"Error WeightsConfig(TString wopt) need 3 arguments "<<endl;
	fSpecies=opts->At(0)->GetName();
	fFile=opts->At(1)->GetName();
	fObjName=opts->At(2)->GetName();
      }

    
    WeightsConfig( TString species, TString weightfile,TString obname)
      :fSpecies(std::move(std::move(species))),fFile(std::move(std::move(weightfile))),fObjName(std::move(std::move(obname))){}

      const TString Species() const {return fSpecies;}
      const TString File()  const {
	if(!fFile.BeginsWith("/"))
	  return  (TString)gSystem->Getenv("PWD")+"/"+fFile;
	return fFile;
      }
      const TString ObjName() const {return fObjName;}

      Bool_t IsValid() const{ if(fSpecies==TString()) return kFALSE;return kTRUE; }
      void Copy(const HS::FIT::WeightsConfig& wcon){
	fSpecies=wcon.Species();
	fFile=wcon.File();
	fObjName=wcon.ObjName();
      }
    private:
      TString fSpecies;
      TString fFile;
      TString fObjName;
    };
    
  }//namespace HS
}

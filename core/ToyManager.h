////////////////////////////////////////////////////////////////
///
///Class:               ToyManager
///Description:
///           

#pragma once

#include <utility>


#include "FitManager.h"
#include "Setup.h"
#include "Weights.h"
 
namespace HS{
  namespace FIT{

    using tree_uptr =std::unique_ptr<TTree>;
    using strings_t = std::vector<TString>;

    class ToyManager  : public FitManager{
      
    public:
      ToyManager()=default;
    ToyManager(Int_t n):fNToys(n){};
      ToyManager(const ToyManager&)=default;
      // ToyManager(Int_t n,const FitManager& fm,TString outDir="",TString resultFile=""):fNToys(n), FitManager(fm),fResultOutDir(std::move(std::move(outDir))),fResultFileName(std::move(std::move(resultFile))){};
      ToyManager(Int_t n,const FitManager& fm,TString outDir="",TString resultFile=""):fNToys(n), FitManager(fm),fResultOutDir(outDir),fResultFileName(resultFile){cout<<"ToyManager() "<<fResultFileName<<endl;};
      ToyManager(ToyManager&&)=delete;
      ~ToyManager() override =default;
      ToyManager& operator=(const ToyManager& other) = default;
      ToyManager& operator=(ToyManager&& other) = delete;
      
      Bool_t Run() override;
      void SaveResults() override;
      Int_t GetN() override { 
      	if(!(Bins().GetSize()))
       	  Bins().InitBins(); 
	if(Bins().GetSize())return Bins().GetSize();
       	return 1;
      } 
      //  Int_t GetCurrToy(){ return GetFiti()%fNToys;}
      Int_t GetCurrToy(){ return fToyi;}
      TString GetCurrTitle() override {return Form("Toy%d",GetCurrToy());}
      TString GetDataTreeName() override{return "ToyData";}
      strings_t GetDataFileNames() override{return fToyFileNames;}

      /* Int_t GetDataBin(Int_t ii) override{ */
      /* 	if(fNToys>0) */
      /* 	  return (int)std::round(ii/fNToys); */
      /* 	return ii; */
      /* } */
      strings_t GetToyFileNames(){return fToyFileNames;}
      
      void  Generate();

      std::unique_ptr<FitManager> Fitter();
 
      static std::shared_ptr<ToyManager> GetFromFit(Int_t N,const TString& filename,const TString& result="");
      static std::shared_ptr<ToyManager> GetFromFit(Int_t N,FitManager& fit,const TString& result="");
      static std::shared_ptr<ToyManager> GetFromFit(Int_t N,const std::shared_ptr<FitManager>& fit,const TString& result="");

      void Summarise();
      void Summarise(Int_t ib);
      void PreRun() override;
      void LoadResult();
      void InitSummary();
      
      static const TString InitialParsName(){return "InitialParameters";}

      void SetResultFileName(TString name){fResultFileName=std::move(name);}
    protected:
    
    private:
      RooDataSet* fGenData=nullptr;//!
      strings_t fToyFileNames;

      TString fResultOutDir;
      TString fResultFileName;
      Double_t fIDval=0;
      Int_t fNToys=1;
      Int_t fToyi=0;

      ClassDefOverride(HS::FIT::ToyManager,1);
    };
    
  }//namespace FIT
}//namespace HS


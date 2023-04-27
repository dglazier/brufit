#pragma once

#include <RooAbsPdf.h>
#include "RooHSEventsPDF.h"
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooAbsCategory.h>
#include <RooFormulaVar.h>
#include <vector>
 
namespace HS{
  namespace FIT{
 
     using std::unique_ptr;

     class RooComponentsPDF : public HS::FIT::RooHSEventsPDF {

      
      using vecUPtrReal = std::vector<unique_ptr<RooRealProxy>>;
      using vecUPtrCat = std::vector<unique_ptr<RooCategoryProxy>>;
      using vecComponents = std::vector<vecUPtrReal>;
 
    public:
      RooComponentsPDF() = default; ; 
      RooComponentsPDF(const char *name, const char *title,Double_t base,const RooArgList& obsList,const vector<RooArgList> compList);
      RooComponentsPDF(const RooComponentsPDF& other, const char* name=nullptr) ;
      TObject* clone(const char* newname) const override { return new RooComponentsPDF(*this,newname); }
       inline ~RooComponentsPDF() override =default;
 
      
      Double_t analyticalIntegral(Int_t code,const char* rangeName) const override;
      Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const override;

      Bool_t SetEvTree(TTree* tree,TString cut,TTree* MCGenTree=nullptr) override;
      void HistIntegrals(const char* rangeName) const override;
      void CalcWeightedBaseLine(const char* rangeName) const;
      void RedirectServersToData();
      void RedirectServersToPdf();
      Bool_t isDirectGenSafe(const RooAbsArg& arg) const override ;
      void initGenerator(Int_t code) override;

    protected:
  
      Double_t evaluateData() const override ;
      Double_t evaluateMC(const vector<Float_t> *vars,const  vector<Int_t> *cats) const override;
      void MakeSets();
      void RecalcComponentIntegrals(Int_t code,const char* rangeName) const;
      Double_t componentIntegral(Int_t icomp) const;
      void initIntegrator() override;
 
       void RecalcComponentIntegralsSampling(Int_t code,const char* rangeName) const;
       Double_t componentVariance(Int_t icomp) const;
       void DoFirstIntegrations(const char* rangeName="") const;

       // Double_t sampleIntegral() const;
       
     private:

       RooListProxy fActualObs;
      RooListProxy fActualCats;
      RooListProxy fActualComps;
      
      vecComponents fComponents;
      vecUPtrReal fObservables;
      vecUPtrCat fCategories;

      vector<vector<RooRealProxy*>> fDependentTermProxy;
      vector<vector<RooRealVar*>> fDependentTermParams;
      vector<vector<RooRealProxy*>> fIndependentTermProxy;

      
      vector<RooRealVar*> fIntegrateObs;
      vector<RooCategory*> fIntegrateCats;
      RooArgSet fIntegrateSet;
      
      mutable vector<Double_t> fCacheCompDepIntegral;
      mutable vector<Double_t> fCacheCompDepSigmaIntegral;
      mutable vector<vector<Double_t>> fPrevParVals;
      mutable vector<UInt_t> fRecalcComponent;
      
      RooArgSet fParameters;
 
      Double_t fBaseLine=0;
      mutable Double_t fWeightedBaseLine=0;
      mutable Double_t fNUsedForIntegral=0;
      UInt_t fNObs=0;
      UInt_t fNCats=0;
      UInt_t fNComps=0;
      mutable Bool_t fFirstCalculation=kTRUE;
       
      ClassDefOverride(HS::FIT::RooComponentsPDF,1);
    };

    template<typename T, typename A>
      bool vecContains( T arg, std::vector<T,A> const& vec ) {
      if(std::find(vec.begin(),vec.end(),arg) == vec.end())
	return false;
      return true;
    }


  }
}

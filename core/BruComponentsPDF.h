/**
 * @file BruComponentsPDF.h
 */

#pragma once

#include "BruEventsPDF.h"
#include <RooAbsPdf.h>
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooAbsCategory.h>
#include <RooFormulaVar.h>
#include <vector>
#include <memory>
 
namespace bru {
 
    using std::unique_ptr;

    class BruComponentsPDF : public BruEventsPDF {
      
        using vecUPtrReal = std::vector<unique_ptr<RooRealProxy>>;
        using vecUPtrCat = std::vector<unique_ptr<RooCategoryProxy>>;
        using vecComponents = std::vector<vecUPtrReal>;
    public:
        BruComponentsPDF() = default; 
        BruComponentsPDF(const char *name, const char *title, Double_t base, const RooArgList& obsList, const std::vector<RooArgList> compList);
        BruComponentsPDF(const BruComponentsPDF& other, const char* name = nullptr);
        TObject* clone(const char* newname) const override { return new BruComponentsPDF(*this, newname); }
        ~BruComponentsPDF() override = default;
 
        Double_t analyticalIntegral(Int_t code, const char* rangeName) const override;
        Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const override;

        void HistIntegrals(const char* rangeName) const override;
        void CalcWeightedBaseLine(const char* rangeName) const;
      
        void RedirectServersToData() const;
        void RedirectServersToPdf() const;
      
        bool isDirectGenSafe(const RooAbsArg& arg) const override;
        void initGenerator(Int_t code) override;

    protected:
        void initIntegrator() const override;

        void InitAssertPositiveCheck() const override {
	  initIntegrator(); 
            RedirectServersToPdf();
            if (_MCAPDepTerm.empty() == kTRUE)
                cacheMCAP(&_AssertPosDataReal, &_AssertPosDataCats);
            _assertPostive = kTRUE;
        };
       
        void FinishAssertPositiveCheck() const override {
            RedirectServersToData(); 
            _assertPostive = kFALSE;
        };
      
        Double_t cacheMCAP(const std::vector<Float_t> *vars, const std::vector<Int_t> *cats) const;
        Double_t evaluateMCAP() const;
	 
        Double_t evaluateData() const override;
        Double_t evaluateMC(const std::vector<Float_t> *vars, const std::vector<Int_t> *cats) const override;
        void MakeSets();
        void RecalcComponentIntegrals(Int_t code, const char* rangeName) const;
        Double_t componentIntegral(Int_t icomp) const;
 
        void RecalcComponentIntegralsSampling(Int_t code, const char* rangeName) const;
        Double_t componentVariance(Int_t icomp) const;
        void DoFirstIntegrations(const char* rangeName = "") const;

        Bool_t CheckChange() const override;
      
    private:
        RooListProxy _ActualObs;
        RooListProxy _ActualCats;
        RooListProxy _ActualComps;

	
        vecComponents _Components;
        vecUPtrReal _Observables;
        vecUPtrCat _Categories;

        mutable std::vector<std::vector<RooRealProxy*>> _DependentTermProxy;
        mutable std::vector<std::vector<RooRealVar*>> _DependentTermParams;
        mutable std::vector<std::vector<RooRealProxy*>> _IndependentTermProxy;

        std::vector<std::unique_ptr<RooRealProxy>> _myVarProxies; 
       
        std::vector<RooRealVar*> _IntegrateObs;
        std::vector<RooCategory*> _IntegrateCats;
        RooArgSet _IntegrateSet;
        RooArgSet _Parameters; 
   
        mutable std::vector<Double_t> _CacheCompDepIntegral;
        mutable std::vector<Double_t> _CacheCompDepSigmaIntegral;
        mutable std::vector<std::vector<Double_t>> _PrevParVals;
        mutable std::vector<UInt_t> _RecalcComponent;

        mutable std::vector<std::vector<Double_t>> _MCAPDepTerm;
       
        Double_t _BaseLine = 0;
        mutable Double_t _WeightedBaseLine = 0;
        mutable Double_t _NUsedForIntegral = 0;
        UInt_t _NObs = 0;
        UInt_t _NCats = 0;
        UInt_t _NComps = 0;
        mutable Bool_t _FirstCalculation = kTRUE;
        mutable Bool_t _once = kTRUE;
        mutable Bool_t _assertPostive = kFALSE;
        mutable Long64_t _NIntegralCalls = 0;
      
        ClassDefOverride(bru::BruComponentsPDF, 1);
    };

    template<typename T, typename A>
    bool vecContains(T arg, std::vector<T, A> const& vec) {
        return std::find(vec.begin(), vec.end(), arg) != vec.end();
    }

} // namespace bru
/* /\** */
/*  * @file BruComponentsPDF.h */
/*  *\/ */

/* #pragma once */

/* #include "BruEventsPDF.h" */
/* #include <RooAbsPdf.h> */
/* #include <RooRealProxy.h> */
/* #include <RooCategoryProxy.h> */
/* #include <RooAbsReal.h> */
/* #include <RooRealVar.h> */
/* #include <RooCategory.h> */
/* #include <RooAbsCategory.h> */
/* #include <RooFormulaVar.h> */
/* #include <vector> */
/* #include <memory> */
 
/* namespace bru { */
 
/*     using std::unique_ptr; */

/*     class BruComponentsPDF : public BruEventsPDF { */
      
/*         using vecUPtrReal = std::vector<unique_ptr<RooRealProxy>>; */
/*         using vecUPtrCat = std::vector<unique_ptr<RooCategoryProxy>>; */
/*         using vecComponents = std::vector<vecUPtrReal>; */
 
/*     public: */
/*         BruComponentsPDF() = default;  */
/*         BruComponentsPDF(const char *name, const char *title, Double_t base, const RooArgList& obsList, const std::vector<RooArgList> compList); */
/*         BruComponentsPDF(const BruComponentsPDF& other, const char* name = nullptr); */
/*         TObject* clone(const char* newname) const override { return new BruComponentsPDF(*this, newname); } */
/*         ~BruComponentsPDF() override = default; */
 
/*         Double_t analyticalIntegral(Int_t code, const char* rangeName) const override; */
/*         Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const override; */

/*         void HistIntegrals(const char* rangeName) const override; */
/*         void CalcWeightedBaseLine(const char* rangeName) const; */
      
/*         void RedirectServersToData() const; */
/*         void RedirectServersToPdf() const; */
      
/*         bool isDirectGenSafe(const RooAbsArg& arg) const override; */
/*         void initGenerator(Int_t code) override; */

/*     protected: */
/*         // BUG FIX: Lazy integration initialization */
/*         mutable Bool_t _IntegratorInit = kFALSE; */
/*         void SetupIntegrator() const; */

/*         void InitAssertPositiveCheck() const override { */
/*             SetupIntegrator(); // Ensure safe array sizing before looping */
/*             RedirectServersToPdf(); */
/*             if (_MCAPDepTerm.empty() == kTRUE) */
/*                 cacheMCAP(&_AssertPosDataReal, &_AssertPosDataCats); */
/*             _assertPostive = kTRUE; */
/*         }; */
       
/*         void FinishAssertPositiveCheck() const override { */
/*             RedirectServersToData();  */
/*             _assertPostive = kFALSE; */
/*         }; */
      
/*         Double_t cacheMCAP(const std::vector<Double_t> *vars, const std::vector<Int_t> *cats) const; */
/*         Double_t evaluateMCAP() const; */
	 
/*         Double_t evaluateData() const override; */
/*         Double_t evaluateMC(const std::vector<Double_t> *vars, const std::vector<Int_t> *cats) const override; */
/*         void MakeSets(); */
/*         void RecalcComponentIntegrals(Int_t code, const char* rangeName) const; */
/*         Double_t componentIntegral(Int_t icomp) const; */
 
/*         void RecalcComponentIntegralsSampling(Int_t code, const char* rangeName) const; */
/*         Double_t componentVariance(Int_t icomp) const; */
/*         void DoFirstIntegrations(const char* rangeName = "") const; */

/*         Bool_t CheckChange() const override; */
      
/*     private: */
/*         RooListProxy _ActualObs; */
/*         RooListProxy _ActualCats; */
/*         RooListProxy _ActualComps; */
      
/*         vecComponents _Components; */
/*         vecUPtrReal _Observables; */
/*         vecUPtrCat _Categories; */

/*         mutable std::vector<std::vector<RooRealProxy*>> _DependentTermProxy; */
/*         mutable std::vector<std::vector<RooRealVar*>> _DependentTermParams; */
/*         mutable std::vector<std::vector<RooRealProxy*>> _IndependentTermProxy; */

/*         std::vector<std::unique_ptr<RooRealProxy>> _myVarProxies;  */
       
/*         std::vector<RooRealVar*> _IntegrateObs; */
/*         std::vector<RooCategory*> _IntegrateCats; */
/*         RooArgSet _IntegrateSet; */
/*         RooArgSet _Parameters; //! */
   
/*         mutable std::vector<Double_t> _CacheCompDepIntegral; */
/*         mutable std::vector<Double_t> _CacheCompDepSigmaIntegral; */
/*         mutable std::vector<std::vector<Double_t>> _PrevParVals; */
/*         mutable std::vector<UInt_t> _RecalcComponent; */

/*         mutable std::vector<std::vector<Double_t>> _MCAPDepTerm; */
       
/*         Double_t _BaseLine = 0; */
/*         mutable Double_t _WeightedBaseLine = 0; */
/*         mutable Double_t _NUsedForIntegral = 0; */
/*         UInt_t _NObs = 0; */
/*         UInt_t _NCats = 0; */
/*         UInt_t _NComps = 0; */
/*         mutable Bool_t _FirstCalculation = kTRUE; */
/*         mutable Bool_t _once = kTRUE; */
/*         mutable Bool_t _assertPostive = kFALSE; */
/*         mutable Long64_t _NIntegralCalls = 0; */
      
/*         ClassDefOverride(bru::BruComponentsPDF, 1); */
/*     }; */

/*     template<typename T, typename A> */
/*     bool vecContains(T arg, std::vector<T, A> const& vec) { */
/*         return std::find(vec.begin(), vec.end(), arg) != vec.end(); */
/*     } */

/* } // namespace bru */

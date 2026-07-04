/**
 * @file BruEventsPDF.h
 */

#pragma once

#include "MCEventCache.h"
#include <RooAbsPdf.h>
#include <RooArgSet.h>
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooRandom.h>
#include <TString.h>
#include <TH1F.h>

#include <memory>
#include <vector>

namespace bru {
    
    class BruEventsPDF : public RooAbsPdf {
      
    public:
        static Bool_t BruEventsPDF_IsPlotting;
        static void SetIsPlotting(Bool_t is);

        BruEventsPDF(const char *name, const char *title) : RooAbsPdf(name, title) {};
        BruEventsPDF(const BruEventsPDF& other, const char* name = nullptr);
        BruEventsPDF() = default; 
        ~BruEventsPDF() override;
  
    protected:
        std::shared_ptr<const MCEventCache> _DataCache; 
        
        BruEventsPDF* _Parent = nullptr; 
        Double_t _ConstInt = 1;

        mutable Double_t _SigmaIntegral = 0;
        mutable std::vector<Double_t> _Last;  
        mutable std::vector<TH1F> _HistIntegrals;
        
        Int_t _LastLength{0};
        Long64_t _NInt = -1;
        Long64_t _IntRangeLow = 0;
        Long64_t _IntRangeHigh = 0;
        mutable Long64_t _TreeEntry = 0;
        Int_t _NRanges = 1;
        Int_t _CheckInt = 0;
        Int_t _Npars = 0;
        Int_t _Nvars = 0;
        Int_t _Ncats = 0;
        
        Bool_t _IsIntegrating = kFALSE;
        Bool_t _IsClone = kFALSE;
        
        // Fix: Defaults to FALSE
        Bool_t _ForceConstInt = kFALSE; 
        Bool_t _ForceNumInt = kFALSE;
        Bool_t _UseWeightsGen = kFALSE;

        std::vector<RooArgSet*> _VarSet;
        std::vector<RooRealProxy*> _ProxSet; 
        std::vector<RooCategoryProxy*> _CatSet;
        std::vector<RooRealProxy*> _ParSet;
  
        std::vector<Float_t> _AssertPosDataReal;
        std::vector<Int_t> _AssertPosDataCats;
        Long64_t _Napd = 10000;

        Double_t _MaxValue = 0; 
        Long64_t _Geni = 0; 
        mutable Int_t _IntCounter = 0;
        Bool_t _IsPlotting = kFALSE;
        Bool_t _UseSamplingIntegral = kFALSE;
        
        TString _TruthPrefix = "gen";
        HS::FIT::WeightsConfig _WgtsConf;
        TString _Cut; 
        TString _InWeightCut; 
        Bool_t _IsValid = kTRUE;
        Bool_t _BranchStatus = kTRUE;

        std::vector<TString> _ProtoRealVars;
        std::vector<TString> _ProtoCatVars;
        
        std::vector<Long64_t> _GeneratedIndices; 
        std::vector<Double_t> _GeneratedWeights; 
      
        void InitSets();
        RooArgSet VarSet(Int_t iset) const;
        virtual void HistIntegrals(const char* rangeName) const;
        void SetLowHighVals(Long64_t& ilow, Long64_t& ihigh) const;

        virtual Double_t evaluateData() const { return 0; }
        virtual void initIntegrator() const ;

        Double_t sampleIntegral(Double_t integral, Double_t sigma) const {
            return RooRandom::gaussian() * sigma + integral;
        }

        // Helper to match legacy GetIntegralWeight(ie)
        inline Double_t GetIntegralWeight(Long64_t ie) const {
            return (_DataCache && _DataCache->_UseEvWeights) ? _DataCache->_EvWeights[ie] : 1.0;
        }
      
    public:
        void SetTruthPrefix(const TString& pre) { _TruthPrefix = pre; }
        Bool_t IsValid() const { return _IsValid; }
        void SetInWeights(const HS::FIT::WeightsConfig& wcon) { _WgtsConf.Copy(wcon); }
        void SetInWeights(const TString& wst) {
            if(wst == TString()) return;
            if(!wst.Contains(",")) { _InWeightCut = wst; return; }
            HS::FIT::WeightsConfig wcon(wst);
            _WgtsConf.Copy(wcon);
        }
        virtual Bool_t SetEvTree(TTree* tree, TString cut, TTree* MCGenTree = nullptr);
	virtual void ResetTree(){
	  if (_DataCache) {
            _DataCache.reset();
	  }
	}
	
        void MakeAssertPostiveData();
        Bool_t AssertPositivePDF() const;
        virtual void InitAssertPositiveCheck() const {};
        virtual void FinishAssertPositiveCheck() const {};
      
        Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const override;
        Double_t analyticalIntegral(Int_t code, const char* rangeName) const override;
        Double_t unnormalisedIntegral(Int_t code, const char* rangeName) const;
        Double_t analyticalIntegralForSampling(const char* rangeName) const;
        void generateEvent(Int_t code) override;
        Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const override;
        void initGenerator(Int_t code) override;
      
        bool forceAnalyticalInt(const RooAbsArg& arg) const override {
            for(const auto var : _ProxSet) {
                if(TString(arg.GetName()) == var->GetName()) return true;
            }
            return false;
        }

        virtual Double_t evaluateMC(const std::vector<Float_t> *vars, const std::vector<Int_t> *cats) const { return 0.; }

        Double_t evaluate() const override {
            if(!BruEventsPDF_IsPlotting) return evaluateData();
            if(!_HistIntegrals.empty()) {
                if(_ProxSet.size() == 1) return _HistIntegrals[0].Interpolate(*_ProxSet[0]); 
            }
            return evaluateData();
        }
       
        virtual Bool_t CheckChange() const; 
        Bool_t CheckRange(const char* rangeName) const; 

        void SetCache(std::shared_ptr<const MCEventCache> cache);
        Bool_t AddProtoData(const RooDataSet* data); 
        
        void SetNInt(Long64_t n) { _NInt = n; }
        void SetUseWeightsGen(Bool_t use=kTRUE) { _UseWeightsGen = use; }
        Bool_t UseWeightsGen() const { return _UseWeightsGen; }
        
        const std::vector<Long64_t>& GetGeneratedIndices() const { return _GeneratedIndices; }
        const std::vector<Double_t>& GetGeneratedWeights() const { return _GeneratedWeights; }
        Long64_t GetNMCGenEntries() const { 
            if (_DataCache && _DataCache->_HasMCGenTree) {
                return _DataCache->_NMCGenTreeEntries;
            }
            return 0;
        }
        void SetGeni(Long64_t gi) { _Geni = gi; }
        Long64_t IncrementGeni() {
            if(_IsClone) _Parent->SetGeni(_Geni);
            ++_Geni;
            return _Geni;
        }
        Long64_t GetGeni() const { return _Geni; }
      
        void SetConstInt(Bool_t force=kTRUE) { _ForceConstInt = force; }
        void SetNumInt(Bool_t force=kTRUE) { _ForceNumInt = force; }
        void CheckIntegralParDep(Int_t Ntests);
        
        Double_t GetMaxValue() const { return _MaxValue; }
        void SetMaxValue(Double_t val) { _MaxValue = val; }
        void SetIntRange(Long64_t low, Long64_t high) { _IntRangeLow = low; _IntRangeHigh = high; }
        Long64_t GetIntRangeLow() const { return _IntRangeLow; }
        Long64_t GetIntRangeHigh() const { return _IntRangeHigh; }
        void SetNRanges(Int_t nr) { _NRanges = nr; }
        void SetNextRange(Int_t ir);
        BruEventsPDF* GetParent() const { return _Parent; }
        Bool_t HasMCGenTree() const { return _DataCache && _DataCache->_HasMCGenTree; }
        void Plotting(Bool_t plotting=kTRUE) { _IsPlotting = plotting; }
        void SetHistIntegrals(std::vector<TH1F> &hists) { _HistIntegrals = hists; }
        void ResetHistIntegrals() { _HistIntegrals.clear(); }

      virtual Double_t GetRelativeVariance() const;
        
      // Setting this forces the integration loops to track sum-of-squares
      // and triggers the cross-matrix sampling engine in sub-classes
      void TrackMCVariance(Bool_t track = kTRUE) { 
	_TrackMCVariance = track; 
	if (track) {
	  _UseSamplingIntegral = kTRUE; 
	} else {
	  _UseSamplingIntegral = kFALSE; 
	}
      }
      
    protected:
      Bool_t _TrackMCVariance = kFALSE;
      mutable Double_t _SumSquares = 0;         // NEW: Tracks E[F^2]
      mutable Long64_t _NUsedForIntegral = 0;   // NEW: Tracks N
      
        ClassDefOverride(bru::BruEventsPDF, 1);
    };

} // namespace bru


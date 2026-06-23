/**
 * @file BruEventsHistPDF.h
 */

#pragma once

#include "BruEventsPDF.h"
#include <RooRealProxy.h>
#include <RooHistPdf.h>
#include <RooCategoryProxy.h>
#include <RooAbsReal.h>
#include <RooAbsCategory.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooConstVar.h>
#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <TH2.h>
#include <vector>

namespace bru {
    
    class BruEventsHistPDF : public BruEventsPDF {
    public:
      
        BruEventsHistPDF() = default; 
       
        BruEventsHistPDF(const char *name, const char *title, RooAbsReal& in_x, RooAbsReal& in_alpha, RooAbsReal& in_offset, RooAbsReal& in_scale, Int_t applySmooth=1, Int_t interp=1, Int_t xbins=100, Int_t nsamp=1000, Int_t abins=200);
      
        BruEventsHistPDF(const BruEventsHistPDF& other, const char* name=nullptr);
        TObject* clone(const char* newname) const override { return new BruEventsHistPDF(*this, newname); }
        ~BruEventsHistPDF() override;

        // Modern RooFit Batch Evaluation
        void doEval(RooFit::EvalContext & ctx) const override;
        
        // Custom MC Batch Evaluation
        virtual void evaluateMCBatch(const std::vector<Double_t>& mcx_array, std::vector<Double_t>& output) const;

    protected:
        Double_t _MCx{};

        RooRealProxy x;
        RooRealProxy offset;
        RooRealProxy scale;
        RooRealProxy alpha;

        TH1D _GenHist; 
      
        Int_t _applySmooth = 1; 
        Int_t _Interpolate = 1; 
        Int_t _NAlphaBins = 200;
        Int_t _NXBins0 = 100;
        Int_t _NIntSamples = 1000;
        Bool_t _UseHistGenerator = kTRUE;
      
        Double_t evaluate() const override;
        Double_t evaluateMC(const std::vector<Float_t> *vars, const std::vector<Int_t> *cats) const override;
        Double_t evaluateMC(Double_t mcx) const;
        void MakeSets();
  
        RooDataHist* _Hist = nullptr;
        TH2D* _RHist = nullptr;
        Double_t _VarMax{};
  
        RooRealVar* _x_off = nullptr; 
        RooRealVar* _alphaVar = nullptr;

	// --- AVX2 FAST-MATH CACHE SYSTEM ---
        std::vector<double> _cachedBins;
        double _xMin = 0.0;
        double _xWidth = 0.0;
        int _nx = 0;
        double _yMin = 0.0;
        double _yWidth = 0.0;
        int _ny = 0;
        
        virtual void initializeCache(); // Eagerly called during PDF setup/cloning

    private:
        RooGaussian *_AlphaConstr = nullptr;
        RooGaussian *_OffConstr = nullptr;
        RooGaussian *_ScaleConstr = nullptr;
  
      
    public:
     
        Bool_t SetEvTree(TTree* tree, TString cut, TTree* MCGenTree = nullptr) override;
        virtual void CreateHistPdf();
        virtual void FillBase1DHist(TH1D& his1);
        void CheckForNegativeBins(TH1D& his1);
        void Construct2DHist(TAxis xaxis, TAxis AlphAxis, Bool_t isAlphaConst);
      
        virtual void ResetTree() override;
  
        Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const override;
        Double_t analyticalIntegral(Int_t code, const char* rangeName) const override;
        void generateEvent(Int_t code) override;
        void initGenerator(Int_t code) override;
        void UseHistGenerator() { _UseHistGenerator = kTRUE; }
        Bool_t UsingHistGenerator() { return _UseHistGenerator; }
        
        RooGaussian* AlphaConstraint() { return _AlphaConstr; }
        RooGaussian* OffConstraint() { return _OffConstr; }
        RooGaussian* ScaleConstraint() { return _ScaleConstr; }

        TH2* GetRootHist() { return _RHist; }

        std::vector<Double_t> GetBinVector(const TAxis& ax) {
            std::vector<Double_t> xedges(ax.GetNbins());
            ax.GetLowEdge(xedges.data());
            xedges.push_back(ax.GetXmax());
            return xedges;
        }
      
        ClassDefOverride(bru::BruEventsHistPDF, 1);
    };

} // namespace bru

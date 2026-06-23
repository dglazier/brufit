/**
 * @file BruEventsHistPeakPDF.h
 */

#pragma once

#include "BruEventsHistPDF.h"
#include <RooRealProxy.h>
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
    
    class BruEventsHistPeakPDF : public BruEventsHistPDF {
    public:
      
        BruEventsHistPeakPDF() = default; 
       
        BruEventsHistPeakPDF(const char *name, const char *title, RooAbsReal& in_x, RooAbsReal& in_alpha, RooAbsReal& in_offset, RooAbsReal& in_scale, Int_t applySmooth=1, Int_t interp=1, Int_t xbins=100, Int_t nsamp=1000, Int_t abins=200);
      
        BruEventsHistPeakPDF(const BruEventsHistPeakPDF& other, const char* name = nullptr);
        TObject* clone(const char* newname) const override { return new BruEventsHistPeakPDF(*this, newname); }
        ~BruEventsHistPeakPDF() override = default;

        // Vectorized Evaluation Override
        void doEval(RooFit::EvalContext & ctx) const override;
      void evaluateMCBatch(const std::vector<Double_t>& mcx_array, std::vector<Double_t>& output) const override;
      
    protected:
        void FillBase1DHist(TH1D& his1) override;
        
        // --- VARIABLE BIN FAST-LOOKUP CACHE ---
        std::vector<int> _lut;           // O(1) Index Lookup Table
        std::vector<double> _varCenters; // Actual centers of variable bins
        
        int _nLut = 0;               // Resolution of the lookup table
        double _lutWidth = 0.0;          // Uniform width of the lookup slices
        
        void initializeCache() override;
      
        ClassDefOverride(bru::BruEventsHistPeakPDF, 1); 
    };

} // namespace bru

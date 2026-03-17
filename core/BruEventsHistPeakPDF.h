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

namespace bru {
    
    class BruEventsHistPeakPDF : public BruEventsHistPDF {
    public:
      
        BruEventsHistPeakPDF() = default; 
       
        BruEventsHistPeakPDF(const char *name, const char *title, RooAbsReal& in_x, RooAbsReal& in_alpha, RooAbsReal& in_offset, RooAbsReal& in_scale, Int_t applySmooth=1, Int_t interp=1, Int_t xbins=100, Int_t nsamp=1000, Int_t abins=200);
      
        BruEventsHistPeakPDF(const BruEventsHistPeakPDF& other, const char* name = nullptr);
        TObject* clone(const char* newname) const override { return new BruEventsHistPeakPDF(*this, newname); }
        ~BruEventsHistPeakPDF() override = default;

    protected:
        void FillBase1DHist(TH1D& his1) override;
      
        ClassDefOverride(bru::BruEventsHistPeakPDF, 1); 
    };

} // namespace bru

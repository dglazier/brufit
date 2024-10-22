#pragma once

#include "RooHSEventsHistPDF.h"
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


namespace HS{
  namespace FIT{
    
    class BruEventsHistPeakPDF : public RooHSEventsHistPDF {
    public:
      
      BruEventsHistPeakPDF() =default ; 
       
      BruEventsHistPeakPDF(const char *name, const char *title, RooAbsReal& _x,RooAbsReal& _alpha,RooAbsReal& _offset, RooAbsReal& _scale, Int_t applySmooth=1, Int_t interp=1, Int_t xbins=100, Int_t nsamp=1000, Int_t abins=200);
      
      BruEventsHistPeakPDF(const BruEventsHistPeakPDF& other, const char* name=nullptr) ;
      TObject* clone(const char* newname) const override { return new BruEventsHistPeakPDF(*this,newname); }
      ~BruEventsHistPeakPDF() =default;

    protected:
      void  FillBase1DHist(TH1D& his1) override;
      
      ClassDefOverride(HS::FIT::BruEventsHistPeakPDF,1); 
    };//Class

    
  }//namespace FIT
}//namespace HS

 

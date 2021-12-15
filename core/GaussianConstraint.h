#pragma once

#include <RooGaussian.h>
#include <RooRandom.h>
#include <RooRealVar.h>
#include <TString.h>


namespace HS{
  namespace FIT{
    
    class GaussianConstraint {

    public:

    GaussianConstraint(TString name):
      fPDF(name.Data(),name.Data(),fRooX,fRooMean,fRooSigma),
	fRooX((name+"X").Data(),(name+"X").Data(),fX,0.8,1.2),
 	fRooMean((name+"Mean").Data(),(name+"Mean").Data(),fMean),
 	fRooSigma((name+"Sigma").Data(),(name+"Sigma").Data(),fSigma)
	{
	  fRooX.setConstant(kTRUE);
	}
      
      void setXMeanSigma(Double_t var,Double_t mean,Double_t sigma){
	fRooX.setVal(fX=var);
	fRooMean.setVal(fMean=mean);
	fRooSigma.setVal(fSigma=sigma);
      }
      void setMeanSigma(Double_t mean,Double_t sigma){
	fRooMean.setVal(fMean=mean);
	fRooSigma.setVal(fSigma=sigma);
      }

      RooGaussian* getPDF() {return &fPDF;}

      Double_t sample() {
	if(gRandom->Uniform()>0.95){
	  fRand=RooRandom::gaussian();
	  // std::cout<<"Chaning integral "<<std::endl;
	}
	fX=fRand*fSigma + fMean;
	fRooX.setVal(fX);
	return fX;//* fSigma + fMean;
      }
      Double_t getX() const { return fRooX.getVal();}//*fSigma + fMean; }
      RooRealVar& getRooVar(){return fRooX;}
      
      Double_t getVal() const {return fPDF.getVal(); } //change from normal dist.
    private:

      Double_t fMean=1;
      Double_t fSigma=1;
      Double_t fX=1;
      Double_t fRand=0.;
      RooRealVar fRooMean;
      RooRealVar fRooSigma;
      RooRealVar fRooX;
      RooGaussian fPDF; //must be initialised last

      ClassDef(HS::FIT::GaussianConstraint,1);
    };

  }
}

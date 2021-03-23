#pragma once

#include "RooHSComplex.h"
#include "RooHSDWigner.h"
#include <RooAbsReal.h>
#include <RooRealProxy.h>
#include <TMath.h>
//#include "arr.h"


namespace HS{
  namespace FIT{

    class RooHSDWignerProduct : public RooHSComplex {

    public:
      RooHSDWignerProduct()=default;
      RooHSDWignerProduct(const char *name, const char *title,  RooAbsReal& dwigner1, RooAbsReal& dwigner2, double factor=1);
      RooHSDWignerProduct(const RooHSDWignerProduct& other, const char* name=nullptr);
      TObject* clone(const char* newname) const override{ return new RooHSDWignerProduct(*this, newname);}
      ~RooHSDWignerProduct() override=default;


      TString FactoryReal() const final{
	return TString("RooHSDWignerProductRe::")
	  +"Re"+GetName()+Form("(%s,%s,%s)",PhiName1(),PhiName2(),GetName());
      }

      TString FactoryImag() const final{
	return TString("RooHSDWignerProductIm::")
	  +"Im"+GetName()+Form("(%s,%s,%s)",PhiName1(),PhiName2(),GetName());
      }

      TString FactoryImagConj() const final {
	return TString("RooHSDWignerProductIm::")
	  +"CoIm"+GetName()+Form("(%s,%s,%s,-1)",PhiName1(),PhiName2(),GetName());
      }


      Int_t M1() const {return _M1;}
      Int_t M2() const {return _M2;}
      const char* PhiName1()const {return _phiName1.Data();}
      const char* PhiName2()const {return _phiName2.Data();}
      
    protected:

      Double_t evaluate() const final;
      
   
      private:
      RooRealProxy _D1;
      RooRealProxy _D2;

      Int_t _M1=0;
      Int_t _M2=0;
      
      TString _phiName1;
      TString _phiName2;
      
      ClassDefOverride(HS::FIT::RooHSDWignerProduct,1);
    };

 
   inline Double_t RooHSDWignerProduct::evaluate() const
    {
    
      return _D1*_D2;
   
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////
    class RooHSDWignerProductIm : public RooAbsReal {
    public:
      RooHSDWignerProductIm() =default;
   
      RooHSDWignerProductIm(const char *name, const char *title,  RooAbsReal& phi1, RooAbsReal& phi2, RooAbsReal& dproduct, Int_t conj=1 );

      RooHSDWignerProductIm(const RooHSDWignerProductIm& other, const char* name = nullptr);
      TObject* clone(const char* newname) const override { return new RooHSDWignerProductIm(*this, newname); }
      ~RooHSDWignerProductIm() override =default;

    protected:
      Double_t evaluate() const final;
  
    private:

      RooRealProxy _phi1;
      RooRealProxy _phi2;
      RooRealProxy _mag;

      Int_t _M1=0;
      Int_t _M2=0;

      Short_t _conj=1;

      ClassDefOverride(HS::FIT::RooHSDWignerProductIm,1);
    };

    ////////////////////////////////////////////////////////////////////////////////
    inline Double_t RooHSDWignerProductIm::evaluate() const{
      if(! (_M1+_M2) )return 0; //M=0 =>real
    
      Double_t angle=_M1 * _phi1 + _M2*_phi2; 
      return _conj*_mag*TMath::Sin(angle);
    }

 
    ///////////////////////////////////////////////////////////////////////////////////////////////
    class RooHSDWignerProductRe : public RooAbsReal {
    public:
      RooHSDWignerProductRe() =default;
      RooHSDWignerProductRe(const char *name, const char *title,  RooAbsReal& phi1,  RooAbsReal& phi2, RooAbsReal& dproduct );

      RooHSDWignerProductRe(const RooHSDWignerProductRe& other, const char* name = nullptr);
      TObject* clone(const char* newname) const override { return new RooHSDWignerProductRe(*this, newname); }
      ~RooHSDWignerProductRe() override =default;

    protected:
      Double_t evaluate() const final;
    private:

      RooRealProxy _phi1;
      RooRealProxy _phi2;
      RooRealProxy _mag;

      Int_t _M1=0;
      Int_t _M2=0;
 
   
      ClassDefOverride(HS::FIT::RooHSDWignerProductRe,1);
    };

    ////////////////////////////////////////////////////////////////////////////////
    inline Double_t RooHSDWignerProductRe::evaluate() const{
      if( _M1+_M2 == 0 ) return _mag; //purely real
      Double_t angle=_M1 * _phi1 + _M2*_phi2;
      
      return _mag*TMath::Cos(angle);
    }
  

  }//FIT
}//HS 

#pragma once

#include "RooHSComplex.h"
#include <RooAbsReal.h>
#include <RooRealProxy.h>
#include <Math/SpecFunc.h>
#include <TMath.h>
//#include "arr.h"


namespace HS{
  namespace FIT{

    class RooHSDWigner : public RooHSComplex {

    public:
      RooHSDWigner()=default;
      RooHSDWigner(const char *name, const char *title, RooAbsReal& theta, RooAbsReal& phi, int l, int m=0,int s=0, double factor=1);
      RooHSDWigner(const RooHSDWigner& other, const char* name=nullptr);
      TObject* clone(const char* newname) const override{ return new RooHSDWigner(*this, newname);}
      ~RooHSDWigner() override=default;


      TString FactoryReal() const override{
	return TString("RooHSDWignerRe::")
	  +"Re"+GetName()+Form("(%s,%s,%s)",_theta.arg().GetName(),_phi.arg().GetName(),GetName());
      }

      TString FactoryImag() const override{
	return TString("RooHSDWignerIm::")
	  +"Im"+GetName()+Form("(%s,%s,%s)",_theta.arg().GetName(),_phi.arg().GetName(),GetName());
      }

      TString FactoryImagConj() const override {
	return TString("RooHSDWignerIm::")
	  +"CoIm"+GetName()+Form("(%s,%s,%s,-1)",_theta.arg().GetName(),_phi.arg().GetName(),GetName());
      }


      Int_t L() const{return _L;}
      Int_t M() const{return _M;}
      Int_t S() const{return _S;}
      Double_t evaluateMC(double th) const;
      const char* PhiName(){return _phi.arg().GetName();}
      const char* ThetaName(){return _theta.arg().GetName();}

    protected:

      Double_t evaluate() const final;

      double SmallWignerD( int aj, int am, int an, double beta ) const;
      
      Bool_t CheckClean() const;
 
    private:

      RooRealProxy _theta;
      RooRealProxy _phi;
      mutable Double_t _lastTheta=1E10;
      mutable Double_t _lastVal=0;
      Int_t _L=0;
      Int_t _M=0;
      Int_t _S=0;
  
      ClassDefOverride(HS::FIT::RooHSDWigner,1);
    };

    ////////////////////////////////////////////////////////////////////////////////
    
    inline Double_t RooHSDWigner::evaluateMC(double th) const
    {
      //     if(CheckClean()) return _lastVal;
     if(th==_lastTheta ) return _lastVal;
     _lastTheta=th;
 
     return  _lastVal = SmallWignerD( _L, _M, _S, th );
      
    
    }

   inline Double_t RooHSDWigner::evaluate() const
    {

      if(CheckClean()) return _lastVal;
      
      return  _lastVal = SmallWignerD( _L, _M, _S, _theta );
   
    }

    ////////////////////////////////////////////////////////////////////////////////

    inline Bool_t RooHSDWigner::CheckClean() const
    {
      Double_t th=_theta;
      if(th==_lastTheta ) return kTRUE;
      _lastTheta=th;
      return kFALSE;
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////
    class RooHSDWignerIm : public RooAbsReal {
    public:
      RooHSDWignerIm() =default;
   
      RooHSDWignerIm(const char *name, const char *title, RooAbsReal& theta, RooAbsReal& phi, RooAbsReal& dwigner, Int_t conj=1 );

      RooHSDWignerIm(const RooHSDWignerIm& other, const char* name = nullptr);
      TObject* clone(const char* newname) const override { return new RooHSDWignerIm(*this, newname); }
      ~RooHSDWignerIm() override =default;

       
    protected:
      Double_t evaluate() const final;
      Bool_t CheckClean() const;

    private:

      RooRealProxy _theta;
      RooRealProxy _phi;
      RooRealProxy _mag;
      mutable Double_t _lastTheta=-1E10;
      mutable Double_t _lastPhi=-1E10;
      mutable Double_t _lastVal=0;
      Int_t _M=0;
      Short_t _conj=1;

      ClassDefOverride(HS::FIT::RooHSDWignerIm,1);
    };

    ////////////////////////////////////////////////////////////////////////////////
    inline Double_t RooHSDWignerIm::evaluate() const{
      if(!_M)return 0; //M=0 =>real
      if(CheckClean()) return _lastVal;
   
      Double_t angle=_M*_lastPhi; 
      return _lastVal=_conj*_mag*TMath::Sin(angle);
    }

    ////////////////////////////////////////////////////////////////////////////////
    
    inline Bool_t RooHSDWignerIm::CheckClean() const
    {
      Double_t th=_theta;
      Double_t ph=_phi;
      if(th==_lastTheta&&ph==_lastPhi ) return kTRUE;
      _lastTheta=th;
      _lastPhi=ph;
      return kFALSE;
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////
    class RooHSDWignerRe : public RooAbsReal {
    public:
      RooHSDWignerRe() =default;
      RooHSDWignerRe(const char *name, const char *title, RooAbsReal& theta, RooAbsReal& phi, RooAbsReal& dwigner );

      RooHSDWignerRe(const RooHSDWignerRe& other, const char* name = nullptr);
      TObject* clone(const char* newname) const override { return new RooHSDWignerRe(*this, newname); }
      ~RooHSDWignerRe() override =default;

    protected:
      Double_t evaluate() const final;
      Bool_t CheckClean() const;
    private:

      RooRealProxy _theta;
      RooRealProxy _phi;
      RooRealProxy _mag;
      mutable Double_t _lastTheta=-1E10;
      mutable Double_t _lastPhi=-1E10;
      mutable Double_t _lastVal=0;
      Int_t _M=0;
 
      ClassDefOverride(HS::FIT::RooHSDWignerRe,1);
    };

    ////////////////////////////////////////////////////////////////////////////////
    inline Double_t RooHSDWignerRe::evaluate() const{
      if(CheckClean()) return _lastVal;
      ////CHANGE FOR WIGNER FUNC
      Double_t angle=_M*_lastPhi;
      return _lastVal=_mag*TMath::Cos(angle);
    }
    ////////////////////////////////////////////////////////////////////////////////
    // inline Double_t RooHSDWignerRe::evaluateMC(Double_t th,Double_t ph) const{
    //    if(th==_lastTheta&&ph==_lastPhi ) return _lastVal;
    //   _lastTheta=th;
    //   _lastPhi=ph;
  
    //   Double_t angle=_M*_lastPhi;
    //   return _lastVal=SmallWignerD( _L, _M, _S, th )*TMath::Cos(angle);
    // }

    ////////////////////////////////////////////////////////////////////////////////
    
    inline Bool_t RooHSDWignerRe::CheckClean() const
    {
      Double_t th=_theta;
      Double_t ph=_phi;
      if(th==_lastTheta&&ph==_lastPhi ) return kTRUE;
      _lastTheta=th;
      _lastPhi=ph;
      return kFALSE; 
    }

  }//FIT
}//HS 

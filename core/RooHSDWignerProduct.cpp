#include "RooHSDWignerProduct.h"

#include "TError.h"

namespace HS{
  namespace FIT{


    RooHSDWignerProduct::RooHSDWignerProduct(const char *name, const char *title,  RooAbsReal& dwigner1, RooAbsReal& dwigner2, Double_t factor)
      : RooHSComplex(name, title)
      ,_D1("dwigner1", "dwigner1", this, dwigner1)
      ,_D2("dwigner2", "dwigner2", this, dwigner2)
    {
      auto mag1=dynamic_cast<RooHSDWigner*>(&dwigner1);
      auto mag2=dynamic_cast<RooHSDWigner*>(&dwigner2);
     _M1=mag1->M();
     _M2=mag2->M();
     _phiName1=mag1->PhiName();
     _phiName2=mag2->PhiName();
    }


 ////////////////////////////////////////////////////////////////////////////////

    RooHSDWignerProduct::RooHSDWignerProduct(const RooHSDWignerProduct& other, const char* name)
      : RooHSComplex(other, name),
	_D1("dwigner1", this, other._D1),
	_D2("dwigner2", this, other._D2),
	_M1(other._M1),
	_M2(other._M2),
	_phiName1(other._phiName1),
	_phiName2(other._phiName2)
    {
   
    }



  ////////////////////////////////////////////////////////////////////////////////

    RooHSDWignerProductRe::RooHSDWignerProductRe(const char *name, const char *title,  RooAbsReal& phi1, RooAbsReal& phi2, RooAbsReal& dproduct )
      :  RooAbsReal(name, title),
	 _phi1("Phi1","Phi1",this,phi1),
	 _phi2("Phi2","Phi2",this,phi2),
	 _mag("Mag","Mag",this,dproduct)
    {
      auto mag=dynamic_cast<RooHSDWignerProduct*>(&dproduct);
      _M1=mag->M1();
      _M2=mag->M2();
    }

  ////////////////////////////////////////////////////////////////////////////////
    RooHSDWignerProductRe::RooHSDWignerProductRe(const RooHSDWignerProductRe& other, const char* name)
      : RooAbsReal(other, name),
	_phi1("Phi1", this, other._phi1),
	_phi2("Phi2", this, other._phi2),
	_mag("Mag", this, other._mag),
	_M1(other._M1),
	_M2(other._M2)
    {
   
    }


 ////////////////////////////////////////////////////////////////////////////////

    RooHSDWignerProductIm::RooHSDWignerProductIm(const char *name, const char *title,  RooAbsReal& phi1,  RooAbsReal& phi2, RooAbsReal& dproduct,Int_t conj)
      :  RooAbsReal(name, title),
	 _phi1("Phi1","Phi1",this,phi1),
	 _phi2("Phi2","Phi2",this,phi2),
	 _mag("Mag","Mag",this,dproduct),
 	 _conj(conj)
    {
      auto mag=dynamic_cast<RooHSDWignerProduct*>(&dproduct);
      _M1=mag->M1();
      _M2=mag->M2();
    }

  ////////////////////////////////////////////////////////////////////////////////
    RooHSDWignerProductIm::RooHSDWignerProductIm(const RooHSDWignerProductIm& other, const char* name)
      : RooAbsReal(other, name),
	_phi1("Phi1", this, other._phi1),
	_phi2("Phi2", this, other._phi2),
	_mag("Mag", this, other._mag),
	_M1(other._M1),
	_M2(other._M2),
	_conj(other._conj)
    {
   
    }
  
  }//FIT
}//HS

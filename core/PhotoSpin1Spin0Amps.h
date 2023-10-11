#pragma once

////////////////////////////////////////////////////////////////
///
///Class:               PhotoSpin1Spin0Amps
///Description: Formalism for photoproduction of gp->Xp with X->yz
///            where y and z are spin 1 and spin 0 respectively and g is polarised
#include <TString.h>
#include "FitManager.h"
#include "AmpConfigure.h"
#include "ComponentsPdfParser.h"

namespace HS{
  namespace FIT{

    class PhotoSpin1Spin0Amps : public AmpConfigure{

    public:
      PhotoSpin1Spin0Amps()=default;
      PhotoSpin1Spin0Amps(const std::string& name):_Name{name} {};
      PhotoSpin1Spin0Amps(const PhotoSpin1Spin0Amps&)=default;
      PhotoSpin1Spin0Amps(PhotoSpin1Spin0Amps&&)=default;
      ~PhotoSpin1Spin0Amps() =default;
      PhotoSpin1Spin0Amps& operator=(const PhotoSpin1Spin0Amps& other)=default;
      PhotoSpin1Spin0Amps& operator=(PhotoSpin1Spin0Amps&& other) = default;

     void SetManager(FitManager* manager){
	_Setup=&(manager->SetUp());
      }
     void SetSetup(Setup* setup){
       _Setup=setup;
      }
      void SetDecayAngleCosThGJ(const TString& var){
	_Setup->LoadVariable(var);
	_DecayAngleCosThGJ=var(0,var.First('['));
      }
      void SetDecayAnglePhiGJ(const TString& var){
	_Setup->LoadVariable(var);
	_DecayAnglePhiGJ=var(0,var.First('['));
      }
      void SetDecayAngleCosThHF(const TString& var){
	_Setup->LoadVariable(var);
	_DecayAngleCosThHF=var(0,var.First('['));
      }
      void SetDecayAnglePhiHF(const TString& var){
	_Setup->LoadVariable(var);
	_DecayAnglePhiHF=var(0,var.First('['));
      }
      void SetPolPhi(const TString& var){
	_Setup->LoadVariable(var);
	_PolPhi=var(0,var.First('['));
      }
      void SetPolarisation(const TString& var){
	_Setup->LoadVariable(var);
	_Polarisation=var(0,var.First('['));
      }
      void SetConstPolarisation(const TString& var){ //e.g. Pol[0.5]
	_Setup->LoadConstant(var);
	//	_Polarisation=var(0,var.First('['));
	//extract numerical value for polairsation
	_Polarisation=var(var.First('[')+1,var.First(']')-var.First('[')-1);
	_constPol=kTRUE;
      }


 void SetLmax(Int_t lm){
	_Lmax = lm;
	_Mmax = lm; //default = Lmax
      }
      void SetMmax(Int_t lm){
	_Mmax = lm;
      }
      void Setlmax(Int_t lm){
	_lmax = lm;
      }
      void Setmmax(Int_t lm){
	_mmax = lm;
      }
      void SetJmax(Int_t lm){
	_Jmax = lm;
      }
      void SetNrefl(Int_t r){
	_Nref = r;
      }
      void PositiveMOnly(){
	_NegativeM=kFALSE;
      }
      void SetOnlyEvenWaves(){
	_OnlyEven=kTRUE;
      }
      std::string  ConfigureMoments() override;
      
      std::string ConfigurePWAs() override;

      void LoadModelPDF(Long64_t Nevents=1) override;

      void PrintModel();
      
      const std::string&  GetName(){return _Name;}

      const std::string&  GetSummation(){return _Sum;}


   protected:

      ComponentsPdfParser PolWignerDFunctionProductMoments();
      
    private:

      std::string _DecayAngleCosThGJ;
      std::string _DecayAnglePhiGJ;
      std::string _DecayAngleCosThHF;
      std::string _DecayAnglePhiHF;
      std::string _PolPhi;
      std::string _Polarisation;
      std::string _Sum;
      std::string _Name;
      
      Int_t _Lmax=0;
      Int_t _Mmax=0;
      Int_t _lmax=0;
      Int_t _Jmax=0;
      Int_t _mmax=0;
      Int_t _Nref=2;
      Bool_t _OnlyEven=kFALSE;
      Bool_t _NegativeM=kTRUE;
      Bool_t _constPol=kFALSE;

      
    };
    
  }
}

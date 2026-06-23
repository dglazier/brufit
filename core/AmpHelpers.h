
////////////////////////////////////////////////////////////////
///
///Class:               AmpHelpers
///Description:
///           


#pragma once

#include <RooRandom.h>
#include <TMath.h>
#include "AmpConfigure.h"

namespace HS{
  namespace FIT{


    class AmpHelpers {

    public:
      
    AmpHelpers(AmpConfigure* configure):_configure{configure}{}

      void  ConfigAmps(Setup* setup);
      void  CopyToMomentPars();
      void  CopyToAmpPars();
      void  RandomiseFitParameters();
      void  RandomiseAmps();
      
      void SetSetup(Setup* setup){
	_fitSetup=setup;
      }
      
      void SetParRange(TString name, Double_t amin, Double_t amax){
	dynamic_cast<RooRealVar*>(_ampSetup.Parameters().find(name))->setRange(amin,amax);
    }

    private:
      
      AmpConfigure* _configure={nullptr}; // just for confuguring
      Setup _ampSetup; // just for configuring
      Setup* _fitSetup={nullptr}; //actual fit setup
      Bool_t _IsAmplitudes=kTRUE;
      Bool_t _IsConfigured=kFALSE;
      
    };

  
   inline  void  AmpHelpers::ConfigAmps(Setup* setup){
     if(_IsConfigured){
       _fitSetup=setup;
       return;//already done
     }
     _fitSetup=setup;

     
     _IsAmplitudes = _configure->IsAmplitudes(); //config should be configured for whatever we are fitting (moments or amps)
      //Copy parameters and formulas to _ampSetup
      auto oldSetup = _configure->GetSetup();
      _configure->SetSetup(&_ampSetup);
      
      //configure my setup
      _configure->ConfigurePWAs();

      //put back old setup
      _configure->SetSetup(oldSetup);
 
      _IsConfigured=kTRUE;
  
    }
    
    inline  void  AmpHelpers::CopyToMomentPars(){
      auto& ampMoments = _ampSetup.Formulas();
      auto& parMoments = _fitSetup->Parameters();
      // static_range_cast does not work until 6.28
      for(auto par:static_range_cast<RooRealVar *>(parMoments)){
    	auto* mom=dynamic_cast<RooFormulaVar*>(ampMoments.find(par->GetName()));
      	if((mom)!=nullptr){
	  // cout<<"copy formula value for "<<mom->GetName()<<" "<<mom->getVal()<<std::endl;
      	  par->setVal(mom->getVal());
      	}
      }

       
      return;
    }


    inline  void  AmpHelpers::CopyToAmpPars(){

      auto& ampPars = _ampSetup.Parameters();
      auto& myPars = _fitSetup->Parameters();
      // static_range_cast does not work until 6.28
      for(auto par:static_range_cast<RooRealVar *>(myPars)){
	if(par->isConstant())
	  continue;
      	auto* mom=dynamic_cast<RooRealVar*>(ampPars.find(par->GetName()));
    	if((mom)!=nullptr){
      	  //cout<<"copy formula value for "<<mom->GetName()<<" "<<mom->getVal()<<std::endl;
      	  par->setVal(mom->getVal());
      	}
      }

      
        return;
    }

    ///////////////////////////////////////////////////////
    inline  void AmpHelpers::RandomiseFitParameters(){
      RandomiseAmps();
      //std::cout<<"AmpHelpers::RandomiseFitParameters() "<<_IsAmplitudes<<std::endl;
      if(_IsAmplitudes)
	CopyToAmpPars();
      else
	CopyToMomentPars();
    }
    ////////////////////////////////////////////////////////
    
    inline  void AmpHelpers::RandomiseAmps(){
      //Take amplitude of form a_x_y... b_x_y... aphi_x_y... bphi_x_y...
      auto& pars = _ampSetup.Parameters();
      Double_t ampNorm=0;

      for(auto par:pars){
	auto realPar = dynamic_cast<RooRealVar*>(par);
	if(realPar==nullptr)
	  continue;
	if(realPar->isConstant())
	  continue;
			       
	realPar->randomize();

	if(TString(realPar->GetName()).Contains("phi")==kFALSE){//track magnitude so < 1
	  ampNorm+=realPar->getVal()*realPar->getVal();
	  realPar->setVal(TMath::Abs(realPar->getVal()));
	}
	else{//make sure intial phases within 2pi
	  realPar->setVal(std::remainder(realPar->getVal(),2*TMath::Pi()));
	}
   
      }
      //don't forget final normalisation parameter = sqrt(1 - others^2)
      auto normalisepar = RooRandom::uniform();
      ampNorm+=normalisepar*normalisepar;
 
      for(auto par:pars){
	auto realPar = dynamic_cast<RooRealVar*>(par);
	if(realPar==nullptr)
	  continue;
	if(TString(realPar->GetName()).Contains("phi")==kFALSE)
	  realPar->setVal(realPar->getVal()/TMath::Sqrt(ampNorm));
      }
      
    
    }//RandomAmps

 
  
    
  }
}

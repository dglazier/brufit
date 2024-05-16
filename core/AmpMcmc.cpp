#include "AmpMcmc.h"
#include "AmpConfigure.h"
#include <TMath.h>

namespace HS{
  namespace FIT{
    
    AmpMcmc::AmpMcmc(AmpConfigure* configure,Int_t Niter,Int_t Nburn, Float_t norm,UInt_t nrefits,Bool_t nozeroinit) : BruMcmcCovariance(Niter,Nburn,norm),fNFits{nrefits}, fNoZeroInitialVal{nozeroinit},_ampHelper{configure} {
      SetNameTitle("HSAmpMcmc","Mcmc multi fit  for amplitudes");
    }

    void AmpMcmc::Run(Setup &setup,RooAbsData &fitdata){
      cout<<"AmpMcmc::Run "<<GetName()<<" "<<fNFits<<endl;
      fSetup=&setup;
      fData=&fitdata;

      _ampHelper.ConfigAmps(fSetup);
 
     //loop over refits
      //BruMcmcSeqHelper::Run(setup,fitdata);
      UInt_t nrefit = 0;
      while(nrefit++<fNFits){
	SetName(Form("HSAmpMcmc_%d_",nrefit));
	cout<<"AmpMcmc::Run "<<nrefit<<" "<<GetName()<<endl;
	
	AmpMcmc::RandomiseParameters();
	
	BruMcmcCovariance::Run(setup,fitdata);

	if(Success()==kFALSE){
	  //failed so try again
	  nrefit--;
	  continue;
	}
	if(nrefit<fNFits) SaveInfo();
      }
    }

    void AmpMcmc::RandomiseParameters(){
      cout<<"AmpMcmc::RandomiseParameters() "<<fSetup<<endl;//exit(0);
      Double_t intensity0=0;
      //auto& pars = _ampSetup.Parameters();
 
      //Try 10000 times to get a physical starting value
      auto Nrand=10000;
      for(Int_t irand=0;irand<Nrand;++irand){

	_ampHelper.RandomiseFitParameters();
	
	if(fNoZeroInitialVal==kFALSE) //if don't care about 0 iniatal intensity break
	  break;
	
	// if we do care keep trying until we get non-zero intensity
	intensity0 = fSetup->Model()->getVal();
	std::cout<<"AmpMcmc::RandomiseParameters() rand "<<irand<<" "<<intensity0<<endl;
	if (intensity0!=0) break;
	if(irand==(Nrand-1)){
	  std::cerr<<"AmpMcmc::Run FATAL : error tried to randomise parameters 9999 times without a non-zero intensity and you specified not have a non zero inital intensity. You may try removing the true from the Minuit constructor if you do not need this. Note on irand = "<<irand<<std::endl;
	  exit(0);
	}
      }
      //fSetup->Parameters().Print("v");
       // std::cout<<"AmpMcmc::RandomiseParameters() for the "<<fIFit++<< " time, current intensity "<<intensity0<<std::endl;
       //exit(0);
    }

    ///////////////////////////////////////////////////////////////
    // void  AmpMcmc::ConfigAmps(AmpConfigure* config){

      
    //   _IsAmplitudes = config->IsAmplitudes(); //config should be configured for whatever we are fitting (moments or amps)

    //   _ampHelper.ConfigAmps(config);
      
    //   // auto copyconfig=config;
    //   // //Copy parameters and formulas to _ampSetup
    //   // copyconfig.SetSetup(&_ampSetup);
    //   // copyconfig.ConfigurePWAs();
    // }

    
    // void  AmpMcmc::CopyToMomentPars(){

    //   auto& ampMoments = _ampSetup.Formulas();
    //   auto& parMoments = fSetup->Parameters();
    //   // static_range_cast does not work until 6.28
    //   for(auto par:static_range_cast<RooRealVar *>(parMoments)){
    //   	auto* mom=dynamic_cast<RooFormulaVar*>(ampMoments.find(par->GetName()));
    //   	if((mom)!=nullptr){
    //   	  cout<<"copy formula value for "<<mom->GetName()<<" "<<mom->getVal()<<std::endl;
    //   	  par->setVal(mom->getVal());
    //   	}
    //   }

      
    //   return;
    // }


    // void  AmpMcmc::CopyToAmpPars(){

    //   auto& ampPars = _ampSetup.Parameters();
    //   auto& myPars = fSetup->Parameters();

    //   // static_range_cast does not work until 6.28
    //   for(auto par:static_range_cast<RooRealVar *>(myPars)){
    //   	auto* mom=dynamic_cast<RooFormulaVar*>(ampPars.find(par->GetName()));
    //   	if((mom)!=nullptr){
    //   	  cout<<"copy formula value for "<<mom->GetName()<<" "<<mom->getVal()<<std::endl;
    //   	  par->setVal(mom->getVal());
    //   	}
    //   }

      
    //   return;
    // }



    
  }
}

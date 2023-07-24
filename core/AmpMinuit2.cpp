#include "AmpMinuit2.h"
#include <RooRandom.h>
#include <TMath.h>

namespace HS{
  namespace FIT{
    
    AmpMinuit2::AmpMinuit2(AmpConfigure* configure,UInt_t nrefits,Bool_t nozeroinit) : Minuit2(nrefits,nozeroinit),_ampHelper{configure} {
      SetNameTitle("HSAmpMinuit2","Minuit2 minimiser for amplitudes");
    }


    void AmpMinuit2::RandomiseParameters(){
      // cout<<"AmpMinuit2::RandomiseParameters()"<<endl;exit(0);
      Double_t intensity0=0;
      //auto& pars = fSetup->Parameters();

      _ampHelper.SetSetup(nullptr); //make sure reconfigure
      _ampHelper.ConfigAmps(fSetup);
      // _ampHelper.SetSetup(fSetup);

      
      for(Int_t irand=0;irand<10000;++irand){
	//Amp::RandomiseAmps(pars);
	
 	_ampHelper.RandomiseFitParameters();

	if(fNoZeroInitialVal==kFALSE) //if don't care about 0 inital intensity break
	  break;
	// if we do care keep trying until we get non-zero intensity
	intensity0 = fSetup->Model()->getVal();
	if (intensity0!=0) break;
	if(irand==9999){
	  std::cerr<<"AmpMinuit2::Run FATAL : error tried to randomise parameters 9999 times without a non-zero intensity and you specified not have a non zero inital intensity. You may try removing the true from the Minuit constructor if you do not need this. Note on irand = "<<irand<<std::endl;
	  exit(0);
	}
      }
      fSetup->Parameters().Print("v");
      std::cout<<"AmpMinuit2::RandomiseParameters() for the "<<fIFit++<< " time, current intensity"<<intensity0<<std::endl;
      //   exit(0);
    }
  
  
  }
}

#include "AmpMinuit2.h"
#include <RooRandom.h>
#include <TMath.h>

namespace HS{
  namespace FIT{
    
    AmpMinuit2::AmpMinuit2(UInt_t nrefits,Bool_t nozeroinit) : Minuit2(nrefits,nozeroinit) {
      SetNameTitle("HSAmpMinuit2","Minuit2 minimiser for amplitudes");
    }


    void AmpMinuit2::RandomiseParameters(){
      // cout<<"AmpMinuit2::RandomiseParameters()"<<endl;exit(0);
      Double_t intensity0=0;
      auto& pars = fSetup->Parameters();

      auto RandomiseAmps = [&pars](){
			     Double_t ampNorm=0;
			     for(auto par:pars){
			       auto realPar = dynamic_cast<RooRealVar*>(par);
			       if(realPar==nullptr)
				 continue;
			       if(realPar->isConstant())
				 continue;
			       
			       realPar->randomize();

			       if(TString(realPar->GetName()).Contains("phi")==kFALSE)//track magnitude so < 1
				 ampNorm+=realPar->getVal()*realPar->getVal();
   
			     }
			     //don't forget final normalisation parameter = sqrt(1 - others^2)
			     auto normalisepar = RooRandom::uniform();
			     ampNorm+=normalisepar*normalisepar;
 
			     if(ampNorm>1){
			       for(auto par:pars){
				 auto realPar = dynamic_cast<RooRealVar*>(par);
				 if(realPar==nullptr)
				   continue;
				 if(TString(realPar->GetName()).Contains("phi")==kFALSE)
				   realPar->setVal(realPar->getVal()/TMath::Sqrt(ampNorm));
			       }
			     }
			   };


      
      for(Int_t irand=0;irand<10000;++irand){
	RandomiseAmps();
 
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
      pars.Print("v");
      std::cout<<"AmpMinuit2::RandomiseParameters() for the "<<fIFit++<< " time, current intensity"<<intensity0<<std::endl;
      //   exit(0);
    }
  }

}

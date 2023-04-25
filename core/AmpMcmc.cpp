#include "AmpMcmc.h"
#include <RooRandom.h>
#include <TMath.h>

namespace HS{
  namespace FIT{
    
    AmpMcmc::AmpMcmc(Int_t Niter,Int_t Nburn, Float_t norm,UInt_t nrefits,Bool_t nozeroinit) : RooMcmcSeqHelper(Niter,Nburn,norm),fNFits{nrefits}, fNoZeroInitialVal{nozeroinit} {
      SetNameTitle("HSAmpMcmc","Mcmc multi fit  for amplitudes");
    }

    void AmpMcmc::Run(Setup &setup,RooAbsData &fitdata){
      cout<<"AmpMcmc::Run "<<GetName()<<" "<<fNFits<<endl;
      fSetup=&setup;
      fData=&fitdata;
     //loop over refits
      //RooMcmcSeqHelper::Run(setup,fitdata);
      UInt_t nrefit = 0;
      while(nrefit++<fNFits){
	SetName(Form("HSAmpMcmc_%d_",nrefit));
	cout<<"AmpMcmc::Run "<<nrefit<<" "<<GetName()<<endl;
	//	RooMcmc::SetupBasicUsage();
	RandomiseParameters();
	
	AmpMcmc::RandomiseParameters();
	RooMcmcSeqHelper::Run(setup,fitdata);
	
	if(nrefit<fNFits) SaveInfo();
      }
    }

    void AmpMcmc::RandomiseParameters(){
      cout<<"AmpMcmc::RandomiseParameters() "<<fSetup<<endl;//exit(0);
      Double_t intensity0=0;
      auto& pars = fSetup->Parameters();
      pars.Print("v");
      auto RandomiseAmps = [&pars](){
			     cout<<"AmpMcmc::RandomiseAmps()"<<endl;
			     Double_t ampNorm=0;
			     for(auto par:pars){
			       cout<<par->GetName()<<endl;
			       auto realPar = dynamic_cast<RooRealVar*>(par);
			       if(realPar==nullptr)
				 continue;
			       if(realPar->isConstant())
				 continue;
			       
			       realPar->randomize();

			       if(TString(realPar->GetName()).Contains("phi")==kFALSE)//track magnitude so < 1
				 ampNorm+=realPar->getVal()*realPar->getVal();
   
			     }
			     cout<<"AmpMcmc::RandomiseAmps()"<<endl;
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
			     cout<<"AmpMcmc::RandomiseAmps()"<<endl;
			   };


      
      for(Int_t irand=0;irand<10000;++irand){
	RandomiseAmps();
 
	if(fNoZeroInitialVal==kFALSE) //if don't care about 0 inital intensity break
	  break;
	// if we do care keep trying until we get non-zero intensity
	intensity0 = fSetup->Model()->getVal();
	if (intensity0!=0) break;
	if(irand==9999){
	  std::cerr<<"AmpMcmc::Run FATAL : error tried to randomise parameters 9999 times without a non-zero intensity and you specified not have a non zero inital intensity. You may try removing the true from the Minuit constructor if you do not need this. Note on irand = "<<irand<<std::endl;
	  exit(0);
	}
      }
      pars.Print("v");
      std::cout<<"AmpMcmc::RandomiseParameters() for the "<<fIFit++<< " time, current intensity"<<intensity0<<std::endl;
      //   exit(0);
    }
  }

}

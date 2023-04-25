////////////////////////////////////////////////////////////////
///
///Class:               AmpMinuit2
///Description:
///           

#include "Minimiser.h"


#pragma once
namespace HS{
  namespace FIT{


    class AmpMinuit2  : public Minuit2 {
      
    public:

      AmpMinuit2(UInt_t nrefits=0,Bool_t nozeroinit=kFALSE);
      AmpMinuit2(const AmpMinuit2&)=default;
      AmpMinuit2(AmpMinuit2&&)=default;
      ~AmpMinuit2() override =default;
      AmpMinuit2& operator=(const AmpMinuit2& other)=default;
      AmpMinuit2& operator=(AmpMinuit2&& other) = default;  

      /* void FitTo() override { */
      /* 	auto fitOptions=fSetup->FitOptions(); */
      /* 	fitOptions.Add(dynamic_cast<RooCmdArg*>(RooFit::Minimizer("Minuit2").Clone())); */
      /* 	fResult=fSetup->Model()->fitTo(*fData,fitOptions); */
      /* }; */
      
      void RandomiseParameters() override;
      UInt_t fIFit=0;
 
      ClassDefOverride(HS::FIT::AmpMinuit2,1);
      
    };

  }//namespaces
}

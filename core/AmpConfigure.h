////////////////////////////////////////////////////////////////
///
///Class:               AmpConfigure
///Description:  Interfaces to amplitude configures, which can be used
///              in minimisation classes
///           


#pragma once

#include "Setup.h"

namespace HS{
  namespace FIT{


    class AmpConfigure {

    public:
      
      virtual void LoadModelPDF(Long64_t Nevents=1) = 0;
      virtual std::string  ConfigureMoments() = 0;   
      virtual std::string ConfigurePWAs() = 0;

      void SetSetup(Setup* setup){
	_Setup=setup;
      }

      Setup* GetSetup() {return _Setup;}
      
     Bool_t IsAmplitudes() const {return _IsAmplitudes;}

    protected:
      
     Setup* _Setup={nullptr};//!
     Bool_t _IsAmplitudes=kTRUE;

 
    };

  }
}

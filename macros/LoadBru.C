#include <TSystem.h>
#include <TString.h>
#include <TInterpreter.h>
#include <TROOT.h>
#include <iostream>

namespace HS{namespace FIT{namespace PROCESS{}}};
using namespace HS;
using namespace HS::FIT;
using namespace HS::FIT::PROCESS;

void LoadBru(TString Selection=""){
  gSystem->Load("libRooStats");
  gSystem->Load("libProof");
  
  TString BRUCODE=gSystem->Getenv("BRUFIT");
  TString fitpath=BRUCODE+"/core";
  TString macpath=BRUCODE+"/macros";
  TString utpath=BRUCODE+"/utility";
  // TString dmpath="/hsdata";
   
  if(!TString(gInterpreter->GetIncludePath()).Contains(fitpath)){
    gInterpreter->AddIncludePath(fitpath);
    gROOT->SetMacroPath(Form("%s:%s",gROOT->GetMacroPath(),(fitpath).Data()));
    gROOT->SetMacroPath(Form("%s:%s",gROOT->GetMacroPath(),(macpath).Data()));
    gROOT->SetMacroPath(Form("%s:%s",gROOT->GetMacroPath(),(utpath).Data()));
    gSystem->Load(BRUCODE+"/lib/libbrufit.so");
  }
  


}

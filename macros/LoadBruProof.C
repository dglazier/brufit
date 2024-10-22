#include <TSystem.h>
#include <TString.h>
#include <TInterpreter.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TProof.h>
#include <iostream>

namespace HS{namespace FIT{namespace PROCESS{}}};
using namespace HS;
using namespace HS::FIT;
using namespace HS::FIT::PROCESS;

void LoadBruProof(Int_t Nworkers=1,TString Selection=""){

  TString PWD=gSystem->Getenv("PWD");

  gSystem->Load("libRooStats");
  gSystem->Load("libProof");
  gSystem->Load("libMathMore");

  TString BRUCODE=gSystem->Getenv("BRUFIT");
  TString fitpath=BRUCODE+"/core";
  TString libpath=BRUCODE+"/lib";
  TString libext=(TString(gSystem->GetBuildCompilerVersion()).Contains("darwin") or TString(gSystem->GetBuildArch()).Contains("macos"))
                  ? ".dylib" : ".so";

  if(!TString(gInterpreter->GetIncludePath()).Contains(fitpath)){
    gInterpreter->AddIncludePath(fitpath);
    gROOT->SetMacroPath(Form("%s:%s",gROOT->GetMacroPath(),(fitpath).Data()));
    gSystem->Load(BRUCODE+"/lib/libbrufit"+libext);
  }

  gSystem->Load("libProof");
  TProof *proof =nullptr;
  if(!gProof)
    proof = TProof::Open("://lite");
  else
    proof=gProof;

  Int_t NCores=Nworkers;
  proof->SetParallel(NCores);


  proof->AddIncludePath(fitpath);
  proof->AddIncludePath(PWD);

   // get the sandbox directroy
  TString sandbox="~/.proof";
  if(TString(gEnv->GetValue("ProofLite.Sandbox",""))!=TString()){
    sandbox=gEnv->GetValue("ProofLite.Sandbox","");
  }
  gSystem->Exec(Form("cp $BRUFIT/lib/libbrufit_rdict.pcm %s/cache/.",sandbox.Data()));
  gProof->Load(gSystem->DynamicPathName("libRooStats"),kTRUE);
  gProof->Load(BRUCODE+"/lib/libbrufit"+libext,kTRUE);

}

#include "EquationSolver.h"

#include "ToyManager.h"
#include "PhotoTwoSpin0Amps.h"

HS::FIT::Setup& ConfigureAmpsNoValues(UInt_t Lmax=1,UInt_t Mmax=1,UInt_t Ref=2){

  using namespace HS::FIT;
  
  auto Generator =new  ToyManager(1); //number of toys to generate (must have generated enough phase space events!)

  //set the output directory for the file Toy*.root to go
  Generator->SetUp().SetOutDir("genBruAmps/");

  //Use amlitude configue class to define model
  PhotoTwoSpin0Amps config("PWA");
  config.SetManager(Generator);
  //set data variables which must be in the simulated tree (see below)
  config.SetDecayAngleCosTh("CosTh[0.21,-1,1]");
  config.SetDecayAnglePhi("Phi[0.2,-3.14159,3.14159]");
  config.SetPolPhi("PolPhi[0.2,-3.14159,3.14159]");

  //if polarisation is in data tree
  config.SetPolarisation("Pol[0.9,0.5,1]");

 
  //Now set model options
  //Lmax
  config.SetLmax(Lmax);
  
  //Mmax = Lmax if not set
  config.SetMmax(Mmax);
  
  //number of reflectivities = 1 or 2
  config.SetNrefl(Ref);

  //Only use even , S,D,... waves
  // config.SetOnlyEvenWaves();

  //Load required functions
  config.ConfigurePWAs();

  //Load fit PDF to generate 1E4 events
  config.LoadModelPDF(1E4);

 
 
  return Generator->SetUp();  
 
}

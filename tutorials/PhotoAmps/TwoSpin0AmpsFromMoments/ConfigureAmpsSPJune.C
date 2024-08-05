#include "EquationSolver.h"

#include "ToyManager.h"
#include "PhotoTwoSpin0Amps.h"

HS::FIT::Setup& ConfigureAmpsSPJune(){

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
  config.SetLmax(1);
  
  //Mmax = Lmax if not set
  config.SetMmax(1);
  
  //number of reflectivities = 1 or 2
  config.SetNrefl(2);

  //Only use even , S,D,... waves
  // config.SetOnlyEvenWaves();

  //Load required functions
  config.ConfigurePWAs();

  //Load fit PDF to generate 1E4 events
  config.LoadModelPDF(1E4);

 
  /******* Amp values **********/
  // S+0+ polar 0.36231 0.0 real  = a_0_0 aphi_0_0
  // S+0- polar 0.20329 0.0 real  = b_0_0 bphi_0_0
  // P+1+ polar 0.52213 2.71276   = a_1_1 aphi_1_1
  // P+1- polar 0.24958 -0.77666  = b_1_1 bphi_1_1
  // P+0+ polar 0.09433 -1.77094  = a_1_0 aphi_1_0
  // P+0- polar 0.23762 0.53551   = b_1_0 bphi_1_0
  // P-1+ polar 0.37141 -2.38804  = a_1_-1 aphi_1_-1
  // P-1- polar 0.53776 -0.62676  = b_1_-1 bphi_1_-1
  Generator->SetUp().SetParVal("aphi_0_0",0,kTRUE); //fix S real
  Generator->SetUp().SetParVal("aphi_1_1",2.71276,kTRUE); //P1 real
  Generator->SetUp().SetParVal("aphi_1_0",-1.77094,kFALSE); //P0
  Generator->SetUp().SetParVal("aphi_1_-1",-2.38804,kFALSE); //P-1

  Generator->SetUp().SetParVal("bphi_0_0",0,kFALSE); //fix S real
  Generator->SetUp().SetParVal("bphi_1_1",-0.77666,kTRUE); //P1 real
  Generator->SetUp().SetParVal("bphi_1_0",0.53551,kFALSE); //P0
  Generator->SetUp().SetParVal("bphi_1_-1",-0.62676,kFALSE); //P-1

  Generator->SetUp().SetParVal("a_0_0",0.36231,kFALSE); //S
  Generator->SetUp().SetParVal("a_1_-1",0.37141,kFALSE); //P-1
  Generator->SetUp().SetParVal("a_1_0",0.09433,kFALSE); //P0
  Generator->SetUp().SetParVal("b_0_0",0.20329,kFALSE); //S
  Generator->SetUp().SetParVal("b_1_-1",0.53776,kFALSE); //P-1
  Generator->SetUp().SetParVal("b_1_0",0.23762,kFALSE); //P0
  Generator->SetUp().SetParVal("b_1_1",0.24958,kFALSE); //P1
  //Note a_1_1 (the last +ve ref amp) is fixed to 1 - sum of squares of others
  // Generator->SetUp().SetParVal("a_1_1",0.39582,kFALSE); //P1

 
  return Generator->SetUp();  
 
}

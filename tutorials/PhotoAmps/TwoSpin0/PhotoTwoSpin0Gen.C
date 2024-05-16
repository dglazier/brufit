{
  ToyManager Generator(1); //number of toys to generate (must have generated enough phase space events!)

  //set the output directory for the file Toy*.root to go
  Generator.SetUp().SetOutDir("genBruAmps/");

  //Use amlitude configue class to define model
  PhotoTwoSpin0Amps config("PWA");
  config.SetManager(&Generator);
  //set data variables which must be in the simulated tree (see below)
  config.SetDecayAngleCosTh("CosTh[0.21,-1,1]");
  config.SetDecayAnglePhi("Phi[0.2,-3.14159,3.14159]");
  config.SetPolPhi("PolPhi[0.2,-3.14159,3.14159]");

  //if polarisation is in data tree
  config.SetPolarisation("Pol[0.9,0.5,1]");

  //if polarisation has fixed value
  //config.SetConstPolarisation("Pol[0.75]");
 
  //In case using weights etc.
  Generator.SetUp().SetIDBranchName("UID"); 

  //give the phase space file to generate from
  //You may run > brufit GeneratePhaseSpace.C , to create flat/Toy0.root
  Generator.LoadSimulated("ToyData","flat/Toy0.root",config.GetName());

  //Now set model options
  //Lmax
  config.SetLmax(2);
  
  //Mmax = Lmax if not set
  config.SetMmax(1);
  
  //number of reflectivities = 1 or 2
  config.SetNrefl(1);

  //Only use even , S,D,... waves
  config.SetOnlyEvenWaves();

  //Load required functions
  config.ConfigurePWAs();

  //Load fit PDF to generate 1E4 events
  config.LoadModelPDF(1E4);

  //set truth amplitudes
  //Notation : refl_l_m , where reflevtivity = 'a'(+) or 'b'(-)
  //phases
  Generator.SetUp().SetParVal("aphi_0_0",0,kTRUE); //fix S real
  Generator.SetUp().SetParVal("aphi_2_-1",15.4*TMath::DegToRad(),kFALSE); //D-1
  Generator.SetUp().SetParVal("aphi_2_0",174*TMath::DegToRad(),kFALSE); //D0
  Generator.SetUp().SetParVal("aphi_2_1",-81.6*TMath::DegToRad(),kFALSE); //D+1
  //magnitudes
  Generator.SetUp().SetParVal("a_0_0",0.499,kFALSE); //S
  Generator.SetUp().SetParVal("a_2_-1",0.201,kFALSE); //D-1
  Generator.SetUp().SetParVal("a_2_0",0.567,kFALSE); //D0
  //!!Note we do not set a_2_1 as it is fixed by the condition
  // a_0_0^2 +a_2_-1^2 + a_2_0^2 + a_2_1^2 =1
  //It is always the highest l, m +ve reflectivity which is constrained
  //Generator.SetUp().SetParVal("a_2_1",0.624,kFALSE); //D+1

  //if using -ve reflectivity also fix its S wave phase
  //Fitter.SetUp().SetParVal("bphi_0_0",0,kTRUE);
  //etc

  std::cout<<"DONE CONFIGURING "<<endl;
  //Generate the events
  Here::Go(&Generator);

 
}

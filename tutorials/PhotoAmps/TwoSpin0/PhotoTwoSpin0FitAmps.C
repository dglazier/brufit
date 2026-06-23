{
  
  FitManager Fitter;// manage the fitting
  //set the output directory for the fit results files Results*.root
  Fitter.SetUp().SetOutDir("fitBruAmps/");

  //Use amlitude configue class to define model
  PhotoTwoSpin0Amps config("PWA");
  config.SetManager(&Fitter);
   //set data variables which must be in the input tree
  config.SetDecayAngleCosTh("CosTh[0.21,-1,1]");
  config.SetDecayAnglePhi("Phi[0.2,-3.14159,3.14159]");
  config.SetPolPhi("PolPhi[0.2,-3.14159,3.14159]");
  config.SetPolarisation("Pol[0.3,0.7]");
  //config.SetConstPolarisation("Pol[0.75]"); //alternative
 
  //In case using weights etc.
  Fitter.SetUp().SetIDBranchName("UID"); 

  //load simulated data for normalisation integral
  //treename, filename, PDF name
  Fitter.LoadSimulated("ToyData","flat/Toy0.root",config.GetName());

  //load data to be fit (this was created by PhotoTwoSpin0Gen.C)
  Fitter.LoadData("ToyData","genBruAmps/Toy0.root");
 
  //Now set model options
  //Lmax
  config.SetLmax(2);
  //Mmax = Lmax if not set
  config.SetMmax(1);
  //config.SetLmax(3);
  //Mmax = Lmax if not set
  //config.SetMmax(3);
  //number of reflectivities = 1 or 2
  config.SetNrefl(1);
 //Only use even , S,D,... waves
  //config.SetOnlyEvenWaves();
   
  //Load required functions
  config.ConfigurePWAs();

  //Load fit PDF
  config.LoadModelPDF();

  //fix S real
  Fitter.SetUp().SetParVal("aphi_0_0",0,kTRUE);//kTRUE=>constant

  //////////////////**could start fir with true values
  //if using -ve reflectivity also fix its S wave phase
  //Fitter.SetUp().SetParVal("bphi_0_0",0,kTRUE);
  //set truth amplitudes
  //Notation : refl_l_m , where reflevtivity = 'a'(+) or 'b'(-)
  //phases
  // Fitter.SetUp().SetParVal("aphi_0_0",0,kTRUE); //fix S real
  // Fitter.SetUp().SetParVal("aphi_2_-1",15.4*TMath::DegToRad(),kFALSE); //D-1
  // Fitter.SetUp().SetParVal("aphi_2_0",174*TMath::DegToRad(),kFALSE); //D0
  // Fitter.SetUp().SetParVal("aphi_2_1",-81.6*TMath::DegToRad(),kFALSE); //D+1
  // //magnitudes
  // Fitter.SetUp().SetParVal("a_0_0",0.499,kFALSE); //S
  // Fitter.SetUp().SetParVal("a_2_-1",0.201,kFALSE); //D-1
  // Fitter.SetUp().SetParVal("a_2_0",0.567,kFALSE); //D0
  

 
  //some plotter options
  // Fitter.TurnOffPlotting();
  // Fitter.SetPlotOptions("MCMC"); //Make MCMC related plots
  //Fitter.SetPlotOptions("goff"); //save plots but do not show (batch)
 
  // (dynamic_cast<BruEventsPDF*>(&Fitter.SetUp().PDFs()[0]))->SetConstInt();
  // (dynamic_cast<RooHSEventsPDF*>(&Fitter.SetUp().PDFs()[0]))->SetConstInt();
 //********************************************
  //Perform single fit with default Minuit2 minimiser
  //Here::Go(&Fitter);
 
  //********************************************
  //Perform fit 20 times Minuit2 minimiser
  //All results are saved in same Results file in the TTree ResultTreeBru
  //Fitter.SetMinimiser(new AmpMinuit2(&config,5));
  //Here::Go(&Fitter);
  
  //********************************************
  //Perform "fit" with an MCMC sampler
  // a tree MCMCTree is included in the Results*.root file
  //most basic sequential proposal (Nsamples,burnin,step size, desired acceptance, min acceptance, max acceptance)
  //auto mcmc=new BruMcmcSeqHelper(2000,1000,0.1,0.23,0.16,0.3);
  //brufit covariance matric based proposal
  // auto mcmc=new BruMcmcCovariance({2000,10000,5000},100,1,0.23,0.16,0.3);
  // mcmc->SetCyclicParameters(Fitter.SetUp().FilterParameters("phi"));
  // Fitter.SetMinimiser(mcmc);
  // Here::Go(&Fitter);
 
  //********************************************
  //Perform "fit" with an MCMC sampler with multiple chains
  //Nsamples,burnin,step size,NChains
  auto mcmc=new AmpMcmc(&config,{5000,10000,5000},2,100,0.1,0.23,0.16,0.3,false);
  mcmc->SetCyclicParameters(Fitter.SetUp().FilterParameters("phi"));
  Fitter.SetMinimiser(mcmc);
  Here::Go(&Fitter);

  
 }

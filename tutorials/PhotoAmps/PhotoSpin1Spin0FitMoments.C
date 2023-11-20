{
  FitManager Fitter;// manage the fitting
  //set the output directory for the fit results files Results*.root
  Fitter.SetUp().SetOutDir("fitBruMoments/");

  //Use amlitude configue class to define model
  PhotoTwoSpin0Amps config("PWA");
  config.SetManager(&Fitter);
  //set data variables which must be in the input tree
  config.SetDecayAngleCosThGJ("CosThGJ[0.21,-1,1]");
  config.SetDecayAnglePhiGJ("PhiGJ[0.2,-3.14159,3.14159]");
  config.SetDecayAngleCosThHF("CosThHF[0.21,-1,1]");
  config.SetDecayAnglePhiHF("PhiHF[0.2,-3.14159,3.14159]");
  config.SetPolPhi("PolPhi[0.2,-3.14159,3.14159]");
  config.SetPolarisation("Pol[0.9,0.5,1]");

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
  config.SetMmax(2);
  //Jmax 
  config.SetJmax(2);
  //lmax
  config.Setlmax(1);//Must be 1 for Spin1Spin0
  //mmax
  config.Setmmax(1);//Must be 1 for Spin1Spin0
  //number of reflectivities = 1 or 2
  config.SetNrefl(1);
  //Only use even , S,D,... waves
  // config.SetOnlyEvenWaves();


  //Load required functions
  config.ConfigureMoments();

  //Load fit PDF
  config.LoadModelPDF();
  

  //set some fitter options
  Fitter.SetUp().AddFitOption(RooFit::PrintEvalErrors(-1));//suppress error messaages
  Fitter.SetUp().AddFitOption(RooFit::NumCPU(6)); //number of CPUs to split likelihood calc.
  
  //default error strategy for Minuit fits is asymptotically correct approach
  //https://arxiv.org/abs/1911.01303, but this may be slow
  Fitter.SetUp().ErrorsWrong();//"naive" error calculation, much faster
  //Fitter.SetUp().ErrorSumW2();//sumW2 correction if using weights

  //some plotter options
  // Fitter.TurnOffPlotting();
  // Fitter.SetPlotOptions("MCMC"); //Make MCMC related plots
  // Fitter.SetPlotOptions("goff"); //save plots but do not show (batch)

  //********************************************
  //Perform fit with default Minuit2 minimiser
  // Here::Go(&Fitter);
 
  //********************************************
  //Perform fit 20 times Minuit2 minimiser
  //All results are saved in same Results file in the TTree ResultTreeBru
  //Fitter.SetMinimiser(new AmpMinuit2(&config,10));
  //Here::Go(&Fitter);
  
  //********************************************
  //Perform "fit" with an MCMC sampler
  // a tree MCMCTree is included in the Results*.root file
  //most basic sequential proposal (Nsamples,burnin,step size, desired acceptance, min acceptance, max acceptance)
  //auto mcmc=new BruMcmcSeqHelper(2000,1000,0.1,0.23,0.16,0.3);
  //brufit covariance matric based proposal
  auto mcmc=new BruMcmcCovariance(10000,1000,0.1,0.23,0.16,0.3);
  ////mcmc->TurnOffCovariance();//BruMcmcCovariance only, do not proceed with covariance based sampling, just perform basic stepping
  Fitter.SetMinimiser(mcmc);
  Here::Go(&Fitter);
 
  //********************************************
  //Perform "fit" with an MCMC sampler with multiple chains
  //Nsamples,burnin,step size,NChains
  //  auto mcmc=new AmpMcmc(&config,10000,1000,0.01,5);
  //mcmc->TurnOffCovariance();//BruMcmcCovariance only, do not proceed with covariance based sampling, just perform basic stepping
  //Fitter.SetMinimiser(mcmc);
  //Here::Go(&Fitter);

  
}

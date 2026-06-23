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
  config.SetPolarisation("Pol[0.9,0.5,1]");
   
  //config.SetConstPolarisation("Pol[0.75]"); //alternative
 
  //In case using weights etc.
  Fitter.SetUp().SetIDBranchName("UID"); 

  //load simulated data for normalisation integral
  //treename, filename, PDF name
  Fitter.LoadSimulated("ToyData","/home/dglazier/Dropbox/HaSpect/dev/brufit/tutorials/PhotoAmps/TwoSpin0/flat/Toy0.root",config.GetName());

  //load data to be fit (this was created by PhotoTwoSpin0Gen.C)
  Fitter.LoadData("ToyData","/home/dglazier/Dropbox/HaSpect/dev/brufit/tutorials/PhotoAmps/TwoSpin0/genBruAmps/Toy0.root");
 
  //Now set model options
  //Lmax
  // config.SetLmax(8);
  config.SetLmax(2);
  //Mmax = Lmax if not set
  //config.SetMmax(4);
  config.SetMmax(2);
  //number of reflectivities = 1 or 2
  config.SetNrefl(2);
 //Only use even , S,D,... waves
 // config.SetOnlyEvenWaves();
   
  //Load required functions
  config.ConfigureMoments();

  //Load fit PDF
  config.LoadModelPDF();

  
  //some plotter options
  // Fitter.TurnOffPlotting();
  //Fitter.SetPlotOptions("MCMC:AUTOCORR"); //Make MCMC related plots
  //Fitter.SetPlotOptions("goff"); //save plots but do not show (batch)

  //********************************************
  //Perform fit with default Minuit2 minimiser
  // Here::Go(&Fitter);
 
  //********************************************
  //Perform fit 10 times Minuit2 minimiser
  //All results are saved in same Results file in the TTree ResultTreeBru
  //Fitter.SetMinimiser(new AmpMinuit2(&config,10));
  //Here::Go(&Fitter);
  //Proof::Go(&Fitter,1);

  //********************************************
  //Perform "fit" with an MCMC sampler
  // a tree MCMCTree is included in the Results*.root file
  //most basic sequential proposal (Nsamples,burnin,step size, desired acceptance, min acceptance, max acceptance)
  //auto mcmc=new BruMcmcSeqHelper(2000,1000,0.1,0.23,0.16,0.3);
  //brufit covariance matric based proposal
  // give a vector of number of iterations for each phase {}
  // auto mcmc=new BruMcmcCovariance({200,2000,500},50,0.1,0.23,0.16,0.3);
  ////mcmc->TurnOffCovariance();//BruMcmcCovariance only, do not proceed with covariance based sampling, just perform basic stepping
  //Fitter.SetMinimiser(mcmc);
  //Here::Go(&Fitter);
  //std::vector<Int_t> Niters,Int_t Nburn=10, Float_t norm=0.01,float target=0.234,float accmin=0.15,float accmax=0.35)
  //  auto mcmc=new BruMcmcCovariance({10000,100000,200000},100,1,0.23,0.16,0.3);
  auto mcmc=new BruMcmcCovariance({5000,20000,10000},100,1,0.23,0.16,0.3);
  Fitter.SetMinimiser(mcmc);
  Here::Go(&Fitter);
 
  //********************************************
  //Perform "fit" with an MCMC sampler with multiple chains
  //Nsamples,burnin,step size,NChains
  //  auto mcmc=new AmpMcmc(&config,{5000,10000,5000},100,10);
  //mcmc->TurnOffCovariance();//BruMcmcCovariance only, do not proceed with covariance based sampling, just perform basic stepping
  //Fitter.SetMinimiser(mcmc);
  //Here::Go(&Fitter);

  
}

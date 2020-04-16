{
  /////////////////////////////////////////////////////////////////
  ////Fit the data
  FitManager Fitter;
  Fitter.SetUp().SetOutDir("outSphHarmonic/");

  Fitter.SetUp().LoadVariable("CosTh[0,-1,1]");
  Fitter.SetUp().LoadVariable("Phi[-3.14159,3.14159]");

   auto configFitPDF=HS::FIT::EXPAND::ComponentsRealSphHarmonic(Fitter.SetUp(),"Moments","CosTh","Phi",3,2);


  Fitter.SetUp().FactoryPDF(configFitPDF);
   
  
  Fitter.SetUp().LoadSpeciesPDF("Moments",1); //2000 events


  //Get the generated data
  Fitter.LoadData("ToyData","genSphHarmonic/Toy0.root");
  //flat data for mc integration
  Fitter.LoadSimulated("ToyData","flatSphHarmonic/Toy0.root","Moments");


  //could try MCMC
  gBenchmark->Start("fit ");
  Fitter.SetMinimiser(new RooMcmcSeq(2000,1000,50));
  Fitter.SetUp().AddFitOption(RooFit::Optimize(1));
  //Here::Go(&Fitter); //try MCMC comment in here
  gBenchmark->Stop("fit ");
  gBenchmark->Print("fit ");


  // Or just Minuit
  gBenchmark->Start("minuit ");
  Fitter.SetMinimiser(new Minuit2());
  Fitter.SetUp().AddFitOption(RooFit::Optimize(1));
  Here::Go(&Fitter);
  gBenchmark->Stop("minuit ");
  gBenchmark->Print("minuit ");
}

{

  //We want to generate toy datasets based on previous fit results
  //Give previous fir directory
  TString fitDir("/hdd/Dropbox/HaSpect/dev/brufit/tutorials/AnalyseToyMC/outFitToData/");
  //Give minimiser (we may have done fits with different minimsers)
  TString fitMinimiser("HSMinuit2");
  //  TString fitMinimiser("RooMcmcSeqThenCov");

  //Open file with saved brufit fitmanager
  //auto sdir=gDirectory;
  std::unique_ptr<TFile> fitFile(new TFile(fitDir+"/HSFit.root"));

  //Get fit manager (this is configured with your PDF and bining etc.)
  auto prevFit= *fitFile->Get<FitManager>("HSFit");
  prevFit.InitPrevResult(fitDir,fitMinimiser);

  //Create toy manager to generate 10 datasets
  auto toy=ToyManager::GetFromFit(10,prevFit,prevFit.MinimiserFileName());
  toy->SetUp().SetOutDir("outCompToy"); //give different output directory
  //Load simulated data for filtering events from
  toy->LoadData("MyModel","DataSignal.root");//need data for protodata PolState and Pol
  toy->LoadSimulated("MyModel","MC.root","SigAsym");
  
 
  
  Here::Go(toy); //generate toys

  //git a new fit mamanger configured from the ToyManager
  //(you could just create a new fit manager yourself)
  auto toyfitter= toy->Fitter();
  toyfitter->IgnorePrevResult(); //make sure initial fit parameters are not equal to previous fit results

  
  Here::Go(toyfitter.get());
  // Proof::Go(toyfitter.get(),4);
  //create pull distributions etc
  toy->Summarise();
}

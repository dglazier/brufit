//brufit  FitHSMCModelBins.C
{
  //PROOF needs full paths!
  
  sPlot RF;
  RF.SetUp().SetOutDir("outPeakHistBins/");
  ///////////////////////////////Load Variables
  RF.SetUp().LoadVariable("Mmiss[0,9.5]");//should be same name as variable in tree
  RF.SetUp().LoadAuxVar("Eg[0,10]");
  // RF.SetUp().AddCut("Eg<3.2");

  RF.SetUp().SetIDBranchName("fgID");

  //Fit  Signal with smoothed and interpolated histogram pdf with 1000 bins, and 5000 for normalisation integral
  Int_t smooth = 1;
  Int_t interpolate = 1;
  Int_t NPdfBins = 100;
  Int_t NNormBins = 50000;
 
  //////////////////////////////Make signal PDF
  auto pdfString = Form("BruEventsHistPeakPDF::Signal(Mmiss, smear_Sig[0,0,0.1], off[0,-0.05,0.05], scale[1,0.95,1.05],%d,%d,%d,%d)",smooth,interpolate,NPdfBins,NNormBins);
  RF.SetUp().FactoryPDF(pdfString);
  
  RF.SetUp().LoadSpeciesPDF("Signal",1); 

  //////////////////////////////Make background PDF
  RF.SetUp().FactoryPDF("RooHSEventsHistPDF::BG(Mmiss,alphaB[0,0,5],offB[0,0,0],scaleB[1.0,0.8,1.2])");
  RF.SetUp().LoadSpeciesPDF("BG",1);

  ////////////////////////////Make Bins
  RF.Bins().LoadBinVar("Eg",5,3,4);

  ///////////////////////////Load Data
  RF.LoadData("MyModel","Data.root");
  RF.LoadSimulated("MyModel","SigData.root", "Signal");
  RF.LoadSimulated("MyModel","BGData.root", "BG");
  //RF.ReloadData("Data.root");
  //RF.ReloadSimulated("SigData.root", "Signal");
  //RF.ReloadSimulated("BGData.root", "BG");

  gBenchmark->Start("timer");
  //Or try an mcmc minimser 1000-># of points, 200->burnin 10 ~ 1/step size
  //auto mcmc=new BruMcmcCovariance(200,100,0.1,0.23,0.16,0.3);
  //mcmc->TurnOffCovariance();//BruMcmcCovariance only, do not proceed with covariance based sampling, just perform basic stepping
  //RF.SetMinimiser(mcmc);

   Here::Go(&RF);
  //Proof::Go(&RF,4); //run proof with 4 workers
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");

  gBenchmark->Start("timer2");
  new TCanvas;
  RF.DrawWeighted("Mmiss>>(100,0,10)","Signal");
  new TCanvas;
  RF.DrawWeighted("M1>>(100,0,10)","Signal");
  gBenchmark->Stop("timer2");
  gBenchmark->Print("timer2");
  //compare to true signal
  FiledTree::Read("MyModel","Data.root")->Tree()->Draw("M1","Sig==1","same");

}

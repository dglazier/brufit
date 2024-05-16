//brufit  FitHSMCModel.C
{

  sPlot RF;
  RF.SetUp().SetOutDir("out/");
  ///////////////////////////////Load Variables
  RF.SetUp().LoadVariable("Mmiss[0,9.5]");//should be same name as variable in tree
  RF.SetUp().LoadAuxVar("Eg[0,10]");
  // RF.SetUp().AddCut("Eg<3.2");

  RF.SetUp().SetIDBranchName("fgID");


  //////////////////////////////Make signal PDF
  RF.SetUp().FactoryPDF("RooHSEventsHistPDF::Signal(Mmiss,alpha[0,0,20],off[0,-2,2],scale[1,0.8,1.2])");
  RF.SetUp().LoadSpeciesPDF("Signal",1); 

  //////////////////////////////Make background PDF
  RF.SetUp().FactoryPDF("RooHSEventsHistPDF::BG(Mmiss,alphaB[0,0,5],offB[0,0,0],scaleB[1.0,0.8,1.2])");
  RF.SetUp().LoadSpeciesPDF("BG",1);

  
  ///////////////////////////Load Data
  RF.LoadData("MyModel","Data.root");
  RF.LoadSimulated("MyModel","SigData.root", "Signal");
  RF.LoadSimulated("MyModel","BGData.root", "BG");

  RF.SetUp().ErrorsWrong();//"naive" error calculation, much faster
  // RF.SetUp().AddFitOption(RooFit::NumCPU(4));

 //or try MCMC algorithm
  //auto mcmc=new BruMcmcCovariance(200,100,0.1,0.23,0.16,0.3);
  //mcmc->TurnOffCovariance();//BruMcmcCovariance only, do not proceed with covariance based sampling, just perform basic stepping
  //RF.SetMinimiser(mcmc);

  Here::Go(&RF);

  new TCanvas;
  RF.DrawWeighted("M1>>(100,0,10)","Signal");
  //compare to true signal
  FiledTree::Read("MyModel","Data.root")->Tree()->Draw("M1","Sig==1","same");

  
}

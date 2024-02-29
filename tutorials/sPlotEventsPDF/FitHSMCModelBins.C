//Run with
//brufit  FitHSMCModelBins.C
{
  //PROOF needs full paths!
  sPlot RF;
  RF.SetUp().SetOutDir("outBins/");
  ///////////////////////////////Load Variables
  RF.SetUp().LoadVariable("Mmiss[0,9.5]");//should be same name as variable in tree
  RF.SetUp().LoadAuxVar("Eg[0,10]");
  // RF.SetUp().AddCut("Eg<3.2");

  RF.SetUp().SetIDBranchName("fgID");


  //////////////////////////////Make signal PDF
  RF.SetUp().FactoryPDF("RooHSEventsHistPDF::Signal(Mmiss,smear_Sig[0,0,20],off_Sig[0,-2,2],scale_Sig[1,0.8,1.2])");
  RF.SetUp().LoadSpeciesPDF("Signal",1);

  //////////////////////////////Make background PDF
  RF.SetUp().FactoryPDF("RooHSEventsHistPDF::BG(Mmiss,smear_Bkg[0,0,5],off_Bkg[0,0,0],scale_Bkg[1.0,0.8,1.2])");
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
  //Run the fit here
  //Or try an mcmc minimser 1000-># of points, 200->burnin 10 ~ 1/step size
  //auto mcmc=new BruMcmcCovariance(200,100,0.1,0.23,0.16,0.3);
  //mcmc->TurnOffCovariance();//BruMcmcCovariance only, do not proceed with covariance-based sampling, just perform basic stepping
  //RF.SetMinimiser(mcmc);

  // Here::Go(&RF);
  Proof::Go(&RF,5); //run proof with 5 workers
  gBenchmark->Show("timer");

  gBenchmark->Start("timer2");
  auto filedTree = FiledTree::Read("MyModel","Data.root");
  auto trueTree = filedTree->Tree();
  TString MmissCut = "((0<=Mmiss)&&(Mmiss<=9.5))";  //Restrict distributions to the same Mmiss range that was used in the fit
  TCanvas* canv = new TCanvas();
  canv->Divide(2,2);
  //Draw signal distributions
  canv->cd(1);
  RF.DrawWeighted("M1>>hM1_Sig(100,0,10)","Signal",MmissCut);
  ((TH1*)gDirectory->Get("hM1_Sig"))->SetLineColor(kRed+1);
  trueTree->Draw("M1","(Sig==1)&&"+MmissCut,"same hist");
  canv->cd(2);
  RF.DrawWeighted("M2>>hM2_Sig(100,0,10)","Signal", MmissCut);
  ((TH1*)gDirectory->Get("hM2_Sig"))->SetLineColor(kRed+1);
  trueTree->Draw("M2","(Sig==1)&&"+MmissCut,"same hist");
  //Draw background distributions
  canv->cd(3);
  RF.DrawWeighted("M1>>hM1_Bkg(100,0,10)","BG",MmissCut);
  ((TH1*)gDirectory->Get("hM1_Bkg"))->SetLineColor(kRed+1);
  trueTree->Draw("M1","(Sig==-1)&&"+MmissCut,"same hist");
  canv->cd(4);
  RF.DrawWeighted("M2>>hM2_Bkg(100,0,10)","BG",MmissCut);
  ((TH1*)gDirectory->Get("hM2_Bkg"))->SetLineColor(kRed+1);
  trueTree->Draw("M2","(Sig==-1)&&"+MmissCut,"same hist");
  gBenchmark->Show("timer2");

  //Make sure weighted tree is written properly
  RF.DeleteWeightedTree();
}

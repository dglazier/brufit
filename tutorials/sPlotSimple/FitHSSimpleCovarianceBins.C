//Run with 
//brufit FitHSSimpleCovarianceBins.C
{
  sPlot RF;
  RF.SetUp().SetOutDir("outCov/");
  ///////////////////////////////Load Variables
  RF.SetUp().LoadVariable("Mmiss[0,10]");//should be same name as variable in tree
  // RF.SetUp().WS().var("Mmiss")->setBins(10);
 
  RF.SetUp().SetIDBranchName("fgID");

  /////////////////////////////Make Model Signal
  RF.SetUp().FactoryPDF("Gaussian::Signal( Mmiss, mean[10,4,7], sigma[0.1,0.0001,3] )");
  RF.SetUp().LoadSpeciesPDF("Signal",1);


  ////////////////////////////////Additional background
  RF.SetUp().FactoryPDF("Chebychev::BG(Mmiss,{a0[-0.1,-1,1],a1[0.1,-1,1],a2[-0.1,-1,1]})");
  //RF.SetUp().FactoryPDF("Chebychev::BG(Mmiss,{a0[-0.1,-1,1],a1[0.1,-1,1]})");
  RF.SetUp().LoadSpeciesPDF("BG",1);

  ///////////////////////////Load Data
  // RF.Data().BootStrap(4);//split the data in 4 seperate fits
  
  RF.Bins().LoadBinVar("Eg",1,3,4);
  RF.LoadData("MyModel","Data.root");

  /////////Set plot options for MCMC
  RF.SetPlotOptions("MCMC:CORNERFULL:CORNERZOOM:AUTOCORR");

  //Do we want to try many fits and use the best?
  //This will randomise the parameters for each fit
  //  RF.SetRefit(2);
  //Run the minuit fit here
  /**/
  gBenchmark->Start("minuit");
  Here::Go(&RF);
  gBenchmark->Stop("minuit");
   /**/
  
  //Choose Non Minuit mimimiser and run the fit
  /**/
  gBenchmark->Start("seq");
  auto mcmc=new RooMcmcSeqHelper(1000,500,1);
  mcmc->SetDesiredAcceptance(0.15,0.25,0.23);
  RF.SetMinimiser(mcmc);
  Here::Go(&RF);
  gBenchmark->Stop("seq");
  /**/
  
   //Choose a proposal using a MH covariance matrix 
   //BOTH mcmc to generate covMat and final chain called here
   /**/
  gBenchmark->Start("seqThencov");
  //(Nsteps (chain1), Nburn(chain1), Norm(chain1), Nsteps(chain2),Nburn(chain2),  Norm(chain2))
  auto mcmcCov=new RooMcmcSeqThenCov(1000,800,10,500,100,1);
  mcmcCov->SetDesiredAcceptance(0.15,0.3,0.23);
  // mcmc->SetUncorrelateYields(0);
  RF.SetMinimiser(mcmcCov);
  // 
  Here::Go(&RF);

  gBenchmark->Stop("seqThencov");
   /**/ 
  gBenchmark->Print("minuit");
  gBenchmark->Print("seq");
  gBenchmark->Print("seqThencov");

  new TCanvas;
  //RF.DrawWeighted("M1","Signal");
  RF.DrawWeighted("M1>>(100,0,10)","Signal");
  //compare to true signal
  FiledTree::Read("MyModel","Data.root")->Tree()->Draw("M1","Sig==1","same");

  //make sure weighted tree is written properly
  RF.DeleteWeightedTree();
}

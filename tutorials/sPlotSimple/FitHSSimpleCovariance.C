//Run with 
//brufit FitHSSimpleCovariance.C
{
  sPlot RF;
  RF.SetUp().SetOutDir("out/");
  ///////////////////////////////Load Variables
  RF.SetUp().LoadVariable("Mmiss[0,10]");//should be same name as variable in tree
  // RF.SetUp().WS().var("Mmiss")->setBins(10);
 
  RF.SetUp().SetIDBranchName("fgID");

  /////////////////////////////Make Model Signal
  RF.SetUp().FactoryPDF("Gaussian::Signal( Mmiss, SIMm[6,4,7], SIMw[0.2,0.0001,3] )");
  RF.SetUp().LoadSpeciesPDF("Signal",1);


  ////////////////////////////////Additional background
  RF.SetUp().FactoryPDF("Chebychev::BG(Mmiss,{a0[-0.1,-1,1],a1[0.1,-1,1]})");
  RF.SetUp().LoadSpeciesPDF("BG",1);

  ///////////////////////////Load Data
  RF.Data().BootStrap(4);//split the data in 4 seperate fits
  
  // RF.Bins().LoadBinVar("Eg",4,3,4);
  RF.LoadData("MyModel","Data.root");

  //Do we want to try many fits and use the best?
  //This will randomise the parameters for each fit
  //  RF.SetRefit(2);

  //Run the minuit fit here
  /**/
  Here::Go(&RF);
   /**/
  
  //Choose Non Minuit mimimiser and run the fit
  /**/
   RF.SetMinimiser(new RooMcmcSeq(1000,500,100));
   Here::Go(&RF);
    /**/

  //Choose a step proposal using minuit covariance matrix
  /**/
   RF.SetMinimiser(new RooMcmcMinuitCov(1000,500,100));
   Here::Go(&RF);
   /**/

  //Choose a proposal using MH covariance matrix
  /**/
   RF.SetMinimiser(new RooMcmcSeqCov(1000,500,20,500));
   //Argument #4 is the burn in for covariance matrix calc
   Here::Go(&RF);
   /**/
   
   //Choose a proposal using a MH covariance matrix 
   //BOTH mcmc to generate covMat and final chain called here
   /**/
   RF.SetMinimiser(new RooMcmcSeqThenCov(1000,2000,500,20,1));
   //(Nsteps (chain1), Nsteps(chain2), Nburn, Norm(chain1), Norm(chain2))
   Here::Go(&RF);
   /**/ 
   
  new TCanvas;
  //RF.DrawWeighted("M1","Signal");
  RF.DrawWeighted("M1>>(100,0,10)","Signal");
  //compare to true signal
  FiledTree::Read("MyModel","Data.root")->Tree()->Draw("M1","Sig==1","same");

  //make sure weighted tree is written properly
  RF.DeleteWeightedTree();
}

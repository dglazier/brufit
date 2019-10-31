//Run with 
//brufit FitObsBins.C
{
 
  Loader::Compile("PhiAsymmetry.cxx");
  
   
  FitManager RF;
  RF.SetUp().SetOutDir("outObsBins/");
  ///////////////////////////////Load Variables
  RF.SetUp().LoadVariable("Phi[-180,180]"); 
  RF.SetUp().LoadVariable("Pol[0,1]"); 
  RF.SetUp().LoadCategory("PolState[Polp=1,Polm=-1]"); 
  
  RF.SetUp().SetIDBranchName("fgID");

  ///////////////////////////////Make additional cut on an AuxVar
  //RF.SetUp().AddCut("AUX>2"); //Additional cut based on vars or aux vars
 
  /////////////////////////////Make Model Signal
  RF.SetUp().FactoryPDF("PhiAsymmetry::SigAsym( Phi,Pol,PolState,A[0,-1,1],B[0,-1,1] )");
  RF.SetUp().LoadSpeciesPDF("SigAsym",1);

  ////////////////////////////Make Bootstrap
  // RF.Data().BootStrap(400);
  ////////////////////////////Make Bins
  RF.Bins().LoadBinVar("Eg",4,3,4);
   
  ///////////////////////////Load Data
  //RF.Data().BootStrap(2);
  RF.LoadData("MyModel","Data.root");
  RF.LoadSimulated("MyModel","MC.root","SigAsym");
  
  //////////////////////////Load Weight
  RF.Data().LoadWeights("Signal","outsPlotBins/Tweights.root");

 
  //Or try an mcmc minimser 1000-># of points, 200->burnin 200 ~ 1/step size
  //RF.SetMinimiser(new RooMcmcSeq(1000,200,200));
  

  Here::Go(&RF);
  //OR run with PROOF-LITE on N=4 cores (you can change the 4)
  // Proof::Go(&RF,4);
  //OR run with FARM
  // Farm::Go(&RF,false);
}

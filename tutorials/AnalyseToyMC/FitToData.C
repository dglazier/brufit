{
  FitManager RF;
  RF.SetUp().SetOutDir("outFitToData/");
  ///////////////////////////////Load Variables
  RF.SetUp().LoadVariable("Phi[-180,180]"); 
  RF.SetUp().LoadVariable("Pol[0,1]"); 
  RF.SetUp().LoadCategory("PolState[Polp=1]");//,Pol0=0,Polm=-1]"); 
  RF.SetUp().SetIDBranchName("fgID");

  //Make Components PDF of 1 + A*Pol*PolState*cos(2Phi)
  //                         + B*Pol*PolState*sin(2Phi)
  RF.SetUp().LoadFormula("COS2=@Pol[]*@PolState[]*cos(2*@Phi[]*TMath::DegToRad())");
  RF.SetUp().LoadFormula("SIN2=@Pol[]*@PolState[]*sin(2*@Phi[]*TMath::DegToRad())");
  RF.SetUp().LoadParameter("A[1,-1,1]");
  RF.SetUp().LoadParameter("B[0.0,-1,1]");

  //////////////////////////////////////////////=> 1+                    A*cos(2Phi)+B*sin(2phi)
  RF.SetUp().FactoryPDF("RooComponentsPDF::SigAsym(1,{Phi,Pol,PolState},=A;COS2:B;SIN2)");
  RF.SetUp().LoadSpeciesPDF("SigAsym",1);

  ////////////////////////////Make Bins
  RF.Bins().LoadBinVar("Eg",10,3,4);

  ///////////////////////////Load Data
  RF.LoadData("MyModel","DataSignal.root");
  RF.LoadSimulated("MyModel","MC.root","SigAsym");

  /* //if want to use MCMC
    auto mcmc=new RooMcmcSeqThenCov(2000,1500,1,2000,1000,1);
    mcmc->SetDesiredAcceptance(0.15,0.5,0.23);
    mcmc->SetUncorrelateYields(0);
    RF.SetMinimiser(mcmc);
  */
  
  //Do the fit
  Here::Go(&RF);
  //  Proof::Go(&RF,4); //4 workers, make sure you have enough CPUs
  RF.WriteThis();
}

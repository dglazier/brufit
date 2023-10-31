{
  //Start from large file of "flat/phase space" events
  TString sigfile  = "flat/Toy0.root";
  
  TString outdir = "toys_test/";
  /****************************************/
  /************Create ToyManager***********/    
  /****************************************/
  ToyManager toyman(1);
  toyman.SetUp().SetOutDir(outdir);
  toyman.SetUp().SetIDBranchName("UID");
  /****************************************/
  /***********Apply Dilepton model**********/    
  /****************************************/
 
  Model(toyman,10000); //create toy data with 100 events

  /**************************************************/
  /**************Load PhaseSpace data****************/
  /**************************************************/
  //We will sample from these events with accept/reject
  //LoadSimulated(Treename,Filename,PDFname)
  toyman.LoadSimulated("ToyData",sigfile,"Dilepton");

  /**************************************************/
  /************Set Toy Parameter Values**************/
  /**************************************************/
  toyman.SetUp().SetParVal("BH",0.6); //Bethe Heitler contribution
  toyman.SetUp().SetParVal("TCS",0.1); //TCS contribution
  toyman.SetUp().SetParVal("ImM",0.7); //Imaginery part of M (ReM^2 =1-ImM^2)
 
  /**************************************************/
  /**************Create Toy  data****************/
  /**************************************************/
  Here::Go(&toyman);


  //do a sanity check fit
  auto toyfitter= toyman.Fitter();
 //********************************************
 // Perform "fit" with an MCMC sampler
  //a tree MCMCTree is included in the Results*.root file
  auto mcmc=new BruMcmcCovariance(20000,1000,1,0.23,0.16,0.3);
  ////mcmc->TurnOffCovariance();//BruMcmcCovariance only, do not proceed with covariance based sampling, just perform basic stepping
  //toyfitter->SetMinimiser(mcmc);
  
    /**************************************************/
  /************Set Toy Parameter Values**************/
  /**************************************************/
  toyfitter->SetUp().SetParVal("BH",0.6); //Bethe Heitler contribution
  toyfitter->SetUp().SetParVal("TCS",0.1); //TCS contribution
  toyfitter->SetUp().SetParVal("ImM",0.8); //Imaginery part of M (ReM^2 =1-ImM^2)

  Here::Go(toyfitter.get());
 
}

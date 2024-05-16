{
  //Start from large file of "flat/phase space" events
  TString sigfile  = "eepData.root";
  
  TString outdir = "sim_toys2/";
  /****************************************/
  /************Create ToyManager***********/    
  /****************************************/
  ToyManager toyman(1);
  toyman.SetUp().SetOutDir(outdir);
  toyman.SetUp().SetIDBranchName("UID");
  toyman.SetTruthPrefix("tru");
  /****************************************/
  /***********Apply Dilepton model**********/    
  /****************************************/
 
  Model(toyman,20000); //create toy data with 50000 events

  /**************************************************/
  /**************Load PhaseSpace data****************/
  /**************************************************/
  //We will sample from these events with accept/reject
  //LoadSimulated(Treename,Filename,PDFname)
  toyman.LoadSimulated("fittree",sigfile,"Dilepton");

  /**************************************************/
  /************Set Toy Parameter Values**************/
  /**************************************************/
  toyman.SetUp().SetParVal("BH",0.6); //Bethe Heitler contribution
  toyman.SetUp().SetParVal("TCS",0.2); //TCS contribution
  toyman.SetUp().SetParVal("ImM",0.8); //Imaginery part of M (ReM^2 =1-ImM^2)
 
  /**************************************************/
  /**************Create Toy  data****************/
  /**************************************************/
  Here::Go(&toyman);


  //do a sanity check fit
  auto toyfitter= toyman.Fitter();
 //********************************************
 // Perform "fit" with an MCMC sampler
  //a tree MCMCTree is included in the Results*.root file
  auto mcmc=new BruMcmcCovariance(5000,1000,1,0.23,0.16,0.3);
  ////mcmc->TurnOffCovariance();//BruMcmcCovariance only, do not proceed with covariance based sampling, just perform basic stepping
  // toyfitter->SetMinimiser(mcmc);
  
   /**************************************************/
   /************Set Toy Parameter Values**************/
   /**************************************************/
  toyfitter->SetUp().SetParVal("BH",0.6); //Bethe Heitler contribution
  toyfitter->SetUp().SetParVal("TCS",0.2); //TCS contribution
  toyfitter->SetUp().SetParVal("ImM",0.5); //Imaginery part of M (ReM^2 =1-ImM^2)

  Here::Go(toyfitter.get());
 
}

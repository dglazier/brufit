{

  TString datafile  = "DATAFILE";
  TString sigfile  = "SIMFILE";
  TString treename  = "TREENAME";
  
  TString outdir = "OUTPUT_DIR/";
  
  /****************************************/
  /************Create FitManager***********/    
  /****************************************/
  FitManager fm;
  fm.SetUp().SetOutDir(outdir);
  fm.SetUp().SetIDBranchName("UID");
  
  /****************************************/
  /***********Apply Dilepton model**********/    
  /****************************************/
  Model(fm); 

  /**************************************************/  
  /********************Make bins*********************/ 
  /**************************************************/
  
  // Double_t tbinLimits[] = {0,0.5,1,1.5,2}; // tbins
  // fm.Bins().LoadBinVar("t",tbinLimits);			  

  /**************************************************/
  /****************Load data and MC******************/
  /**************************************************/

  fm.LoadData(treename,datafile);
  fm.LoadSimulated(treename,sigfile,"Dilepton"); //"Dilepton" is given in Model.C as name of the RooComponentsPDF
  

  /**************************************************/
  /***********Choose minimiser and run***************/ 
  /**************************************************/
 //********************************************
 // Perform "fit" with an MCMC sampler
  //a tree MCMCTree is included in the Results*.root file
  //auto mcmc=new BruMcmcCovariance(5000,1000,1,0.23,0.16,0.3);
  ////mcmc->TurnOffCovariance();//BruMcmcCovariance only, do not proceed with covariance based sampling, just perform basic stepping
  // toyfitter->SetMinimiser(mcmc);
  
   /**************************************************/
   /************Set Toy Parameter Values**************/
   /**************************************************/
  toyfitter->SetUp().SetParVal("BH",0.); //Bethe Heitler contribution
  toyfitter->SetUp().SetParVal("TCS",0.); //TCS contribution
  toyfitter->SetUp().SetParVal("ImM",0.); //Imaginery part of M (ReM^2 =1-ImM^2)

  Here::Go(&fm);
 
  
}

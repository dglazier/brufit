{

  
  /****************************************/
  /***************Filenames****************/    
  /****************************************/
  TString datafile  = "sim_toys2/Toy0.root";
  TString sigfile  = "eepData.root";
  
  TString outdir = "fit_out/";
  
  /****************************************/
  /************Create FitManager***********/    
  /****************************************/
  FitManager fm;
  fm.SetUp().SetOutDir(outdir);
  // fm.SetUp().SetIDBranchName("UID");
  
  /****************************************/
  /***********Apply Dilepton model**********/    
  /****************************************/
  Model(fm); 

  /**************************************************/  
  /********************Make bins*********************/ 
  /**************************************************/
  
  Double_t tbinLimits[] = {0,0.5,1,1.5,2}; // tbins
  fm.Bins().LoadBinVar("t",4,tbinLimits);			  

  /**************************************************/
  /****************Load data and MC******************/
  /**************************************************/

  fm.LoadData("fittree",datafile);
   //"Dilepton" is given in Model.C as name of the RooComponentsPDF
  fm.LoadSimulated("fittree",sigfile,"Dilepton");

  /**************************************************/
  /***********Choose minimiser and run***************/ 
  /**************************************************/
  //number of CPUs to split likelihood calc.
  fm.SetUp().AddFitOption(RooFit::NumCPU(4)); 
  Here::Go(&fm);
 
  
}

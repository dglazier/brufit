void Model(FitManager& fm, Int_t Nevents=1){


  /****************************************/
  /*************Load Variables*************/    
  /****************************************/
  //Watch theta limits!!! cannot go to 0 as discontinuity
  fm.SetUp().LoadVariable(Form("Theta[%lf,%lf]",TMath::Pi()/10,TMath::Pi()*9/10));
  //fm.SetUp().LoadVariable(Form("Theta[%lf,%lf]",0.,TMath::Pi()) );
  fm.SetUp().LoadVariable(Form("Phi[-%lf,%lf]",TMath::Pi(),TMath::Pi()));
  fm.SetUp().LoadVariable(Form("Pol[-1,1]"));

  /****************************************/
  /*************Load Formula*************/    
  /****************************************/
  //TCS only term
  fm.SetUp().LoadFormula("TCS_TH=1+TMath::Cos(@Theta[])*TMath::Cos(@Theta[])");

  //Interference terms
  fm.SetUp().LoadFormula("INT_TH=(1+TMath::Cos(@Theta[])*TMath::Cos(@Theta[]))/(TMath::Sin(@Theta[]))");
  fm.SetUp().LoadFormula("INT_COSPHI=(TMath::Cos(@Phi[]))");
  fm.SetUp().LoadFormula("INT_hSINPHI=(@Pol[]*TMath::Sin(@Phi[]))");

  //Bethe Heitler term
  fm.SetUp().LoadFormula("BH_TH=(1+TMath::Cos(@Theta[])*TMath::Cos(@Theta[]))/(TMath::Sin(@Theta[])*TMath::Sin(@Theta[]))");


  /****************************************/
  /*************Load Parameters************/    
  /****************************************/
  fm.SetUp().LoadParameter("BH[0.5,0.4,1]");
  fm.SetUp().LoadParameter("TCS[0,0,1]");
  //fm.SetUp().LoadParameter("INT[0,0,1]");
  //fm.SetUp().LoadParameter("ReM[0,-1,1]");
  fm.SetUp().LoadParameter("ImM[0,0,1]");

   /****************************************/
  /*************Load Constraints***********/    
  /****************************************/
 //Constrain Total contribution ==1
  fm.SetUp().LoadFormula("INT=(1 - @BH[] - @TCS[])");
  fm.SetUp().LoadFormula("ReM=TMath::Sqrt(1 - @ImM[]*@ImM[])");

  // fm.SetUp().SetIDBranchName("UID");

  /****************************************/
  /*************Make model PDF*************/ //DONE
  /****************************************/
  fm.SetUp().FactoryPDF("RooComponentsPDF::Dilepton(0,{Theta,Phi,Pol},=BH;BH_TH:TCS;TCS_TH:INT;ReM;INT_TH;INT_COSPHI:INT;ImM;INT_TH;INT_hSINPHI)"); 
  fm.SetUp().LoadSpeciesPDF("Dilepton",Nevents);
}

void Model(){
}

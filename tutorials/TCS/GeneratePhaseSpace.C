 {
   //////////////////////////////////////////////////////////////////
   //Create "Simulated data" for generator and normalisation integral
   ToyManager Flat{1};//# number of files to produce
   Flat.SetUp().SetOutDir("flat/");//output directory

   //polarised fit variables
   // Flat.SetUp().LoadVariable(Form("Theta[0,%lf]",TMath::Pi()));
   Flat.SetUp().LoadVariable(Form("Theta[%lf,%lf]",TMath::Pi()/10,TMath::Pi()*9/10));
   // Flat.SetUp().LoadVariable(Form("Theta[%lf,%lf]",TMath::Pi()/4,3*TMath::Pi()/4));
   
   Flat.SetUp().LoadVariable(Form("Phi[-%lf,%lf]",TMath::Pi(),TMath::Pi()));
   Flat.SetUp().LoadVariable(Form("Pol[-1,1]"));
   Flat.SetUp().LoadVariable(Form("t[0,1]"));

   //Flat but give some t-dependence
   Flat.SetUp().LoadFormula("TDEP=TMath::Log(1./@t[])");
   //And phase space theta depedence (==Flat cosTheta)
   Flat.SetUp().LoadFormula("THETADEP=TMath::Sin(@Theta[])");
 
   //Define PDF as constant RooComponentsPDF
   Flat.SetUp().LoadFunctionVar("RooConstVar::Val(1)");
   
   Flat.SetUp().FactoryPDF("RooComponentsPDF::Flat(0,{Theta,Phi,Pol,t},=Val;TDEP;THETADEP)");
   //Load generator to produce 1E6 events
   Flat.SetUp().LoadSpeciesPDF("Flat",1E7); //# events

   // Create a sample of data
   Here::Go(&Flat);
 }

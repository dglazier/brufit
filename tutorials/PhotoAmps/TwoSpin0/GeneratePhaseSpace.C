 {
   //////////////////////////////////////////////////////////////////
   //Create "Simulated data" for generator and normalisation integral
   //Here it is just flat distributions in all variables
   //no acceptance effects are included
   ToyManager Flat{1};//# number of files to produce
   Flat.SetUp().SetOutDir("flat/");//output directory

   //polarised fit variables
   Flat.SetUp().LoadVariable("CosTh[-1,1]"); 
   Flat.SetUp().LoadVariable("Phi[-3.14159,3.14159]");
   Flat.SetUp().LoadVariable("PolPhi[-3.14159,3.14159]");
   Flat.SetUp().LoadVariable("Pol[0.5,1]");


   //Define PDF as constant RooComponentsPDF
   Flat.SetUp().LoadFunctionVar("RooConstVar::Val(1)");
   Flat.SetUp().FactoryPDF("RooComponentsPDF::Flat(1,{CosTh,Phi,PolPhi,Pol},=Val)");
   //Load generator to produce 1E6 events
   Flat.SetUp().LoadSpeciesPDF("Flat",1E6); //# events

   // Create a sample of data
   Here::Go(&Flat);
 }

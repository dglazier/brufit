{
  //////////////////////////////////////////////
  //Create flat data for MC integration
  ToyManager Flat{1};
  Flat.SetUp().SetOutDir("flat/");

  //polarised fit variables
  Flat.SetUp().LoadVariable("Phi[-180,180]");
  Flat.SetUp().LoadVariable("Pol[0.5,1]"); 
  Flat.SetUp().LoadCategory("PolState[Polp=1,Polm=-1]"); 

  
  Flat.SetUp().LoadFunctionVar("RooConstVar::Val(1)");
  Flat.SetUp().FactoryPDF("RooComponentsPDF::Flat(1,{Phi,Pol,PolState},=Val)");

  Flat.SetUp().LoadSpeciesPDF("Flat",1000000); //2000 events

  // Create a sample of data
  Here::Go(&Flat);

  ///////////////////////////////////////////
  //Create Toy data to fit with my PDF
  Loader::Compile("PhiAsymmetry.cxx");
  
  ToyManager Generator{1};
  Generator.SetUp().SetOutDir("gen/");
  Generator.SetUp().LoadVariable("Phi[-3.14159,3.14159]");
  Generator.SetUp().LoadVariable("Pol[0.5,1]"); 
  Generator.SetUp().LoadCategory("PolState[Polp=1,Polm=-1]"); 

  //Give A and B initial values for generating data with 
  Generator.SetUp().FactoryPDF("PhiAsymmetry::SigAsym( Phi,Pol,PolState,A[0.5,-1,1],B[0.1,-1,1] )");
  Generator.SetUp().LoadSpeciesPDF("SigAsym",1000);

  
  Generator.LoadSimulated("ToyData","flat/Toy0.root","SigAsym");

  //Create a sample of data
  Here::Go(&Generator);

  //Now try  
  //   root gen/Toy0.root 
  //   root [4] ToyData->Draw("Phi","PolState==1")
}

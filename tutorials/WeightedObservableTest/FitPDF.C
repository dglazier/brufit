 {
   ///////////////////////////////////////////
   //Fit Toy data with my PDF
  Loader::Compile("PhiAsymmetry.cxx");
  
  FitManager Fitter;
  Fitter.SetUp().SetOutDir("fit/");
  Fitter.SetUp().LoadVariable("Phi[-180,180]");
  Fitter.SetUp().LoadVariable("Pol[0.5,1]"); 
  Fitter.SetUp().LoadCategory("PolState[Polp=1,Polm=-1]"); 

  Fitter.SetUp().FactoryPDF("PhiAsymmetry::SigAsym( Phi,Pol,PolState,A[0.,-1,1],B[0.,-1,1] )");
  Fitter.SetUp().LoadSpeciesPDF("SigAsym",1);

  
  Fitter.LoadData("ToyData","gen/Toy0.root");
  Fitter.LoadSimulated("ToyData","flat/Toy0.root","SigAsym");

  //Fit toy data
  Here::Go(&Fitter);
}

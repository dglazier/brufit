{
  //////////////////////////////////////////////////////////////////
  //Create "Simulated data for generator and normalisation integral
  ToyManager Flat{1};
  Flat.SetUp().SetOutDir("flatSphHar/");

  //polarised fit variables
  Flat.SetUp().LoadVariable("CosTh[-1,1]"); 
  Flat.SetUp().LoadVariable("Phi[-3.14159,3.14159]");
 
  
  Flat.SetUp().LoadFunctionVar("RooConstVar::Val(1)");
  Flat.SetUp().FactoryPDF("RooComponentsPDF::Flat(1,{CosTh,Phi},=Val)");

  Flat.SetUp().LoadSpeciesPDF("Flat",1000000); //2000 events

  // Create a sample of data
  Here::Go(&Flat);
  
  cout<<" GENERATED pseudo simulated events "<<endl<<endl;
  //////////////////////////////////////////
  //Creaet pseuso data for fitting
  ToyManager Generator{1};
  Generator.SetUp().SetOutDir("genSphHar/");
  Generator.SetUp().LoadVariable("CosTh[-1,1]");
  Generator.SetUp().LoadVariable("Phi[-3.14159,3.14159]");

  //use parser to expand spherical harmonics with  0<L<4 0<M<L
  //and create all constants
  //SphHarmonicMoments(TString name,TString cth,TString phi,Int_t Lmax,Int_t Mmin,Int_t Mmax)
  ComponentsPdfParser  parser=HS::FIT::SphHarmonicMoments("SphHarmonicMoments","CosTh","Phi",4,0,4);

  //Now expansion in spherical harmonic moments
  //H_Alpha_L_M are fit parameters => the moments
  //K_L = =TMath::Sqrt(2*L.+1.)/TMath::Sqrt(4*TMath::Pi())
  //Y_L_M are spherical harmonic functions depending on CosTh and Phi
  string sum = "H_0_0_0[1]"; //constant == 1;
  sum +=       "+ SUM(L[1|4],M[0|4<L+1]){H_0_L_M[0,-1,1]*K_L*ReY_L_M(CosTh,Phi,Y_L_M)}";
  
  Generator.SetUp().ParserPDF(sum,parser);
  Generator.SetUp().LoadSpeciesPDF("SphHarmonicMoments",10000); //2000 events

  //Give moments some initial values
  Generator.SetUp().WS().var("H_0_4_4")->setVal(0.5);
  Generator.SetUp().WS().var("H_0_2_1")->setVal(0.2);
  
  //"Simulated" data to project spherical harmonic distributions onto
  Generator.LoadSimulated("ToyData","flatSphHar/Toy0.root","SphHarmonicMoments");

  //Create a sample of data
  gBenchmark->Start("gen");
  Here::Go(&Generator);
  gBenchmark->Stop("gen");
  gBenchmark->Print("gen");
  cout<<" GENERATED pseudo data events "<<endl<<endl;


  //now create a fitter from the toymanager
  auto Fitter=Generator.Fitter();
  
  //And fit the sample data
  Here::Go(Fitter);

  cout<<" FITTED pseudo data events "<<endl<<endl;


}

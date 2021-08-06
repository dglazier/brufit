{
  //////////////////////////////////////////////////////////////////
  //Create "Simulated data for generator and normalisation integral
  ToyManager Flat{1};
  Flat.SetUp().SetOutDir("flatDWigner/");

  //polarised fit variables
  Flat.SetUp().LoadVariable("CosTh[-1,1]"); 
  Flat.SetUp().LoadVariable("Phi[-3.14159,3.14159]");
 
  
  Flat.SetUp().LoadFunctionVar("RooConstVar::Val(1)");
  Flat.SetUp().FactoryPDF("RooComponentsPDF::Flat(1,{CosTh,Phi},=Val)");

  Flat.SetUp().LoadSpeciesPDF("Flat",1000000); //2000 events

  // Create a sample of data
  Here::Go(&Flat);
  
  //////////////////////////////////////////
  //Create pseuso data for fitting
  ToyManager Generator{1};
  Generator.SetUp().SetOutDir("genDWigner/");
  Generator.SetUp().LoadVariable("CosTh[-1,1]");
  Generator.SetUp().LoadVariable("Phi[-3.14159,3.14159]");

  Generator.SetUp().LoadFormula("Theta=TMath::ACos(@CosTh[])");

//use parser to expand Wigner D functions with  0<L<4 0<M<L 0<N<L (for >=0 factorials)
  //WignerDFunctionMoments(TString name,TString theta,TString cth,TString phi,Int_t Lmax,Int_t Mmin,Int_t Mmax,Int_t Nmin,Int_t Nmax)
  //Note if cth="" then Theta assumed to be actual data not formula,
  //if cth has a string then this is assumed to be the tree data and Theta caclulated via formula
  ComponentsPdfParser  parser=HS::FIT::WignerDFunctionMoments("Moments","Theta","CosTh","Phi",4,0,4,0,4);

  //Now expansion in Wigner D function  moments
  //G_L_M_N are fit parameters => the moments
  //D_L_M_N are the Wigner D functions, depending on Theta and Phi
  //(note formula calculation of Theta from CosTh)
  //K_L =TMath::Sqrt(2*L.+1.)/TMath::Sqrt(4*TMath::Pi())

   string sum = "G_0_0_0[1]"; //constant == 1;
  sum +=       "+ SUM(L[1|4],M[0|1<L+1],N[0|1<L+1]){K_L*G_L_M_N[0,-1,1]*ReD_L_M_N(Theta,Phi,D_L_M_N)}";
  
  Generator.SetUp().ParserPDF(sum,parser);
  Generator.SetUp().LoadSpeciesPDF("Moments",40000); //2000 events

  //Give some dummy values need to be quite small <1/(sqrt(4pi))
  Generator.SetUp().WS().var("G_4_1_0")->setVal(0.5);
   
  //"Simulated" data to project spherical harmonic distributions onto
  Generator.LoadSimulated("ToyData","flatDWigner/Toy0.root","Moments");

  //Create a sample of data
  gBenchmark->Start("gen");
  Here::Go(&Generator);
  gBenchmark->Stop("gen");
  gBenchmark->Print("gen");

  
  //now create a fitter from the toymanager
  auto Fitter=Generator.Fitter();
 //Fitter->SetMinimiser(new RooMcmcSeq(5000,500,50));
 
  //And fit the sample data
  gBenchmark->Start("fit");
  Here::Go(Fitter);
  gBenchmark->Stop("fit");
  gBenchmark->Print("fit");
  

}

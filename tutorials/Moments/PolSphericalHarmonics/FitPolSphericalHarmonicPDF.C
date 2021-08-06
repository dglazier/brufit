{
  FitManager Fitter{};
  //output directory for results, split files etc.
  Fitter.SetUp().SetOutDir("fitSphHar/");
  //variables corresponding to branch names in input Tree (should all be double)
  Fitter.SetUp().LoadVariable("CosTh[-1,1]");
  Fitter.SetUp().LoadVariable("Phi[-3.14159,3.14159]");
  Fitter.SetUp().LoadVariable("PolPhi[-3.14159,3.14159]");
  Fitter.SetUp().LoadVariable("Pol[0.7,0.8]");

  //use parser to expand spherical harmonics with  0<L<4 0<M<L
  //and create all constants
  //PolarisedSphHarmonicMoments(TString name,TString cth,TString phi,TString phiPol,TString Pol,Int_t Lmax,Int_t Mmin,Int_t Mmax)
  ComponentsPdfParser  parser=HS::FIT::PolarisedSphHarmonicMoments("SphHarmonicMoments","CosTh","Phi","PolPhi","Pol",4,0,4);

  //Now expansion in spherical harmonic moments
  //H_Alpha_L_M are fit parameters => the moments
  //K_L = =TMath::Sqrt(2*L.+1.)/TMath::Sqrt(4*TMath::Pi())
  //Y_L_M are spherical harmonic functions depending on CosTh and Phi
  //Note the factor tau from JPAC paper is cared for in
  //PolarisedSphHarmonicMoments function

  string sum = "H_0_0_0[1]"; //constant == 1;
  sum +=       "+ SUM(L[1|4],M[0|4<L+1]){H_0_L_M[0,-1,1]*K_L*ReY_L_M(CosTh,Phi,Y_L_M)}";
  sum +=       "+ SUM(L[0|4],M[0|4<L+1]){H_1_L_M[0,-1,1]*K_L*ReY_L_M(CosTh,Phi,Y_L_M)*COS2PHI}";
  sum+=        "+ SUM(L[1|4],M[1|4<L+1]){H_2_L_M[0,-1,1]*K_L*ImY_L_M(CosTh,Phi,Y_L_M)*SIN2PHI}";

  
  Fitter.SetUp().ParserPDF(sum,parser);
  Fitter.SetUp().LoadSpeciesPDF("SphHarmonicMoments",100000); //100000 events

  //Give moments some initial values
  // Fitter.SetUp().WS().var("H_0_4_4")->setVal(0.5);
  // Fitter.SetUp().WS().var("H_0_2_1")->setVal(0.2);


  Fitter.LoadData("ToyData","genSphHar/Toy0.root");
  Fitter.LoadSimulated("ToyData","flatSphHar/Toy0.root","SphHarmonicMoments");

  //And fit the sample data
  // Fitter->SetUp().AddFitOption(RooFit::NumCPU(4));
  gBenchmark->Start("fit");
  //try MCMC algorithms
  //  Fitter->SetMinimiser(new RooMcmcSeq(1000,500,100));
  Here::Go(&Fitter);
  //Proof::Go(&Fitter,4);
  gBenchmark->Stop("fit");
  gBenchmark->Print("fit");

}

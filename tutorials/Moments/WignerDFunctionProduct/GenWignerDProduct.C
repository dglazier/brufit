{
  //////////////////////////////////////////////////////////////////
  //Create "Simulated data for generator and normalisation integral
  ToyManager Flat{1};
  Flat.SetUp().SetOutDir("flatDWignerProduct/");

  //polarised fit variables
  Flat.SetUp().LoadVariable("CosThA1[-1,1]"); 
  Flat.SetUp().LoadVariable("PhiA1[-3.14159,3.14159]");
  Flat.SetUp().LoadVariable("CosThB1[-1,1]"); 
  Flat.SetUp().LoadVariable("PhiB1[-3.14159,3.14159]");
 
  
  Flat.SetUp().LoadFunctionVar("RooConstVar::Val(1)");
  Flat.SetUp().FactoryPDF("RooComponentsPDF::Flat(1,{CosThA1,PhiA1,CosThB1,PhiB1},=Val)");

  Flat.SetUp().LoadSpeciesPDF("Flat",1000000); //2000 events

  // Create a sample of data
  Here::Go(&Flat);
  
  //////////////////////////////////////////
  //Create pseuso data for fitting
  ToyManager Generator{1};
  Generator.SetUp().SetOutDir("genDWignerProduct/");
  Generator.SetUp().LoadVariable("CosThA1[-1,1]");
  Generator.SetUp().LoadVariable("PhiA1[-3.14159,3.14159]");
  Generator.SetUp().LoadVariable("CosThB1[-1,1]");
  Generator.SetUp().LoadVariable("PhiB1[-3.14159,3.14159]");

  Generator.SetUp().LoadFormula("ThetaA1=TMath::ACos(@CosThA1[])");
  Generator.SetUp().LoadFormula("ThetaB1=TMath::ACos(@CosThB1[])");
 
  //use parser to expand Wigner D functions with  0<L<4 0<M<4 0<N<4 0<L2<4(for >=0 factorials)
  //ComponentsPdfParser WignerDFunctionProductMoments(TString name,TString theta1,TString cth1,TString phi1,TString theta2,TString cth2,TString phi2,Int_t Lmax1,Int_t Mmin1,Int_t Mmax1,Int_t Nmin1,Int_t Nmax1,Int_t Lmax2)
   //Note if cth1="" then Theta assumed to be actual data not formula,
  //if cth1 has a string then this is assumed to be the tree data and Theta1 calculated via formula

  ComponentsPdfParser  parser=HS::FIT::WignerDFunctionProductMoments("Moments","ThetaA1","CosThA1","PhiA1","ThetaB1","CosThB1","PhiB1",2,0,2,0,2,2);
 
  //Now expansion in Wigner D function product (D.D)  moments
  //H_l_m_L_M are fit parameters => the moments
  //Dab_l_m_L_M are the products of Wigner D functions D^L_Mm.D^l_m0
  //See Chung Spin Formalisms eqn 7.32
  //depending on ThetaA1,ThetaB1 and PhiA1,PhiB1
  //(note formula calculation of Theta from CosTh)
  //K_L =TMath::Sqrt(2*L.+1.)/TMath::Sqrt(4*TMath::Pi())
  Generator.SetUp().LoadConstant("H_0_0_0_0[1]");
  string sum;//="K_0*K_0*h_0_0_0_0[1]";
  sum +=     "SUM(L[0|2],l[0|2],M[0|2<L+1],m[0|2<L+1<l+1]){K_l*K_L*H_l_m_L_M[0,-1,1]*ReDab_l_m_L_M}";

  Generator.SetUp().ParserPDF(sum,parser);
  //Generator.SetUp().SetParVal("H_0_0_0_0",1,kTRUE);//and make const
  Generator.SetUp().LoadSpeciesPDF("Moments",40000); //2000 events

  //Set some moment values to generate toy data with
  Generator.SetUp().WS().var("H_2_2_2_2")->setVal(0.6);
  Generator.SetUp().WS().var("H_1_0_1_0")->setVal(0.2);
  Generator.SetUp().WS().var("H_2_0_0_0")->setVal(-0.4);
  
  //"Simulated" data to project spherical harmonic distributions onto
  Generator.LoadSimulated("ToyData","flatDWignerProduct/Toy0.root","Moments");

  //Create a sample of data
  gBenchmark->Start("gen");
  Here::Go(&Generator);
  gBenchmark->Stop("gen");
  gBenchmark->Print("gen");
  
  
  //now create a fitter from the toymanager
  auto Fitter=Generator.Fitter();
  Fitter->SetMinimiser(new RooMcmcSeq(1000,500,50));
  //Fitter->SetUp().SetParVal("H_0_0_0_0",1,kTRUE);//and make const

  //And fit the sample data
  gBenchmark->Start("fit");
  Here::Go(Fitter);
  gBenchmark->Stop("fit");
  gBenchmark->Print("fit");
  

}

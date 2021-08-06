{
  //////////////////////////////////////////////////////////////////
  //Create "Simulated data for generator and normalisation integral
  ToyManager Flat{1};
  Flat.SetUp().SetOutDir("flatDWignerProductSym/");

  //polarised fit variables
  Flat.SetUp().LoadVariable("CosThA1[-1,1]"); 
  Flat.SetUp().LoadVariable("PhiA1[-3.14159,3.14159]");
  Flat.SetUp().LoadVariable("CosThB1[-1,1]"); 
  Flat.SetUp().LoadVariable("PhiB1[-3.14159,3.14159]");
  Flat.SetUp().LoadVariable("CosThA2[-1,1]"); 
  Flat.SetUp().LoadVariable("PhiA2[-3.14159,3.14159]");
  Flat.SetUp().LoadVariable("CosThB2[-1,1]"); 
  Flat.SetUp().LoadVariable("PhiB2[-3.14159,3.14159]");
 
  
  Flat.SetUp().LoadFunctionVar("RooConstVar::Val(1)");
  Flat.SetUp().FactoryPDF("RooComponentsPDF::Flat(1,{CosThA1,PhiA1,CosThB1,PhiB1,CosThA2,PhiA2,CosThB2,PhiB2},=Val)");

  Flat.SetUp().LoadSpeciesPDF("Flat",1000000); //2000 events

  // Create a sample of data
  Here::Go(&Flat);
  
  //////////////////////////////////////////
  //Create pseuso data for fitting
  ToyManager Generator{1};
  Generator.SetUp().SetOutDir("genDWignerProductSym/");
  Generator.SetUp().LoadVariable("CosThA1[-1,1]");
  Generator.SetUp().LoadVariable("PhiA1[-3.14159,3.14159]");
  Generator.SetUp().LoadVariable("CosThB1[-1,1]");
  Generator.SetUp().LoadVariable("PhiB1[-3.14159,3.14159]");
  Generator.SetUp().LoadVariable("CosThA2[-1,1]");
  Generator.SetUp().LoadVariable("PhiA2[-3.14159,3.14159]");
  Generator.SetUp().LoadVariable("CosThB2[-1,1]");
  Generator.SetUp().LoadVariable("PhiB2[-3.14159,3.14159]");

  Generator.SetUp().LoadFormula("ThetaA1=TMath::ACos(@CosThA1[])");
  Generator.SetUp().LoadFormula("ThetaB1=TMath::ACos(@CosThB1[])");
  Generator.SetUp().LoadFormula("ThetaA2=TMath::ACos(@CosThA2[])");
  Generator.SetUp().LoadFormula("ThetaB2=TMath::ACos(@CosThB2[])");
 
  //use parser to expand Wigner D functions with  0<L<4 0<M<4 0<N<4 0<L2<4(for >=0 factorials)
  //ComponentsPdfParser WignerDFunctionProductMoments(TString name,TString theta1,TString cth1,TString phi1,TString theta2,TString cth2,TString phi2,Int_t Lmax1,Int_t Mmin1,Int_t Mmax1,Int_t Nmin1,Int_t Nmax1,Int_t Lmax2)
   //Note if cth1="" then Theta assumed to be actual data not formula,
  //if cth1 has a string then this is assumed to be the tree data and Theta1 calculated via formula

  //1,2 => symmeterised angles
  ComponentsPdfParser  parser=HS::FIT::WignerDFunctionProductMomentsSym("Moments","ThetaA1","CosThA1","PhiA1","ThetaB1","CosThB1","PhiB1","ThetaA2","CosThA2","PhiA2","ThetaB2","CosThB2","PhiB2",2,0,2,0,2,2);
 
 
  //Now expansion in Wigner D function product (D.D)  moments
  //H_l_m_L_M are fit parameters => the moments
  //Dab_l_m_L_M are the products of Wigner D functions D^L_Mm.D^l_m0
  //See Chung Spin Formalisms eqn 7.32
  //depending on ThetaA1,ThetaB1 and PhiA1,PhiB1
  //(note formula calculation of Theta from CosTh)
  //K_L =TMath::Sqrt(2*L.+1.)/TMath::Sqrt(4*TMath::Pi())
  Generator.SetUp().LoadConstant("H_0_0_0_0[2]");//Normalisation constant
  string sum;
  sum +="SUM(L[0|2],l[0|2],M[0|2<L+1],m[0|2<L+1<l+1]){K_l*K_L*H_l_m_L_M[0,-1,1]*ReDab_l_m_L_M}";
  sum +="+ SUM(L[0|2],l[0|2],M[0|2<L+1],m[0|2<L+1<l+1]){K_l*K_L*H_l_m_L_M*ResymDab_l_m_L_M}";

  Generator.SetUp().ParserPDF(sum,parser);
 
  Generator.SetUp().LoadSpeciesPDF("Moments",40000); //2000 events

  //Set some moment values to generate toy data with
  Generator.SetUp().WS().var("H_2_2_2_2")->setVal(0.3);
  Generator.SetUp().WS().var("H_1_0_1_0")->setVal(0.2);
  Generator.SetUp().WS().var("H_2_0_0_0")->setVal(-0.2);
  cout<<" value "<<Generator.SetUp().WS().function("ResymDab_0_0_0_0")->getVal()<<endl;exit(0);
  //"Simulated" data to project spherical harmonic distributions onto
  Generator.LoadSimulated("ToyData","flatDWignerProductSym/Toy0.root","Moments");

  //Create a sample of data
  gBenchmark->Start("gen");
  Here::Go(&Generator);
  gBenchmark->Stop("gen");
  gBenchmark->Print("gen");
  
  
  //now create a fitter from the toymanager
  auto Fitter=Generator.Fitter();
  //  Fitter->SetMinimiser(new RooMcmcSeq(5000,500,100));
 
  //And fit the sample data
  gBenchmark->Start("fit");
  //Proof::Go(Fitter,1);
  Here::Go(Fitter);
  gBenchmark->Stop("fit");
  gBenchmark->Print("fit");
  

}

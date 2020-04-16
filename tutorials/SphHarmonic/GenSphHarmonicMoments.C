//run with brufit GenSphHarmonicMoments.C
{


  //Create flat data for MC integration
  ToyManager Flat{1};
  Flat.SetUp().SetOutDir("flatSphHarmonic/");

  //polarised fit variables
  Flat.SetUp().LoadVariable("CosTh[0,-1,1]");
  Flat.SetUp().LoadVariable("Phi[-3.14159,3.14159]");

  
  Flat.SetUp().LoadFunctionVar("RooConstVar::Val(1)");
  Flat.SetUp().FactoryPDF("RooComponentsPDF::Flat(1,{CosTh,Phi},=Val)");

  Flat.SetUp().LoadSpeciesPDF("Flat",1000000); //2000 events

  // Create a sample of data
  gBenchmark->Start("flat");
  Here::Go(&Flat);
  gBenchmark->Stop("flat");
  gBenchmark->Print("flat");
  
  //////////////////////////////////////////////////////////////////////
  //////Create pseudo data
  ////// You can try different moments below
  ToyManager Generator{1};
  Generator.SetUp().SetOutDir("genSphHarmonic/");

  Generator.SetUp().LoadVariable("CosTh[0,-1,1]");
  Generator.SetUp().LoadVariable("Phi[-3.14159,3.14159]");

  //Configure the polarised moments, the last 2 numbers give the Lmax and Mmax
  //Moments is the name the PDF will take
  //All Spherical harmonics and Cos2PolPhi functions are created here
  auto configGenPDF=HS::FIT::EXPAND::ComponentsRealSphHarmonic(Generator.SetUp(),"Moments","CosTh","Phi",3,2);
  cout<<endl<<endl<<endl<<configGenPDF<<endl;

  //Give some dummy values need to be quite small <1/(sqrt(4pi))
  Generator.SetUp().WS().var("H0_2_1")->setVal(0.6);
  Generator.SetUp().WS().var("H0_2_2")->setVal(-0.5);
  Generator.SetUp().WS().var("H0_1_1")->setVal(-0.8);
  Generator.SetUp().WS().var("H0_2_0")->setVal(-0.1);

  //create and load PDF into the intensity function
  Generator.SetUp().FactoryPDF(configGenPDF);
  Generator.SetUp().LoadSpeciesPDF("Moments",10000); //10000 events

  //Create a sample of data
  gBenchmark->Start("gen");
  Here::Go(&Generator);
  gBenchmark->Stop("gen");
  gBenchmark->Print("gen");


}

#include "ConfigureAmpsNoValues.C"

void RunGivenBruMoments(){
  
  gRandom->SetSeed(0);
  
  //Define amplitudes (in additional file)
  auto& setup = ConfigureAmpsNoValues(2,1,1);//(Lmax,MMax,Nref)
  setup.SetParVal("aphi_0_0",0,kTRUE); //fix S real

  //Fix with moments from ConfigureAmpsSPJune.C
  MomentHelper moments;
  moments.Set("../TwoSpin0/fitBruMoments/ResultsHSAmpMinuit2.root");
 
  
  //setup the solver, arguments :
  //  setup = BruFit setup
  //  resolution = smear moments by resolution
  //  ignore_observables = do not include the following polarised moments
  //                       i.e. alpha = 0,1,2, or 3 => H_0,H_1,H_2,H_3
  m2pw::EquationSolver solver{setup,0.0,{"H_3"}}; //ignore H_3
  //  m2pw::EquationSolver solver{setup,0.05,{"H_3","H_2","H_1"}};
  solver.SetEquationValues(moments);
  solver.Print("v");

  //create output tree
  solver.MakeResultTree("resultsGivenBruMomentsNoH3.root");

  gBenchmark->Start("solver");

  //loop and perform 10,000 minimisations with random starting amplitudes
  for(int i = 0; i<10000;i++){
    if(i%100==0) cout<<i<<" "<<endl;
    solver.GetPars().Randomise();
    solver.Solve();
    solver.FillTree();
  }

  gBenchmark->Stop("solver");
  gBenchmark->Print("solver");

  //Save results tree
  solver.GetPars().CloseTree();

  solver.PrintResult();
 
}

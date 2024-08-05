#include "ConfigureAmpsNoValues.C"

void RunGivenMoments(){
  
  gRandom->SetSeed(0);
  
  //Define amplitudes (in additional file)
  auto& setup = ConfigureAmpsNoValues();
  setup.Formulas().Print();

  //Fix with moments from ConfigureAmpsSPJune.C
  MomentHelper moments;
  moments.Set("H_0_0_0", 2.00000);
  moments.Set("H_1_0_0", 0.321818);
  moments.Set("H_0_1_0", 0.0802487);
  moments.Set("H_1_1_0", -0.111632);
  moments.Set("H_0_1_1", -0.145798);
  moments.Set("H_1_1_1", -0.0248914);
  moments.Set("H_2_1_1", 0.456002);
  moments.Set("H_0_2_0", -0.252528);
  moments.Set("H_1_2_0", -0.085465);
  moments.Set("H_0_2_1", -0.0521929);
  moments.Set("H_1_2_1", -0.0028709);
  moments.Set("H_2_2_1", 0.0336173);
  moments.Set("H_0_2_2", -0.201987);
  moments.Set("H_1_2_2", 0.0289487);
  moments.Set("H_2_2_2", -0.177133);
  
  
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
  solver.MakeResultTree("resultsGivenMomentsNoH3.root");

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

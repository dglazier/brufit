#include "ConfigureAmpsSPJune.C"

void RunMoments(){
  
  gRandom->SetSeed(0);
  
  //Define amplitudes (in additional file)
  auto& setup = ConfigureAmpsSPJune();
  setup.Formulas().Print();

  
  //setup the solver, arguments :
  //  setup = BruFit setup
  //  resolution = smear moments by resolution
  //  ignore_observables = do not include the following polarised moments
  //                       i.e. alpha = 0,1,2, or 3 => H_0,H_1,H_2,H_3
  m2pw::EquationSolver solver{setup,0.0,{"H_3"}}; //ignore H_3
  //  m2pw::EquationSolver solver{setup,0.05,{"H_3","H_2","H_1"}};
  solver.Print("v");

  //create output tree
  solver.MakeResultTree("resultsMomentsNoH3.root");

  gBenchmark->Start("solver");

  //loop and perform 50,000 minimisations with random starting amplitudes
  for(int i = 0; i<50000;i++){
    if(i%100==0) cout<<i<<" "<<endl;
    solver.GetPars().Randomise();
    solver.Solve();
    solver.FillTree();
  }

  gBenchmark->Stop("solver");
  gBenchmark->Print("solver");

  //Save results tree
  solver.GetPars().CloseTree();
 
}

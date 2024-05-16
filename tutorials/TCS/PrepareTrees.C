{

  //we need out input trees to have theta, phi, polarisation, t, Egamma branches
  //use this script to adapt any trees from data or simulation
  ROOT::RDataFrame df("particle", "/hdd/jlab/chanser_pass2_eep/eep_eep_none_all_pd__/particleData/ParticleVariables_0.root");
  auto df2 = df.Filter("Truth==1&&eepQ2<1&&eepIMep>1.5&&EBCuts==1&&eepEgamma>7&&trueepPol>0&&trueepPol<1&&eept<2").Alias("Phi", "eepMesonPhi").Alias("t", "eept").Define("Theta", "TMath::ACos(eepMesonCosTh)").Define("Hel","-1+2*((Int_t)gRandom->Integer(2))").Define("Pol", "sqrt(1-eepPol*eepPol)*0.8*(Hel)");
  auto dftru = df2.Alias("truPhi", "trueepMesonPhi").Alias("trut", "trueept").Define("truTheta", "TMath::ACos(trueepMesonCosTh)").Alias("truHel", "Hel").Define("truPol", "sqrt(1-trueepPol*trueepPol)*0.8*Hel");
  

  //add UID branch
  int counter=0;
  auto uid = [&counter]() { return counter++; };
  auto  dfuid = dftru.Define("UID", uid, {});

  // write out new dataset. this triggers the event loop
  dfuid.Snapshot("fittree", "eepData.root",{"Theta","Phi","Pol","t","eepEgamma","truTheta","truPhi","truPol","trut","trueepEgamma","UID"});


}

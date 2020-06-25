////Usage: root 'macros/AddIDBranch.C("TreeName","FileName.root")'

void AddIDBranch(TString treename,TString filename, Long64_t initVal=0){

  cout<<"AddIDBranch : this is going to add a branch to "<<treename<<" in "<<filename<<endl;
  cout<<"type q now to quit, any other key to continue."<<endl;
  char choice;
  cin>>choice;
  if(choice=='q') return;

  //safe to continue
  auto treefile=TFile::Open(filename,"update");
  auto tree =static_cast<TTree* >(treefile->Get(treename));

  if(tree->GetBranch("UID")){
    cout<<"Branch UID already exists, exiting "<<endl;
    return; //aready exists
  }

  Double_t UID{static_cast<Double_t>(initVal)};//make it a double so can read in RooFit
  TBranch* branch=tree->Branch("UID",&UID,"UID/D");

  Long64_t id=(Long64_t)UID;
  
  for(Long64_t i=0;i<tree->GetEntries();i++){
    branch->Fill();
    id++;
    UID=id;   
  }

  tree->Write();
  
  delete treefile;
 
}

//Run with 
//root CreateData.C()

void CreateData(){

	TFile* fileData=new TFile("data.root","recreate");

	Double_t Eg, fgID;
	TTree* treeData=new TTree("data","data");
	treeData->Branch("Eg",&Eg,"Eg/D");
	treeData->Branch("fgID",&fgID,"fgID/D");

	Int_t Nev=1E5;
	
	for(Int_t i=0;i<Nev;i++){
		fgID=i;
		Eg=gRandom->Gaus(9,10);
		treeData->Fill();
	}

	treeData->Write();
	fileData->Close();

	
	TFile* fileSim=new TFile("sim.root","recreate");

	TTree* treeSim=new TTree("sim","sim");
	treeSim->Branch("Eg",&Eg,"Eg/D");
	treeSim->Branch("fgID",&fgID,"fgID/D");

	for(Int_t i=0;i<Nev;i++){
		fgID=i;
		Eg=gRandom->Gaus(9,7);
		treeSim->Fill();
	}

	treeSim->Write();
	fileSim->Close();
}

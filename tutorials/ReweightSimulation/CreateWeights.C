TH1* weightHist = new TH1F("impWeights","impWeights",50,5,10);

void CreateWeights() {
	
	TFile *fileData=new TFile("data.root");
	TTree *treeData= (TTree*) fileData->Get("data");
	
	TFile *fileMC=new TFile("sim.root");
	TTree *treeMC= (TTree*) fileMC->Get("sim");
	
	
	Weights* impWeights = new Weights("MCWeights");
	impWeights->SetSpecies("LikeData");
	impWeights->SetIDName("fgID");
	impWeights->SetFile("Weights.root");
	
	impWeights->ImportanceSampling(treeMC, treeData, weightHist, "Eg");
	impWeights->PrintWeight();
	
	impWeights->Draw1DWithWeights(treeMC,weightHist,"Eg","LikeData");
	
	weightHist->SetLineColor(2);
	weightHist->Draw("hist");
	treeData->Draw("Eg","","same");
	impWeights->Save();
}

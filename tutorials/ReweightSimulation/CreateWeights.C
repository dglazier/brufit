TH1* weightHist = new TH1F("impWeights","impWeights",50,5,10);

void CreateWeights() {
	//Here we want the Eg distribution in MC to match that from data
	TFile *fileData=new TFile("data.root");
	TTree *treeData= (TTree*) fileData->Get("data");
	
	TFile *fileMC=new TFile("sim.root");
	TTree *treeMC= (TTree*) fileMC->Get("sim");
	
	
	Weights* impWeights = new Weights("MCWeights");
	impWeights->SetSpecies("LikeData");
	impWeights->SetIDName("fgID");
	impWeights->SetFile("Weights.root");
	
	TH1* weightHist = new TH1F("impWeights","impWeights",100,-50,50);
	impWeights->ImportanceSampling(treeMC, treeData, weightHist, "Eg");
	impWeights->PrintWeight();
	
	TH1* weightedHist = new TH1F("WeightedMC","WeightedMC",50,-10,30);
	
	//Draw MC distribution with weights
	impWeights->Draw1DWithWeights(treeMC,weightedHist,"Eg","LikeData");
	
	weightedHist->SetLineColor(2);
	weightedHist->SetMinimum(0);
	weightedHist->Draw("hist");

	//Draw target data distribution
	treeData->Draw("Eg","","same");

	//Draw unweighted MC distribution
	treeMC->Draw("Eg","","same");

	impWeights->Save();
}

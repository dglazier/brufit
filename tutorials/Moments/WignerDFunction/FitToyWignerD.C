//Just an example how to run weighted simulation for mc integration
//when using PdfParser.
//must first run
//brufit GenWignerD.C
//then
//brufit FitToyWignerD.C
//See WEIGHTS@ addition to parser string
//Weights matching simulation to data are created in CreateWeights
//We weight CosTh which means the simulation will match the data
//with parameters==0, so we should get all parameters conssitent with 0
//rather than their true value.


void CreateWeights();

void FitToyWignerD(){
  CreateWeights();

   
  FitManager Fitter;
  Fitter.SetUp().SetOutDir("weightFitDWigner/");
  Fitter.SetUp().LoadVariable("CosTh[-1,1]");
  Fitter.SetUp().LoadVariable("Phi[-3.14159,3.14159]");

  Fitter.SetUp().LoadFormula("Theta=TMath::ACos(@CosTh[])");

//use parser to expand Wigner D functions with  0<L<4 0<M<L 0<N<L (for >=0 factorials)
  //WignerDFunctionMoments(TString name,TString theta,TString cth,TString phi,Int_t Lmax,Int_t Mmin,Int_t Mmax,Int_t Nmin,Int_t Nmax)
  //Note if cth="" then Theta assumed to be actual data not formula,
  //if cth has a string then this is assumed to be the tree data and Theta caclulated via formula
  ComponentsPdfParser  parser=HS::FIT::WignerDFunctionMoments("Moments","Theta","CosTh","Phi",4,0,4,0,4);

  //Now expansion in Wigner D function  moments
  //G_L_M_N are fit parameters => the moments
  //D_L_M_N are the Wigner D functions, depending on Theta and Phi
  //(note formula calculation of Theta from CosTh)
  //K_L =TMath::Sqrt(2*L.+1.)/TMath::Sqrt(4*TMath::Pi())

   string sum = "G_0_0_0[1]"; //constant == 1;
  sum +=       "+ SUM(L[1|4],M[0|1<L+1],N[0|1<L+1]){K_L*G_L_M_N[0,-1,1]*ReD_L_M_N(Theta,Phi,D_L_M_N)}";

  // sum+="WEIGHTS@WeightSpeciesName,WeightFileName,WeightObjectName";
  sum+="WEIGHTS@LikeData,flatDWigner/Weights.root,MCWeights";
  
  Fitter.SetUp().ParserPDF(sum,parser);
  Fitter.SetUp().LoadSpeciesPDF("Moments",40000); //2000 events
  
  Fitter.LoadData("ToyData","genDWigner/Toy0.root");
  Fitter.LoadSimulated("ToyData","flatDWigner/Toy0.root","Moments");


  
  Here::Go(&Fitter); //fit toy

  

}
#include "../macros/AddIDBranch.C"
void CreateWeights() {
  //Make sure trees have ID branch "UID"
  AddIDBranch("ToyData","genDWigner/Toy0.root");
  AddIDBranch("ToyData","flatDWigner/Toy0.root");

  //Here we want the CosTh  distribution in MC to match that from data
  unique_ptr<TFile>  fileData{new TFile("genDWigner/Toy0.root")};
  TTree *treeData= fileData->Get<TTree>("ToyData");
  
  unique_ptr<TFile>  fileMC{new TFile("flatDWigner/Toy0.root")};
  TTree *treeMC= fileMC->Get<TTree>("ToyData");
	

  //configure a weights object name MCWeights
  Weights impWeights("MCWeights");
  impWeights.SetSpecies("LikeData");
  impWeights.SetIDName("UID");
  impWeights.SetFile("flatDWigner/Weights.root");

  //create a histogram as basis for reweighting
  TH1F weightHist("impWeights","impWeights",100,-1,1);
  impWeights.ImportanceSampling(treeMC, treeData, &weightHist, "CosTh");
  impWeights.PrintWeight();
  
  TH1* weightedHist = new TH1F("WeightedMC","WeightedMC",50,-1,1);
	
  //Draw MC distribution with weights
  impWeights.Draw1DWithWeights(treeMC,weightedHist,"CosTh","LikeData");
	
  weightedHist->SetLineColor(2);
  weightedHist->SetMinimum(0);
  weightedHist->DrawCopy("hist");

  //Draw target data distribution
  treeData->Draw("CosTh","","same");

  //Draw unweighted MC distribution
  treeMC->Draw("CosTh","","same");
  
  impWeights.Save();
  
}

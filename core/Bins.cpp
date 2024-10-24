//#include "Riostream.h"
#include "Bins.h"
#include "TROOT.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TFileMerger.h"
#include "TObjectTable.h"
#include "TRandom3.h"
#include "TSystem.h"
#include <algorithm>
#include <utility>
//#include "ProcInfo_t.h"

namespace HS{
  namespace FIT{

    Bins::Bins(const TString& name) :TNamed(name,name){
    }

    Bins::Bins(const TString& name,const TString& filename):TNamed(name,name){
      TDirectory *saveDir=gDirectory;
      fFile=new TFile(filename);
      if(!fFile->IsOpen()) {Error("Bins::Bins(TString name,TString filename)"," File does not exist %s",filename.Data());return;}
      auto* filebins=dynamic_cast<Bins*>(fFile->Get(name));
      if(!filebins) {Error("Bins::Bins(TString name,TString filename)","did not find bins %s in file %s",name.Data(),filename.Data());return;}
      fBinNames=filebins->GetBinNames();
      fNbins=filebins->GetN();
      fNaxis=filebins->GetNAxis();
      fVarAxis=filebins->GetVarAxis();
      fBinnedTreeName = filebins->GetBinnedTreeName();
      fFileNames = filebins->GetFileNames();
  
      //tree is not written to file as data member
      saveDir->cd();
    }
    Bins::Bins(const Bins& other, const char* name): TNamed(name,name){
      fBinNames=other.fBinNames;
      fNbins=other.fNbins;
      fNaxis=other.fNaxis;
      fVarAxis=other.fVarAxis;
      fFile=nullptr;
      fBinnedTreeName = other.fBinnedTreeName;
      fMAXFILES=other.fMAXFILES;
      fMaxEntries=other.fMaxEntries;
    }
    Bins::~Bins(){
      if(fFile){fFile->Close(); delete fFile;}
    }
    void Bins::AddAxis(TString name,Int_t nbins,Double_t min,Double_t max){
      //Add a new axis for a given variable, name should be tree name
      //Want to make an array for use with FindBin(). must call the other contructor
      TArrayD bins(nbins+1);
      Double_t binwidth = (max - min) / Double_t(nbins);
      for(Int_t i=0;i<nbins+1;i++)
	bins[i]= min + i * binwidth;
      AddAxis(std::move(name),nbins,bins.GetArray());
  
    }
    void Bins::AddAxis(const TString& name,Int_t nbins,Double_t* xbins){
      //Add a new axis for a given variable, name should be tree name
      TAxis axis(nbins,xbins);
      axis.SetName(name);
      fVarAxis.push_back(axis);
      fNaxis++;
    }
    void Bins::InitialiseBins(){
      if(fNaxis==0) return;
      //Make bin names for every individual bin
      TString binName;
      fNbins=0;
      fBinNames.clear();
      IterateAxis(0,binName);

    }
    void Bins::IterateAxis(Int_t iA,const TString& binName) {
      //iterate through all bins possible with given axis and constuct
      //unique names for them
      if (iA >= fNaxis){ //stop clause
	fBinNames.push_back(binName);
	fNbins++;
	//Info("Bins::IterateAxis"," %s",binName.Data());
	return;
      }
      VecString_t part;
      for (int iB = 1; iB <= fVarAxis[iA].GetNbins(); iB++) { 
	fVarAxis[iA].SetBinLabel(iB,Form("%1.6f_",fVarAxis[iA].GetBinCenter(iB)));
	part.push_back(TString(fVarAxis[iA].GetName())+fVarAxis[iA].GetBinLabel(iB));
	IterateAxis(iA+1,binName+fVarAxis[iA].GetName()+fVarAxis[iA].GetBinLabel(iB));
      }
      if(std::find(fPartName.begin(),fPartName.end(),part)==fPartName.end()) {
	fPartName.insert(fPartName.begin(),part);//for correct ordering with fVar.Axis vector
      }
    }
    void Bins::RunBinTree(TTree* tree,TString selection){
      fSelection=std::move(selection);
      fFileNames.clear();

      if(fNbins==0) InitialiseBins();//1 time initialisation
  
      if(fNbins<fMAXFILES){
	RunBinTree(tree,0,fNbins);
	return;
      }
      auto Nlots=(Int_t)(fNbins/fMAXFILES);
      Int_t Nrem=fNbins%fMAXFILES;
      for(Int_t i=0;i<Nlots;i++)
	RunBinTree(tree,fMAXFILES*i,fMAXFILES*(i+1));
      //and remainder
      RunBinTree(tree,fMAXFILES*(Nlots),fMAXFILES*(Nlots)+Nrem);
  
    }

    void Bins::RunBinTree(TTree* tree,Int_t BMin,Int_t BMax){
      //Create all sub trees
      std::cout<<"Bins::RunBinTree Running bins from "<<BMin<<" to "<<BMax<<std::endl;
      //create entry lists for tree
      TDirectory *saveDir=gDirectory;
      //  fFile->cd();

      Bool_t GotAnInt=kFALSE;
      TVectorD vVal(fNaxis);//values of variables for given entry
      vector<Int_t> vValI(fNaxis);//int values of variables for given entry
      vector<Int_t> vIntIndex;
      for(Int_t j=0;j<fNaxis;j++){
	tree->SetBranchStatus(fVarAxis[j].GetName(),"1");//STATUS must be called before ADDRESS!! see Important remarkse in TChain SetBranchStatus!
	tree->GetLeaf(fVarAxis[j].GetName())->GetTypeName();//this isneeded to check type properly when setting branch address
	tree->GetBranch(fVarAxis[j].GetName())->SetAutoDelete();
	if(tree->SetBranchAddress(fVarAxis[j].GetName(),&vVal[j])==-2){
	  //In case we have Int_t branches
	  if(tree->SetBranchAddress(fVarAxis[j].GetName(),&vValI[j])==0){
	    vIntIndex.push_back(j);
	    GotAnInt=kTRUE;
	  }
	}
      }

      saveDir->cd();

      //prepare the binned trees
      Int_t Nhere=BMax-BMin;
      fTrees.reserve(Nhere);
      fTrees.resize(Nhere);
      //make output directory if not existing
      gSystem->MakeDirectory(fOutDir+"/");
      fBinnedTreeName=tree->GetName();
      for(Int_t ib=BMin;ib<BMax;ib++){
	gSystem->MakeDirectory(fOutDir+"/"+GetBinName(ib));
	fFileNames.push_back(fOutDir+"/"+GetBinName(ib)+"/Tree"+fDataName+".root");
	fTrees[ib-BMin]=new BinTree(Nhere,fOutDir+"/"+GetBinName(ib)+"/Tree"+fDataName,tree,fOmitBranches);	
      }

      saveDir->cd();

      Int_t totalBytes=0;

 
      //Turn on any branches needed to evaluate selection
      auto leaves=tree->GetListOfLeaves();
      auto branches=tree->GetListOfBranches();
      vector<TString> on_branches;
      for(Int_t ib=0;ib<branches->GetEntries();ib++){
	if(tree->GetBranchStatus(branches->At(ib)->GetName())){
	  on_branches.emplace_back(branches->At(ib)->GetName());
	}
      }
      //turn on all branches so can make formula
      tree->SetBranchStatus("*",true);
  
      //create selection cut
      if(fSelection==TString()) fSelection="1";
      TTreeFormula treeCut("selection",fSelection,tree);
      vector<TString> cut_branches;
      for(Int_t jl=0;jl<treeCut.GetNcodes();jl++){
	if(treeCut.GetLeaf(jl)->GetBranch()){
	  cut_branches.emplace_back(treeCut.GetLeaf(jl)->GetBranch()->GetName());
	}
      }
      //Now only turn on required branches
      tree->SetBranchStatus("*",false);
      for(const auto& brname: on_branches)
	tree->SetBranchStatus(brname,true);
      for(const auto& brname: cut_branches)
	tree->SetBranchStatus(brname,true);
  
      for(Long64_t i=0;i<tree->GetEntries();i++){//loop over events
	tree->GetEntry(i);
	if(!static_cast<Bool_t>(treeCut.EvalInstance()))
	  continue;
	if(GotAnInt){//put the integer value in the double array
	  for(int iv : vIntIndex){
	    vVal[iv]=vValI[iv];
	  }
	}
	if(i%100000==0){
	  std::cout<<"On event "<<i<<" = "<<100.*i/tree->GetEntries()<<"%"<<std::endl;
	}
	fBin=FindBin(vVal);
	//check if bin in current range
	if(fBin>=BMax||fBin<BMin) continue;
	Int_t aBin=fBin-BMin;
	//Fill the tree associated with this bin
	Int_t evSize=fTrees[aBin]->ReadEvent();
	totalBytes+=evSize;
	if(fTrees[aBin]->GetEntries()==(Long64_t)fMaxEntries/fNbins/evSize) {
	  fTrees[aBin]->Reset();
	}
      }
      for(const auto& brname: cut_branches)
	if(std::find(on_branches.begin(),on_branches.end(),brname)==on_branches.end()){
	  tree->SetBranchStatus(brname,false);
	}
  
      tree->ResetBranchAddresses();
      saveDir->cd();
      //cleanup
      for(Int_t ib=BMin;ib<BMax;ib++){
	Int_t aBin=ib-BMin;
	delete  fTrees[aBin];
	fTrees[aBin]=nullptr;
      }
      fTrees.clear();
    }
    void Bins::Save(const TString& filename){
      Info(" Bins::Save()"," Saving %s to %s",GetName(),filename.Data());
      fFile=new TFile(filename,"recreate");
      Write();
      if(fFile){
	delete fFile;
	fFile=nullptr;}
    }
    void Bins::PrintAxis(){

      for(Int_t iA=0;iA<fNaxis;iA++)
	std::cout<<fVarAxis[iA].GetName()<<" "<<fVarAxis[iA].GetNbins()<<" "<<fVarAxis[iA].GetXmin()<<" "<<fVarAxis[iA].GetXmax()<<" "<<std::endl;
    }

    void Bins::MakeDirectories(){
      if(fNbins==0)
	InitialiseBins();
      if(fNbins==0)
	return;
      cout<<"Make dirs "<< fNbins<<" "<<fOutDir<<" "<<GetBinName(0)<<endl;
      gSystem->MakeDirectory(fOutDir+"/");
      for(Int_t ib=0;ib<fNbins;ib++){
	gSystem->MakeDirectory(fOutDir+"/"+GetBinName(ib));
      }
    }

    Int_t Bins::FindBin(Double_t v0){
      TVectorD vals(1);vals[0]=v0;
      return FindBin(vals);
    }
    Int_t Bins::FindBin(Double_t v0,Double_t v1){
      TVectorD vals(2);vals[0]=v0;vals[1]=v1;
      return FindBin(vals);
    }
    Int_t Bins::FindBin(Double_t v0,Double_t v1,Double_t v2){
      TVectorD vals(3);vals[0]=v0;vals[1]=v1;vals[2]=v2;
      return FindBin(vals);
    }
    Int_t Bins::FindBin(Double_t v0,Double_t v1,Double_t v2,Double_t v3,Double_t v4,Double_t v5){
      TVectorD vals(6);vals[0]=v0;vals[1]=v1;vals[2]=v2;vals[3]=v3;vals[4]=v4;vals[5]=v5;
      return FindBin(vals);
    }
    Int_t Bins::FindBin(TVectorD vals){
      //Loop over each axis and find bin for each
      Bool_t InLimits=kTRUE;
      for(Int_t iA=0;iA<fNaxis;iA++)//first check var is within variable ranges
	if(vals[iA]<fVarAxis[iA].GetXmin()||vals[iA]>fVarAxis[iA].GetXmax()) InLimits=kFALSE;
      if(!InLimits) {return -1;}
      //now find bin for each axis
      vector<Int_t> vBin(fNaxis); //store for the bin number of each axis
      for(Int_t iA=0;iA<fNaxis;iA++){//loop over vars/axis
	vBin[iA]=1+TMath::BinarySearch(fVarAxis[iA].GetNbins(),fVarAxis[iA].GetXbins()->GetArray(),vals[iA]);
      }
      //now have the bin for each axis, find single bin
      Int_t theBin=-1;
      for(Int_t iA1=0;iA1<fNaxis-1;iA1++){
	Int_t tbin=vBin[iA1]-1;//-1 as bin indexing starts at 1 with 0 underflow
	for(Int_t iA2=iA1+1;iA2<fNaxis;iA2++)
	  tbin*=fVarAxis[iA2].GetNbins();
	theBin+=tbin;
      }
      theBin+=vBin[fNaxis-1];
      return theBin;
    }

    ////////////////////////////////////////////////////////////////
    ///BinTree utility class
    ///Duplicates a tree but keeps its branches/memory etc seperate
    ///This allows us to make many copies without memory issues
    ///CloneTree give trouble with memory, when lots of copies
    BinTree::BinTree(Int_t nbins,const TString& name,TTree* tree0,vector<TString> omit){
      std::cout<<"Constructing Bin Tree "<<name<<std::endl;
      fName=name;
      fFile=TFile::Open(fName+".root","recreate");
      vector<TString> turnOn;
      for( auto& bname : omit ){//turn off omitted branches
	if(tree0->GetBranchStatus(bname)) turnOn.push_back(bname);
	tree0->SetBranchStatus(bname,false);
      }
      fTree=tree0->CloneTree(0);	
      for(auto&  bname : turnOn ){ //turn back on omitted branches now tree is made
	tree0->SetBranchStatus(bname,true);
      }
  
      fTree->SetName(tree0->GetName());
      fTree->SetDirectory(fFile);
      fTree->SetAutoSave(1E12); //We do our won autosave as this one changes basket size greatly increasing memory when large number of bins
      fTree->SetBasketSize("*",64000); //cloned trees have the parent basket size which can be very large and use large amount of memeory when we great many bins
      fTree->SetAutoFlush(1E12); //Don't let root flush or it will make basket sizes
    }
    BinTree::~BinTree(){
      if(fTree&&fFile)
	Save();
 
    }
    void BinTree::Save(){
      std::cout<<"BinTree::Save() "<<fName<<std::endl;
      if(!fTree) return;
      fFile->cd();
      /// fTree->FlushBaskets();
      fTree->Write();
      fTree->SetDirectory(nullptr);
      fTree->ResetBranchAddresses();
      delete fTree;
      // fFile->Close();
      delete fFile;fFile=nullptr;
      fTree=nullptr;
    }
    void BinTree::Reset(){
      std::cout<<"reset "<<fName<<std::endl;
      fTree->AutoSave("FlushBaskets");
      fTree->SetBasketSize("*",16000); //just in case...
      return;
 
    }
  }
}

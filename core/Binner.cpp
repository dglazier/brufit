#include "Binner.h"
#include "FiledTree.h"
#include <algorithm> 
#include <utility>

namespace HS{
  namespace FIT{
    
    Binner::Binner(Setup &setup):fBins("HSBins") {
      LoadSetup(setup);
    }
    
    void Binner::LoadSetup(Setup &setup) {
      fIsSetup=kTRUE;
      
      fBins.SetName("HSBins");
      
      fOutDir=setup.GetOutDir();
      
      auto vars=setup.DataVars();//get dataset variables
      auto it=vars.iterator(); //iterate over them

      TObject * var=nullptr;
      while ((var=it.Next())){
	fVarNames.push_back(var->GetName()); //store name of dataset varible for filtering tree branches
      }
 
    }

    void  Binner::LoadBinVar(TString opt,Int_t nbins,Double_t min,Double_t max){
      fBins.AddAxis(std::move(opt),nbins,min,max);     
    }
    void  Binner::LoadBinVar(const TString& opt,Int_t nbins,Double_t* xbins){
      fBins.AddAxis(opt,nbins,xbins);     
    }

    void Binner::ReloadData(const TString& fname,const TString& name){
      Bins bins("HSBins",Form("%s/%sBinsConfig.root",fOutDir.Data(),name.Data()));
      //fBins=bins;
      fNameToFiles[name]=bins.GetFileNames();
      fNameToTree[name]=bins.GetBinnedTreeName();
      fBinNames=bins.GetBinNames();
      
      }
    void Binner::SplitData(const TString& tname,TString fname,const TString& name){
      if(fBins.GetNAxis()==0){ //no splits required

	//PROOF needs full paths so enforce it here
	if(!fname.BeginsWith("/"))
	  fname = TString(gSystem->Getenv("PWD"))+"/"+fname;

	fNameToFiles[name]={{fname}};
	fNameToTree[name]=tname;
	fBinNames={{""}};

	return;
      }
      if(!fIsSetup) {
	cout<<"Binner::SplitData ERROR not setup yet!"<<endl;
	exit(0);
      }
      auto filetree=FiledTree::Read(tname,fname);
      SplitData(filetree->Tree().get(), name);
    }

    void Binner::SplitData(TTree* tree,const TString& name){
      //turn off all branches we do not require to save space and memory
      tree->SetBranchStatus("*",false);
      
      for(auto &vname : fVarNames) {//only copy variable branches for speed
	tree->SetBranchStatus(vname,true);
      }
      for(auto &vname : fKeepBranches) {//plus any others 
	if(tree->GetBranch(vname))tree->SetBranchStatus(vname,1);
      }
       cout<<"SplitData "<<fOutDir<<" "<<fBins.GetN()<<endl;
      //now split the tree into bins and save in subdirs of fOutDir
      fBins.SetOutDir(fOutDir);
      fBins.SetDataName(name);
      fBins.RunBinTree(tree,fSelection);
      fBins.Save(fOutDir+name+"BinsConfig.root");
      fNameToFiles[name]=fBins.GetFileNames();
      fNameToTree[name]=tree->GetName();
      fBinNames=fBins.GetBinNames();
  }
    const TString Binner::BinName(UInt_t i)  {

      InitBins();
      
      if(!GetSize())
	return TString(); //no bins
      
      if(i>=GetSize()){
	cout<<"Error Binner::BinName not enough bins "<<i <<GetSize()<<endl;
	return TString();
      }
      return fBinNames[i];
    }
    void Binner::InitBins(){
      if(!GetSize()){
	fBins.InitialiseBins();
	fBinNames=fBins.GetBinNames();
      }
    }
 
  
  }//namepsace FIT
}//namespace HS

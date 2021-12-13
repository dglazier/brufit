#include "Data.h"

#include <utility>

#include <utility>

#include <memory>

namespace HS{
  namespace FIT{

  

    DataEvents::DataEvents(Setup &setup,TString tname,const strings_t& files) :
      fSetup(&setup),fTreeName(std::move(std::move(tname))),fFileNames(files),
      fFiledTrees(files.size())
    {
	return;

    }
    TString DataEvents::GetItemName(Int_t ii){
      TString itemName;
      
      if(fNBoots>0&&!fBootStrap.get())
	BootStrap(fNBoots);//.recreate bootstrapper
      
      
      if(fBootStrap.get())
	itemName+=Form("Boot%d",fBootStrap->GetBootID(ii));
      
      else if(fNToys>0) //only toys if no bootstrap
	itemName+=Form("Toy%d",(Int_t) ii%fNToys);
      
      return itemName;
    }

    void DataEvents::Load(Setup &setup,const TString& tname,const strings_t& files)
    {
      fSetup=&setup;
      cout<<"DataEvents::Load "<<tname<<" with "<<files.size()<<" files"<<endl;
      //check if bootstrapping
      if(fBootStrap.get()){
	LoadBootStrap(tname,files);
	return;
      }
      //just load give files
      fTreeName=tname;
      fFileNames=files;
      fFiledTrees.resize(files.size());
 
    }
    
    void DataEvents::LoadBootStrap(const TString& tname,strings_t files)
    {
      cout<<"DataEvents::LoadBootStrap("<<tname <<" "<<files.size()<<endl;
      fBootStrap->SetOutDir(fSetup->GetOutDir());
      fFileNames.clear();
      //Loop over all the filenames (e.g different bins) and split the data
      for(auto &filename : files){
	fBootStrap->DivideData(tname,filename);
      }
      auto newFiles=fBootStrap->GetFileNames();
      fFileNames.insert( fFileNames.end(), newFiles.begin(), newFiles.end() );
      fTreeName=tname;
      fFiledTrees.resize(fFileNames.size());
    }
    void  DataEvents::LoadWeights(TString wname,TString fname,TString wobj){
      fWgtsConf=WeightsConfig{std::move(wname),std::move(fname),std::move(wobj)};
      LoadWeights();
    }
    void  DataEvents::LoadWeights(){
      fInWeights = std::unique_ptr<Weights>(new Weights{});
      //fInWeights->LoadSaved(fWgtsConf.File(),fWgtsConf.ObjName());
      fInWeights->LoadSavedDisc(fWgtsConf.File(),fWgtsConf.ObjName());
      fInWeights->PrintWeight();
      fInWeightName=fWgtsConf.Species().Data();
      fInWeightFile=fWgtsConf.File().Data();
      fInWeightObjName=fWgtsConf.ObjName().Data();
      cout<<"  DataEvents::LoadWeights using "<<fInWeightName<<" weights "<<fWgtsConf.File()<<" "<<fWgtsConf.ObjName()<<endl;
    }
 
    dset_uptr DataEvents::Get(const UInt_t iset) {

      if(fFileNames.size()<=iset)
	return dset_uptr();
      
      cout<<" RooAbsData& DataEvents::Get "<<" "<<fFileNames[iset]<<" tree "<<fTreeName<<" weights "<<fInWeights.get()<<" "<<fInWeightName<<endl;
      
      fFiledTrees[iset]=FiledTree::Read(fTreeName,fFileNames[iset]); //will be delted at end of function
  
     auto rawtree= fFiledTrees[iset]->Tree().get() ;
     auto vars = fSetup->DataVars();
     
     if(!fInWeights.get()&&fInWeightName!=TString()){ //if Data object read from root file
       LoadWeights();
     }

     const char* useWeightName=nullptr;
     if(fInWeights.get()){//if weights add branches and vars
       //create a copy in a new file to append the weights to
       //Keep it in file as large trees can use too much memory
       auto weightedFileTree=FiledTree::CloneFull(rawtree,fSetup->GetOutDir()+Form("/DataInWeightedTree%d.root",iset));
       //auto weightedFileTree=FiledTree::RecreateCopyFull(rawtree,fSetup->GetOutDir()+"DataInWeightedTree.root");
       fFiledTrees[iset].reset();
       fFiledTrees[iset]=std::move(weightedFileTree);

       rawtree= fFiledTrees[iset]->Tree().get() ;	
       //Add weights to tree
       fInWeights->AddToTree(rawtree);	
      //fInWeights->AddToTreeDisc(rawtree,fSetup->GetOutDir()+"DataInWeights.root");	
       fWeightVar = std::unique_ptr<RooRealVar>(new RooRealVar{fInWeightName,fInWeightName,0});
       vars.add(*fWeightVar.get());
       useWeightName=fInWeightName.Data(); //get char*
     }
     
     //only let datset clone active branches
     TIter iter=vars.createIterator();
     rawtree->SetBranchStatus("*",false);
     while(auto* arg=dynamic_cast<RooAbsArg*>(iter()))	
       rawtree->SetBranchStatus(arg->GetName(),true);	

     auto ds=std::unique_ptr<RooDataSet>(new RooDataSet{"DataEvents","DataEvents", rawtree,vars, fSetup->DataCut(),useWeightName});

     fFiledTrees[iset].reset(); //delete rawtree 
     if(fInWeights.get()){
       fInWeights.reset();
     }		
     
     
     ds->Print();
     return std::move(ds); 
    }
   
  }//namespace FIT
}//namespace HS

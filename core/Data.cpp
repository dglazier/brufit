/**
 * @file Data.cpp
 * @brief Implementation of the DataEvents class.
 */

#include "Data.h"

#include <utility>
#include <memory>
#include <iostream>

namespace HS {
namespace FIT {

    // ========================================================================
    // Initialization & Configuration
    // ========================================================================

    DataEvents::DataEvents(Setup &setup, TString tname, const strings_t& files) :
        fSetup(&setup), fTreeName(std::move(tname)), fFileNames(files) 
    {
    }
    
    TString DataEvents::GetItemName(Int_t ii) {
        TString itemName;
        
        if (fNBoots > 0 && !fBootStrap.get()) BootStrap(fNBoots);
        
        if (fBootStrap.get()) {
            itemName += Form("Boot%d", fBootStrap->GetBootID(ii));
        } else if (fNToys > 0) { 
            itemName += Form("Toy%d", (Int_t)ii % fNToys);
        }
        
        return itemName;
    }

    void DataEvents::Load(Setup &setup, const TString& tname, const strings_t& files) {
        fSetup = &setup;
        std::cout << "DataEvents::Load " << tname << " with " << files.size() 
                  << " files, boot : " << fBootStrap.get() << " " << fNBoots << std::endl;
        
        if (fBootStrap.get()) {
            LoadBootStrap(tname, files);
            return;
        }
        
        fTreeName = tname;
        fFileNames = files;
    }
    
    void DataEvents::LoadBootStrap(const TString& tname, strings_t files) {
        fBootStrap->SetOutDir(fSetup->GetOutDir());
        fFileNames.clear();
        
        // Loop over all the filenames (e.g. different bins) and split the data
        for (auto& filename : files) {
            fBootStrap->DivideData(tname, filename);
        }
        
        auto newFiles = fBootStrap->GetFileNames();
        fFileNames.insert(fFileNames.end(), newFiles.begin(), newFiles.end());
        fTreeName = tname;
    }

    // ========================================================================
    // Weight Handling
    // ========================================================================

    void DataEvents::LoadWeights(TString wname, TString fname, TString wobj) {
        fWgtsConf = WeightsConfig{std::move(wname), std::move(fname), std::move(wobj)};
        
        fInWeightName = fWgtsConf.Species().Data();
        fInWeightFile = fWgtsConf.File().Data();
        fInWeightObjName = fWgtsConf.ObjName().Data();
        
        // Construct the dummy variable ONCE here, thread-safely. 
        // It will be shared across workers, but RooDataSet takes its own copy upon import.
        fWeightVar = std::make_shared<RooRealVar>(fInWeightName, fInWeightName, 0);
        
        std::cout << "  DataEvents::LoadWeights configured using " << fInWeightName 
                  << " weights from " << fWgtsConf.File() << " " << fWgtsConf.ObjName() << std::endl;
    }
 
    // ========================================================================
    // Core Dataset Extraction
    // ========================================================================

    dset_uptr DataEvents::Get(const UInt_t iset) {
        if (fFileNames.size() <= iset) return dset_uptr();
      
        std::cout << " RooAbsData& DataEvents::Get " << fFileNames[iset] 
                  << " tree " << fTreeName << " applying weights: " 
                  << (fInWeightName != TString() ? fInWeightName.Data() : "None") << std::endl;

        // --- TRANSIENT FILE HANDLING ---
        auto filetree = FiledTree::Read(fTreeName, fFileNames[iset]); 
        auto rawtree = filetree->Tree().get();
        auto vars = fSetup->DataVars();
     
        const char* useWeightName = nullptr;
        
        // By instantiating a local 'Weights' object, we completely avoid thread collisions 
        if (fInWeightName != TString()) { 
            Weights localWeights;
            localWeights.LoadSavedDisc(fInWeightFile, fInWeightObjName);
            
            // Note: Added an extra random memory identifier to the file name to prevent 
            // any fringe chance of two threads processing the same iset and colliding on disk.
            TString tempFileName = fSetup->GetOutDir() + Form("/DataInWeightedTree_%d_%p.root", iset, (void*)rawtree);
            auto weightedFileTree = FiledTree::CloneFull(rawtree, tempFileName);
            
            // Safely transfer ownership to our local RAII wrapper
            filetree = std::move(weightedFileTree);
            rawtree = filetree->Tree().get();	
            
            localWeights.AddToTree(rawtree);	
            
            // Add the pre-configured weight variable
            vars.add(*fWeightVar.get());
            useWeightName = fInWeightName.Data(); 
        }
     
        // Only activate branches explicitly required by the PDF model
        rawtree->SetBranchStatus("*", false);
        for (auto* arg : vars) {
            rawtree->SetBranchStatus(arg->GetName(), true);	
        }
     
        // Import the unbinned TTree into RooFit's memory space
        auto ds = std::make_unique<RooDataSet>(
            "DataEvents", "DataEvents", vars, 
            RooFit::Import(*rawtree), 
            RooFit::Cut(fSetup->DataCut()), 
            RooFit::WeightVar(useWeightName)
        );

        ds->Print();
        return ds; 
    }
   
} // namespace FIT
} // namespace HS

// #include "Data.h"

// #include <utility>

// #include <utility>

// #include <memory>

// namespace HS{
//   namespace FIT{

  

//     DataEvents::DataEvents(Setup &setup,TString tname,const strings_t& files) :
//       fSetup(&setup),fTreeName(std::move(std::move(tname))),fFileNames(files),
//       fFiledTrees(files.size())
//     {
// 	return;

//     }
//     TString DataEvents::GetItemName(Int_t ii){
//       TString itemName;
      
//       if(fNBoots>0&&!fBootStrap.get())
// 	BootStrap(fNBoots);//.recreate bootstrapper
      
      
//       if(fBootStrap.get())
// 	itemName+=Form("Boot%d",fBootStrap->GetBootID(ii));
      
//       else if(fNToys>0) //only toys if no bootstrap
// 	itemName+=Form("Toy%d",(Int_t) ii%fNToys);
      
//       return itemName;
//     }

//     void DataEvents::Load(Setup &setup,const TString& tname,const strings_t& files)
//     {
//       fSetup=&setup;
//       cout<<"DataEvents::Load "<<tname<<" with "<<files.size()<<" files"<<" boot : "<<fBootStrap.get()<<" "<<fNBoots<<endl;
//       //check if bootstrapping
//       if(fBootStrap.get()){
// 	LoadBootStrap(tname,files);
// 	return;
//       }
//       //just load give files
//       fTreeName=tname;
//       fFileNames=files;
//       fFiledTrees.resize(files.size());
 
//     }
    
//     void DataEvents::LoadBootStrap(const TString& tname,strings_t files)
//     {
//       fBootStrap->SetOutDir(fSetup->GetOutDir());
//       fFileNames.clear();
//       //Loop over all the filenames (e.g different bins) and split the data
//       for(auto &filename : files){
// 	fBootStrap->DivideData(tname,filename);
//       }
//       auto newFiles=fBootStrap->GetFileNames();
//       fFileNames.insert( fFileNames.end(), newFiles.begin(), newFiles.end() );
//       fTreeName=tname;
//       fFiledTrees.resize(fFileNames.size());
//     }
//     void  DataEvents::LoadWeights(TString wname,TString fname,TString wobj){
//       fWgtsConf=WeightsConfig{std::move(wname),std::move(fname),std::move(wobj)};
//       LoadWeights();
//     }
//     void  DataEvents::LoadWeights(){
//       fInWeights = std::unique_ptr<Weights>(new Weights{});
//       //fInWeights->LoadSaved(fWgtsConf.File(),fWgtsConf.ObjName());
//       fInWeights->LoadSavedDisc(fWgtsConf.File(),fWgtsConf.ObjName());
//       fInWeights->PrintWeight();
//       fInWeightName=fWgtsConf.Species().Data();
//       fInWeightFile=fWgtsConf.File().Data();
//       fInWeightObjName=fWgtsConf.ObjName().Data();
//       cout<<"  DataEvents::LoadWeights using "<<fInWeightName<<" weights "<<fWgtsConf.File()<<" "<<fWgtsConf.ObjName()<<endl;
//     }
 
//     dset_uptr DataEvents::Get(const UInt_t iset) {

//       if(fFileNames.size()<=iset)
// 	return dset_uptr();
      
//       cout<<" RooAbsData& DataEvents::Get "<<" "<<fFileNames[iset]<<" tree "<<fTreeName<<" weights "<<fInWeights.get()<<" "<<fInWeightName<<endl;

      
//       fFiledTrees[iset]=FiledTree::Read(fTreeName,fFileNames[iset]); //will be delted at end of function
  
//      auto rawtree= fFiledTrees[iset]->Tree().get() ;
//      auto vars = fSetup->DataVars();
     
//      if(!fInWeights.get()&&fInWeightName!=TString()){ //if Data object read from root file
//        LoadWeights();
//      }

//      const char* useWeightName=nullptr;
//      if(fInWeights.get()){//if weights add branches and vars
//        //create a copy in a new file to append the weights to
//        //Keep it in file as large trees can use too much memory
//        auto weightedFileTree=FiledTree::CloneFull(rawtree,fSetup->GetOutDir()+Form("/DataInWeightedTree%d.root",iset));
//        //auto weightedFileTree=FiledTree::RecreateCopyFull(rawtree,fSetup->GetOutDir()+"DataInWeightedTree.root");
//        fFiledTrees[iset].reset();
//        fFiledTrees[iset]=std::move(weightedFileTree);

//        rawtree= fFiledTrees[iset]->Tree().get() ;	
//        //Add weights to tree
//        fInWeights->AddToTree(rawtree);	
//       //fInWeights->AddToTreeDisc(rawtree,fSetup->GetOutDir()+"DataInWeights.root");	
//        fWeightVar = std::unique_ptr<RooRealVar>(new RooRealVar{fInWeightName,fInWeightName,0});
//        vars.add(*fWeightVar.get());
//        useWeightName=fInWeightName.Data(); //get char*
//      }
     
//      //only let datset clone active branches
//      //TIter iter=vars.createIterator();
//      rawtree->SetBranchStatus("*",false);
//      //while(auto* arg=dynamic_cast<RooAbsArg*>(iter()))
//      for(auto* arg:vars){
//        rawtree->SetBranchStatus(arg->GetName(),true);	
//      }
     
//      //auto ds=std::unique_ptr<RooDataSet>(new RooDataSet{"DataEvents","DataEvents", rawtree,vars, fSetup->DataCut(),useWeightName});
//      auto ds=std::unique_ptr<RooDataSet>(new RooDataSet{"DataEvents","DataEvents", vars,RooFit::Import(*rawtree), RooFit::Cut(fSetup->DataCut()),RooFit::WeightVar(useWeightName)});

//      fFiledTrees[iset].reset(); //delete rawtree 
//      if(fInWeights.get()){
//        fInWeights.reset();
//      }		
     
     
//      ds->Print();
//      return std::move(ds); 
//     }
   
//   }//namespace FIT
// }//namespace HS

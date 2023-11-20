#pragma once

#include <TTree.h>
#include <RooDataSet.h>
#include <TString.h>
#include <iostream>
//#include <>
#include "Setup.h"
//#include ""

namespace HS{
  namespace FIT{
    void AddFormulaToTree(TTree* outTree, const RooDataSet& forTree, HS::FIT::Setup& setup){
  
      auto& pars = setup.Parameters();
      auto formulas=setup.ParameterFormulas(); //formulas that just depend on parameters, not variables/observables
      if(!formulas.getSize()) return;
 
      auto formVals = std::vector<Double_t>(formulas.getSize());
      auto formBranches=std::vector<TBranch*>(formulas.getSize());
 
      TIter iter=formulas.createIterator();
      Int_t iform=0;

      //create branches and make formVals their references
      while(auto* formu=dynamic_cast<RooFormulaVar*>(iter())){
	TString formuName=formu->GetName();
	formVals[iform]=0;
	formBranches[iform]=nullptr;
	formBranches[iform]=outTree->Branch(formuName,&formVals[iform],formuName+"/D");
	iform++;
      }

      //Loop over events in dataset
      for(int ids=0;ids<forTree.numEntries();++ids){
	auto eventPars = forTree.get(ids);
	pars.assignFast(*eventPars);
 
	//now calculate value of formula for these parameters
	iter.Reset();
	iform=0;
	while(auto* formu=dynamic_cast<RooFormulaVar*>(iter())){
   
	  formVals[iform]=formu->getValV();
	  formBranches[iform]->Fill();
	  //cout<< formBranches[iform]->GetName()<<" "<<formVals[iform]<<endl;
	  iform++;
   
	}
	if(ids%1000==0)std::cout<<"created event # "<<ids<<" out of "<<forTree.numEntries()<<std::endl;
      }
    }






  }
}

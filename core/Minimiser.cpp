#include "Minimiser.h"
#include <RooStats/RooStatsUtils.h>
#include <RooDataSet.h>
#include <TTree.h>
#include <TMath.h>

namespace HS{
  namespace FIT{
    
    Minuit::Minuit(UInt_t nrefits,Bool_t nozeroinit) : fNRefits(nrefits),fNoZeroInitialVal(nozeroinit),_saveDataSet{} {
      SetNameTitle("HSMinuit","Minuit minimiser");
    }
 
    void Minuit::Run(Setup &setup,RooAbsData &fitdata){
      fSetup=&setup;
      fData=&fitdata;
      
      //original fit using intial parameters
      FitTo();
      if(!fNRefits){ //One fit is enough...
	fResult->Print();
	
	return;
      }
      AppendFitToTree();
      
      //Perform many fits with differnt initial parameters
      UInt_t nrefit = 0;
      fSetup->SaveSnapShot(Form("Refit_%d",nrefit)); //original

      fLikelies.clear();
      using result_uptr=std::unique_ptr<RooFitResult>;
	
      vector<result_uptr> results;

      //save original likelihood and results
      StoreLikelihood(fLikelies);
      results.push_back(result_uptr{dynamic_cast<RooFitResult*>(fResult->clone())});

      //loop over refits
      while(nrefit++<fNRefits){

	RandomiseParameters();
	
	FitTo();

	AppendFitToTree();
	
	fSetup->SaveSnapShot(Form("Refit_%d",nrefit));
	
	StoreLikelihood(fLikelies);
	
	results.push_back(result_uptr{dynamic_cast<RooFitResult*>(fResult->clone())});
      }
	
      //Find the best fit result and save it
      Int_t best=std::distance(fLikelies.begin(), std::min_element(fLikelies.begin(), fLikelies.end()));
	
      cout<<"Minuit::FitTo() best fit likelihood "<<fLikelies[best]<<" from fit "<<best<<" all likelihoods "<<endl;

      for(auto& lh:fLikelies)
	cout<<lh<<" ";
      cout<<endl;

      fSetup->LoadSnapShot(Form("Refit_%d",best));
      fResult=dynamic_cast<RooFitResult*>(results[best]->clone());    
      fResult->Print();

      //cleanup
      fLikelies.clear();
      _Statuses.clear();
      
      return;
    }
    ////////////////////////////////////////////////////////////////
    file_uptr Minuit::SaveInfo(){
      
      TString fileName=fSetup->GetOutDir()+fSetup->GetName()+"/Results"+fSetup->GetTitle()+GetName()+".root";
      //TString fileName=fSetup->GetOutDir()+fSetup->GetName()+"/"+FileName();

      file_uptr file(TFile::Open(fileName,"recreate"));
      
      fSetup->Parameters().Print();
      fSetup->Yields().Print();
      //save paramters and chi2s in  dataset (for easy merging)
      RooArgSet saveArgs(fSetup->Parameters());
      saveArgs.add(fSetup->Yields());
     
      RooRealVar Nllval("NLL","NLL",fResult->minNll());
      saveArgs.add(Nllval);
      RooRealVar state("status","status",fResult->status());
      saveArgs.add(state);
	
      RooDataSet saveDS(FinalParName(),TString(GetName())+"Results",saveArgs);
      saveDS.add(saveArgs);
      saveDS.Write();
      
      TTree* treeDS=RooStats::GetAsTTree(ResultTreeName(),ResultTreeName(),saveDS);
      treeDS->Write();
      fResult->Write("MinuitResult");

      if(_treeDSbru!=nullptr){
	_treeDSbru->SetDirectory(file.get());
	_treeDSbru->Write();
      }
      // SaveRefits(saveArgs);

      //cleanup
      if(treeDS){treeDS=nullptr;}
      if(_treeDSbru){_treeDSbru=nullptr;}
      
      return std::move(file);
 
    }
    void Minuit::SaveRefits(RooArgSet& saveArgs){
      //save all refit parameters  into tree
      RooDataSet saveDS(FinalParName(),TString(GetName())+"Results",saveArgs);
      for(Int_t irefit=0;irefit<fNRefits;++irefit){

     
	fSetup->LoadSnapShot(Form("Refit_%d",irefit));
	//fResult=dynamic_cast<RooFitResult*>(results[irefit]->clone());    
	RooArgSet refitArgs(fSetup->Parameters());
	refitArgs.add(fSetup->Yields());
	
	RooRealVar Nllval("NLL","NLL",fLikelies[irefit]);
	refitArgs.add(Nllval);
	
	RooRealVar state("status","status",_Statuses[irefit]);
	refitArgs.add(state);
	
	saveDS.add(refitArgs);

      }
      saveDS.Print("v");
      TTree* treeDS=RooStats::GetAsTTree(ResultTreeName()+"_Refits",ResultTreeName()+"_Refits",saveDS);
      treeDS->Write();
  
    }
    
    void Minuit::StoreLikelihood(vector<Double_t> &fLikelies){
      if(!fResult) return;
      Bool_t nan=TMath::IsNaN(fResult->minNll());
      //check covariance OK or externally provide (SumW2Error) =-1
      Bool_t edm=(fResult->covQual()>1)||(fResult->covQual()==-1);
      Bool_t fail=(fResult->minNll()!=-1e+30);
      if((!nan)&&edm&&fail)fLikelies.push_back(fResult->minNll());
      else fLikelies.push_back(1E300);

      _Statuses.push_back(fResult->status());
    }
    
    void Minuit::RandomiseParameters(){
      cout<<"Minuit::RandomiseParameters()"<<endl;exit(0);
      Double_t intensity0=0;

	for(Int_t irand=0;irand<10000;++irand){
	  fSetup->RandomisePars(); 
	  if(fNoZeroInitialVal==kFALSE) //if don't care about 0 inital intensity break
	    break;
	  // if we do care keep trying until we get non-zero intensity
	  intensity0 = fSetup->Model()->getVal();
	  if (intensity0!=0) break;
	  if(irand==9999){
	    std::cerr<<"Minuit::Run FATAL : error tried to randomise parameters 9999 times without a non-zero intensity and you specified not have a non zero inital intensity. You may try removing the true from the Minuit constructor if you do not need this. Note on irand = "<<irand<<std::endl;
	    exit(0);
	  }
	}
     
	

    }
    //////////////////////////////////////////////////////////////

    void Minuit::AppendFitToTree(){
      

	RooArgSet isaveArgs(fSetup->Parameters());
	isaveArgs.add(fSetup->Yields());
     
	RooRealVar iNllval("NLL","NLL",fResult->minNll());
	isaveArgs.add(iNllval);
	RooRealVar istate("status","status",fResult->status());
	isaveArgs.add(istate);
 
   	RooDataSet saveDS(FinalParName(),TString(GetName())+"Results",isaveArgs);
 	saveDS.add(isaveArgs);

	auto treeDSbru=std::unique_ptr<TTree>(RooStats::GetAsTTree(ResultTreeName()+"Bru",ResultTreeName(),saveDS));

	TList list;
	if(_treeDSbru) list.Add(_treeDSbru);
	list.Add(treeDSbru.get());
	_treeDSbru= TTree::MergeTrees(&list);
    }
    
    //////////////////////////////////////////////////////////////
  Minuit2::Minuit2(UInt_t nrefits,Bool_t nozeroinit) : Minuit(nrefits,nozeroinit) {
      SetNameTitle("HSMinuit2","Minuit2 minimiser");
    }

    
  }

}

#pragma once

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooArgList.h>
#include <TString.h>
#include <TFile.h>
#include <map>

class MomentHelper{

 public:
  
  MomentHelper()=default;


  void Set(TString mom,Double_t val){
    _moments[mom]= val;
  }

  void Set(const TString& filename){
    auto file = std::unique_ptr<TFile>(TFile::Open(filename));
    auto fileMoments = file->Get<RooDataSet>("FinalParameters")->get();
    Set(fileMoments);
    Set("H_0_0_0", 2.00000);

  }
  void Set(const RooArgSet* pars){
    for(auto* mom:*pars){
      RooRealVar* rmom= dynamic_cast<RooRealVar*>(mom);
      if(rmom==nullptr) continue;
      if(mom->GetName()[0]!='H') continue;
      _moments[TString(mom->GetName())]=rmom->getVal();
    }
  }

  Double_t GetVal(const TString& name){
    if(_moments.find(name)==_moments.end()) return 0.0;
    return _moments[name];
  }
  
 private:
  
  std::map<TString,Double_t> _moments;
  
};

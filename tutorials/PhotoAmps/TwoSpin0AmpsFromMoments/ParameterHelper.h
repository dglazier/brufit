#pragma once

#include "Setup.h"
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <RooRealVar.h>
#include <map>
#include <vector>

enum class ParType{EParameter=-1,EDependancy=-2,EConstant=-3};

namespace m2pw{

  using HS::FIT::Setup;
  
  class  ParameterHelper{

  public:

    ParameterHelper() = default;
    
    ParameterHelper(const Setup& setup){
      auto& formulas=setup.ParameterFormulas();
      for(auto formp:formulas){
	auto var  = dynamic_cast<RooFormulaVar*>(formp);
	//valid formula contain H_ or a_ where the a_ is the last wave used for normalisation
	if(TString(formp->GetName()).Contains("H_") //a moment
	   ||TString(formp->GetName()).Contains("a_")) //or constraint amp
	   {
	     UseRooFormula(var);
	   }
      }//formulas
      
    }
    
  public :

    Bool_t AddParameter(const RooRealVar* var);

    Bool_t AddConstant(const RooRealVar* cvar);
    
    //    void UseRooFormula(const RooFormula* form);
    void UseRooFormula(const RooFormulaVar* form);

    std::vector<Double_t > ConstantValues(const RooFormulaVar* form) const;
    
    std::vector<Int_t > Indices(const RooFormulaVar* form) const;

    std::vector<TString > Dependencies(const RooFormulaVar* form) const;


    Bool_t IsConst(const TString& name)const {
      return (_constToIndex.find(name) == _constToIndex.end()) ? kFALSE :kTRUE;
    }
    Double_t GetConstVal(const TString& name) const {
      return IsConst(name)? _constVals[_constToIndex.at(name)] : 0;
      //      return IsConst(name)? _constVals[_constToIndex[name]] : 0;
    }
    Bool_t ParExists(const TString& name)const {
      return ( _nameToIndex.find(name) == _nameToIndex.end() ) ?false:true;
    }
    Int_t Index(const TString& name) const{
      if(IsConst(name)==true) return static_cast<Int_t>(ParType::EConstant);
      if(ParExists(name)==false) return static_cast<Int_t>(ParType::EDependancy);
      //not a dependency or a const so return variable index
      return _nameToIndex.at(name);
    }
 
    const TString& GetParName(Int_t ipar)const {return _names[ipar];}

    void Print(const TString opt="") const;

    Double_t Nvars() const {return _nextIndex;}
    Double_t Nconst() const {return _nextConstIndex;}

    Bool_t IsMagnitude(Int_t ipar)const {return _isMagnitude[ipar];}

    const Double_t* CurrentVals() const { return _currentVals.data();}
    std::vector<Double_t> CurrentValsCopy() const { return _currentVals;}

    void SetCurrentVal(const TString& name,Double_t val){
      SetCurrentVal(_nameToIndex[name],val);
    }
    void SetCurrentVal(UInt_t i,Double_t val){
      _currentVals[i]=val;
    }
    void SetCurrentVals(const double* pars){
      for(UInt_t i=0;i<_currentVals.size();++i)
	_currentVals[i]=pars[i];
    }
    
    void Randomise();
    void TransformConstrained(double* unitPars) const;
    void ZeroSmallAmps() ;
    void MakePertubation(Double_t psize);
    Bool_t CheckInRange(const double* pars) const ;

    Double_t SumMags() const;
    
    Double_t Min(UInt_t i) const {return _minVals[i];}
    Double_t Max(UInt_t i) const {return _maxVals[i];}

    const std::vector<Double_t>& StepSizes() const {return _stepSize;}
    void ScaleSteps(Double_t factor){
      for(auto& step:_stepSize)
	step*=factor;
    };

    const TString& GetName(Int_t ipar)const {return _names[ipar];}

    void MakeParTree(const TString& filename);
 
    void RerangeParameters();
 
    void FillTree(){
      RerangeParameters();
      _tree->Fill();
    }
    
    TTree* GetTree() const{
      return _tree;
    }
    
    void CloseTree(){
      _file->cd();
      _tree->Write();
      _file.reset(); //deletes _tree
      _tree=nullptr;
    }

  private :

    
    std::map<TString,Int_t> _nameToIndex;
    std::map<TString,Int_t> _constToIndex; //parameters that have been fixed constant
    //parameter info
    std::vector<TString> _names;  //order of parameters
    std::vector<Double_t> _minVals;
    std::vector<Double_t> _maxVals;
    std::vector<Double_t> _stepSize;
    std::vector<Bool_t > _isMagnitude;
    mutable Double_t _sumMags=0;
    //constant info
    std::vector<TString> _constNames;  //order of constants

    //current set of parameter values
    std::vector<Double_t> _currentVals;
    mutable std::vector<Double_t> _cachedVals;
    std::vector<Double_t> _constVals;
    
    UInt_t _nextIndex=0;
    UInt_t _nextConstIndex=0;

    std::unique_ptr<TFile> _file;
    TTree* _tree=nullptr; //only option for tree is regular pointer

  };//ParameterHelper


}//namepace m2pw

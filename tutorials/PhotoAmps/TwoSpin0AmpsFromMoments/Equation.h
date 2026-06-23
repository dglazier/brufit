#pragma once

#include "ParameterHelper.h"
#include <TString.h>
#include <TMath.h>
#include <TFormula.h>
#include <RooRealVar.h>
#include <RooFormulaVar.h>

namespace m2pw{

  static Bool_t IsConstraint(const TString& name){return (name[0] == 'H');}

  class Equation {

  public:
    Equation() = default;
    Equation(const Equation& other)= default;
    Equation(RooFormulaVar* var,ParameterHelper* pars,Double_t noise=0.E-20);
    
  
    const TString& GetName() const {return _name;}
    const TString& GetOrigFormula() const {return _origFormula;}
    
    Double_t EqnValue() const {return _eqnValue;}
    
    unsigned int LocalNDim() const {
      return _localNdim;
    }
    Double_t DoEval(const double* x) const ;
    Double_t DoEvalSq(const double* x) const { 
      auto result=DoEval(x);
      return result*result;
    }
    Bool_t NeedsRecalc() const{
      return _needsRecalc;
    }
    void SetNeedsRecalc() const {_needsRecalc=kTRUE;}
    void SetNoRecalc() const{_needsRecalc=kFALSE;}

    void FindDependencies(const std::vector<Equation >& eqns);

    void AddValBranch(TTree* tree) const{
      tree->Branch(GetName(),&_formuVal,GetName()+"/D");
    }
    void CalcFormuVal(const double* x) const{
      _formuVal = DoEval(x) - _eqnValue; //remove constraint term
    }
    void FindL(){
      auto  name = GetName();
      _L=TString(name(4,1)).Atoi();
      _LWeight = (2*_L+1);
    }
    void Print(const TString opt="") const;
    void SetEquationValue(Double_t val);

  protected:
    void SetName(const TString& name){_name=name;}
    void SetOrigFormula(const TString& name){_origFormula=name;}
 
    void OrganiseVariables(const double* x) const;

  private:
    
    mutable Double_t _cachedVal=0;
    Double_t _eqnValue=0;
    std::vector<Double_t> _constantVals;
    mutable std::vector<Double_t> _localX;
    mutable Double_t _formuVal=0;
    std::vector<const Equation*> _dependencies;
    Double_t _LWeight=1;
    
    std::vector<Int_t> _indices; //position of variables in _pars
  
    UInt_t _localNdim=0;
    UInt_t _L=0;
    
    RooFormulaVar _rooFormulaVar; //from RooFormulaVar
    TFormula _formula; //constraint equation (origformula - eqnValue)
    TString _name;
    TString _origFormula; //original formula string

    ParameterHelper* _pars=nullptr;
    // RooFormulaVar* _rooVar=nullptr;
  
    mutable Bool_t _needsRecalc=kTRUE;

  };//Equation


}

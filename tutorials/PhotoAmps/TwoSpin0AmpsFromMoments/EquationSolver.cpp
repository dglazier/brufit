#include "EquationSolver.h"
#include <TBenchmark.h>


namespace m2pw{
  EquationSolver *geqn_solver=nullptr;

  EquationSolver::EquationSolver(const Setup& setup, Double_t noise,std::vector<TString> noUse):
    _pars{setup}
  {
    //create a equation for each formula in setup
    //Currently rely on moments starting with H
    //been given constraints, other formulas left as are.
    auto& formulas=setup.ParameterFormulas();

    const auto validEqn = [&noUse](const TString& name) {
      for(const auto& match:noUse){
	if(name.Contains(match))
	  return kFALSE;
      }
      return kTRUE;
    };

    for(auto form:formulas){
      RooFormulaVar* var  = dynamic_cast<RooFormulaVar*>(form);
      //check if we are using it
      if(validEqn(var->GetName())==kFALSE) continue; 

     if( IsConstraint(var->GetName()) ) //e.g. Moments
	_eqns.push_back(Equation(var,&_pars,noise));
      else //e.g. normalisation equation
	_others.push_back(Equation(var,&_pars));
    }
    //Check if any equations depend on other formula
    for(auto& eqn:_eqns){
      eqn.FindDependencies(_others);
    }
    //initalise with current values
    for(auto& eqn:_eqns){
      eqn.FindL();
      eqn.DoEval(_pars.CurrentVals());
    }
    for(auto& eqn:_others){
      eqn.DoEval(_pars.CurrentVals());
    }

    
    //configure fitter
    _fitter.Config().SetParamsSettings(NDim(),_pars.CurrentVals());
    for(UInt_t i=0;i<NDim();i++){
      _fitter.Config().ParSettings(i).SetLimits(_pars.Min(i),_pars.Max(i));
    }
 
    if(geqn_solver==nullptr)geqn_solver=this;
    
  }
  //////////////////////////////////////////////////////////
  ///Override true equation values with measured etc
  void EquationSolver::SetEquationValues(MomentHelper& moms){
    for(auto& eqn:_eqns){
      eqn.FindL();
      //Set the moment to the value given in MomentsHelper
      eqn.SetEquationValue(moms.GetVal(eqn.GetName()));
      eqn.DoEval(_pars.CurrentVals());
    }
    // for(auto& eqn:_others){
    //   eqn.DoEval(_pars.CurrentVals());
    // }

  }
  ////////////////////////////////////////////////////////
  ///calculate all relevent variables and fill ParameterHelper Tree
  void EquationSolver::FillTree(){
    //my own values
    _val = DoEval(_pars.CurrentVals());
    _log_val=TMath::Log(_val);

    //equation values
    for(const auto& eqn:_eqns){
      eqn.CalcFormuVal(_pars.CurrentVals());
    }
    for(const auto& eqn:_others){
      eqn.CalcFormuVal(_pars.CurrentVals());
    }
    _pars.FillTree();

    //keep lowest values parameters
    if(_min_val>_val){
      _min_val = _val;
      _min_pars = _pars.CurrentValsCopy();
    }
  }
  ////////////////////////////////////////////////////////
  /// Make tree and branches
  void EquationSolver::MakeResultTree(const TString& name){
      _pars.MakeParTree(name);
      AddValBranch(_pars.GetTree());
      //equation values
      for(const auto& eqn:_eqns){
	eqn.AddValBranch(_pars.GetTree());
      }
      for(const auto& eqn:_others){
	eqn.AddValBranch(_pars.GetTree());
      }
   }
  /////////////////////////////////////////////////////////
  ///cache all equations together for dependecies
  void EquationSolver::CacheVals (const double *pars) const{
    for(auto& eqn:_eqns){//recalculate all
      eqn.SetNeedsRecalc();
    }
    for(auto& eqn:_others){//recalculate all
      eqn.SetNeedsRecalc();
    }
    
    for(auto& eqn:_others){//recalculate all
      eqn.DoEval(pars);
    }
    for(auto& eqn:_eqns){//recalculate all
      eqn.DoEval(pars);
    }
  
    // for(auto& eqn:_eqns){//recalculate all
    //   eqn.SetNoRecalc();
    // }
    // for(auto& eqn:_others){//recalculate all
    //   eqn.SetNoRecalc();
    // }
  }
  ///////////////////////////////////////////////////
  ///remove equations that match part of string
  // void EquationSolver::RemoveEquations(const TString& match){
  //   const auto pred = [&match](const Equation& eq) {
  //     return eq.GetName().Contains(match);
  //   };
  //   // _eqns.erase(std::remove_if(_eqns.begin(), _eqns.end(), pred), _eqns.end());
  //   Size_t ipos=0;
  //   for(auto& eqn:_eqns){//recalculate all
  //     if(pred(eqn)){
  // 	cout<<"Match "<<match<<endl;
  // 	_eqns.erase(_eqns.begin()+ipos);
  //     }
  //     ipos++;
  //   }
 
  // }


  ////////////////////////////////////////////////////////
  ///Solve the system of equations with minuit
  ROOT::Fit::FitResult EquationSolver::Solve(){
     _fitter.Config().SetMinimizer("Minuit2","Migrad");
    _fitter.Config().SetParamsSettings(NDim(),_pars.CurrentVals(),_pars.StepSizes().data());
 
     // gBenchmark->Start("solver");
     
     _fitter.FitFCN(NDim(),*this);
     // gBenchmark->Stop("solver");
     // gBenchmark->Print("solver");
     auto resultpars=_fitter.Result().GetParams();
     _pars.SetCurrentVals(resultpars);
     //     _val = DoEval(_pars.CurrentVals());
     // _log_val=TMath::Log(_val);
     return _fitter.Result();
    }
   

}

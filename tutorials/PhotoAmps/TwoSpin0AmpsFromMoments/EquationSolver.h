#pragma once

#include "ParameterHelper.h"
#include "MomentHelper.h"
#include "Equation.h"
#include "Setup.h"
#include <TString.h>
#include <Math/IParamFunction.h>
#include <Fit/Fitter.h>


namespace m2pw{

  using equations_t = std::vector<Equation>;
  using HS::FIT::Setup;


  class EquationSolver : public ROOT::Math::IMultiGradFunction {

  public:

    EquationSolver(const Setup& setup,Double_t noise=0.E-20,std::vector<TString> noUse={});
    ~EquationSolver() = default;
     
    void SetEquationValues(MomentHelper& moms);

      //virtual inherited functions
    ROOT::Math::IMultiGradFunction *Clone() const
    {
      // auto* that=new EquationSolver(*this);
      //return that;
      return nullptr;
    }
    
    double DoEval (const double *par) const {
      Double_t result=0;
      CacheVals(par);
      for(auto& eqn:_eqns){
	Double_t val = eqn.DoEvalSq(par);
	result+=(val);
      }
      _val = result;
      _log_val=TMath::Log(_val);
 
      return  result;
    }
  
    double DoDerivative(const double * x, unsigned int ipar) const {
      Double_t result=0;
      return result;
    }

    unsigned int NDim() const {return _pars.Nvars();}

    double operator() (const double *par) const {  
      return DoEval(par);
    }
    
    //////////Others
    
    ROOT::Fit::FitResult Solve();
    
    void CacheVals (const double *pars) const;

    void Print(const TString opt="") const{
      cout<<"EquationSolver : "<<endl;
      _pars.Print(opt);
      for(const auto& eqn:_eqns){
	eqn.Print(opt);
      }
      for(const auto& eqn:_others){
	eqn.Print(opt);
      }
      cout<<"\t Current Value = "<<_val<<endl;
    }

    void PrintResult(){
      _pars.SetCurrentVals(_min_pars.data());
      _pars.Print();
    }
    ParameterHelper& GetPars() { return _pars;}
    const equations_t& GetEqns() const { return _eqns;}
    
    ROOT::Fit::Fitter& Fitter(){return _fitter;}

    void MakeResultTree(const TString& name);
    
    void AddValBranch(TTree* tree) const{
      tree->Branch("val",&_val,"val/D");
      tree->Branch("log_val",&_log_val,"val/D");
    }
    void FillTree();
    Double_t GetVal() const {return _val;}

    // void RemoveEquations(const TString& match);
    
  private:
    
    ParameterHelper _pars;
    equations_t _eqns;
    equations_t _others;
    ROOT::Fit::Fitter _fitter;
    mutable Double_t _val=0;
    mutable Double_t _log_val=0;
    mutable Double_t _min_val=1E12;
    mutable std::vector<Double_t>_min_pars;
    
  };//EquationSolver

  extern EquationSolver *geqn_solver;

extern "C" {
  double cpp_eval(const double * params, size_t d=0) {
    Double_t result= geqn_solver->DoEval(params);
    return result;
  }
};

  
}

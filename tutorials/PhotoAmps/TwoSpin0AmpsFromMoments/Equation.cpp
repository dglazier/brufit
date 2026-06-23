#include "Equation.h"
#include <TRandom.h>

namespace m2pw{
  
  ///////////////////////////////////////////////////////////////
  Equation::Equation(RooFormulaVar* var,ParameterHelper* pars, Double_t noise):
    //_rooVar{var},
    _pars{pars},
    //_rooFormula{var->formula()}
    // _rooFormulaVar{var->GetName(),var->expression(),var->dependents()}
    _rooFormulaVar{*var}
  {
    SetName(var->GetName());
    SetOrigFormula(var->GetTitle());
    SetEquationValue(var->getVal()+gRandom->Gaus(0,noise));
    cout<<"\""<<GetName()<<"\""<<", "<<_eqnValue<<endl;
    
    //get indices so that _indices[i] gives local
    //position of parameter in this formula
    _indices = _pars->Indices(&_rooFormulaVar);
    _constantVals = _pars->ConstantValues(&_rooFormulaVar);

    _localNdim = _indices.size();
    _localX.resize(_localNdim);

  }
  //////////////////////////////////////////////////
  //Set the value for the equation formula
   void Equation::SetEquationValue(Double_t val){
     
    _eqnValue=0.0;

    //and now construct equation formula
    TString equation(_rooFormulaVar.expression());
    //In case of a moment
    //subtract val off from the formula, to get something =0
    if(IsConstraint(GetName())){
      _eqnValue=val;
      equation.Append(Form(" - %0.16f",_eqnValue));
     }
    //else leave alone as dependent
    
    //now make the formula to be minimised
    _formula=TFormula(GetName()+"_Eqn",equation);
  }
  
  ////////////////////////////////////////////////////
  ///calculate current value of cosntraint equation
  double Equation::DoEval (const double* x) const {

    if(NeedsRecalc()==kFALSE) return _cachedVal;
  
    OrganiseVariables(x);
    
    _cachedVal = _formula.EvalPar( _localX.data());
    
    SetNoRecalc();
    
    return  _cachedVal; 
  }
  //////////////////////////////////////////////////////////////////////////
  ///Find my dependencies so can use their current value
  void Equation::FindDependencies(const std::vector<Equation >& eqns){
     auto deps = _pars->Dependencies(&_rooFormulaVar);
     for(auto& depname:deps){
       bool got_it=false;
       //search for depname in the equations
       for(const auto& eqn:eqns){
	 if(eqn.GetName()==depname){
	   _dependencies.push_back(dynamic_cast<const Equation*>(&eqn));
	   got_it=true;
	   break;
	 }
       }
       if(got_it==false){
	 std::cout<<"TakeDependencies didn't find dependency equation  "<<depname<<" in " <<GetOrigFormula()<<std::endl;
       }
     }
     
   }
  //////////////////////////////////////////////////////////////////////////
  ///Get all the parameter values for this evaluation
  void Equation::OrganiseVariables(const double* x) const{
    UInt_t ilocal=0;
    UInt_t idep=0;
    UInt_t iconst=0;

    //local parameters can be real parameters, dependencies or constants
    for(auto& index:_indices){
      if(index>=0){
	_localX[ilocal]=x[index];
      }
      else if(index==static_cast<Int_t>(ParType::EDependancy)){//dependency
	_localX[ilocal]=_dependencies[idep]->DoEval(x);
	++idep;
      }
      else if(index==static_cast<Int_t>(ParType::EConstant)){//constant
	_localX[ilocal]=_constantVals[iconst];
	++iconst;
      }
      else cout<<"DoEval should not be here "<<endl;
      ++ilocal;
    }
  }
 
  /////////////////////////////////////////////////////////
  void Equation::Print(const TString opt) const{
    std::cout<<"\t Equation::Print() :" <<GetName()<<std::endl;
    std::cout<<"\t\t Formula = "<<_origFormula<<endl;
    std::cout<<"\t\t Constrained Value = "<< _cachedVal<<endl;
    std::cout<<"\t\t Equation Value = "<< _eqnValue<<endl;
    std::cout<<"\t\t _L = "<< _L <<" "<<_LWeight<<endl;
    //_formula.Print();
    if(opt!= "v") return;
    
    UInt_t isynch=0;
    std::cout<<"\t\t Dependencies : "<<endl;
    for(auto dep :_dependencies){
      std::cout<<"\t\t"<<dep->GetName()<<std::endl;
      isynch++;
    }
    isynch=0;
    
  }//Print

}

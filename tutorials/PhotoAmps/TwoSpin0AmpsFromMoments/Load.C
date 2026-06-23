#include <TROOT.h>

void Load(){
  gROOT->ProcessLine(".L ParameterHelper.cpp+");
  gROOT->ProcessLine(".L Equation.cpp+");
  gROOT->ProcessLine(".L EquationSolver.cpp+");
  gROOT->ProcessLine(".L MomentHelper.h+");
  
}

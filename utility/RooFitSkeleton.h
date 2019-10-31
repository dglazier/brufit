#pragma once

#include <TString.h>
#include <TMacro.h>

namespace HS{
  namespace FIT{
  
    class RooFitSkeleton   {
    
    
    public :
    
      RooFitSkeleton()=default;
      virtual ~RooFitSkeleton()=default;
    
    public :
    
      void MakeCode();

      void CreateRooFitEventsPDF(TString pdfName,TString obsNames,TString parNames);
      TString FindNextLineLike(TString linelike);
      void AddLineAfter(TString line0,TString line1,Int_t off=0);
      void ContinueLineAfter(TString line1,Int_t off=0);
      void MoveToLine(TString line0);
      void ReplaceMacroText(TString text0,TString text1);
      void ReplaceAllMacroText(TString text0,TString text1);
      void ReplaceInCurrLine(TString text0,TString text1);

    protected:
      
      TMacro fCurMacro;
      Int_t fPlace=0;
      TString fOption;
 
 
    };
  }//namespace FIT
}//namespace HS

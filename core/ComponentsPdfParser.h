////////////////////////////////////////////////////////////////
///
///Class:               ComponentsPdfParser
///Description:   Construct ComponentsPdf string
///

#pragma once
#include "PdfParser.h"

#include <utility>

namespace HS{
  namespace FIT{

   class ComponentsPdfParser : public PdfParser {
      
      
    public:
    ComponentsPdfParser(string name):PdfParser(std::move(name)){};
     ComponentsPdfParser() =default;
   
      //Deal with HS::FIT
      string ConstructPDF(string str) override;
      string NextComponentsTerm(string str, size_t& pos, size_t& tpos);
       
    private :

       
      
    };//class ComponentsPdfParser

  
  }
}

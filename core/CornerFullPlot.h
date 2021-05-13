#ifndef HS_FIT_CORNERFULLPLOT_h
#define HS_FIT_CORNERFULLPLOT_h

#include "Setup.h"
#include "PlotResults.h"
#include "RooMcmc.h"
#include <RooDataSet.h>
#include <RooHist.h>
#include <TList.h>
#include <TH2.h>

namespace HS{
  namespace FIT{

    class CornerFullPlot : public PlotResults {

    public:
      CornerFullPlot(Setup *setup, RooMcmc* mcmc, TList* canvases);
      CornerFullPlot()=default;
      CornerFullPlot(const CornerFullPlot& )=default;
      CornerFullPlot(CornerFullPlot&& )=default;
      virtual ~CornerFullPlot()=default;
      CornerFullPlot& operator = (const CornerFullPlot& other)=default;
      CornerFullPlot& operator = (CornerFullPlot&& other) = default;

    protected:
       
    private:

    };

  }//FIT
}//HS

#endif

#ifndef HS_FIT_CORNERPLOT_h
#define HS_FIT_CORNERPLOT_h

#include "Setup.h"
#include "PlotResults.h"
#include "RooMcmc.h"
#include <RooDataSet.h>
#include <RooHist.h>
#include <TList.h>
#include <TH2.h>

namespace HS{
  namespace FIT{

    class CornerPlot : public PlotResults {

    public:
      CornerPlot(Setup *setup, RooMcmc* mcmc,TList* canvases);
      CornerPlot()=default;
      CornerPlot(const CornerPlot& )=default;
      CornerPlot(CornerPlot&& )=default;
      virtual ~CornerPlot()=default;
      CornerPlot& operator = (const CornerPlot& other)=default;
      CornerPlot& operator = (CornerPlot&& other) = default;

    protected:

    private:

    };

  }//FIT
}//HS

#endif

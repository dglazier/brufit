#ifndef HS_FIT_AUTOCORRPLOT_h
#define HS_FIT_AUTOCORRPLOT_h

#include "Setup.h"
#include "PlotResults.h"
#include "RooMcmc.h"
#include <RooDataSet.h>
#include <RooHist.h>
#include <TList.h>

namespace HS{
  namespace FIT{

    class AutocorrPlot : public PlotResults {

    public:
      AutocorrPlot(Setup *setup,  RooMcmc* mcmc, TList* canvases);
      AutocorrPlot()=default;
      AutocorrPlot(const AutocorrPlot& )=default;
      AutocorrPlot(AutocorrPlot&& )=default;
      virtual ~AutocorrPlot()=default;
      AutocorrPlot& operator = (const AutocorrPlot& other)=default;
      AutocorrPlot& operator = (AutocorrPlot&& other) = default;

    protected:

    private:

    };

  }//Fit
}//HS

#endif

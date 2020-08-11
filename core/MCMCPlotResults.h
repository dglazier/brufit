#ifndef HS_FIT_MCMCPLOTRESULTS_h
#define HS_FIT_MCMCPLOTRESULTS_h

#include "Setup.h"
#include "PlotResults.h"
#include <RooDataSet.h>
#include <RooHist.h>
#include <TList.h>

namespace HS{
  namespace FIT{

    using roohist_uptr=std::unique_ptr<RooHist>;

    class MCMCPlotResults : public PlotResults {

    public:
      MCMCPlotResults(Setup *setup, const RooDataSet* data, const TString& tag, const TString& mcmcFile, const Int_t& NthDraw);
      MCMCPlotResults()=default;
      MCMCPlotResults(const MCMCPlotResults&)=default;
      MCMCPlotResults(MCMCPlotResults&&)=default;
      virtual ~MCMCPlotResults()=default;
      MCMCPlotResults& operator = (const MCMCPlotResults& other)=default;
      MCMCPlotResults& operator = (MCMCPlotResults&& other) = default;

      //   void PlotMCMCDataModel();

    protected:

    private:
      std::unique_ptr<TList> fCanvases{new TList()};
      std::vector<roohist_uptr> fRooHists;
      TString fFileName;
    };

  }//Fit
}//HS

#endif

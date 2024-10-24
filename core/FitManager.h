////////////////////////////////////////////////////////////////
///
///Class:               FitManager
///Description:
///           

#pragma once

#include "Setup.h"
#include "PlotResults.h"
#include "MCMCPlotResults.h"
#include "AutocorrPlot.h"
#include "CornerPlot.h"
#include "CornerFullPlot.h"
#include "RooMcmc.h"
#include "BruMcmc.h"
#include "Data.h"
#include "Binner.h"
#include "Minimiser.h"
#include <TNamed.h>
#include <RooMinimizer.h>
#include <RooAbsData.h>
#include <RooFitResult.h>

#include <utility>

#include <memory>


namespace HS{
  namespace FIT{

    using dataevs_ptr=std::shared_ptr<HS::FIT::DataEvents>;
    using strings_t = std::vector<TString>;
    using plotresult_uptr=std::unique_ptr<PlotResults>;
    
    
    class FitManager  : public TNamed{
      
    public:
      FitManager()=default;
      FitManager(const FitManager& other);
      FitManager(FitManager&&)=delete;
      ~FitManager() override =default;
      FitManager& operator=(const FitManager& other);
      FitManager& operator=(FitManager&& other) = delete;

      Setup *PointerSetUp() {return &fSetup;};
      Setup &SetUp() {return fSetup;};
      const Setup &ConstSetUp() {return fSetup;};
      Setup *CurrSetUp()  {return fCurrSetup.get();};

      //Note the default name and title are given by the bin and bootstrap
      //combination, Data GetGroup and GetItemName are BootStrap related
      //Default name= binname
      virtual TString GetCurrName(){return Bins().BinName(GetDataBin(fFiti));}
      //Default title data item (bootstrap number)
      virtual TString GetCurrTitle(){return Data().GetItemName(fFiti);}
      virtual Int_t GetDataBin(Int_t ii){ return Data().GetDataBin(ii);}
      virtual TString GetDataTreeName() {return fData.ParentTreeName();}
      virtual  strings_t GetDataFileNames() {return fData.FileNames();}
      
      void CopySetup(TObject* obj){fSetup=*(dynamic_cast<Setup*>(obj));}
      void CopyBinner(const Binner* obj){fBinner=*obj;}
      void SaveSetup();
      void CreateCurrSetup();
      void LoadSetup(const TString& dir);

      virtual void WriteThis();
      virtual void PreRun(){}//WriteThis();}
      
      Binner &Bins(){
	if(!fBinner.IsSetup())
	  fBinner.LoadSetup(fSetup);
	return fBinner;
      }
      const Binner *PointerBinner() const{return &fBinner;};

      virtual Int_t GetN(){return fData.GetN();}
      virtual Int_t GetFiti(){return fFiti;}
      
      virtual Bool_t Run();
      virtual void RunAll();
      virtual void RunOne(Int_t ifit);
      virtual void FitTo();
      
      virtual void Reset(){
	fFiledTrees.clear();

	//Keep last results
	fLastPars.reset(dynamic_cast<RooArgSet*>(fCurrSetup->ParsAndYields().snapshot()));
	fLastForms.reset(dynamic_cast<RooArgList*>(fCurrSetup->Formulas().snapshot()));
	
	fCurrSetup.reset();
	fCurrDataSet.reset();
      }
      
      void InitPrevResult(const TString& resultDir="",const TString& resultMinimiser="");
      void LoadPrevResult(const TString& resultDir,const TString& resultMinimiser);
      void IgnorePrevResult(){fUsePrevResult=kFALSE;}
      
      void LoadData(const TString& tname,const strings_t& fnames){
	 fData.Load(fSetup,tname,fnames);
      }
      //  void LoadData(const TString& tname,const TString& fname,const TString& name="Data"){
      void LoadData(const TString& tname,const TString& fname){
	const TString name="Data";
	fBinner.SplitData(tname,fname,name);
	LoadData(fBinner.TreeName(name),fBinner.FileNames(name));
	fData.SetParentName(fname);
 	fData.SetParentTreeName(tname);
      }
      void ReloadData(const TString& fname,const TString& name="Data"){
	fBinner.ReloadData(fname,name);
  	LoadData(fBinner.TreeName(name),fBinner.FileNames(name));
 	fData.SetParentName(fname);
 	fData.SetParentTreeName(fBinner.TreeName(name));
     }
      void ReloadData(const TString& tname,const TString& fname,const TString& name){
	ReloadData(fname,name);
      }
      
      void LoadSimulated(const TString& tname,const TString& fname,const TString& name){
	fBinner.SplitData(tname,fname,name);
      }
     void LoadSimulatedWithoutBinning(const TString& tname,const TString& fname,const TString& name){
       fBinner.SplitData(tname,fname,name); //create maps etc
       fBinner.SetAllFileNamesTo(fname,name); //but point all bins to same file
      }
      
      void ReloadSimulated(const TString& fname,const TString& name){
	fBinner.ReloadData(fname,name);
      }
      void ReloadSimulated(const TString& tname,const TString& fname,const TString& name){
	fBinner.ReloadData(fname,name);
      }
      
      void LoadGenerated(const TString& tname,TString fname,const TString& name, Bool_t ignoreCuts=kFALSE){
            if(ignoreCuts)
                  fBinner.SplitData(tname,std::move(fname),name+"__MCGen");
            else{
                  TString buffer = fBinner.GetCut();
                  fBinner.RemoveAllCuts();
                  fBinner.SplitData(tname,std::move(fname),name+"__MCGen");
                  fBinner.AddCut(buffer);
            }
      }
      void ReloadGenerated(const TString& fname,const TString& name){
	fBinner.ReloadData(fname,name+"__MCGen");
      }
      void ReloadGenerated(const TString& tname,const TString& fname,const TString& name){
	fBinner.ReloadData(fname,name+"__MCGen");
      }

      // dataevs_ptr& Data() {return fData;}
      DataEvents& Data() {return fData;}
      
      void SetMinimiser(Minimiser* mi){
	fMinimiser.reset(mi);
	SetMinimiserType(fMinimiser->GetName());
      }
      void SetMinimiserType(const TString& mtype){fMinimiserType=(mtype);}
      TString GetMinimiserType() const {return fMinimiserType;}
      //    Minimiser* GetMinimiser() const {return fMinimiser;}
      TString MinimiserFileName(){return TString("Results")+fMinimiserType+".root";}

      virtual void FillEventsPDFs();

      void PlotDataModel()
      {
	if(dynamic_cast<RooMcmc*>(fMinimiser.get()))
	  { 
	    fPlots.push_back((std::unique_ptr<MCMCPlotResults>(new MCMCPlotResults{fCurrSetup.get(),fCurrDataSet.get(),GetCurrName()+GetCurrTitle(),dynamic_cast<RooMcmc*>(fMinimiser.get()),fPlotOptions})));
	  }
	else if(dynamic_cast<BruMcmc*>(fMinimiser.get())){
	  fPlots.push_back((std::unique_ptr<MCMCPlotResults>(new MCMCPlotResults{fCurrSetup.get(),fCurrDataSet.get(),GetCurrName()+GetCurrTitle(),dynamic_cast<BruMcmc*>(fMinimiser.get()),fPlotOptions})));

	}
	else
	  fPlots.push_back((std::unique_ptr<PlotResults>(new PlotResults{fCurrSetup.get(),fCurrDataSet.get(),GetCurrName()+GetCurrTitle(),fPlotOptions})));
      }
      
      void RedirectOutput(const TString& log="");
      void SetRedirectOutput(){fRedirect=kTRUE;}

      void SetCompiledMacros(strings_t macs){
	fCompiledMacros=std::move(macs);
      }
      strings_t GetCompiledMacros(){return fCompiledMacros;}

      void SetPlotOptions(const TString& opt){fPlotOptions=opt;}
      void SetYieldMaxFactor(Double_t factor){fYldMaxFactor=factor;}
      //void SetIsSamplingIntegrals(){fIsSamplingIntegrals=kTRUE;}

      const RooArgSet* GetFitParameters() {
	if(fLastPars.get()==nullptr)//not fit, just use setup values
	  fLastPars.reset(dynamic_cast<RooArgSet*>(fSetup.ParsAndYields().snapshot()));
	return fLastPars.get();
     }
      const RooArgList* GetFitFormulas() {
	if(fLastForms.get()==nullptr)//not fit, just use setup values
	  fLastForms.reset(dynamic_cast<RooArgList*>(fSetup.Formulas().snapshot()));
	return fLastForms.get();
      }

      void TurnOffPlotting(){
	fDoPlotting = kFALSE;
      }

      void DoBinnedFits(Bool_t dbf=kTRUE){fuseBinnedFit=dbf;}

      void SetTruthPrefix(const TString& pre){fTruthPrefix=pre;}
 
     protected:
      std::unique_ptr<Setup> fCurrSetup={}; //!
      std::unique_ptr<RooDataSet> fCurrDataSet={}; //!
      
      virtual void SaveResults();
      //virtual void PlotSavedResults();
      
    private:
      
      Setup fSetup;
      
      DataEvents fData; 
      
      Binner fBinner;

      minimiser_uptr fMinimiser; //!
      TString fMinimiserType;
      
      std::vector<filed_uptr> fFiledTrees;//!
      std::vector<plotresult_uptr> fPlots;//!
      RooFitResult* fResult=nullptr;//!
      
      std::unique_ptr<RooArgSet> fLastPars=nullptr; //!
      std::unique_ptr<RooArgList> fLastForms=nullptr; //!

      strings_t fCompiledMacros; //!
  
      Bool_t fRedirect=kFALSE;

      UInt_t fFiti=0;
      Double_t fYldMaxFactor=2.;

      
      Bool_t fUsePrevResult=kFALSE;
      TString fPrevResultDir;
      TString fPrevResultMini;
      
      TString fTruthPrefix="xxxxxx";

      TString fPlotOptions;
      Bool_t fDoPlotting = kTRUE;
      Bool_t fuseBinnedFit = kFALSE;

      //Bool_t fIsSamplingIntegrals=kFALSE;
      
      ClassDefOverride(HS::FIT::FitManager,1);
     };

  }//namespace FIT
}//namespace HS


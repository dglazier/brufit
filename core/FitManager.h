/**
 * @file FitManager.h
 * @brief Main manager class for configuring and executing BruFit operations.
 * @author D. Glazier (Original) & Refactored for C++17
 */

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
#include <vector>

namespace HS {
namespace FIT {

    // Type aliases for cleaner code
    using dataevs_ptr     = std::shared_ptr<HS::FIT::DataEvents>;
    using strings_t       = std::vector<TString>;
    using plotresult_uptr = std::unique_ptr<PlotResults>;
    using minimiser_uptr  = std::unique_ptr<Minimiser>;

    /**
     * @class FitManager
     * @brief High-level controller for RooFit/BruFit execution.
     * * The FitManager is responsible for bringing together the model setup,
     * the data processing (binning), and the minimizer execution. It handles
     * the loop over data bins and manages the transient state for each fit.
     */
    class FitManager : public TNamed {
      
    public:
        /** @brief Default constructor. */
        FitManager() = default;
        
        /** @brief Copy constructor. Safely duplicates the fit configuration. */
        FitManager(const FitManager& other);
        
        /** @brief Move constructor deleted due to complex ROOT object ownership. */
        FitManager(FitManager&&) = delete;
        
        /** @brief Default destructor. Smart pointers handle cleanup. */
        ~FitManager() override = default;
        
        /** @brief Copy assignment operator. */
        FitManager& operator=(const FitManager& other);
        
        /** @brief Move assignment deleted due to complex ROOT object ownership. */
        FitManager& operator=(FitManager&& other) = delete;

        /** @name Setup Accessors
         * Methods to access the underlying physics model configuration.
         */
        ///@{
        Setup* PointerSetUp() { return &fSetup; }
        Setup& SetUp() { return fSetup; }
        const Setup& ConstSetUp() const { return fSetup; }
        Setup* CurrSetUp() { return fCurrSetup.get(); }
        ///@}

        /** @brief Gets the name of the current bin being processed. */
        virtual TString GetCurrName() { return Bins().BinName(GetDataBin(fFiti)); }
        
        /** @brief Gets the title (usually the bootstrap item name) of the current fit. */
        virtual TString GetCurrTitle() { return Data().GetItemName(fFiti); }
        
        virtual Int_t GetDataBin(Int_t ii) { return Data().GetDataBin(ii); }
        virtual TString GetDataTreeName() { return fData.ParentTreeName(); }
        virtual strings_t GetDataFileNames() { return fData.FileNames(); }
      
        /** @brief Copies a Setup configuration from a generic ROOT TObject. */
        void CopySetup(TObject* obj) { fSetup = *(dynamic_cast<Setup*>(obj)); }
        
        /** @brief Copies a Binner configuration. */
        void CopyBinner(const Binner* obj) { fBinner = *obj; }
        
        /** @brief Saves the current Setup object to a ROOT file. */
        void SaveSetup();
        
        /** @brief Creates a transient deep copy of the Setup for the current bin execution. */
        void CreateCurrSetup();
        
        /** @brief Loads a Setup configuration from a specified directory. */
        void LoadSetup(const TString& dir);

        /** @brief Serializes the FitManager state to a ROOT file. */
        virtual void WriteThis();
        virtual void PreRun() {}
      
        /** @brief Gets the Binner object, initializing it from Setup if necessary. */
        Binner& Bins() {
            if (!fBinner.IsSetup()) fBinner.LoadSetup(fSetup);
            return fBinner;
        }
        const Binner* PointerBinner() const { return &fBinner; }

        /** @brief Returns the total number of fits/bins to process. */
        virtual Int_t GetN() { return fData.GetN(); }
        
        /** @brief Returns the index of the current fit/bin being processed. */
        virtual Int_t GetFiti() { return fFiti; }
      
        /** @name Execution Controls */
        ///@{
        /** @brief Executes the fit for the current loaded dataset. 
         * @return True if successful, False if skipped (e.g., empty bin). */
        virtual Bool_t Run();
        
        /** @brief Loops over all data bins and executes fits sequentially. */
        virtual void RunAll();
        
        /** @brief Executes the fit for a specific bin index. 
         * @param ifit The index of the bin/dataset to fit. */
        virtual void RunOne(Int_t ifit);
        
        /** @brief Core minimization dispatch routine. Passes data to the configured Minimiser. */
        virtual void FitTo();
        ///@}
      
        /** @brief Resets transient state (current setup, datasets, files) after a fit completes. */
        virtual void Reset() {
            fFiledTrees.clear();
            // Keep last results securely
            fLastPars.reset(dynamic_cast<RooArgSet*>(fCurrSetup->ParsAndYields().snapshot()));
            fLastForms.reset(dynamic_cast<RooArgList*>(fCurrSetup->Formulas().snapshot()));
            
            // Release current run memory
            fCurrSetup.reset();
            fCurrDataSet.reset();
        }
      
        /** @brief Configures the manager to use parameters from a previous fit result. */
        void InitPrevResult(const TString& resultDir = "", const TString& resultMinimiser = "");
        void LoadPrevResult(const TString& resultDir, const TString& resultMinimiser);
        void IgnorePrevResult() { fUsePrevResult = kFALSE; }
      
        /** @name Data Ingestion APIs */
        ///@{
        void LoadData(const TString& tname, const strings_t& fnames) {
            fData.Load(fSetup, tname, fnames);
        }

        void LoadData(const TString& tname, const TString& fname) {
            const TString name = "Data";
            fBinner.SplitData(tname, fname, name);
            LoadData(fBinner.TreeName(name), fBinner.FileNames(name));
            fData.SetParentName(fname);
            fData.SetParentTreeName(tname);
        }

        void ReloadData(const TString& fname, const TString& name = "Data") {
            fBinner.ReloadData(fname, name);
            LoadData(fBinner.TreeName(name), fBinner.FileNames(name));
            fData.SetParentName(fname);
            fData.SetParentTreeName(fBinner.TreeName(name));
        }

        void ReloadData(const TString& tname, const TString& fname, const TString& name) {
            ReloadData(fname, name);
        }
      
        void LoadSimulated(const TString& tname, const TString& fname, const TString& name) {
            fBinner.SplitData(tname, fname, name);
        }

        void LoadSimulatedWithoutBinning(const TString& tname, const TString& fname, const TString& name) {
            fBinner.SplitData(tname, fname, name); 
            fBinner.SetAllFileNamesTo(fname, name); 
        }
      
        void ReloadSimulated(const TString& fname, const TString& name) {
            fBinner.ReloadData(fname, name);
        }

        void ReloadSimulated(const TString& tname, const TString& fname, const TString& name) {
            fBinner.ReloadData(fname, name);
        }
      
        void LoadGenerated(const TString& tname, TString fname, const TString& name, Bool_t ignoreCuts = kFALSE) {
            if (ignoreCuts) {
                fBinner.SplitData(tname, std::move(fname), name + "__MCGen");
            } else {
                TString buffer = fBinner.GetCut();
                fBinner.RemoveAllCuts();
                fBinner.SplitData(tname, std::move(fname), name + "__MCGen");
                fBinner.AddCut(buffer);
            }
        }

        void ReloadGenerated(const TString& fname, const TString& name) {
            fBinner.ReloadData(fname, name + "__MCGen");
        }

        void ReloadGenerated(const TString& tname, const TString& fname, const TString& name) {
            fBinner.ReloadData(fname, name + "__MCGen");
        }

        DataEvents& Data() { return fData; }
        ///@}
      
        /** @name Minimizer Configuration */
        ///@{
        /**
         * @brief Sets the minimizer engine (Legacy API).
         * @param mi Raw pointer to a Minimiser. FitManager takes ownership.
         */
        void SetMinimiser(Minimiser* mi) {
            fMinimiser.reset(mi);
            if (fMinimiser) SetMinimiserType(fMinimiser->GetName());
        }

        /**
         * @brief Sets the minimizer engine (Modern C++ API).
         * @param mi Unique pointer to a Minimiser.
         */
        void SetMinimiser(minimiser_uptr mi) {
            fMinimiser = std::move(mi);
            if (fMinimiser) SetMinimiserType(fMinimiser->GetName());
        }

        void SetMinimiserType(const TString& mtype) { fMinimiserType = mtype; }
        TString GetMinimiserType() const { return fMinimiserType; }
        TString MinimiserFileName() { return TString("Results") + fMinimiserType + ".root"; }
        ///@}

        virtual void FillEventsPDFs();
        
        /** @brief Generates visualizations of the data overlaid with the fit model. */
        void PlotDataModel();
      
        /** @brief Redirects standard output to a log file during fitting. */
        void RedirectOutput(const TString& log = "");
        void SetRedirectOutput() { fRedirect = kTRUE; }

        void SetCompiledMacros(strings_t macs) { fCompiledMacros = std::move(macs); }
        strings_t GetCompiledMacros() { return fCompiledMacros; }

        void SetPlotOptions(const TString& opt) { fPlotOptions = opt; }
        void SetYieldMaxFactor(Double_t factor) { fYldMaxFactor = factor; }

        const RooArgSet* GetFitParameters() {
            if (!fLastPars) // not fit, just use setup values
                fLastPars.reset(dynamic_cast<RooArgSet*>(fSetup.ParsAndYields().snapshot()));
            return fLastPars.get();
        }

        const RooArgList* GetFitFormulas() {
            if (!fLastForms) // not fit, just use setup values
                fLastForms.reset(dynamic_cast<RooArgList*>(fSetup.Formulas().snapshot()));
            return fLastForms.get();
        }

        /** @brief Disables automatic plotting after fits. */
        void TurnOffPlotting() { fDoPlotting = kFALSE; }
        
        /** @brief Forces RooFit to perform a binned fit instead of unbinned. */
        void DoBinnedFits(Bool_t dbf = kTRUE) { fuseBinnedFit = dbf; }
        
        void SetTruthPrefix(const TString& pre) { fTruthPrefix = pre; }
 
    protected:
        std::unique_ptr<Setup> fCurrSetup;        ///< Transient active physics model
        std::unique_ptr<RooAbsData> fCurrDataSet; ///< Transient active dataset being fitted
      
        virtual void SaveResults();
      
    private:
      
        Setup fSetup;                             ///< Master configuration template
        DataEvents fData;                         ///< Data ingestion manager
        Binner fBinner;                           ///< Kinematic binning manager

        minimiser_uptr fMinimiser;                ///<! Do not serialize active minimizer engine
        TString fMinimiserType;                   ///< String identifier for the minimizer
      
        std::vector<filed_uptr> fFiledTrees;      ///<! Transient file handlers keeping TTrees open
        std::vector<plotresult_uptr> fPlots;      ///<! Transient plotting objects
        RooFitResult* fResult = nullptr;          ///<! Active fit result pointer
      
        std::unique_ptr<RooArgSet> fLastPars;     ///<! Cache of parameters from the previous bin
        std::unique_ptr<RooArgList> fLastForms;   ///<! Cache of formulas from the previous bin

        strings_t fCompiledMacros;                ///<! Helper macros compiled via ACLiC
  
        Bool_t fRedirect = kFALSE;                ///< Flag indicating if stdout is redirected
        UInt_t fFiti = 0;                         ///< Current bin index executing
        Double_t fYldMaxFactor = 2.;              ///< Constraint factor for maximum species yields
      
        Bool_t fUsePrevResult = kFALSE;           ///< Flag to seed current fit with previous bin's results
        TString fPrevResultDir;                   ///< Directory containing previous results
        TString fPrevResultMini;                  ///< Minimizer used in previous results
      
        TString fTruthPrefix = "xxxxxx";          ///< Prefix for truth-level branch matching

        TString fPlotOptions;                     ///< Formatting options for generated plots
        Bool_t fDoPlotting = kTRUE;               ///< Master flag for plot generation
        Bool_t fuseBinnedFit = kFALSE;            ///< Master flag for forced binned execution

        ClassDefOverride(HS::FIT::FitManager, 1);
    };

} // namespace FIT
} // namespace HS



////////////////////////////////////////////////////////////////
///
///Class:               FitManager
///Description:
///           

// #pragma once

// #include "Setup.h"
// #include "PlotResults.h"
// #include "MCMCPlotResults.h"
// #include "AutocorrPlot.h"
// #include "CornerPlot.h"
// #include "CornerFullPlot.h"
// #include "RooMcmc.h"
// #include "BruMcmc.h"
// #include "Data.h"
// #include "Binner.h"
// #include "Minimiser.h"
// #include <TNamed.h>
// #include <RooMinimizer.h>
// #include <RooAbsData.h>
// #include <RooFitResult.h>

// #include <utility>

// #include <memory>


// namespace HS{
//   namespace FIT{

//     using dataevs_ptr=std::shared_ptr<HS::FIT::DataEvents>;
//     using strings_t = std::vector<TString>;
//     using plotresult_uptr=std::unique_ptr<PlotResults>;
    
    
//     class FitManager  : public TNamed{
      
//     public:
//       FitManager()=default;
//       FitManager(const FitManager& other);
//       FitManager(FitManager&&)=delete;
//       ~FitManager() override =default;
//       FitManager& operator=(const FitManager& other);
//       FitManager& operator=(FitManager&& other) = delete;

//       Setup *PointerSetUp() {return &fSetup;};
//       Setup &SetUp() {return fSetup;};
//       const Setup &ConstSetUp() {return fSetup;};
//       Setup *CurrSetUp()  {return fCurrSetup.get();};

//       //Note the default name and title are given by the bin and bootstrap
//       //combination, Data GetGroup and GetItemName are BootStrap related
//       //Default name= binname
//       virtual TString GetCurrName(){return Bins().BinName(GetDataBin(fFiti));}
//       //Default title data item (bootstrap number)
//       virtual TString GetCurrTitle(){return Data().GetItemName(fFiti);}
//       virtual Int_t GetDataBin(Int_t ii){ return Data().GetDataBin(ii);}
//       virtual TString GetDataTreeName() {return fData.ParentTreeName();}
//       virtual  strings_t GetDataFileNames() {return fData.FileNames();}
      
//       void CopySetup(TObject* obj){fSetup=*(dynamic_cast<Setup*>(obj));}
//       void CopyBinner(const Binner* obj){fBinner=*obj;}
//       void SaveSetup();
//       void CreateCurrSetup();
//       void LoadSetup(const TString& dir);

//       virtual void WriteThis();
//       virtual void PreRun(){}//WriteThis();}
      
//       Binner &Bins(){
// 	if(!fBinner.IsSetup())
// 	  fBinner.LoadSetup(fSetup);
// 	return fBinner;
//       }
//       const Binner *PointerBinner() const{return &fBinner;};

//       virtual Int_t GetN(){return fData.GetN();}
//       virtual Int_t GetFiti(){return fFiti;}
      
//       virtual Bool_t Run();
//       virtual void RunAll();
//       virtual void RunOne(Int_t ifit);
//       virtual void FitTo();
      
//       virtual void Reset(){
// 	fFiledTrees.clear();

// 	//Keep last results
// 	fLastPars.reset(dynamic_cast<RooArgSet*>(fCurrSetup->ParsAndYields().snapshot()));
// 	fLastForms.reset(dynamic_cast<RooArgList*>(fCurrSetup->Formulas().snapshot()));
	
// 	fCurrSetup.reset();
// 	fCurrDataSet.reset();
//       }
      
//       void InitPrevResult(const TString& resultDir="",const TString& resultMinimiser="");
//       void LoadPrevResult(const TString& resultDir,const TString& resultMinimiser);
//       void IgnorePrevResult(){fUsePrevResult=kFALSE;}
      
//       void LoadData(const TString& tname,const strings_t& fnames){
// 	 fData.Load(fSetup,tname,fnames);
//       }
//       //  void LoadData(const TString& tname,const TString& fname,const TString& name="Data"){
//       void LoadData(const TString& tname,const TString& fname){
// 	const TString name="Data";
// 	fBinner.SplitData(tname,fname,name);
// 	LoadData(fBinner.TreeName(name),fBinner.FileNames(name));
// 	fData.SetParentName(fname);
//  	fData.SetParentTreeName(tname);
//       }
//       void ReloadData(const TString& fname,const TString& name="Data"){
// 	fBinner.ReloadData(fname,name);
//   	LoadData(fBinner.TreeName(name),fBinner.FileNames(name));
//  	fData.SetParentName(fname);
//  	fData.SetParentTreeName(fBinner.TreeName(name));
//      }
//       void ReloadData(const TString& tname,const TString& fname,const TString& name){
// 	ReloadData(fname,name);
//       }
      
//       void LoadSimulated(const TString& tname,const TString& fname,const TString& name){
// 	fBinner.SplitData(tname,fname,name);
//       }
//      void LoadSimulatedWithoutBinning(const TString& tname,const TString& fname,const TString& name){
//        fBinner.SplitData(tname,fname,name); //create maps etc
//        fBinner.SetAllFileNamesTo(fname,name); //but point all bins to same file
//       }
      
//       void ReloadSimulated(const TString& fname,const TString& name){
// 	fBinner.ReloadData(fname,name);
//       }
//       void ReloadSimulated(const TString& tname,const TString& fname,const TString& name){
// 	fBinner.ReloadData(fname,name);
//       }
      
//       void LoadGenerated(const TString& tname,TString fname,const TString& name, Bool_t ignoreCuts=kFALSE){
//             if(ignoreCuts)
//                   fBinner.SplitData(tname,std::move(fname),name+"__MCGen");
//             else{
//                   TString buffer = fBinner.GetCut();
//                   fBinner.RemoveAllCuts();
//                   fBinner.SplitData(tname,std::move(fname),name+"__MCGen");
//                   fBinner.AddCut(buffer);
//             }
//       }
//       void ReloadGenerated(const TString& fname,const TString& name){
// 	fBinner.ReloadData(fname,name+"__MCGen");
//       }
//       void ReloadGenerated(const TString& tname,const TString& fname,const TString& name){
// 	fBinner.ReloadData(fname,name+"__MCGen");
//       }

//       // dataevs_ptr& Data() {return fData;}
//       DataEvents& Data() {return fData;}
      
//       void SetMinimiser(Minimiser* mi){
// 	fMinimiser.reset(mi);
// 	SetMinimiserType(fMinimiser->GetName());
//       }
//       void SetMinimiserType(const TString& mtype){fMinimiserType=(mtype);}
//       TString GetMinimiserType() const {return fMinimiserType;}
//       //    Minimiser* GetMinimiser() const {return fMinimiser;}
//       TString MinimiserFileName(){return TString("Results")+fMinimiserType+".root";}

//       virtual void FillEventsPDFs();

//       void PlotDataModel()
//       {
// 	if(dynamic_cast<RooMcmc*>(fMinimiser.get()))
// 	  { 
// 	    fPlots.push_back((std::unique_ptr<MCMCPlotResults>(new MCMCPlotResults{fCurrSetup.get(),fCurrDataSet.get(),GetCurrName()+GetCurrTitle(),dynamic_cast<RooMcmc*>(fMinimiser.get()),fPlotOptions})));
// 	  }
// 	else if(dynamic_cast<BruMcmc*>(fMinimiser.get())){
// 	  fPlots.push_back((std::unique_ptr<MCMCPlotResults>(new MCMCPlotResults{fCurrSetup.get(),fCurrDataSet.get(),GetCurrName()+GetCurrTitle(),dynamic_cast<BruMcmc*>(fMinimiser.get()),fPlotOptions})));

// 	}
// 	else
// 	  fPlots.push_back((std::unique_ptr<PlotResults>(new PlotResults{fCurrSetup.get(),fCurrDataSet.get(),GetCurrName()+GetCurrTitle(),fPlotOptions})));
//       }
      
//       void RedirectOutput(const TString& log="");
//       void SetRedirectOutput(){fRedirect=kTRUE;}

//       void SetCompiledMacros(strings_t macs){
// 	fCompiledMacros=std::move(macs);
//       }
//       strings_t GetCompiledMacros(){return fCompiledMacros;}

//       void SetPlotOptions(const TString& opt){fPlotOptions=opt;}
//       void SetYieldMaxFactor(Double_t factor){fYldMaxFactor=factor;}
//       //void SetIsSamplingIntegrals(){fIsSamplingIntegrals=kTRUE;}

//       const RooArgSet* GetFitParameters() {
// 	if(fLastPars.get()==nullptr)//not fit, just use setup values
// 	  fLastPars.reset(dynamic_cast<RooArgSet*>(fSetup.ParsAndYields().snapshot()));
// 	return fLastPars.get();
//      }
//       const RooArgList* GetFitFormulas() {
// 	if(fLastForms.get()==nullptr)//not fit, just use setup values
// 	  fLastForms.reset(dynamic_cast<RooArgList*>(fSetup.Formulas().snapshot()));
// 	return fLastForms.get();
//       }

//       void TurnOffPlotting(){
// 	fDoPlotting = kFALSE;
//       }

//       void DoBinnedFits(Bool_t dbf=kTRUE){fuseBinnedFit=dbf;}

//       void SetTruthPrefix(const TString& pre){fTruthPrefix=pre;}
 
//      protected:
//       std::unique_ptr<Setup> fCurrSetup={}; //!
//       std::unique_ptr<RooDataSet> fCurrDataSet={}; //!
      
//       virtual void SaveResults();
//       //virtual void PlotSavedResults();
      
//     private:
      
//       Setup fSetup;
      
//       DataEvents fData; 
      
//       Binner fBinner;

//       minimiser_uptr fMinimiser; //!
//       TString fMinimiserType;
      
//       std::vector<filed_uptr> fFiledTrees;//!
//       std::vector<plotresult_uptr> fPlots;//!
//       RooFitResult* fResult=nullptr;//!
      
//       std::unique_ptr<RooArgSet> fLastPars=nullptr; //!
//       std::unique_ptr<RooArgList> fLastForms=nullptr; //!

//       strings_t fCompiledMacros; //!
  
//       Bool_t fRedirect=kFALSE;

//       UInt_t fFiti=0;
//       Double_t fYldMaxFactor=2.;

      
//       Bool_t fUsePrevResult=kFALSE;
//       TString fPrevResultDir;
//       TString fPrevResultMini;
      
//       TString fTruthPrefix="xxxxxx";

//       TString fPlotOptions;
//       Bool_t fDoPlotting = kTRUE;
//       Bool_t fuseBinnedFit = kFALSE;

//       //Bool_t fIsSamplingIntegrals=kFALSE;
      
//       ClassDefOverride(HS::FIT::FitManager,1);
//      };

//   }//namespace FIT
// }//namespace HS


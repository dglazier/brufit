/**
 * @file Setup.h
 * @brief Class for configuring and managing RooFit Workspaces and PDFs.
 * @details The Setup class acts as a high-level wrapper around a RooWorkspace. 
 * It parses string commands to build complex PDFs, manages variables, and handles 
 * constraint organization. It implements a string-replay mechanism for safe 
 * copying across threads (e.g., for PROOF or RDataFrame).
 */

#pragma once

#include "PdfParser.h"

#include <RooStats/ModelConfig.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooDataSet.h>
#include <TNamed.h>
#include <TString.h>
#include <TSystem.h>

#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <map>

namespace HS {
namespace FIT {

    using std::vector;
    using std::string;
    using std::cout;
    using std::endl;

    // Modern C++ Type Aliases
    using realvar_ptr = std::shared_ptr<RooRealVar>;
    using realvars_t  = std::vector<RooRealVar*>;
    using catvar_ptr  = std::shared_ptr<RooCategory>;
    using catvars_t   = std::vector<RooCategory*>;
    using strings_t   = std::vector<TString>;

    // Helper functions
    RooArgSet MakeArgSet(realvars_t vars);
    RooArgSet MakeArgSet(catvars_t cats);
    void SetAllValLimits(RooArgList&, Double_t val, Double_t low = 0, Double_t high = 0);
    void ReadFormula(TString forma, strings_t& svars, strings_t& sranges);
    Bool_t ArgListContainsName(RooArgList& items, TString name);

    class RandomConstrained;
    
    /**
     * @class Setup
     * @brief High-level controller for RooFit PDF assembly.
     */
    class Setup : public TNamed {
      
    public:
        /** @brief Default constructor. */
        Setup();
        
        /** @brief Construct with a specific name. */
        explicit Setup(const TString& name);
        
        /** @brief Thread-safe copy constructor (rebuilds Workspace via string replay). */
        Setup(const Setup& other);
        
        /** @brief Move semantics deleted due to strict ROOT workspace ownership rules. */
        Setup(Setup&&) = delete;
        
        /** @brief Destructor. Cleans up any manually owned RooFit models. */
        ~Setup() override {
            if (fModel) delete fModel;
            fModel = nullptr;
        }
        
        /** @brief Thread-safe assignment operator (rebuilds Workspace via string replay). */
        Setup& operator=(const Setup& other);
        Setup& operator=(Setup&& other) = delete;

        /** @name PDF and Variable Parsers */
        ///@{
        /** @brief Passes a factory string directly to the underlying RooWorkspace. */
        void FactoryPDF(TString opt);
        
        /** @brief Parses custom BruFit string syntax to assemble complex PDFs. */
        void ParserPDF(TString str, PdfParser& parse);
        
        void LoadVariable(const TString& opt);
        void LoadCategory(const TString& opt);
        void LoadAuxVar(const TString& opt);
        void LoadFormula(TString formu);
        void LoadParameter(const TString& opt);
        void LoadConstant(const TString& opt);
        void LoadFunctionVar(const TString& opt);
        
        /** @brief Loads a PDF as an individual species component for an extended ML fit. */
        void LoadSpeciesPDF(TString opt, Float_t Scale0 = 1);
        
        /** @brief Sums all loaded species PDFs into a single RooAddPdf model. */
        void TotalPDF();
        ///@}
      
        /** @name Accessors */
        ///@{
        const RooWorkspace& WS() { return fWS; }
        RooAbsPdf* Model() const { return fModel; }
        RooAddPdf ExtendModel() const { return RooAddPdf(*dynamic_cast<RooAddPdf*>(fModel)); }

        RooArgSet& DataVars();
        RooArgSet& Cats();
        RooArgSet& FitVarsAndCats();
        RooArgSet& ParsAndYields();
        RooArgSet& NonConstParsAndYields();
        void CopyRealProperties(RooArgSet& change, const RooArgSet& tothis);
        
        const realvars_t& AuxVars() const { return fAuxVars; }
        ///@}

        /** @name Cut and Directory Management */
        ///@{
        /** @brief Returns the combined string of logical cuts applied to the dataset. */
        const TString Cut() const {
            if (fAddCut.Sizeof() > 1 && fVarCut.Sizeof() > 1) return fAddCut + "&&" + fVarCut;
            else if (fAddCut.Sizeof() > 1) return fAddCut;
            else if (fVarCut.Sizeof() > 1) return fVarCut;
            else return TString();
        }	
        
        void AddCut(const TString& cut) {
            if (fAddCut.Sizeof() > 1) { fAddCut += "&&"; }
            fAddCut += cut;
        }
     
        void SetDataOnlyCut(TString cut) { fDataOnlyCut = std::move(cut); }
        const TString DataCut() const {
            auto cut = Cut();
            if (cut.Sizeof() > 1 && fDataOnlyCut.Sizeof() > 1) return cut + "&&" + fDataOnlyCut;
            else if (cut.Sizeof() > 1) return cut;
            else if (fDataOnlyCut.Sizeof() > 1) return fDataOnlyCut;
            else return TString();
        }
 
        const TString GetIDBranchName() const { return fIDBranchName; }
        void SetIDBranchName(TString name) { fIDBranchName = std::move(name); }
        
        const TString GetOutDir() const {
            if (fOutDir == TString()) return "./";
            return fOutDir + "/";
        }
        
        void SetOutDir(TString name) {
            if (!name.BeginsWith("/")) name = TString(gSystem->WorkingDirectory()) + "/" + name;
            fOutDir = name;
            gSystem->Exec(Form("mkdir -p %s", fOutDir.Data()));
        }
        ///@}

        /** @name Parameter Iterators and Utilities */
        ///@{
        const realvars_t& FitVars() const { return fFitVars; }
        const catvars_t& FitCats() const { return fFitCats; }
        
        RooArgList& Yields() { return fYields; }
        RooArgList& Parameters() { return fParameters; }
        RooArgList& Constants() { return fConstants; }
        const RooArgList& Formulas() const { return fFormulas; }
        const RooArgList& ParameterFormulas() const { return fParameterFormulas; }
        RooArgList& PDFs() { return fPDFs; }
        const RooArgList& constPDFs() const { return fPDFs; }
        RooArgList& Constraints() { return fConstraints; }

        Double_t SumOfYields();
        
        void AddGausConstraint(RooGaussian *pdf);
        void AddFormulaConstraint(RooFormulaVar *formu);
        void AddPdfConstraint(RooAbsPdf *pdf);
        
        void AddFitOption(const RooCmdArg& cmd) { fFitOptions.Add(dynamic_cast<RooCmdArg*>(cmd.Clone())); }
        RooLinkedList FitOptions() { return fFitOptions; }
        
        void DefaultFitOptions() { AddFitOption(RooFit::Warnings(kFALSE)); }
        void RequiredFitOptions();

        void ErrorsSumW2() { fErrorsSumW2 = kTRUE; fErrorsAsym = kFALSE; }
        void ErrorsAsymp() { fErrorsSumW2 = kFALSE; fErrorsAsym = kTRUE; }
        void ErrorsWrong() { fErrorsSumW2 = kFALSE; fErrorsAsym = kFALSE; }
        
        void RandomisePars();
        void OrganiseConstraints();

        RooRealVar* GetVar(const TString& name) { return (dynamic_cast<RooRealVar*>(DataVars().find(name))); }
        RooRealVar* GetPar(const TString& name) { return (dynamic_cast<RooRealVar*>(fParameters.find(name))); }
        
        void SetParVal(const TString& par, Double_t val, Bool_t co = kFALSE);
        void SetYldVal(const TString& yld, Double_t val, Bool_t co = kFALSE);
        void SetConstPar(const TString& par, Bool_t co = kTRUE);
        void SetConstPDFPars(const TString& pdf, Bool_t co = kTRUE);
        Bool_t IsParSetConst(const TString& name) const;
        
        void SaveSnapShot(const TString& name) { fWS.saveSnapshot(name, RooArgSet(fYields, fParameters), kTRUE); }
        void LoadSnapShot(const TString& name) { fWS.loadSnapshot(name); }

        RooStats::ModelConfig* GetModelConfig();
        TString GetPDFInWeights(const TString& name) { return fPDFInWeights[name]; }

        RooAbsPdf* ComponentsPDF(TString opt);
        ///@}
      
    protected:
        void LoadParameterOnTheFly(const TString& opt);

    private:
        realvars_t fFitVars; 
        realvars_t fAuxVars;
        catvars_t  fFitCats;
        RooArgSet fVars;
        RooArgSet fCats; 
        RooArgSet fPars;
        RooArgSet fFuncVars; 
        RooArgList fFormulas;            ///<! Cannot be serialized directly
        RooArgList fParameterFormulas;   ///<! Cannot be serialized directly
        RooArgSet fVarsAndCats;
        RooArgSet fParsAndYields;
        RooArgSet fNCParsAndYields;
        RooArgList fYields; 
        RooArgList fPDFs; 
        RooArgList fParameters; 
        RooArgList fConstants;
        RooArgList fConstraints; 
        RooLinkedList fFitOptions; 
        vector<std::unique_ptr<RandomConstrained>> _parConstraints;
        
        RooAbsPdf* fModel = nullptr;     ///<! Owned by Setup manually (deleted in destructor)
 
        RooWorkspace fWS;                ///< Master ROOT Workspace
        TString fVarCut; 
        TString fAddCut; 
        TString fIDBranchName = "UID";
        TString fOutDir;
        TString fDataOnlyCut;
        TList fNeedToDeleteThis;         ///< Legacy garbage collector for custom objects
        
        // --- String Replay Cache ---
        // These vectors cache configuration commands so the Workspace can be
        // identically rebuilt when copying Setup across threads.
        strings_t fVarString;
        strings_t fCatString;
        strings_t fParString;
        strings_t fConstString;
        strings_t fFormString;
        strings_t fAuxVarString;
        strings_t fPDFString;
        strings_t fFuncVarString;        ///<!
        TString fParserPDFString;
        PdfParser* fParserPDFparser = nullptr;

        Bool_t fErrorsSumW2 = kTRUE; 
        Bool_t fErrorsAsym = kFALSE;

        std::map<TString, Bool_t> fConstPars; 
        std::map<TString, Bool_t> fConstPDFPars; 

        std::vector<std::pair<TString, Float_t>> fSpecString;
        std::map<TString, TString> fPDFInWeights;
        TString fYld = "Yld_";

        ClassDefOverride(HS::FIT::Setup, 1);
    };
    
    /**
     * @class RandomConstrained
     * @brief Helper class to draw randomized starting parameter values based on a constraint PDF.
     */
    class RandomConstrained {
    public:
        RandomConstrained(RooAbsPdf* con, RooRealVar* par, Int_t N) :
            _constraint{con}, _par{*par}, _nCache{N}, _entry{N}, _parName{par->GetName()} {}
      
        Double_t get() {
            if (_entry >= _nCache) {
                _entry = 0;
                if (_cache) {
                    _cache->reset();
                    delete _cache;
                }
                _cache = _constraint->generate(_par, _nCache);
            }
            return _cache->get(_entry++)->getRealValue(_parName);
        }
      
    private:
        RooDataSet* _cache = nullptr;
        RooAbsPdf* _constraint = nullptr;
        RooArgSet _par; 
        TString _parName;
        Int_t _nCache = 1;
        Int_t _entry = 1;
    };
    
} // namespace FIT
} // namespace HS

// ////////////////////////////////////////////////////////////////
// ///
// ///Class:               Setup
// ///Description:
// ///           

// #pragma once

// #include "PdfParser.h"

// #include <RooStats/ModelConfig.h>
// #include <RooWorkspace.h>
// #include <RooRealVar.h>
// #include <RooCategory.h>
// #include <RooAddPdf.h>
// #include <RooStats/ModelConfig.h>
// #include <RooGaussian.h>
// #include <RooFormulaVar.h>
// #include <RooDataSet.h>
// #include <TNamed.h>
// #include <TString.h>
// #include <TSystem.h>

// #include <utility>
// #include <vector>
// #include <string>

// //#pragma link C++ class vector<std::pair<TString,Float_t> >+;

// //#pragma link C++ class std::map<TString, TString>+;

// namespace HS{
//   namespace FIT{

//     using std::vector;
//     using std::string;
//     using std::cout;
//     using std::endl;

//     using realvar_ptr = std::shared_ptr<RooRealVar>;
//     using realvars_t    = std::vector<RooRealVar*>;
//     using catvar_ptr  = std::shared_ptr<RooCategory>;
//     using catvars_t     = std::vector<RooCategory*>;
//     using strings_t = std::vector<TString>;

//     //Helper functions
//     RooArgSet MakeArgSet(realvars_t vars);
//     RooArgSet MakeArgSet(catvars_t cats);
//     void SetAllValLimits(RooArgList&,Double_t val,Double_t low=0,Double_t high=0);
//     void ReadFormula(TString forma, strings_t& svars,strings_t& sranges);
//     Bool_t ArgListContainsName(RooArgList& items,TString name);

//     class RandomConstrained;
    
//     class Setup : public TNamed {
      
//     public:
//       Setup(const TString& name);
//       Setup();
//       Setup(const Setup& other);
//       Setup(Setup&&)=delete;
//       ~Setup() override{
// 	if(fModel) delete fModel;
// 	fModel=nullptr;
//       }
//       Setup& operator=(const Setup& other);
//       Setup& operator=(Setup&& other) = delete;//because RooWorkSpace


//       void FactoryPDF(TString opt);
//       void ParserPDF(TString str, PdfParser& parse);
//       void LoadVariable(const TString& opt);
//       void LoadCategory(const TString& opt);
//       void LoadAuxVar(const TString& opt);
//       void LoadFormula(TString formu);
//       void LoadParameter(const TString& opt);
//       void LoadConstant(const TString& opt);
//       void LoadFunctionVar(const TString& opt);
//       void LoadSpeciesPDF(TString opt,Float_t Scale0=1);
//       void TotalPDF();
      
//       const RooWorkspace& WS(){return fWS;}

//       RooArgSet& DataVars();
//       RooArgSet& Cats();
//       RooArgSet& FitVarsAndCats();
//       RooArgSet& ParsAndYields();
//       RooArgSet& NonConstParsAndYields();
//       void CopyRealProperties(RooArgSet& change, const RooArgSet& tothis);
      
//       const realvars_t& AuxVars()const {return fAuxVars;}
	
//       RooAbsPdf* Model()  const {return fModel;}
//       RooAddPdf ExtendModel() const{return RooAddPdf(*dynamic_cast<RooAddPdf*>(fModel));}


      
//       const TString Cut() const {
// 	if(fAddCut.Sizeof()>1&&fVarCut.Sizeof()>1)
// 	  return fAddCut+"&&"+fVarCut;
// 	else 	if(fAddCut.Sizeof()>1) return fAddCut;
// 	else 	if(fVarCut.Sizeof()>1) return fVarCut;
// 	else return TString();
//       }	
//       void AddCut(const TString& cut){if(fAddCut.Sizeof()>1){fAddCut+="&&";}fAddCut+=cut;};
     
//       void SetDataOnlyCut(TString cut) {fDataOnlyCut=std::move(cut);}
//       const TString DataCut() const {
// 	//in case you want to cut on variable not in simulated data
// 	auto cut=Cut();
// 	if(cut.Sizeof()>1&&fDataOnlyCut.Sizeof()>1)
// 	  return cut+"&&"+fDataOnlyCut;
// 	else 	if(cut.Sizeof()>1) return cut;
// 	else 	if(fDataOnlyCut.Sizeof()>1) return fDataOnlyCut;
// 	else return TString();
//       }

 
//       const TString GetIDBranchName() const {return fIDBranchName;}
//       void SetIDBranchName(TString name){fIDBranchName=std::move(name);}
//       const TString GetOutDir() const {
// 	if(fOutDir==TString())
// 	  return "./";
// 	return fOutDir+"/";
//       }
//       void SetOutDir(TString name){
// 	if(!name.BeginsWith("/"))
// 	  name = TString(gSystem->WorkingDirectory())+"/"+name;
// 	fOutDir=name;
// 	gSystem->Exec(Form("mkdir -p %s",fOutDir.Data()));
//       }

//       const realvars_t &FitVars() const {return fFitVars;}
//       const catvars_t &FitCats()const {return fFitCats;}
      
//       RooArgList& Yields()  {return fYields;}
//       RooArgList& Parameters() {return fParameters;}
//       RooArgList& Constants() {return fConstants;}
//       const RooArgList& Formulas() const {return fFormulas;}
//       const RooArgList& ParameterFormulas() const {return fParameterFormulas;}
//       RooArgList& PDFs()   {return fPDFs;}
//       const RooArgList& constPDFs()   const {return fPDFs;}
//       RooArgList& Constraints(){return fConstraints;}

//       Double_t SumOfYields();
      
//       void AddGausConstraint(RooGaussian *pdf){
// 	if(!pdf) return;
// 	fConstraints.add(*(pdf));
//       }
//       void AddFormulaConstraint(RooFormulaVar *formu){
// 	if(!formu) return;
// 	fConstraints.add(*(dynamic_cast<RooAbsArg*>(formu)));
//       }
//       void AddPdfConstraint(RooAbsPdf *pdf){
// 	if(!pdf) return;
// 	fConstraints.add(*pdf);
//       }
      
//       void AddFitOption(const RooCmdArg& cmd){fFitOptions.Add(dynamic_cast<RooCmdArg*>(cmd.Clone()));}
//       RooLinkedList FitOptions(){return fFitOptions;}
      
//       void DefaultFitOptions(){
// 	AddFitOption(RooFit::Warnings(kFALSE));
// 	//AddFitOption(RooFit::Minos(kFALSE));
// 	//AddFitOption(RooFit::Minimizer("Minuit2"));
//       }
//       void RequiredFitOptions(){
// 	if(fFitOptions.find("Save")!=nullptr) return;
	
// 	AddFitOption(RooFit::Save(kTRUE));
// 	if(fErrorsSumW2) AddFitOption(RooFit::SumW2Error(kTRUE));
// 	else if(fErrorsAsym) AddFitOption(RooFit::AsymptoticError(true));
	
//       }

//       void ErrorsSumW2(){
// 	fErrorsSumW2=kTRUE;
// 	fErrorsAsym=kFALSE;
//       }
//       void ErrorsAsymp(){
// 	fErrorsSumW2=kFALSE;
// 	fErrorsAsym=kTRUE;
//       }
//       void ErrorsWrong(){
// 	fErrorsSumW2=kFALSE;
// 	fErrorsAsym=kFALSE;
//       }
      
//       void RandomisePars();
//       void OrganiseConstraints();

//       RooRealVar* GetVar(const TString& name){
// 	return (dynamic_cast<RooRealVar*>(DataVars().find(name)));
//       }
//       RooRealVar* GetPar(const TString& name){
// 	return (dynamic_cast<RooRealVar*>(fParameters.find(name)));
//       }
//       void SetParVal(const TString& par,Double_t val,Bool_t co=kFALSE){
// 	if(dynamic_cast<RooRealVar*>(fParameters.find(par))==nullptr){
// 	  std::cerr<<"Error Setup::SetParVal invalid parameter"<<par<<std::endl;
// 	  exit(0);
// 	}
      
// 	(dynamic_cast<RooRealVar*>(fParameters.find(par)))->setVal(val);
// 	(dynamic_cast<RooRealVar*>(fParameters.find(par)))->setConstant(co);
// 	fConstPars[par]=co;
//       }
//       void SetYldVal(const TString& yld,Double_t val,Bool_t co=kFALSE){
// 	(dynamic_cast<RooRealVar*>(fYields.find(yld)))->setVal(val);
// 	(dynamic_cast<RooRealVar*>(fYields.find(yld)))->setConstant(co);
//       }
//       void SetConstPar(const TString& par,Bool_t co=kTRUE){
// 	(dynamic_cast<RooRealVar*>(fParameters.find(par)))->setConstant(co);
// 	fConstPars[par]=co;
//       }
//       void SetConstPDFPars(const TString& pdf,Bool_t co=kTRUE){
// 	(dynamic_cast<RooAbsPdf*>(fPDFs.find(pdf)))->getParameters(DataVars())->setAttribAll("Constant",co);
// 	fConstPDFPars[pdf]=co;
//       }
//       Bool_t IsParSetConst(const TString& name) const{
// 	//std::cout<<"IsParSetConst "<<name <<" "<<(fConstPars.find(name) != fConstPars.end())<<" "<< fConstPars.at(name) <<endl;
// 	return fConstPars.find(name) != fConstPars.end()? fConstPars.at(name) : kFALSE;}
      
//       void SaveSnapShot(const TString& name){fWS.saveSnapshot(name,RooArgSet(fYields,fParameters),kTRUE);};
//       void LoadSnapShot(const TString& name){fWS.loadSnapshot(name);}

//       RooStats::ModelConfig*  GetModelConfig();
//       TString GetPDFInWeights(const TString& name) {return fPDFInWeights[name];}


//       RooAbsPdf* ComponentsPDF(TString opt);
//     protected:
//       void LoadParameterOnTheFly(const TString& opt);

//     private:
 
//       //note fWS owns all of these vector pointers
//       realvars_t fFitVars; 
//       realvars_t fAuxVars;
//       catvars_t  fFitCats;
//       RooArgSet fVars;
//       RooArgSet fCats; //only categories
//       RooArgSet fPars;
//       RooArgSet fFuncVars; 
//       RooArgList fFormulas;//! CANT WRITE formulas ArgSet!
//       RooArgList fParameterFormulas;//! CANT WRITE formulas ArgSet!
//       RooArgSet fVarsAndCats;
//       RooArgSet fParsAndYields;
//       RooArgSet fNCParsAndYields;//Non constant parameters and yields
//       RooArgList fYields; //species yields
//       RooArgList fPDFs; //species pdfs
//       RooArgList fParameters; //model parameters
//       RooArgList fConstants;//model constants
//       RooArgList fConstraints; //constraints on  parameters
//       RooLinkedList fFitOptions; //
//       vector< std::unique_ptr<RandomConstrained> >_parConstraints;// pdf constraints on parameters
      
//       RooAbsPdf* fModel=nullptr; //!owned by workspace
 
//       RooWorkspace fWS; //
//       TString fVarCut; //
//       TString fAddCut; //
//       TString fIDBranchName="UID";//
//       TString fOutDir;//
//       TString fDataOnlyCut;//
//       TList fNeedToDeleteThis;//
      
//       strings_t fVarString;//
//       strings_t fCatString;//
//       strings_t fParString;//
//       strings_t fConstString;//
//       strings_t fFormString;//
//       strings_t fAuxVarString;//
//       strings_t fPDFString;//
//       strings_t fFuncVarString;//!
//       TString fParserPDFString;//
//       PdfParser* fParserPDFparser=nullptr;//

//       Bool_t	fErrorsSumW2=kTRUE; //default to this
//       Bool_t	fErrorsAsym=kFALSE;

//       std::map<TString,Bool_t> fConstPars; //
//       std::map<TString,Bool_t> fConstPDFPars; //

//       std::vector<std::pair<TString,Float_t> > fSpecString;//
//       std::map<TString,TString> fPDFInWeights;//
//       TString fYld="Yld_";//yield variable prepend

//       ClassDefOverride(HS::FIT::Setup,1);
//     };
    
//     class RandomConstrained {
//       //For randoming parameter values, which are constrained
//     public:

//       RandomConstrained(RooAbsPdf* con,RooRealVar* par,Int_t N):
// 	_constraint{con},_par{*par},_nCache{N},_entry{N},_parName{par->GetName()}
//       {

//       }
      
//      Double_t get(){
// 	if(_entry>=_nCache){
// 	  _entry=0;
// 	  if(_cache){
// 	    _cache->reset();
// 	    delete _cache;
// 	  }
// 	  _cache = _constraint->generate(_par,_nCache);
// 	}
// 	return _cache->get(_entry++)->getRealValue(_parName);
//       }
      
//     private:
//       RooDataSet *_cache=nullptr;
//       RooAbsPdf* _constraint=nullptr;
//       RooArgSet _par; //make an argset from this 1 par as needed for..
//       TString _parName;
      
//       Int_t _nCache=1;
//       Int_t _entry=1;

 
//     };
    
//   }//namespace FIT
// }//namespace HS

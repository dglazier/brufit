////////////////////////////////////////////////////////////////
///
///Class:               Setup
///Description:
///           

#pragma once

#include "PdfParser.h"

#include <RooStats/ModelConfig.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooAddPdf.h>
#include <RooStats/ModelConfig.h>
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooDataSet.h>
#include <TNamed.h>
#include <TString.h>
#include <TSystem.h>

#include <utility>
#include <vector>
#include <string>

//#pragma link C++ class vector<std::pair<TString,Float_t> >+;

//#pragma link C++ class std::map<TString, TString>+;

namespace HS{
  namespace FIT{

    using std::vector;
    using std::string;
    using std::cout;
    using std::endl;

    using realvar_ptr = std::shared_ptr<RooRealVar>;
    using realvars_t    = std::vector<RooRealVar*>;
    using catvar_ptr  = std::shared_ptr<RooCategory>;
    using catvars_t     = std::vector<RooCategory*>;
    using strings_t = std::vector<TString>;

    //Helper functions
    RooArgSet MakeArgSet(realvars_t vars);
    RooArgSet MakeArgSet(catvars_t cats);
    void SetAllValLimits(RooArgList&,Double_t val,Double_t low=0,Double_t high=0);
    void ReadFormula(TString forma, strings_t& svars,strings_t& sranges);
    Bool_t ArgListContainsName(RooArgList& items,TString name);

    class RandomConstrained;
    
    class Setup : public TNamed {
      
    public:
      Setup(const TString& name);
      Setup();
      Setup(const Setup& other);
      Setup(Setup&&)=default;
      ~Setup() override{if(fModel) delete fModel;fModel=nullptr;}
      Setup& operator=(const Setup& other);
      Setup& operator=(Setup&& other) = delete;//because RooWorkSpace


      void FactoryPDF(TString opt);
      void ParserPDF(const TString& str, PdfParser& parse);
      void LoadVariable(const TString& opt);
      void LoadCategory(const TString& opt);
      void LoadAuxVar(const TString& opt);
      void LoadFormula(TString formu);
      void LoadParameter(const TString& opt);
      void LoadConstant(const TString& opt);
      void LoadFunctionVar(const TString& opt);
      void LoadSpeciesPDF(TString opt,Float_t Scale0=1);
      void TotalPDF();
      
      const RooWorkspace& WS(){return fWS;}

      RooArgSet& DataVars();
      RooArgSet& Cats();
      RooArgSet& FitVarsAndCats();
      RooArgSet& ParsAndYields();
      RooArgSet& NonConstParsAndYields();
      const realvars_t& AuxVars()const {return fAuxVars;}
	
      RooAbsPdf* Model()  const {return fModel;}
      RooAddPdf ExtendModel() const{return RooAddPdf(*dynamic_cast<RooAddPdf*>(fModel));}


      
      const TString Cut() const {
	if(fAddCut.Sizeof()>1&&fVarCut.Sizeof()>1)
	  return fAddCut+"&&"+fVarCut;
	else 	if(fAddCut.Sizeof()>1) return fAddCut;
	else 	if(fVarCut.Sizeof()>1) return fVarCut;
	else return TString();
      }	
      void AddCut(const TString& cut){if(fAddCut.Sizeof()>1){fAddCut+="&&";}fAddCut+=cut;};
     
      void SetDataOnlyCut(TString cut) {fDataOnlyCut=std::move(cut);}
      const TString DataCut() const {
	//in case you want to cut on variable not in simulated data
	auto cut=Cut();
	if(cut.Sizeof()>1&&fDataOnlyCut.Sizeof()>1)
	  return cut+"&&"+fDataOnlyCut;
	else 	if(cut.Sizeof()>1) return cut;
	else 	if(fDataOnlyCut.Sizeof()>1) return fDataOnlyCut;
	else return TString();
      }

 
      const TString GetIDBranchName() const {return fIDBranchName;}
      void SetIDBranchName(TString name){fIDBranchName=std::move(name);}
      const TString GetOutDir() const {
	if(fOutDir==TString())
	  return "./";
	return fOutDir+"/";
      }
      void SetOutDir(TString name){
	if(!name.BeginsWith("/"))
	  name = TString(gSystem->Getenv("PWD"))+"/"+name;
	fOutDir=name;
	gSystem->Exec(Form("mkdir -p %s",fOutDir.Data()));
      }

      const realvars_t &FitVars() const {return fFitVars;}
      const catvars_t &FitCats()const {return fFitCats;}
      
      RooArgList& Yields()  {return fYields;}
      RooArgList& Parameters() {return fParameters;}
      RooArgList& Constants() {return fConstants;}
      RooArgList& Formulas() {return fFormulas;}
      RooArgList& ParameterFormulas() {return fParameterFormulas;}
      RooArgList& PDFs()   {return fPDFs;}
      const RooArgList& constPDFs()   const {return fPDFs;}
      RooArgList& Constraints(){return fConstraints;}

      Double_t SumOfYields();
      
      void AddGausConstraint(RooGaussian *pdf){
	if(!pdf) return;
	fConstraints.add(*(pdf));
      }
      void AddFormulaConstraint(RooFormulaVar *formu){
	if(!formu) return;
	fConstraints.add(*(dynamic_cast<RooAbsArg*>(formu)));
      }
      void AddPdfConstraint(RooAbsPdf *pdf){
	if(!pdf) return;
	fConstraints.add(*pdf);
      }
      
      void AddFitOption(const RooCmdArg& cmd){fFitOptions.Add(dynamic_cast<RooCmdArg*>(cmd.Clone()));}
      RooLinkedList FitOptions(){return fFitOptions;}
      
      void DefaultFitOptions(){
	AddFitOption(RooFit::SumW2Error(kTRUE));
	//AddFitOption(RooFit::NumCPU(4));
	AddFitOption(RooFit::Save(kTRUE));
	AddFitOption(RooFit::Warnings(kFALSE));
	//AddFitOption(RooFit::Minos(kFALSE));
	//AddFitOption(RooFit::Minimizer("Minuit2"));
      }
      void RandomisePars();
      void OrganiseConstraints();
      
      void SetParVal(const TString& par,Double_t val,Bool_t co=kFALSE){
	(dynamic_cast<RooRealVar*>(fParameters.find(par)))->setVal(val);
	(dynamic_cast<RooRealVar*>(fParameters.find(par)))->setConstant(co);
	fConstPars[par]=co;
      }
      void SetYldVal(const TString& yld,Double_t val,Bool_t co=kFALSE){
	(dynamic_cast<RooRealVar*>(fYields.find(yld)))->setVal(val);
	(dynamic_cast<RooRealVar*>(fYields.find(yld)))->setConstant(co);
      }
      void SetConstPar(const TString& par,Bool_t co=kTRUE){
	(dynamic_cast<RooRealVar*>(fParameters.find(par)))->setConstant(co);
	fConstPars[par]=co;
      }
      void SetConstPDFPars(const TString& pdf,Bool_t co=kTRUE){
	(dynamic_cast<RooAbsPdf*>(fPDFs.find(pdf)))->getParameters(DataVars())->setAttribAll("Constant",co);
	fConstPDFPars[pdf]=co;
      }
      Bool_t IsParSetConst(const TString& name) const{
	//std::cout<<"IsParSetConst "<<name <<" "<<(fConstPars.find(name) != fConstPars.end())<<" "<< fConstPars.at(name) <<endl;
	return fConstPars.find(name) != fConstPars.end()? fConstPars.at(name) : kFALSE;}
      
      void SaveSnapShot(const TString& name){fWS.saveSnapshot(name,RooArgSet(fYields,fParameters),kTRUE);};
      void LoadSnapShot(const TString& name){fWS.loadSnapshot(name);}

      RooStats::ModelConfig*  GetModelConfig();
      TString GetPDFInWeights(const TString& name) {return fPDFInWeights[name];}


      RooAbsPdf* ComponentsPDF(TString opt);
    protected:
      void LoadParameterOnTheFly(const TString& opt);

    private:
 
      //note fWS owns all of these vector pointers
      realvars_t fFitVars;      
      realvars_t fAuxVars;      
      catvars_t  fFitCats;
      RooArgSet fVars;
      RooArgSet fCats; //only categories
      RooArgSet fPars;//!
      RooArgSet fFuncVars;//!
      RooArgList fFormulas;//! CANT WRITE formulas ArgSet!
      RooArgList fParameterFormulas;//! CANT WRITE formulas ArgSet!
      RooArgSet fVarsAndCats;
      RooArgSet fParsAndYields;
      RooArgSet fNCParsAndYields; //Non constant parameters and yields
      RooArgList fYields;//species yields
      RooArgList fPDFs;//species pdfs
      RooArgList fParameters;//model parameters
      RooArgList fConstants;//model constants
      RooArgList fConstraints;//constraints on  parameters
      RooLinkedList fFitOptions;//
      vector< std::unique_ptr<RandomConstrained> >_parConstraints;//! pdf constraints on parameters
      
      RooAbsPdf* fModel=nullptr; //!owned by workspace
 
      RooWorkspace fWS;
      TString fVarCut;
      TString fAddCut;
      TString fIDBranchName="UID";
      TString fOutDir;
      TString fDataOnlyCut;

      strings_t fVarString;
      strings_t fCatString;
      strings_t fParString;
      strings_t fConstString;
      strings_t fFormString;
      strings_t fAuxVarString;
      strings_t fPDFString;
      strings_t fFuncVarString;

      std::map<TString,Bool_t> fConstPars;
      std::map<TString,Bool_t> fConstPDFPars;

      std::vector<std::pair<TString,Float_t> > fSpecString;
      std::map<TString,TString> fPDFInWeights;
      TString fYld="Yld_";//yield variable prepend

      ClassDefOverride(HS::FIT::Setup,1);
    };
    
    class RandomConstrained {

    public:

      RandomConstrained(RooAbsPdf* con,RooRealVar* par,Int_t N):
	_constraint{con},_par{*par},_nCache{N},_entry{N},_parName{par->GetName()}
      {

      }
      
     Double_t get(){
	if(_entry>=_nCache){
	  _entry=0;
	  if(_cache){
	    _cache->reset();
	    delete _cache;
	  }
	  _cache = _constraint->generate(_par,_nCache);
	}
	return _cache->get(_entry++)->getRealValue(_parName);
      }
      
    private:
      RooDataSet *_cache=nullptr;
      RooAbsPdf* _constraint=nullptr;
      RooArgSet _par; //make an argset from this 1 par as needed for..
      TString _parName;
      
      Int_t _nCache=1;
      Int_t _entry=1;

 
    };
    
  }//namespace FIT
}//namespace HS


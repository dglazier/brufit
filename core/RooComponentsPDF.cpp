#include "RooComponentsPDF.h" 
#include <RooAbsReal.h> 
#include <RooAbsArg.h>
#include <RooAbsCategory.h> 
#include <cmath> 
#include "TMath.h" 

namespace HS{
  namespace FIT{
    
    RooComponentsPDF::RooComponentsPDF(const char *name, const char *title,Double_t base,const RooArgList& obsList,const vector<RooArgList> compList)
      :  HS::FIT::RooHSEventsPDF(name,title), fBaseLine(base)
      
    {

      fNObs=0;
      fNCats=0;
      for(Int_t i=0;i<obsList.getSize();i++){
	//case real variable
	if(dynamic_cast<RooRealVar*>(&obsList[i])){
	  std::unique_ptr<RooRealProxy> tempR{new RooRealProxy(obsList[i].GetName(),obsList[i].GetName(),this,dynamic_cast<RooAbsReal&>(obsList[i]))};
	//cout<<"adding obs "<<obsList[i].GetName()<<" "<<dynamic_cast<RooRealVar*>(&obsList[i])<<" "<<dynamic_cast<RooCategory*>(&obsList[i])<<endl;
	  fNObs++;
	  fActualObs.add(dynamic_cast<RooAbsReal&>(obsList[i]));
	  fObservables.push_back(std::move(tempR));
	  continue;
	}
	//case category
	if(dynamic_cast<RooCategory*>(&obsList[i])){
	  std::unique_ptr<RooCategoryProxy> tempC{new RooCategoryProxy(obsList[i].GetName(),obsList[i].GetName(),this,dynamic_cast<RooAbsCategory&>(obsList[i]))};
	//cout<<"adding cat "<<obsList[i].GetName()<<endl;
	  fNCats++;
	  fActualCats.add(dynamic_cast<RooAbsCategory&>(obsList[i]));
	  fCategories.push_back(std::move(tempC));
	  continue;
	}
	
	
      }

      fNComps=compList.size();
      //Arrange the components into a std::vec in order
      //This seems to speed up calcuations
      for(auto& comp: compList){
	vecUPtrReal vterms;
	for(Int_t i=0;i<comp.getSize();i++){
	  //add this component term to the list of all components
	  fActualComps.add(dynamic_cast<RooAbsReal&>(comp[i]));
	  //now create a proxy with it
	  std::unique_ptr<RooRealProxy> temp{new RooRealProxy(comp[i].GetName(),comp[i].GetName(),this,dynamic_cast<RooAbsReal&>(comp[i]))};
	  vterms.push_back(std::move(temp));
	}
	fComponents.push_back(std::move(vterms));
      }
      MakeSets();

   
      initIntegrator();

    }

    RooComponentsPDF::RooComponentsPDF(const RooComponentsPDF& other, const char* name) :
      HS::FIT::RooHSEventsPDF(other,name),
      fBaseLine(other.fBaseLine),
      fActualComps("AllComponents",this,other.fActualComps),
      fActualCats("AllCategories",this,other.fActualCats),
      fActualObs("AllObservables",this,other.fActualObs)
    {

      fWeightedBaseLine=other.fWeightedBaseLine;

      Int_t counter=0;
      
      for(const auto & fObservable : other.fObservables){
	unique_ptr<RooRealProxy> temp{new RooRealProxy(fObservable->GetName(),this,*fObservable)};
	fObservables.push_back(std::move(temp));
      }
      for(const auto & fCategorie : other.fCategories){
	unique_ptr<RooCategoryProxy> temp{new RooCategoryProxy(fCategorie->GetName(),this,*fCategorie)};
	fCategories.push_back(std::move(temp));
      }
      
      for(auto& comp: other.fComponents){
	vecUPtrReal vterms;
	for(const auto & i : comp){
	  unique_ptr<RooRealProxy> temp{new RooRealProxy(i->GetName(),this,*i)};
	  vterms.push_back(std::move(temp));
	}
	fComponents.push_back(std::move(vterms));
      }

      fNObs=other.fNObs;
      fNCats=other.fNCats;
      fNComps=other.fNComps;
      
      MakeSets();
    
      
      initIntegrator();

        //get the original cached integrals      
      fCacheCompDepIntegral=other.fCacheCompDepIntegral;

    } 
    void RooComponentsPDF::MakeSets(){
      //roorealvars
      for(auto &obs: fObservables){

      	fProxSet.push_back(obs.get());
	//for linking tree events to integral calculation
	fIntegrateObs.push_back(dynamic_cast<RooRealVar*>(fIntegrateSet.addClone(obs->arg())));
      }
      //roocategories
      for(auto &obs: fCategories){
 	
     	fCatSet.push_back(obs.get());
	//for linking tree events to integral calculation
	fIntegrateCats.push_back(dynamic_cast<RooCategory*>(fIntegrateSet.addClone(obs->arg())));
      }
 
      
      for(auto &comp: fComponents)
	for(auto &term: comp){
	  fParSet.push_back(term.get()); //par set just used to check change
	  //get the RooRealVar for this Proxy term
	  auto  argTerm=fActualComps.find(term->GetName());
	  //get its variables, if any
	  auto vars=argTerm->getVariables();
	  TIter iter=vars->createIterator();
	  while(auto* arg=dynamic_cast<RooAbsArg*>(iter())){
	    //if new variable and not observable
	    //include it as parameter
	    if(!fActualObs.contains(*arg)&&!fActualCats.contains(*arg)&&!fParameters.contains(*arg)){
	      fParameters.add(*arg);
	    }
	  }
	}
      InitSets();
    }
    Int_t RooComponentsPDF::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
    {	
      Info("RooHSEventsPDF::getGenerator","Looking for generator");
      if(!fEvTree) return 0; //no MC events to generate from
      //case generate all variables
      if (matchArgs(directVars,generateVars,VarSet(0))) return 1 ;
      return 0;

    }



    Double_t RooComponentsPDF::evaluateData() const 
    {
      //cout<<"RooComponentsPDF::evaluateData() "<<endl;
      Double_t val=fBaseLine;
       for(auto &comp: fComponents){
	Double_t product=1;
	for(auto &term: comp){
	  // cout<<"term "<<term->GetName()<<" "<< *term.get()<<endl;
	  product*= *term.get(); //take the product of all the terms for this component
	}
	//	cout<<"product "<<product<<endl;
	val+=product; //add them to total
      }

       return val;
    }
    
    void RooComponentsPDF::RedirectServersToPdf(){
      cout<<"                 RooComponentsPDF::RedirectServersToPdf()"<<endl;
      //point the terms to the integral events rather than data events
      for(UInt_t icomp=0;icomp<fNComps;icomp++){
	for(const auto &term:fDependentTermProxy[icomp]){
	  auto unconstTerm=const_cast<RooAbsReal*>(&term->arg());
	  unconstTerm->recursiveRedirectServers(fIntegrateSet);
	}	
      }	
    }
    void RooComponentsPDF::RedirectServersToData(){
      //point the terms back to the data events rather than integral events
      for(UInt_t icomp=0;icomp<fNComps;icomp++){
	for(const auto &term:fDependentTermProxy[icomp]){
	  auto unconstTerm=const_cast<RooAbsReal*>(&term->arg());
	  unconstTerm->recursiveRedirectServers(fActualObs);
	}	
      }	
      
    }

    void RooComponentsPDF::HistIntegrals(const char* rangeName) const{
   //point the terms to the integral events rather than data events
      for(UInt_t icomp=0;icomp<fNComps;icomp++){
	for(const auto &term:fDependentTermProxy[icomp]){
	  auto unconstTerm=const_cast<RooAbsReal*>(&term->arg());
	  unconstTerm->recursiveRedirectServers(fIntegrateSet);
	}
      }
       RooHSEventsPDF::HistIntegrals(rangeName);

      //point the terms back to the data events rather than integral events
      for(UInt_t icomp=0;icomp<fNComps;icomp++){
	for(const auto &term:fDependentTermProxy[icomp]){
	  auto unconstTerm=const_cast<RooAbsReal*>(&term->arg());
	  unconstTerm->recursiveRedirectServers(fActualObs);
	}
      }

    }
    Double_t RooComponentsPDF::evaluateMC(const vector<Float_t> *vars,const  vector<Int_t> *cats) const
    {
      // cout<< "RooComponentsPDF::evaluateMC"<<endl;
      //read in observable value for this event
      for(Int_t ii=0;ii<fNvars;ii++){
	//cout<<vars->at(fTreeEntry*fNvars+ii)<<" ";
	fIntegrateObs[ii]->setVal(vars->at(fTreeEntry*fNvars+ii));
      }
      for(Int_t ii=0;ii<fNcats;ii++){
	fIntegrateCats[ii]->setIndex(cats->at(fTreeEntry*fNcats+ii));
      }
      //cout<<endl;
      return evaluateData();
    }
    Bool_t RooComponentsPDF::isDirectGenSafe(const RooAbsArg& arg) const {
      if(fActualObs.find(arg.GetName())) return kTRUE;
      if(fActualCats.find(arg.GetName())) return kTRUE;
      return kFALSE;
    }

    void RooComponentsPDF::initGenerator(Int_t code)
    {
      RedirectServersToPdf();
      RooHSEventsPDF::initGenerator(code);
    }
    void RooComponentsPDF::initIntegrator()
    {
      //Each Component is arranged in terms which are
      //  Observable independent
      //        Can just use current value without looping over events
      //  Observable dependent, parameter independent
      //        Can just loop over event once and cache result
      //  Observable dependent, parameter dependenent
      //        Must loop over events whenever a parameter value changes
      //The final integral for each component is the product of these
      //three super terms
      
        //Initiliase components cache
      //fDependentTerm.resize(fNComps);
      fDependentTermProxy.resize(fNComps);
      fDependentTermParams.resize(fNComps);
      fPrevParVals.resize(fNComps);
      fCacheCompDepIntegral.resize(fNComps);
      
      //fIndependentTerm.resize(fNComps);
      fIndependentTermProxy.resize(fNComps);
      //fParameterTerm.resize(fNComps);
      UInt_t icomp=0;
      for(auto &comp: fComponents){
	fCacheCompDepIntegral[icomp]=1;
	UInt_t iterm=0;
	Double_t product=1;
	for(auto &term: comp){
	  //Identify which terms are dependent on fit observables and cats(VarSet)
	  auto  arg=fActualComps.find(term->GetName());
	  auto deps=arg->getDependents(VarSet(0));
	
	  if(deps->getSize()){
	   
	    fDependentTermProxy[icomp].push_back(term.get());
	    //Identify which terms are dependent on fit parameters (ParSet)
	    auto parDeps=arg->getDependents(fParameters);
	    if(parDeps->getSize()){
	     
	      TIter iter=parDeps->createIterator();
	      while(auto* arg=dynamic_cast<RooAbsArg*>(iter())){
		auto *rarg=dynamic_cast<RooRealVar*>(arg);	    
		if(!(vecContains(rarg,fDependentTermParams[icomp]))){
		  fDependentTermParams[icomp].push_back(rarg);
		  Double_t initf=-1E6;
		  fPrevParVals[icomp].push_back(initf);
		}
	      }
	    }
	    else{
	    }
	    
	  }
	  else{

	    fIndependentTermProxy[icomp].push_back(term.get());
	  }
	  
	  iterm++;
	}
	
	icomp++;
	
      }	  
     }
    
    
    
    
    
    Double_t RooComponentsPDF::analyticalIntegral(Int_t code,const char* rangeName) const
    {
       if(code!=1) return RooHSEventsPDF::analyticalIntegral(code,rangeName);
       if(code==1&&fForceConstInt&&!fEvTree) {fLast[0]=1;return fLast[0];}

       //Check baseline caclulated
      if(fWeightedBaseLine==0&&fBaseLine!=0&&fUseEvWeights)
	CalcWeightedBaseLine(rangeName);
      else
	fWeightedBaseLine=fBaseLine;

      //Check which dependent terms need recalculation
      //This will be 1) if they are dependent on parameters
      //           2) one or more of the parameters have changed
      Bool_t needRecalc=kFALSE;
      fRecalcComponent.clear();
      for(UInt_t icomp=0;icomp<fNComps;icomp++){
	  if(fDependentTermProxy[icomp].size()) {
	    
	    UInt_t ipar=0;
	    for(auto par:fDependentTermParams[icomp]){
	      Double_t previous=fPrevParVals[icomp][ipar];
	      Double_t pval=par->getVal();
	     
	      if(pval!=previous){ //trigger recalc
		needRecalc=kTRUE;
		fPrevParVals[icomp][ipar]=pval;
		if(!vecContains(icomp,fRecalcComponent)) fRecalcComponent.push_back(icomp);
	      }
	      //Store parameter values to check for a change
	      ipar++;
	    }
	  }
	  
      }
      ///////////////////////////////
      if(needRecalc)RecalcComponentIntegrals(code,rangeName);

      Double_t integral=fWeightedBaseLine;
    
      for(UInt_t icomp=0;icomp<fNComps;icomp++)
	integral+=componentIntegral(icomp);
      // cout<<"                   INTEGRAL "<<integral<<endl;
      //   exit(0);
      return integral;
    }
    
    void RooComponentsPDF::CalcWeightedBaseLine(const char* rangeName) const{
     Long64_t ilow,ihigh=0;
      SetLowHighVals(ilow,ihigh);
         //point the terms to the integral events rather than data events
      for(const auto& icomp:fRecalcComponent){
	for(const auto &term:fDependentTermProxy[icomp]){
	  auto unconstTerm=const_cast<RooAbsReal*>(&term->arg());
	  unconstTerm->recursiveRedirectServers(fIntegrateSet);
	}
      }
      fWeightedBaseLine=0;
      //Loop over events and recalcaulte partial integrals
      //that depend on parameters that have changed
      Long64_t accepted=0;   
      for(Long64_t ie=ilow;ie<ihigh;ie++){
	fTreeEntry=ie;
	if(!CheckRange(TString(rangeName).Data())) continue;
	accepted++;
	fWeightedBaseLine+=GetIntegralWeight(ie);
      }
      fWeightedBaseLine/=accepted;//normalise to number of events to prevent huge integrals
      cout<<"RooComponentsPDF::CalcWeightedBaseLine "<<fWeightedBaseLine<<endl;
    }
    void RooComponentsPDF::RecalcComponentIntegrals(Int_t code,const char* rangeName) const{
      Long64_t ilow,ihigh=0;
      SetLowHighVals(ilow,ihigh);
         //point the terms to the integral events rather than data events
      for(const auto& icomp:fRecalcComponent){
	for(const auto &term:fDependentTermProxy[icomp]){
	  auto unconstTerm=const_cast<RooAbsReal*>(&term->arg());
	  unconstTerm->recursiveRedirectServers(fIntegrateSet);
	}
      }

      //Loop over events and recalcaulte partial integrals
      //that depend on parameters that have changed
      Long64_t accepted=0;
        for(Long64_t ie=ilow;ie<ihigh;ie++){
	fTreeEntry=ie;
	if(!CheckRange(TString(rangeName).Data())) continue;
	accepted++;
	//read in observable value for this event
	for(Int_t ii=0;ii<fNvars;ii++)
	  fIntegrateObs[ii]->setVal(fvecReal[fTreeEntry*fNvars+ii]);
	for(Int_t ii=0;ii<fNcats;ii++)
	  fIntegrateCats[ii]->setIndex(fvecCat[fTreeEntry*fNcats+ii]);
	//calculate the partial integrals
	for(const auto& icomp:fRecalcComponent){
	  Double_t product=1;
	  for(const auto &term:fDependentTermProxy[icomp]){
	    product*= *term;
	  }
	  product*=GetIntegralWeight(ie);
	  
	  fCacheCompDepIntegral[icomp]+=product;
	}
      }
   
      //Normalise to number of events
      for(const auto& icomp:fRecalcComponent){
	fCacheCompDepIntegral[icomp]=fCacheCompDepIntegral[icomp]/accepted;
      }
     //point the terms back to the data events rather than integral events
      for(const auto& icomp:fRecalcComponent){
	for(const auto &term:fDependentTermProxy[icomp]){
	  auto unconstTerm=const_cast<RooAbsReal*>(&term->arg());
	  unconstTerm->recursiveRedirectServers(fActualObs);
	}
      }

      
    }
    Double_t RooComponentsPDF::componentIntegral(Int_t icomp) const{
      //calculate integral of this component
      //First take product of terms independent of observables;
      Double_t product=1;
      product*=fCacheCompDepIntegral[icomp];
      int counter=0;
      for(auto& term:fIndependentTermProxy[icomp]){
	product*= *term;
      }
       return product; 
    }

    Bool_t RooComponentsPDF::SetEvTree(TTree* tree,TString cut,TTree* MCGenTree){
      auto val = RooHSEventsPDF::SetEvTree(tree,cut,MCGenTree);

      //Caclulate current value of component integrals
      for(UInt_t icomp=0;icomp<fNComps;icomp++)
	fRecalcComponent.push_back(icomp);
      RecalcComponentIntegrals(0,"");
      return val;
    }
  }
}

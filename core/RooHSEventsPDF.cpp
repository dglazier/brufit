#include "RooHSEventsPDF.h"
 
#include <RooRealVar.h>
#include <RooCategory.h> 
#include <RooRandom.h>
#include <RooDataSet.h> 
#include <RooMsgService.h> 
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooPullVar.h> 
#include <TMath.h>
#include <TObjArray.h> 
#include <TCanvas.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TSystem.h>
#include <TEntryList.h>
#include <algorithm> 
#include <random>


namespace HS{
  namespace FIT{

    Bool_t RooHSEventsPDF::RooHSEventsPDF_IsPlotting=kFALSE;

    void RooHSEventsPDF::SetIsPlotting(Bool_t is){
      RooHSEventsPDF_IsPlotting=is;	
    }

    RooHSEventsPDF::RooHSEventsPDF(const RooHSEventsPDF& other, const char* name) :  RooAbsPdf(other,name) 
    {
      // cout<<"RooHSEventsPDF::RooHSEventsPDF "<<GetName()<<other.fNTreeEntries<< " "<<other.fvecReal.size()<<" is cloen "<<other.fIsClone<<" "<<&other<<endl;
      fIsClone=kTRUE;
      fParent=const_cast<RooHSEventsPDF*>(&other);

      fvecReal=other.fvecReal;
      fvecCat=other.fvecCat;
      fvecRealGen=other.fvecRealGen;
      fvecCatGen=other.fvecCatGen;
      fNTreeEntries=other.fNTreeEntries;
      fTreeEntryNumber=other.fTreeEntryNumber;
      
      fAssertPosDataReal=other.fAssertPosDataReal;
      fAssertPosDataCats=other.fAssertPosDataCats;
      fNapd =other.fNapd;
      
      if(other.fEvTree)fEvTree=other.fEvTree;
      fNInt=other.fNInt;
      fGeni=other.fGeni;
      fTruthPrefix=other.fTruthPrefix;
      //if(other.fEntryList)fEntryList=(TEntryList*)other.fEntryList->Clone();
      fForceConstInt=other.fForceConstInt;
      fForceNumInt=other.fForceNumInt;
      fConstInt=other.fConstInt;
      fCheckInt=other.fCheckInt;
      fUseWeightsGen=other.fUseWeightsGen;
      fCut=other.fCut;
      fInWeightCut=other.fInWeightCut;
      fIsValid=other.fIsValid;
      fUseEvWeights=other.fUseEvWeights;
      // if(other.fUseSamplingIntegral){
      // 	fUseSamplingIntegral=other.fUseSamplingIntegral;
      // 	fIntegralPDF=other.fIntegralPDF;
      // }
      fWgtsConf=other.fWgtsConf;
									     
      fEvWeights=other.fEvWeights;
      fHistIntegrals=other.fHistIntegrals;
      fMaxValue=other.fMaxValue;
      fIntRangeLow=other.fIntRangeLow;
      fIntRangeHigh=other.fIntRangeHigh;
    }
    RooHSEventsPDF::~RooHSEventsPDF(){

      /*   std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"<<"delete RooHSEventsPDF::~RooHSEventsPDF() "
	<<fIsClone<<" "<<fParent<<" "<<fEntryList<<" "
	<<fMaxValue<<" "
	<< endl;
      */
      //If I made a generator entry list save it
      //so it can be used to filter the original tree
      //object the entrylist if I want to use it!
      if(fEntryList!=nullptr){
	if(fEntryList->GetN()!=0){//has this clone used generator?
	  TFile entryFile(TString("entryFile_")+GetName()+".root","recreate");
	  fEntryList->Write();
	}
      }
      if(fEntryList) delete fEntryList;
      if(fLast) delete fLast;
      
   
      if(fWeights){
	fWeights->Save();
	delete fWeights;
      }
      if(fInWeights){
	delete fInWeights;
      }
      for(auto & i : fVarSet)
      	delete i;
      fVarSet.clear();
      }

    void RooHSEventsPDF::InitSets(){
      fNpars=fParSet.size();
      fNvars=fProxSet.size();
      fNcats=fCatSet.size();
      fLastLength=fNpars+1;
      fLast=new Double_t[fNpars+1]; //Number of fit parameters
      for(Int_t i=0;i<fNpars+1;i++)
	fLast[i]=100.;
    }
    RooArgSet RooHSEventsPDF::VarSet(Int_t iset) const{
      RooArgSet aset(Form("VarSet_%d",iset));
      if(iset==0){
	for(auto i : fProxSet){
	  aset.add(i->arg());
	}
	for(auto i : fCatSet){
	  aset.add(i->arg());
	}
      }
      else{
	for(auto j : fProxSet){//add if not proxy being removed
	  if(fProxSet[iset-1]->GetName()!=j->GetName()) aset.add(j->arg());
	}
	for(auto i : fCatSet){ //Just add all categories
	  aset.add(i->arg());
	}
      }
      return aset;
    }
    Int_t RooHSEventsPDF::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
    {
      Info("RooHSEventsPDF::getGenerator","Looking for generator");
      if(!fEvTree) return 0; //no MC events to generate from
      //case generate all variables
      if (matchArgs(directVars,generateVars,VarSet(0))) return 1 ;
      return 0;

    }
    void RooHSEventsPDF::initIntegrator()
    {

    }
    void RooHSEventsPDF::initGenerator(Int_t code)
    {
       Info("RooHSEventsPDF::initGenerator","Going to generate starting from %lld with weights %d",fGeni,(Int_t)fUseWeightsGen);
      //Calculate the max value for accept reject purposes
      //Note we use parent to make sure this is only done once
      //RooFit creates a clone PDF instance each time it wants to generate
      if(fParent->GetMaxValue()==0||fParent->CheckChange()){	
	Double_t value=0;		
	if(code==1){	
	  //Brute force find maximum value
	  fMaxValue=0;
	  for(Int_t i=0;i<fNTreeEntries;i++){
	    fTreeEntry=i;
	    value=evaluateMC(&fvecRealGen,&fvecCatGen);
	    if(value<0){ std::cout<<" RooHSEventsPDF::initGenerator -ve intensity !!! "<<value<<" while max was "<<fMaxValue<<std::endl;exit(0);}
	    
	    //	    std::cout<<"RooHSEventsPDF::initGenerator "<<value<<" "<<i<<std::endl;
	    if(value>fMaxValue)fMaxValue=value*1.01;//make it a little larger
	  }
 
	  fParent->SetMaxValue(fMaxValue);
	}
      }		
      //construct entry list so can reproduce full tree branches,
      //not jist those loaded as variables
      if(!fEntryList){
	fEntryList=new TEntryList("GenEvents","GenEvents",fEvTree);
      }
      else{
	fEntryList->Reset();
	fEntryList->SetTree(fEvTree);
      }
      if(fUseWeightsGen){
	fWeights=new Weights("genWeights");
	fWeights->SetSpecies(GetName());
	fWeights->SetFile(TString(GetName())+"Weights.root");
      }

        
      Info("RooHSEventsPDF::initGenerator","Max value %lf",fMaxValue);
    }
    void RooHSEventsPDF::generateEvent(Int_t code){
      // Info("RooHSEventsPDF::generateEvent","Going to generate starting from %lld with ",fGeni);
        Double_t value=0;
      if(!fUseWeightsGen){
	while(fGeni<fNTreeEntries){
	  fTreeEntry=IncrementGeni();
	  value=evaluateMC(&fvecRealGen,&fvecCatGen); //evaluate true values
	  if(value>fMaxValue*RooRandom::uniform()){//accept
	    for(Int_t i=0;i<fNvars;i++)
	      (*(fProxSet[i]))=fvecReal[fTreeEntry*fNvars+i]; //write reconstructed
	    for(Int_t i=0;i<fNcats;i++)
	      (*(fCatSet[i]))=fvecCat[fTreeEntry*fNcats+i];

	    //Add actual entry number from original tree
	    //this can then be used to filter original tree
	    //with all branches
	    fEntryList->Enter(fTreeEntryNumber[fTreeEntry]);
	    return;
	  }
	}
      }
      else{
	//using weights
	while(fGeni<fNTreeEntries){
	  fParent->SetGeni(fGeni);
	  //fEvTree->GetEntry(fGeni++);
	  fTreeEntry=fGeni++;
	  if(!CheckRange("")) continue;
	  value=evaluateMC(&fvecRealGen,&fvecCatGen);
	  for(Int_t i=0;i<fNvars;i++)
	    (*(fProxSet[i]))=fvecReal[fTreeEntry*fNvars+i];
	  for(Int_t i=0;i<fNcats;i++)
	    (*(fCatSet[i]))=fvecCat[fTreeEntry*fNcats+i];
	  fWeights->FillWeight(fGeni-1,value); 
	  fEntryList->Enter(fGeni-1);
	  return;
	}
      }
      Fatal("RooHSEventsPDF::generateEvent","Ran out of events at %lld",fGeni);
      //Used up all the events in the tree!
    }
    Int_t RooHSEventsPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,const char* rangeName) const
    {
    
      if(fForceNumInt) return 0; //might be good to check numerical integral sometimes
      if(!fEvTree&&!fForceConstInt) return 0; //no MC events to integrate over

      if(fProxSet.size()==1&&fCatSet.size()==0){//special case 1 variable

	if (matchArgs(allVars,analVars,VarSet(0))){
	  //if now plotting create histograms
	  if(RooHSEventsPDF_IsPlotting&&fHistIntegrals.size()==0)
	    HistIntegrals(rangeName);

	  return 1;
	}
      }
      else{//For variables
	for(UInt_t i=0;i<1+fProxSet.size();i++){
	  if(!fEvTree&&fForceConstInt&&i==0) {return 1;} //no tree, but const int
	  else if(!fEvTree) return 0; //no const integral for projections
	  if (matchArgs(allVars,analVars,VarSet(i))) {return i+1 ;}
	}
      }
      //Note not implemented for cats
      return 0;
    }
    //new function for uncertainty on MC integral
    //Need to copy from ComponentsPDF which has been optimised
    /*******************************************/
    Double_t RooHSEventsPDF::analyticalIntegralForSampling(const char* rangeName) const
    {
  
      Double_t integral=0;
      Long64_t accepted=0;
      Long64_t all=0;
      Long64_t ilow=0;
      Long64_t ihigh=0;

      //Set range of events to integrate over
      SetLowHighVals(ilow,ihigh); 
      //Loop over events and add to integral
      if(CheckChange()){
	std::vector<double> values(ihigh-ilow);
	
	for(Long64_t ie=ilow;ie<ihigh;ie++){
	  fTreeEntry=ie;
	  if(!CheckRange(rangeName)){
	    values[all]=0;
	    ++all;
	    continue;
	  }
	  values[all]=evaluateMC(&fvecReal,&fvecCat)*GetIntegralWeight(ie);
	  integral+=values[all];
	  ++accepted; //actual entries to count
	  ++all; //just for array sizing
	}
		
	//normalise integral by number of events accepted
	integral/=accepted;
	
	double sum_of_diffs = 0.;
	//variance of MC integral/Volume =  variance(Sum(f(x)-<f>))^2 / N
	//https://en.wikipedia.org/wiki/Monte_Carlo_integration
	std::for_each(values.begin(), values.end(),
		      [&sum_of_diffs,&integral] (double n) {
			double term = (n-integral);
			sum_of_diffs += term*term;});
	fSigmaIntegral = TMath::Sqrt(sum_of_diffs/accepted);
	
	fLast[0]= integral;
      }      

      return sampleIntegral(fLast[0],fSigmaIntegral);
    }
    
    Double_t RooHSEventsPDF::analyticalIntegral(Int_t code,const char* rangeName) const
    {
       if(code==1&&fForceConstInt&&!fEvTree) {fLast[0]=1;return fLast[0];}
       Long64_t NEv=0;
       
 
      //In case changed for generation
  
      Double_t integral=0.;
 
      //only recalculate if a par changes when all variables included(ie code=1)
      if(code==1)
	if(!CheckChange()) return fLast[0];
      if(code==1){
	auto check= AssertPositivePDF();
	if(check==kFALSE) return fLast[0]=0;
	//cout<<"RooHSEventsPDF::analyticalIntegral was it OK "<<check<<endl;

	   //	if(fUseSamplingIntegral==kFALSE){
	  Long64_t accepted=0;
	  Long64_t ilow=0;
	  Long64_t ihigh=0;

	  //Set range of events to integrate over
	  SetLowHighVals(ilow,ihigh); 
	  //Loop over events and add to integral
	  for(Long64_t ie=ilow;ie<ihigh;ie++){
	    fTreeEntry=ie;
	    if(!CheckRange(rangeName)) continue;
	    accepted++;
	    integral+=evaluateMC(&fvecReal,&fvecCat)*GetIntegralWeight(ie);
	  }
	
	  //normalise integral by number of events accepted
	  integral/=accepted;
	  
	//Needs fixed to componentsPDF method
	//else{//use sampled method
	  //	  integral = analyticalIntegralForSampling(rangeName);
	//}
      }
      else{
	if(fHistIntegrals.size()==0)
	  HistIntegrals(rangeName);
	
	Int_t vindex=code-2;
	Double_t vval=*(fProxSet[vindex]);
	integral=fHistIntegrals[vindex].Interpolate(vval);
	if(integral<0) integral=0;
	//std::cout<<"DEBUG RooHSEventsPDF::integral "<<integral<<std::endl;
	return integral;
      }
      // Set Last[0] so we can just return that if no parameter changes
      fLast[0]=integral;
      //std::cout<<"DEBUG RooHSEventsPDF::integral "<<fLast[0]<<std::endl;
      return fLast[0];
    }

    Double_t RooHSEventsPDF::unnormalisedIntegral(Int_t code,const char* rangeName) const{
      Double_t integral=0;
      Double_t nev=0;
      Double_t nMC=0;
      if(code==1){
	for(Long64_t ie=0;ie<fNTreeEntries;ie++){
	  fTreeEntry=ie;
	  if(!CheckRange(TString(rangeName).Data())) continue;
          // cout << "evaluateMC: " << evaluateMC(&fvecReal,&fvecCat) << " GetIntegralWeight: " << GetIntegralWeight(ie) << endl;

	  integral+=evaluateMC(&fvecReal,&fvecCat)*GetIntegralWeight(ie);
	  nev++;
	}
	cout << "RooHSEventsPDF::unnormalisedIntegral #MC=" << nev << ". Integral=" << integral << endl;
      }
      else if(code==2 && fHasMCGenTree){
	for(Long64_t ie=0;ie<fNMCGenTreeEntries;ie++){
	  fTreeEntry=ie;
	  integral+=evaluateMC(&fvecRealMCGen,&fvecCatMCGen);
	  nMC++;
	}
	cout << "RooHSEventsPDF::unnormalisedIntegral #GEN= " << nMC << ". Integral=" << integral  << endl;
      }
      else{
	return 0;
      }
      return integral;
    }
    
    void RooHSEventsPDF::HistIntegrals(const char* rangeName) const{
      //create histograms for each observable and fill integrating
      //over other observables, filling with evaluateMC
      Long64_t ilow=0;
      Long64_t ihigh=0;
      SetLowHighVals(ilow,ihigh);
      //  cout<<"RooHSEventsPDF::HistIntegrals"<<endl;
      for(Int_t i=0;i<fNvars;i++){
	auto  arg=dynamic_cast<const RooRealVar*>(&fProxSet[i]->arg());
	if(arg)
	  fHistIntegrals.emplace_back(arg->GetName(),arg->GetName(),arg->getBins(),arg->getMin(),arg->getMax());
      }
      Long64_t accepted=0;
      for(Int_t ie=ilow;ie<ihigh;ie++){
	fTreeEntry=ie;
	if(!CheckRange(TString(rangeName).Data())){continue;}
	accepted++;
	Double_t value=evaluateMC(&fvecReal,&fvecCat)*GetIntegralWeight(ie);
	for(Int_t vindex=0;vindex<fNvars;vindex++){
	  fHistIntegrals[vindex].Fill(fvecReal[fTreeEntry*fNvars+vindex],value/fHistIntegrals[vindex].GetBinWidth(1));
	}
      }
      //normalise to number of accepted events
      for(Int_t vindex=0;vindex<fNvars;vindex++)
	fHistIntegrals[vindex].Scale(1./accepted);
  
      fParent->SetHistIntegrals(fHistIntegrals);
  
    }

    void RooHSEventsPDF::SetLowHighVals(Long64_t& ilow,Long64_t& ihigh) const{
      ilow=0;
      ihigh=0;
  
      if(fParent){
	ilow=fParent->GetIntRangeLow();
	ihigh=fParent->GetIntRangeHigh();
      }
      else{
	ilow=GetIntRangeLow();
	ihigh=GetIntRangeHigh();
      }

      if(ihigh==0&&fNInt>-1) ihigh=fNInt;
  
      else if(ihigh==0) ihigh=fNTreeEntries; 
  
      if(ihigh>fNTreeEntries) ihigh=fNTreeEntries;
    }

    Bool_t RooHSEventsPDF::CheckRange(const char* rangeName) const{
      //bool brange=TString(rangeName)==TString("");
      //if(brange) return kTRUE;
      for(UInt_t i=0;i<fProxSet.size();i++){
	//	RooRealVar* var=(dynamic_cast<RooRealVar*>(&(fProxSet[i]->arg())));
	auto var=(dynamic_cast<const RooRealVar*>(&(fProxSet[i]->arg())));
	if(!var->inRange(fvecReal[fTreeEntry*fNvars+i],TString(rangeName).Data())){return kFALSE;}
      }
      return kTRUE;

    }
    Bool_t RooHSEventsPDF::CheckChange() const{
      //Note analytical integral is const funtion so can only change data members
      //which are pointed to something, thus need Double_t *fLast
      //and construct a N-D array where we can change elements

      // std::cout<<"RooHSEventsPDF::CheckChange() "<<fParSet.size()<<std::endl;
      Bool_t hasChanged=false;
      for(Int_t i=1;i<fNpars+1;i++)
	if(fLast[i]!=(*(fParSet[i-1]))){
	  hasChanged=true;
	  //  std::cout<<"RooHSEventsPDF::CheckChange() "<<fParSet[i-1]->GetName()<<" "<<fLast[i]<<" to "<<(*(fParSet[i-1]))<<std::endl;
	}
      if(hasChanged){
	for(Int_t i=1;i<fNpars+1;i++){
	  fLast[i]=*(fParSet[i-1]);
	 
	}
      }
      return hasChanged;
    }
 
    Bool_t RooHSEventsPDF::SetEvTree(TTree* tree,TString cut,TTree* MCGenTree){


      if(!tree->GetEntries())return kFALSE;
      Info("RooHSEventsPDF::SetEvTree"," with name %s and cut  = %s and %lld events",tree->GetName(),cut.Data(), tree->GetEntries());

      //Set the cut
      //Note weight cut can be set with WEIGHT@expr in factory constructor
      if(cut==TString())
	fCut=fInWeightCut;
      else if (fInWeightCut==TString())	
	fCut=cut;
      else
	fCut=cut+"&&"+fInWeightCut;

      //      ProcInfo_t info;
      fEvTree=tree;
      if(MCGenTree!=nullptr){ // generated events used for acceptance correction, do only if tree is available
	fMCGenTree=MCGenTree;
	fHasMCGenTree=kTRUE;
     }
      
     
      fConstInt=fEvTree->GetEntries();//use if constant integral requested
      fEvTree->ResetBranchAddresses();
      //fEvTree->SetBranchStatus("*",0);
      if(MCGenTree){ // generated events used for acceptance correction, do only if tree is available
	fMCGenTree->ResetBranchAddresses();
      }
      
      fBranchStatus=kTRUE;

      //create arrays to connect to tree branches
      TVectorD MCVar(fProxSet.size());
      TVectorD GenVar(fProxSet.size());
      TVectorD MCGenVar(fProxSet.size()); // generated events used for acceptance correction
      vector<Int_t> MCCat(fCatSet.size());
      vector<Int_t> GenCat(fCatSet.size());
      vector<Int_t> MCGenCat(fCatSet.size()); // generated events used for acceptance correction
      fGotGenVar.resize(fProxSet.size());
      fGotGenCat.resize(fProxSet.size());

      //Set branch addresses of tree to data arrays
      for(UInt_t i=0;i<fProxSet.size();i++){
	fGotGenVar[i]=0;
    
	if(fEvTree->GetBranch(fProxSet[i]->GetName())){
	  fEvTree->SetBranchStatus(fProxSet[i]->GetName(),true);
	  fEvTree->SetBranchAddress(fProxSet[i]->GetName(),&MCVar[i]);
	  if(fEvTree->GetBranch(fTruthPrefix+fProxSet[i]->GetName())){
	    fEvTree->SetBranchStatus(fTruthPrefix+fProxSet[i]->GetName(),true);
	    fEvTree->SetBranchAddress(fTruthPrefix+fProxSet[i]->GetName(),&GenVar[i]);
	    fGotGenVar[i]=1;
	    cout<<"Using Generated branch "<<fTruthPrefix+fProxSet[i]->GetName()<<endl;
	  }
	}
	else{
	  Warning("RooHSEventsPDF::SetEvTree","Branch %s not found",fProxSet[i]->GetName()); //May still get set as prototype data
	  //Check if this branch exists in cut and remove it if it does
	  //It may be added later via AddProtoData which should already
	  //have Cut applied to it
	  if(fCut.Contains(fProxSet[i]->GetName())){
	    TString newcut=fCut;
	    newcut.Replace(newcut.Index(fProxSet[i]->GetName())-2,2,"");
	    newcut.Replace(newcut.Index(TString(fProxSet[i]->GetName())+">"),(newcut.Index(TString(fProxSet[i]->GetName())+"<")-newcut.Index(TString(fProxSet[i]->GetName())+">"))*2-1,"");
	    if(newcut==TString("&&"))newcut="";
	    fCut=newcut;
	    cout<<"RooHSEventsPDF::SetEvTree Ammended cut "<<fCut<<endl;
	  }
	  fBranchStatus=kFALSE;
	}

	if(MCGenTree){// generated events used for acceptance correction, do only if tree is available
	  if(fMCGenTree->GetBranch(fProxSet[i]->GetName())){
	    fMCGenTree->SetBranchStatus(fProxSet[i]->GetName(),true);
	    fMCGenTree->SetBranchAddress(fProxSet[i]->GetName(),&MCGenVar[i]);
	  }
	  else{
	    Warning("RooHSEventsPDF::SetEvTree","Branch %s not found in MCGen tree. Acceptance correction will be wrong!!!",fProxSet[i]->GetName()); 
	  }
	}
      }
      for(UInt_t i=0;i<fCatSet.size();i++){
	fGotGenCat[i]=0;
	if(fEvTree->GetBranch(fCatSet[i]->GetName())){
	  fEvTree->SetBranchStatus(fCatSet[i]->GetName(),true);
	  fEvTree->SetBranchAddress(fCatSet[i]->GetName(),&MCCat[i]);
	  if(fEvTree->GetBranch(fTruthPrefix+fCatSet[i]->GetName())){
	    fEvTree->SetBranchStatus(fTruthPrefix+fCatSet[i]->GetName(),true);
	    fEvTree->SetBranchAddress(fTruthPrefix+fCatSet[i]->GetName(),&GenCat[i]);
	    fGotGenCat[i]=1;
	    cout<<"Using Generated branch "<<fTruthPrefix+fCatSet[i]->GetName()<<endl;
	  }	
	}
	else{
	  Warning("RooHSEventsPDF::SetEvTree","Branch %s not found",fCatSet[i]->GetName()); //May still get set as prototype data
	  //Check if this branch exists in cut and remove it if it does
	  //It may be added later via AddProtoData which should already
	  //have Cut applied to it
	  if(fCut.Contains(fCatSet[i]->GetName())){
	    TString newcut=fCut;
	    TString catCut;
	    auto typeIter=fCatSet[i]->arg().typeIterator();
	    catCut+="(";
	    Bool_t first=kTRUE;
	    while(auto type=dynamic_cast<RooCatType*>(typeIter->Next())){
	      if(first){
		catCut+=Form("%s==%d",fCatSet[i]->GetName(),type->getVal());
		first=kFALSE;
	      }
	      else
		catCut+=Form("||%s==%d",fCatSet[i]->GetName(),type->getVal());
	  
	    }	
	    catCut+=")";
	    newcut.ReplaceAll(catCut,"");
	    if(newcut==TString("&&"))newcut="";
	    fCut=newcut;
	    cout<<"RooHSEventsPDF::SetEvTree Ammended cut "<<fCut<<endl;
	  }
	  fBranchStatus=kFALSE;
	}
	if(MCGenTree){// generated events used for acceptance correction, do only if tree is available
	  if(fMCGenTree->GetBranch(fCatSet[i]->GetName())){
	    fMCGenTree->SetBranchStatus(fCatSet[i]->GetName(),true);
	    fMCGenTree->SetBranchAddress(fCatSet[i]->GetName(),&MCGenCat[i]);
	  }
	  else{
	    Warning("RooHSEventsPDF::SetEvTree","Branch %s not found in MCGen tree. Acceptance correction will be wrong!!!",fCatSet[i]->GetName());
	  }
	}
      }


      //Branches are set now can loop over and extract values
      //fEvTree->GetEntry(0);
 
      //Create arrays to store data
      UInt_t ProxSize=fNvars;
      UInt_t CatSize=fNcats;
      fNTreeEntries=fEvTree->GetEntries();
      fvecReal.resize(fNTreeEntries*ProxSize);
      fvecRealGen.resize(fNTreeEntries*ProxSize);
      fvecCat.resize(fNTreeEntries*CatSize);
      fvecCatGen.resize(fNTreeEntries*CatSize);
      if(MCGenTree){// generated events used for acceptance correction, do only if tree is available
	fNMCGenTreeEntries=fMCGenTree->GetEntries();
	fvecRealMCGen.resize(fNMCGenTreeEntries*ProxSize);
	fvecCatMCGen.resize(fNMCGenTreeEntries*CatSize);
       }
      
      Double_t idVal=0;
      Int_t spId=-1;
      if(fWgtsConf.IsValid()){ //add in ID branch for weighted sim data
	fEvWeights.clear();
	LoadInWeights();
	if(fEvTree->GetBranch(fInWeights->GetIDName())){ //the weight ID branch is in fEvTree
	  fUseEvWeights=kTRUE;
	  fEvTree->SetBranchStatus(fInWeights->GetIDName(),true);
	  fEvTree->SetBranchAddress(fInWeights->GetIDName(),&idVal);
	  fEvWeights.resize(fNTreeEntries);
	  spId=fInWeights->GetSpeciesID(fWgtsConf.Species());
	}
	else cout<<"WARNING RooHSEventsPDF::SetEvTree InWeights ID : "<<fInWeights->GetIDName()<<" does not exist in event tree"<<endl;

      }
  
      //Get entries that pass cut
      //A little subtle but this must be done before SetMakeClass or it
      //doesn't find any entries
      tree->Draw(">>elist", fCut, "entrylist");
      auto *elist = dynamic_cast<TEntryList*>(gDirectory->Get("elist"));
      fEvTree->SetEntryList(elist);
      Long64_t entryNumber=0;
      Long64_t localEntry=0;
      fNTreeEntries=elist->GetN();
	  
      TEntryList* elistMCGen=nullptr;
      if(MCGenTree){// generated events used for acceptance correction, do only if tree is available
	MCGenTree->Draw(">>elistMCGen", "", "entrylistMCGen"); // TODO include fCut???
	elistMCGen = dynamic_cast<TEntryList*>(gDirectory->Get("elistMCGen"));
	fMCGenTree->SetEntryList(elistMCGen);
	fNMCGenTreeEntries=elistMCGen->GetN();
      }
      
  
      //Read weights into fEvWeights
      TBranch* idBranch=nullptr;
      vector<Long64_t> idEntries;


      //Now ready to loop over events and store data
      Long64_t corrEvent=0;
      for(Long64_t iEvent=0;iEvent<fNTreeEntries;iEvent++){
	entryNumber = fEvTree->GetEntryNumber(iEvent);
	if (entryNumber < 0) break;
	localEntry = fEvTree->LoadTree(entryNumber);
	if (localEntry < 0) break;
	fEvTree->GetEntry(localEntry);
	
	Bool_t removeNaNEvent=false;//in case of NaN

	for(UInt_t ip=0;ip<ProxSize;ip++){
	  if( TMath::IsNaN(MCVar[ip]) ){
	    //this event contains a NaN so will ignore
	    cout<<"RooHSEventsPDF::SetEvTree "<<GetName()<<" event with NaN will  remove it "<< localEntry<<endl;
	    removeNaNEvent=true;
	  }
	  else{
	    fvecReal[corrEvent*ProxSize+ip]=MCVar[ip];
	    //Read the generated values if exist if not
	    //use mcvar again, this duplicates data so should
	    //be better optimised
	    if(!fGotGenVar[ip]) fvecRealGen[corrEvent*ProxSize+ip]=MCVar[ip];
	    else fvecRealGen[corrEvent*ProxSize+ip]=GenVar[ip];
	  }
	}
	if(removeNaNEvent) continue;
	//This event has passed all requirements and we are going to keep it
	fTreeEntryNumber.push_back(localEntry);

	//Get weights if used
	if(fUseEvWeights==kTRUE){ 
	  fInWeights->GetEntryBinarySearch(static_cast<Long64_t>(idVal));
	  fEvWeights[corrEvent]=fInWeights->GetWeight(spId);
	}
	//and any categories
	for(UInt_t ip=0;ip<CatSize;ip++){
	  fvecCat[corrEvent*CatSize+ip]=MCCat[ip];
	  if(!fGotGenCat[ip]) fvecCatGen[corrEvent*CatSize+ip]=MCCat[ip];
	  else fvecCatGen[corrEvent*CatSize+ip]=GenCat[ip];
	}
	corrEvent++;
      }
      //finished data loop

      
      if(fNTreeEntries!=corrEvent){
	cout<<"RooHSEventsPDF::SetEvTree note only accepted "<<corrEvent<<" out of "<<fNTreeEntries<<" original events"<<endl;
	fNTreeEntries=corrEvent;
      }
      delete fInWeights;fInWeights=nullptr;
      fEvTree->SetEntryList(nullptr);
      delete elist;elist=nullptr;
 
      entryNumber=0;
      localEntry=0;
      if(MCGenTree){// generated events used for acceptance correction, do only if tree is available
	for(Long64_t iEvent=0;iEvent<fNMCGenTreeEntries;iEvent++){
	  entryNumber = fMCGenTree->GetEntryNumber(iEvent);
	  if (entryNumber < 0)
	    break;
	  localEntry = fMCGenTree->LoadTree(entryNumber);
	  if (localEntry < 0)
	    break;
	  fMCGenTree->GetEntry(localEntry);
	  for(UInt_t ip=0;ip<ProxSize;ip++){
	    //  cout<<iEvent<<" "<<MCVar[ip]<<endl;
	    fvecRealMCGen[iEvent*ProxSize+ip]=MCGenVar[ip];
	  }
	  for(UInt_t ip=0;ip<CatSize;ip++){
	    fvecCatMCGen[iEvent*CatSize+ip]=MCGenCat[ip];
	  }
	}
	fMCGenTree->SetEntryList(nullptr);
	delete elistMCGen;elistMCGen=nullptr;
      }
      
      //reset everything so we don't screw up memory
      fEvTree->ResetBranchAddresses();  
      fEvTree->Reset();  //empty tree to save memory

      if(fMCGenTree){
	fMCGenTree->ResetBranchAddresses();
  	fMCGenTree->Reset();
      }
      Info("RooHSEventsPDF::SetEvTree"," with name %s and cut  = %s and kept %lld events",tree->GetName(),cut.Data(), fNTreeEntries);
       return fBranchStatus;
    }
    
    void  RooHSEventsPDF::LoadInWeights(){
      //GetWeights object 
      cout<<"void RooHSEventsPDF::LoadWeights and use species "<< fWgtsConf.Species()<<" "<<fWgtsConf.File()<<" "<<fWgtsConf.ObjName()<<endl;
      if(fInWeights) delete fInWeights;
      fInWeights=nullptr;
      fInWeights=new Weights();

      fInWeights->LoadSavedDisc(fWgtsConf.File(),fWgtsConf.ObjName());

      if(fInWeights->GetSpeciesID(fWgtsConf.Species())==-1){
	cout<<"ERROR RooHSEventsPDF::LoadInWeights() requested species "<<fWgtsConf.Species()<<" not found in given weights"<<endl;
	fInWeights->Print();
	fIsValid=kFALSE;
      }
    }

    void  RooHSEventsPDF::CheckIntegralParDep(Int_t Ntests){
      fCheckInt=Ntests;
      if(!fEvTree) return; //will check later when tree is set
      if(!fBranchStatus) return; //missing branches may need to add protodata
  
      Long64_t saveNint=fNInt;
      // fNInt=nint;//REOMVE FOR NOW nint is not used
      fNInt=fNTreeEntries;//just use all entries
      //scale Ntests by Ndimensions
      //   Ntests=TMath::Power((Double_t)Ntests,(Double_t)fParSet.size());
      Ntests=(Ntests*fParSet.size());

      Info("RooHSEventsPDF::CheckIntegralParDep","Going to run %d calculations of integral with random parameters",Ntests);
  
      RooRealVar integral("integral","integral",0,0,2);
      integral.setError(sqrt(fNInt)/fNInt); //Error needs to be set before entering in ds
      RooDataSet ds("intds","intds",RooArgSet(integral));
      //want to set random parameter values
      //loop over each parameter and calculate integral
      vector<Double_t> SavedPars;
      for(auto & ip : fParSet){//first save parameters
        auto par=(dynamic_cast<const RooRealVar*>(&(ip->arg())));
	SavedPars.push_back(par->getValV());
      }
      for(Int_t ir=0;ir<Ntests;ir++){ //loop over tests
	for(auto & ip : fParSet){//loop over parameters
	  auto par=(RooRealVar*)(&(ip->arg()));
	  par->setVal((par->getMax("")-par->getMin(""))*RooRandom::uniform()+par->getMin(""));
	}
	//Now calculate integral
	integral.setVal(analyticalIntegral(1,""));
	ds.add(RooArgSet(integral));
      }

 
      Double_t low=0;
      Double_t high=0;
      ds.getRange(integral,low,high);
      integral.setRange(low,high);
      RooPlot *frame=integral.frame();
      ds.plotOn(frame);

      frame->Draw();

      RooRealVar  mean("mean","mean",ds.mean(integral));
      RooRealVar pvar("IntPull","Integral Pull Dist.",-5,5);
      RooPullVar pull("IntPull","Integral Pull Dist.",integral,mean);
      ds.addColumn(pull,kFALSE);

      ds.Print();
      ds.getRange(pvar,low,high);
      pvar.setRange(low,high);
      RooPlot *framePull=pvar.frame();
      ds.plotOn(framePull);
      RooRealVar mp("mp","mp",0,-5,5) ;
      RooRealVar sp("sp","sp",1,0,100);
      RooGaussian gp("gp","gp",pvar,mp,sp);
      gp.fitTo(ds);
      gp.paramOn(framePull);
      gp.plotOn(framePull);
  
      pvar.Print();
      new TCanvas();
      framePull->Draw();
  
      //numerical check gives constant integral, can force const for fit speed
      fConstInt=mean.getVal();
      fNInt=saveNint;
      if(sp.getVal()<2){
	Info("RooHSEventsPDF::CheckIntegralParDep","Numerical integral constant. Will not recalculate during fits to save time, if you want to force calculation do not call this function");
	SetConstInt();
      }
      //reset pars
      for(UInt_t ip=0;ip<fParSet.size();ip++){//first save parameters
	auto par=(RooRealVar*)(&(fParSet[ip]->arg()));
	par->setVal(SavedPars[ip]);
      }
      fCheckInt=kFALSE; //only do once
    }
    Bool_t RooHSEventsPDF::AddProtoData(const RooDataSet* data){
      //merge the current tree with data from another dataset
      //Default it will add any branches in data not in fEvTree
      //      cout<<"RooHSEventsPDF::AddProtoData "<<data<<" "<<fEvTree<<endl;
      if(!fEvTree) return kFALSE;
      if(!fNTreeEntries) return kFALSE;
  
      //Loop over data branches and check if any missing
      const RooArgSet *dataVars=data->get();
      RooArgSet thisVars =VarSet(0);
  
      //Each entry in fEvTree should be given a random value of new variables
      Long64_t Nentries=data->numEntries();
      vector<Long64_t> vrandom(Nentries);
      for(Long64_t ir=0;ir<Nentries;ir++)
	vrandom[ir]=ir;
      std::shuffle(vrandom.begin(),vrandom.end(), std::mt19937(std::random_device()()));

      TIter iter=dataVars->createIterator();
      vector<Double_t> brD;
      vector<Int_t> brI;
      Int_t ND=0;
      Int_t NI=0;
      vector<Bool_t> isD;
      vector<Int_t> index;
      TObjArray branches;
      //Add new branches to EvTree
      TDirectory* saveDir=gDirectory;
      fEvTree->GetDirectory()->cd();

      vector<Short_t> protoDataForVar;
      vector<Short_t> protoDataForCat;
  
      //Look for variables in data that were not in fEvTree
      while(auto* arg=dynamic_cast<RooAbsArg*>(iter())){
	if(TString("UID")==arg->GetName()) continue; //don't replicate ID branch
	if(fEvTree->GetBranch(arg->GetName())) continue; //already exists
	for(Int_t ip=0;ip<fNvars;ip++)
	  if(TString(arg->GetName())==TString(fProxSet[ip]->GetName()))
	    protoDataForVar.push_back(ip);
	for(Int_t ip=0;ip<fNcats;ip++)
	  if(TString(arg->GetName())==TString(fCatSet[ip]->GetName()))
	    protoDataForCat.push_back(ip);
	Info("RooHSEventsPDF::AddProtoData","Added data branch %s",arg->GetName());
      }


      //Loop over data in random order and fill new branches
      Long64_t idata=0;

      UInt_t Nreal=protoDataForVar.size();
      UInt_t Ncat=protoDataForCat.size();
  
      if(!(Nreal+Ncat)) return kTRUE;
  
      //copy proto data to vecs
      for(Long64_t id=0;id<fNTreeEntries;id++){
	dataVars=data->get(vrandom[idata]);
	for(short ip : protoDataForVar){
	  Double_t val=dataVars->getRealValue(fProxSet[ip]->GetName());
	  fvecReal[id*fNvars+ip]=val;
	  fvecRealGen[id*fNvars+ip]=val;
	}  
	for(short ip : protoDataForCat){
	  Int_t val=dataVars->getCatIndex(fCatSet[ip]->GetName());
	  fvecCat[id*fNcats+ip]=val;
	  fvecCatGen[id*fNcats+ip]=val;     
	}
    
	if(idata==(Long64_t)vrandom.size()-1){//Need to reuse data until done all MC events
	  std::shuffle(vrandom.begin(),vrandom.end(), std::mt19937(std::random_device()()));
	  idata=0;
	}
	idata++;
      }
  
      saveDir->cd();
  
      return fBranchStatus;  
    }
    void RooHSEventsPDF::ResetTree(){
  
      if(fEvTree) {delete fEvTree;fEvTree=nullptr;}
    }

    void RooHSEventsPDF::SetNextRange(Int_t ir){
      Long64_t Nentries=fNTreeEntries;
      Int_t range=((Double_t)Nentries)/fNRanges;

      fIntRangeLow=ir*range;
      fIntRangeHigh=(ir+1)*range;

    }

    void RooHSEventsPDF::MakeAssertPostiveData(){
      //   RooArgSet vars = VarSet(0);
      auto saveTreeEntry=fTreeEntry;
      fTreeEntry=0;

      auto NVars=fProxSet.size();
      fAssertPosDataReal.resize(fNapd*NVars);

      auto NCats=fCatSet.size();
      fAssertPosDataCats.resize(fNapd*NCats);

       for(Long64_t iapd = 0 ;iapd<fNapd; ++iapd ){
	
	UInt_t ivar=0;
	for(auto v:fProxSet){
	  //	  cout<<"RooHSEventsPDF::MakeAssertPostiveData() "<<iapd<<" "<<v->GetName()<<std::endl;
	  auto vr = dynamic_cast<const RooRealVar*>(&v->arg());
	  if(vr!=nullptr){
	    fAssertPosDataReal[fTreeEntry*NVars+ivar]=(gRandom->Uniform(vr->getMin(""),vr->getMax("")));
	    ++ivar;
	  }
	}
	UInt_t icat=0;
	for(auto v:fCatSet){
	  auto vc = dynamic_cast<const RooCategory*>(&v->arg());
	
	  if(vc!=nullptr){
	    auto catstate = gRandom->Integer(vc->size());
	    auto val = vc->getOrdinal(catstate).second;
	   
	    fAssertPosDataCats[fTreeEntry*NCats+icat]=val;
	    ++icat;
	  }
	}
	
	fTreeEntry++;
      }
       // cout<<fAssertPosDataReal.size()<<" "<<fAssertPosDataReal[0]<<" "<<fAssertPosDataCats[0]<<std::endl;exit(0);
      fTreeEntry=saveTreeEntry;
    }
   // void RooHSEventsPDF::MakeAssertPostiveData(){
   //    RooArgSet vars = VarSet(0);
   //    auto saveTreeEntry=fTreeEntry;
   //    fTreeEntry=0;
   //    auto NVarsAndCats=vars.getSize();
   //    fAssertPosDataReal.resize(fNapd*NVarsAndCats);
   //    cout<<"RooHSEventsPDF::MakeAssertPostiveData() "<<fAssertPosDataReal.size()<<" "<<fAssertPosDataReal[0]<<std::endl;
   //    for(Long64_t iapd = 0 ;iapd<fNapd; ++iapd ){
	
   // 	//whatabout categories !!!
   // 	UInt_t ivar=0;
   // 	for(auto v:vars){
   // 	  cout<<"RooHSEventsPDF::MakeAssertPostiveData() "<<iapd<<" "<<v->GetName()<<std::endl;
   // 	  auto vr = dynamic_cast<RooRealVar*>(v);
   // 	  if(vr!=nullptr){
   // 	    fAssertPosDataReal[fTreeEntry*NVarsAndCats+ivar]=(gRandom->Uniform(vr->getMin(""),vr->getMax("")));

   // 	  }
   // 	  else{//category
   // 	    auto vc = dynamic_cast<RooCategory*>(v);
   // 	    auto catstate = gRandom->Integer(vc->size());
   // 	    auto val = vc->getOrdinal(catstate).second;
   // 	    fAssertPosDataCat[fTreeEntry*NVarsAndCats+ivar]=val;
   // 	  }
   // 	  ++ivar;
   // 	}
   // 	fTreeEntry++;
   //    }
   //    //cout<<fAssertPosDataReal.size()<<" "<<fAssertPosDataReal[0]<<std::endl;
   //    fTreeEntry=saveTreeEntry;
   //  }
    Bool_t RooHSEventsPDF::AssertPositivePDF() const{
      //make sure this PDF is >=0 for its full allowed range
      // std::cout<<"RooHSEventsPDF::AssertPositivePDF()"<<std::endl;
      InitAssertPositiveCheck() ;
      //std::cout<<"RooHSEventsPDF::AssertPositivePDF() done init"<<std::endl;
      auto saveTreeEntry=fTreeEntry;
      fTreeEntry=0;
      // cout<<"RooHSEventsPDF::AssertPositivePDF()"<<endl;
      for(Long64_t iapd = 0 ;iapd<fNapd; ++iapd ){
       	auto val = evaluateMC(&fAssertPosDataReal,&fAssertPosDataCats);
	++fTreeEntry;
	if(val<-1E-4){ //some tolerance
	  //  cout<<"RooHSEventsPDF::AssertPositivePDF() PDF cannot be -ve. "<<val<<" "<<fTreeEntry<<endl;
	  logEvalError("RooHSEventsPDF::AssertPositivePDF() PDF cannot be -ve...");
	  fTreeEntry=saveTreeEntry;
	  FinishAssertPositiveCheck(); 
	  return kFALSE;
	}
      }
      fTreeEntry=saveTreeEntry;
      //std::cout<<"RooHSEventsPDF::AssertPositivePDF() done "<<std::endl;
       FinishAssertPositiveCheck() ;
     
      return kTRUE;
    }

    
  }
}

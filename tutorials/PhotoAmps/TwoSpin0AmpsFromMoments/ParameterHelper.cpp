#include "ParameterHelper.h"
#include <TMath.h>
#include <TRandom.h>

namespace m2pw{
  
  ///////////////////////////////////////////////////////////////
  ///Record a variable parameter
  Bool_t ParameterHelper::AddParameter(const RooRealVar* var){
    
    if(var->isConstant()){ //need to store constant value
	 return  AddConstant(var);
    }

    //now deal with variables
    TString name=var->GetName();
    //check if already loaded
    if(_nameToIndex.find(name)==_nameToIndex.end()){
      //add it to map and give it an index
      //this is keeping everything in synch
      _nameToIndex[name]=_nextIndex;
      _names.push_back(name);
      //init values, will be overwritten
      _currentVals.push_back(var->getVal());
      _maxVals.push_back(var->getMax());
      _minVals.push_back(var->getMin());
      _stepSize.push_back((var->getMax()-var->getMin())/10);

      //Keep note of magnitudes and phases
      //Assume using PhotoAmps e.g a_Xor aphi_X etc
      _isMagnitude.push_back(name.Contains("phi")==kFALSE);

      //increment
      _nextIndex++;

      
      return true;
    }
    //already have it
    return false;
  }
  //////////////////////////////////////////////////////////
  ///Record a constant parameter
  Bool_t ParameterHelper::AddConstant(const RooRealVar* cvar){

    //Record a constant parameter
    TString name = cvar->GetName();
    if(_constToIndex.find(name)==_constToIndex.end()){
      _constToIndex[name]=_nextConstIndex;
      _constNames.push_back(name);
      _constVals.push_back(cvar->getVal());
     _nextConstIndex++;
      return true;
    }
    return false;
  }
  ///////////////////////////////////////////////////////////
  ///Parse RooFormulaVar to get parameters and Add them
  //void ParameterHelper::UseRooFormula(const RooFormula* form){   
  void ParameterHelper::UseRooFormula(const RooFormulaVar* form){   

    RooAbsArg *par=nullptr;
    //loop over parameters in formula
    Int_t ipar=0;
    while( (par=form->getParameter(ipar++)) != nullptr ){
      RooRealVar* var=dynamic_cast<RooRealVar*>(par);   
      if( var !=nullptr){//only deal with RooRealVars
	AddParameter(var);
      }
    }//next parameter
  }
  ///////////////////////////////////////////////////////////////////////
  ///Return constant parameters in given formula
  std::vector<Double_t > ParameterHelper::ConstantValues(const RooFormulaVar* form) const {
    RooAbsArg *par=nullptr;
    Int_t ipar=0;
    std::vector<Double_t > vals;
    while( (par=form->getParameter(ipar++)) != nullptr ){
      if(IsConst(par->GetName()))
	vals.push_back(GetConstVal(par->GetName()));
    }
      
    return vals;
  }
  ///////////////////////////////////////////////////////////////////
  /// Return indices of parameters used in given formula
  std::vector<Int_t > ParameterHelper::Indices(const RooFormulaVar* form) const{
    Int_t ipar=0;
    std::vector<Int_t > indices;
    RooAbsArg *par=nullptr;
    while( (par=form->getParameter(ipar++)) != nullptr ){
      auto index = Index(par->GetName());
      indices.push_back(index); 
     }
    return indices;
  }
  //////////////////////////////////////////////////////////////////
  ///return indices of parameters which are formula in this equation
  std::vector<TString > ParameterHelper::Dependencies(const RooFormulaVar* form) const{
    Int_t ipar=0;
    std::vector<TString > deps;
    RooAbsArg *par=nullptr;
    while( (par=form->getParameter(ipar++)) != nullptr ){
      auto index = Index(par->GetName());
      if(index==static_cast<Int_t>(ParType::EDependancy)){
	deps.push_back(par->GetName()); 
      }
    }
    return deps;
  }

  /*
  ////////////////////////////////////////////////////////////
  ///Check if values have changed, needs recaching
  Bool_t ParameterHelper::CheckChange(std::vector<Int_t > indices) const{
    for(auto& index:indices){
      if(_currentVals[index] != _cachedVals[index]) return kTRUE;
    }
    return kFALSE;
  }
  ///////////////////////////////////////////////////////////
  ///store variable values used for caching
  void ParameterHelper::StoreCacheVals() const{
    for(UInt_t index = 0; index<NVars();++index){
      _cachedVals[index] = currentVals[index];
    }
  }
  */
  /////////////////////////////////////////////////////////////
  ///Create a tree with parameter branches
  void ParameterHelper::MakeParTree(const TString& filename){
    _file.reset(TFile::Open(filename,"recreate"));
    _tree = new TTree("PartialWaves","pw");
    for(UInt_t i=0;i<_currentVals.size();++i){
      TString brname = GetName(i);
      if(brname.Contains("-")) brname.ReplaceAll("-","m");
      _tree->Branch(brname,&_currentVals[i],brname+"/D");
    }
    _tree->Branch("sumMags",&_sumMags,"sumMags/D");
   
  }
  ////////////////////////////////////////////////////////////
  ///Put parameters in nominal ranges [0,1],[-pi,pi]
  void ParameterHelper::RerangeParameters(){
    for(UInt_t i=0;i<_currentVals.size();++i){
      if(IsMagnitude(i)==kFALSE){
	_currentVals[i]=std::remainder(_currentVals[i],2*TMath::Pi());
      }
      else{
	_currentVals[i]=TMath::Abs( _currentVals[i]);
      }
    }
  }
  ////////////////////////////////////////////////////////////
  ///Randomise parameters by creating and passing unit cube
  void ParameterHelper::Randomise(){
    gRandom->RndmArray(Nvars(),_currentVals.data());
    TransformConstrained(_currentVals.data());
  }
  ////////////////////////////////////////////////////////////
  ///Randomise parameters by passing unit cube
  ///use Dirichlet_distribution for magnitudes^2
  ///https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation
  void ParameterHelper::TransformConstrained(double* transPars) const{
    std::vector<Double_t > gamma_par; //going to renormalise so need to keep individual components
    std::vector<Int_t > gamma_par_i;
    Double_t sum_mag_sq=0;
   //loop over all pars
    for(Int_t i=0;i<Nvars();i++){
      if(IsMagnitude(i)==kFALSE){//phase
	//random in full range
	transPars[i]=Min(i) + (transPars[i])*(Max(i)-Min(i)) ;
	//rerange to -pi to pi
	transPars[i]=std::remainder(transPars[i],2*TMath::Pi());

      }
      else{
	// first inverse transform sample from Gamma(alpha=1,beta=1), which is Exponential(1)
	//For amplitudes we need to square as sum of squares ==1
	//transform first into allowed range
	transPars[i]=Min(i) + (transPars[i])*(Max(i)-Min(i)) ;
	auto sqmag=-TMath::Log(transPars[i]*transPars[i]);
	//if(gRandom->Uniform()>0.3) sqmag=0;
	gamma_par.push_back(sqmag);
	gamma_par_i.push_back(i);
	sum_mag_sq+=sqmag;

      }
    }
    //Now include magnitude of constrained amplitude
    Double_t randMag=gRandom->Uniform(); //usse 0-1
    Double_t conmag_sq=-TMath::Log(randMag*randMag);
    sum_mag_sq+=conmag_sq;
    
    //renormalise so sum of squares == 1
    for(UInt_t imag=0;imag<gamma_par_i.size();++imag){
       transPars[imag]=TMath::Sqrt(gamma_par[imag]/(sum_mag_sq)); 
     }
    

  }
  void ParameterHelper::ZeroSmallAmps() {
    Double_t sumsq=0;
    for(Int_t i=0;i<Nvars();i++){
      if(IsMagnitude(i)==kTRUE){//phase
	if(TMath::Abs(_currentVals[i])<0.05)
	  _currentVals[i]=0;
	sumsq+=_currentVals[i]*_currentVals[i];
      }
    }
    //if constrained amp is <0.05 set it to 0 by increasing first non zero amp
    // Double_t conAmp=sqrt(1-sumsq);
    // if(conAmp<0.05){
    //   for(Int_t i=0;i<Nvars();i++){
    // 	if(IsMagnitude(i)==kTRUE){//phase
    // 	  if(_currentVals[i]>0.05)
    // 	    _currentVals[i]+=(conAmp-1E-6);
    // 	}
    //   }
    // }
  }
      
  void ParameterHelper::MakePertubation(Double_t psize){
    Double_t sum_sq_mags=0;
    //Find the constrained value
    for(UInt_t i=0;i<_currentVals.size();++i){
      if(IsMagnitude(i))
	sum_sq_mags+=_currentVals[i]*_currentVals[i];
    }
    Double_t constrained=(TMath::Sqrt(1-sum_sq_mags));
    for(UInt_t i=0;i<_currentVals.size();++i){
      auto original= _currentVals[i];
     
      if(IsMagnitude(i)==false){
	_currentVals[i]+=gRandom->Uniform()*(_maxVals[i]-_minVals[i])*psize;
      }
      else{
	if(TMath::Abs(_currentVals[i])<0.01){
	  _currentVals[i]=0;
	}
	else{
	  _currentVals[i]+=gRandom->Uniform()*(_maxVals[i]-_minVals[i])*psize*_currentVals[i];
	}
      }
    }

    //MAke sure we stay in range
    sum_sq_mags=0;
   for(UInt_t i=0;i<_currentVals.size();++i){
     if(IsMagnitude(i)){
       sum_sq_mags+=_currentVals[i]*_currentVals[i];
     }
     else {
       if(_currentVals[i]<_minVals[i])_currentVals[i]=_minVals[i]+gRandom->Gaus(TMath::Pi()/2,0.1);// -(_currentVals[i]-_minVals[i]);
       if(_currentVals[i]>_maxVals[i])_currentVals[i]=_maxVals[i]-gRandom->Gaus(TMath::Pi()/2,0.1);// -(_currentVals[i]-_maxVals[i]);
       if(_currentVals[i]<_minVals[i])_currentVals[i]=_maxVals[i]-0.1;//+(_currentVals[i]-_minVals[i]);
       if(_currentVals[i]>_maxVals[i])_currentVals[i]=_minVals[i]+0.1;//+(_currentVals[i]-_maxVals[i]);// -(_currentVals[i]-_maxVals[i]);
       
     }
   }
   sum_sq_mags+=constrained*constrained;
   double check =0;
   for(UInt_t i=0;i<_currentVals.size();++i){
     if(IsMagnitude(i)){
       //	cout<<sum_sq_mags<<endl;
       _currentVals[i]=TMath::Sqrt(_currentVals[i]*_currentVals[i]/sum_sq_mags);
       check+=_currentVals[i]*_currentVals[i];
     }
   }
   // cout<<"check pertubation "<<check<<" "<<constrained<<" "<<sum_sq_mags<<endl;
  }
  Bool_t ParameterHelper::CheckInRange(const double* pars) const {

    for(Int_t ipar=0;ipar<Nvars();ipar++){
      if(pars[ipar]<_minVals[ipar] || pars[ipar]>_maxVals[ipar] )
	return kFALSE;
    }
    return kTRUE;
    
  }
  Double_t ParameterHelper::SumMags() const{
    Double_t sum=0;
    Double_t sumsq=0;
    for(UInt_t i=0;i<_currentVals.size();++i){
      if(IsMagnitude(i)){
	sum += _currentVals[i];
	sumsq += _currentVals[i]*_currentVals[i];
	
      }
    }
    sum+=sqrt(1-sumsq);

    _sumMags=sum;
    return sum;
  }
  void ParameterHelper::Print(const TString opt) const{
    std::cout<<"\t ParameterHelper::Print() variable parameters :" <<std::endl;
    UInt_t isynch=0;
    for(auto name :_names){
      std::cout<<"\t\t"<<name<<" = "<<_currentVals[isynch]
	       <<" > "<<_minVals[isynch]<<" < "<<_maxVals[isynch]
	       <<" is mag "<<_isMagnitude[isynch]<<std::endl;
      isynch++;
    }
    isynch=0;
    std::cout<<"\t ParameterHelper::Print() constant parameters :" <<std::endl;
    for(auto name :_constNames){
      std::cout<<"\t\t"<<name<<" = "<<_constVals[isynch]<<std::endl;
      isynch++;
    }
  }//Print
}

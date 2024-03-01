/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include <Riostream.h> 

#include "RooHSEventsHistPDF.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <cmath> 
#include <TMath.h> 
#include <TF1.h>
#include <TH1.h>
#include <TRandom3.h>

namespace HS{
  namespace FIT{

    RooHSEventsHistPDF::RooHSEventsHistPDF(const char *name, const char *title, RooAbsReal& _x,RooAbsReal& _alpha,RooAbsReal& _offset, RooAbsReal& _scale) :
      RooHSEventsPDF(name,title),
      x("x","x",this,_x),
      offset("offset","offset",this,_offset),
      scale("scale","scale",this,_scale),
      alpha("alpha","alpha",this,_alpha)
    {
  
      MakeSets();
      x.SetName(_x.GetName());
      offset.SetName(_offset.GetName());
      scale.SetName(_scale.GetName());
      alpha.SetName(_alpha.GetName());
      //Additional stuff
      auto *rx=dynamic_cast<RooRealVar*>(&_x);
      auto *ra=dynamic_cast<RooRealVar*>(&_alpha);
      auto *rs=dynamic_cast<RooRealVar*>(&_scale);
      auto *ro=dynamic_cast<RooRealVar*>(&_offset);
      cout<<"RooHSEventsHistPDF::RooHSEventsHistPDF using fudge parameters:"<<endl;
      std::vector<RooRealVar*> vars = {ra, rs, ro};
      for(auto& var : vars){
        if(var->getMax()==var->getMin()){  //turn parameters with zero intervals into constants
          var->setConstant();
        }
        else if(var->isConstant()){  //ensure that constant parameters have zero interval; otherwise RooAbsPdf complains about zero normalization integral
          var->setMin(var->getVal());
          var->setMax(var->getVal());
        }
        var->Print();
      }
      //work out hist bin limits when scale paramter used
      Int_t NBins0=rx->getBins();
      Double_t rsmin=1;
      //fOldScale=rs->getVal();
      if(rs->getMin()) rsmin=rs->getMin();
      else cout<<"RooHSEventsHistPDF::RooHSEventsHistPDF Warning no scale minimum set take = 1"<<endl;
      Double_t mid=(rx->getMax()+rx->getMin())/2;
      Double_t diff=(rx->getMax()-rx->getMin())/2;
      Double_t rMin=mid-diff/rsmin + ro->getMin(); //additional range for possible tranformation or scaling
      Double_t rMax=mid+diff/rsmin + ro->getMax();
      Int_t NbinX=NBins0/rsmin;
      if(NbinX<10) NbinX=10; //force minimum of 10 bins
      cout<<"RooHSEventsHistPDF::RooHSEventsHistPDF using binning ("<<NbinX<<", "<<rMin<<", "<<rMax<<") for x-axis variable '"<<_x.GetName()<<"' of PDF '"<<name<<"'"<<endl;

      Int_t NAlphBins=500;
      TAxis AlphAxis(NAlphBins,ra->getMin(),ra->getMax());
   
     if(ra->isConstant())
	fRHist=new TH2D(TString("hmc_model_")+_x.GetName()+name,TString("MC model for ")+_x.GetName(),NbinX,rMin,rMax,1,ra->getMin()-1,ra->getMin()+1);
      else//Note want to centre first alpha bin on 0
	fRHist=new TH2D(TString("hmc_model_")+_x.GetName()+name,TString("MC model for ")+_x.GetName(),NbinX,rMin,rMax,NAlphBins+1,ra->getMin()-AlphAxis.GetBinCenter(1),ra->getMax()+AlphAxis.GetBinCenter(1));

      fRHist->SetDirectory(nullptr);
   
      fx_off=new RooRealVar(TString("off")+_x.GetName(),"Vx_off",mid,rMin,rMax);
      falpha=new RooRealVar("Valpha","Valpha",0,alpha.min(),alpha.max());

      //Mak Gaussian constraints for parameters
      //alpha = gaussian centre on val, width = max/5
      if(not ra->isConstant())fAlphaConstr=new RooGaussian(TString("AlphaConstr")+GetName(),"AlphaConstr",_alpha,RooFit::RooConst(ra->getVal()),RooFit::RooConst(ra->getMax()/5));
      //off = gaussian centre on val, width = range/5/2
      if(not ro->isConstant())fOffConstr=new RooGaussian(TString("OffConstr")+GetName(),"OffConstr",_offset,RooFit::RooConst(ro->getVal()),RooFit::RooConst((ro->getMax()-ro->getMin())/5/2));
      if(not rs->isConstant())fScaleConstr=new RooGaussian(TString("ScConstr")+GetName(),"ScConstr",_scale,RooFit::RooConst(rs->getVal()),RooFit::RooConst((rs->getMax()-rs->getMin())/5/2));
      // if(fAlphaConstr)fAlphaConstr->Print();
      // if(fOffConstr)fOffConstr->Print();
      // if(fScaleConstr)fScaleConstr->Print();

    } 


    RooHSEventsHistPDF::RooHSEventsHistPDF(const RooHSEventsHistPDF& other, const char* name) :  
      RooHSEventsPDF(other,name),
      x("x",this,other.x),
      offset("offset",this,other.offset),
      scale("scale",this,other.scale),
      alpha("alpha",this,other.alpha)
    {
      // cout<<"COPYNG "<<other.fRHist<<" "<<other.fEvTree<<endl;
      MakeSets();
      x.SetName(other.x.GetName());
      offset.SetName(other.offset.GetName());
      scale.SetName(other.scale.GetName());
      alpha.SetName(other.alpha.GetName());
      //Additional stuff
      if(other.fx_off)fx_off=dynamic_cast<RooRealVar*>(other.fx_off->Clone());
      if(other.falpha)falpha=dynamic_cast<RooRealVar*>(other.falpha->Clone());
      if(other.fHist)fHist=dynamic_cast<RooDataHist*>(other.fHist->Clone(other.fHist->GetName()));
      if(other.fRHist)fRHist=dynamic_cast<TH2D*>(other.fRHist->Clone(other.fRHist->GetName()));
      if(other.fAlphaConstr)fAlphaConstr=dynamic_cast<RooGaussian*>(other.fAlphaConstr->Clone());
      if(other.fOffConstr)fOffConstr=dynamic_cast<RooGaussian*>(other.fOffConstr->Clone());
      if(other.fScaleConstr)fScaleConstr=dynamic_cast<RooGaussian*>(other.fScaleConstr->Clone());
      fVarMax=other.fVarMax;
      // if(fEvTree) SetEvTree(fEvTree,fCut);//Needs fProxSet filled first

      fRHist->SetDirectory(nullptr);
      MakeSets();
    
    }

    RooHSEventsHistPDF::~RooHSEventsHistPDF(){

      if(fHist)delete fHist;
      if(fRHist)delete fRHist; fRHist=nullptr;
      if(fx_off) delete fx_off;
      if(falpha)delete falpha;
      if(fAlphaConstr) delete fAlphaConstr;
      if(fOffConstr) delete fOffConstr;
      if(fScaleConstr) delete fScaleConstr;

    }
    void RooHSEventsHistPDF::MakeSets(){
      fProxSet.clear();
      fVarSet.clear();
      fParSet.clear();
      fProxSet.push_back(&x);
      fParSet.push_back(&offset);
      fParSet.push_back(&scale);
      fParSet.push_back(&alpha);
      InitSets();
    }



    Double_t RooHSEventsHistPDF::evaluate() const 
    {
      if(!fHist) return 1;
      // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
      Double_t arg=(x-fVarMax)*scale+fVarMax;
      arg=arg-offset;
      fx_off->setVal(arg);
      falpha->setVal(Double_t(alpha));
      if(fHist->weight(RooArgSet(*fx_off,*falpha),1,kFALSE)<0)   cout<<fHist->weight(RooArgSet(*fx_off,*falpha),1,kFALSE)<<" "<<arg<<" "<<x<<endl;

      return  fHist->weight(RooArgSet(*fx_off,*falpha),1,kFALSE);
    } 

    Double_t RooHSEventsHistPDF::evaluateMC(const vector<Float_t> *vars,const  vector<Int_t> *cats) const {
      Double_t mcx=(*vars)[0];
 
      return evaluateMC(mcx);  
    }
    Double_t RooHSEventsHistPDF::evaluateMC(Double_t mcx) const {
      Double_t arg=(mcx-fVarMax)*scale+fVarMax;
      //cout<<"DEBUG RooHSEventsHistPDF::evaluateMC "<<fHist<<" "<<arg<<" "<<fx_off<<" "<<falpha<<" "<<fParent<<endl;
      // if(fParent) cout<<dynamic_cast<RooHSEventsHistPDF*>(fParent)->GetRootHist()<<endl;
      arg=arg-offset;
      fx_off->setVal(arg);
      falpha->setVal(Double_t(alpha));
      //cout<<"DEBUG RooHSEventsHistPDF::evaluateMC "<<fHist->weight(RooArgSet(*fx_off,*falpha),1,kFALSE)<<" "<<arg<<" "<<mcx<<endl;
      return  fHist->weight(RooArgSet(*fx_off,*falpha),1,kFALSE);

  
    }

    // Bool_t RooHSEventsHistPDF::SetEvTree(TChain* tree,TString cut,Long64_t ngen){

    //   return RooHSEventsHistPDF::SetEvTree(tree->CloneTree(),cut,ngen); //standard intilisation
    // }
    Bool_t RooHSEventsHistPDF::SetEvTree(TTree* tree,TString cut,TTree* MCGenTree){
      Bool_t OK=RooHSEventsPDF::SetEvTree(tree,cut,MCGenTree); //standard intilisation
      //Now also need to create underlying HistPdf
      if(!fHist)CreateHistPdf();
      if(OK&&fHist) return kTRUE;
      else return kFALSE;
    }
    void RooHSEventsHistPDF::CreateHistPdf(){
      cout<<"      RooHSEventsHistPDF::CreateHistPdf()  "<<endl;
      fConstInt=fNTreeEntries;
      //Essentially cache the PDF into a histogram
      //x-axis = variable
      //y-axis = smearing parameter alpha


      //For each bin make PDF of var convoluted with gaussian
      Long64_t NFT=fNTreeEntries;
      auto* his1=new TH1D("his1D","his1D",fRHist->GetXaxis()->GetNbins(),fRHist->GetXaxis()->GetXmin(),fRHist->GetXaxis()->GetXmax());
      his1->SetDirectory(nullptr);
  
      //create 1D template of x variable
      for(Int_t itr=0;itr<NFT;itr++){//loop over events tree
	// fEvTree->GetEntry(itr);
	// Double_t tvar=fMCVar[0];
	fTreeEntry=itr;
	Double_t tvar=fvecReal[fTreeEntry*fNvars+0];
	his1->Fill(tvar,GetIntegralWeight(itr));
      }

      //Could be problems if weights result in -ve bins...
      for(Int_t ix=0;ix<his1->GetEntries();ix++)
	if(his1->GetBinContent(ix)<0){
	  cout<<" RooHSEventsHistPDF::CreateHistPdf() weights have resulted in some -ve bins. Will try averaging with 2 nearest nieghbours, if still -ve will set to 0. You may want to try using less bins for your histogram PDF via RooRealVar::setBins()"<<endl;
	  Double_t b0,b1,b2;
	  if(ix==1)
	    b0=his1->GetBinContent(ix); //lower edge
	  else
	    b0=his1->GetBinContent(ix-1);
	  b1=his1->GetBinContent(ix);
	  if(ix==his1->GetEntries())
	    b2=his1->GetBinContent(ix);//upper edge
	  else
	    b2=his1->GetBinContent(ix+1);
      
	  Double_t bmean=(b0+b1+b2)/3;
	  if(bmean<0) bmean=0;
	  his1->SetBinContent(ix-1,bmean);
	  his1->SetBinContent(ix,bmean);
	  his1->SetBinContent(ix+1,bmean);
	}
      if(fUseEvWeights) cout<<fEvWeights[0]<<" "<<fEvWeights[1]<<endl;
      his1->Smooth();
      //Fill first y bin of 2D hist (no smearing)
      for(Int_t jtemp=1;jtemp<=fRHist->GetNbinsX();jtemp++)//First alpha bin, no semaring!
	fRHist->Fill(fRHist->GetXaxis()->GetBinCenter(jtemp),fRHist->GetYaxis()->GetBinCenter(1),his1->GetBinContent(jtemp));
  
      //Loop over bins of smearing parameter
      TF1 gausnX("gausnX","gausn(0)",fRHist->GetXaxis()->GetXmin(),fRHist->GetXaxis()->GetXmax());
      for(Int_t ia=2;ia<=fRHist->GetNbinsY();ia++){//alpha bins
	Double_t vAlphb=fRHist->GetYaxis()->GetBinCenter(ia);
	Double_t vAlph=fRHist->GetYaxis()->GetBinLowEdge(ia);
	for(Int_t itemp=1;itemp<=fRHist->GetNbinsX();itemp++){//fill with Gaussian function
	  Double_t vari=fRHist->GetXaxis()->GetBinCenter(itemp);
	  Double_t NX=his1->GetBinContent(itemp);
	  if(!NX) continue;
	  gausnX.SetParameters(NX,vari,vAlphb);
	  //Add gaussian for each non xero bin with width = yaxis value
	  for(Int_t jtemp=1;jtemp<=fRHist->GetNbinsX();jtemp++){//fill with Gaussian function
	    Double_t varj=fRHist->GetXaxis()->GetBinCenter(jtemp);
	    fRHist->Fill(varj,vAlphb,gausnX.Eval(varj));
	  }
	}
      }
      //Cannot calculate covariance if pdf==0 for some events
      //Set every empty bin with very small value to prevent this
      for(Int_t ia=1;ia<=fRHist->GetNbinsY();ia++){//alpha bins
	Double_t max_cont=0;
	for(Int_t jtemp=1;jtemp<=fRHist->GetNbinsX();jtemp++){
	  Double_t cont=fRHist->GetBinContent(jtemp,ia);
	  if(cont>max_cont)max_cont=cont;
	}
	for(Int_t jtemp=1;jtemp<=fRHist->GetNbinsX();jtemp++){
	  Double_t cont=fRHist->GetBinContent(jtemp,ia);
	  if(cont==0) fRHist->SetBinContent(jtemp,ia,1E-10*max_cont);
	}
    
      }
      fRHist->Smooth();//some additional smoothing

      //Store max value of distributions for scaling around
      // fMean=fRHist->GetMean();
      Int_t bx,by,bz;
      fRHist->GetBinXYZ(fRHist->GetMaximumBin(),bx,by,bz);
      fVarMax=fRHist->GetXaxis()->GetBinCenter(bx);

      //import Root TH2 into RooFit 
      fHist = new RooDataHist(fRHist->GetName(),fRHist->GetName(),RooArgSet(*fx_off,*falpha),RooFit::Import(*fRHist));
      //cleanup
      delete his1;
    }
    Int_t RooHSEventsHistPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,const char* rangeName) const
    {
      return RooHSEventsPDF::getAnalyticalIntegral(allVars,analVars,rangeName);
    }
    Double_t RooHSEventsHistPDF::analyticalIntegral(Int_t code,const char* rangeName) const
    {
      if(code==1){//sum over 5000 samples of the PDF
	if(!CheckChange()) return fLast[0];
	Double_t integral=0;
	Int_t Nbins=5000;
	Double_t min=fRHist->GetXaxis()->GetXmin();
	Double_t max=fRHist->GetXaxis()->GetXmax();
	Double_t delta=(max-min)/Nbins;
	auto var=(RooRealVar*)(&(x.arg()));
	for(Int_t ie=1;ie<=Nbins;ie++){
	  Double_t val=min+delta*ie;
	  if(!(var->inRange(val,rangeName))) continue;
	  integral+=evaluateMC(val);
	}
	fLast[0]=integral*delta;
    
	return fLast[0];
      }
  
      return 1; 
    }
    void RooHSEventsHistPDF::ResetTree(){
  
      RooHSEventsPDF::ResetTree();
      if(fHist) {
	delete fHist;
	fHist=nullptr;
	fRHist->Reset();
      }
    }
  }
}

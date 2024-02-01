#include "PhotoTwoSpin0Amps.h"
#include "PhotoTwoSpin0AmpLoader.h"
#include "RooHSEventsPDF.h"

namespace HS{
  namespace FIT{


    std::string PhotoTwoSpin0Amps::ConfigurePWAs(){
      _Sum = "H_0_0_0*K_0*K_0";//constraint == H_0_0_0*Y_0_0*K_0   as K_0=Y_0_0=1/sqrt(4pi);
      //SUM(L[1|4]) => Sum over indice L between 1->4
      //H_0_L_M[0,-1,1] => create parameters e.g. H_0_1_1 initial value 0, range -1 to 1
      //M[0|4<L+1] sum over M up to 4 but <= L+1
      //COS2PHI => formula defined in  PolarisedSphHarmonicMoments
      //ReY_L_M(cosThGJ,Phi,Y_L_M)} => real part of Y^M_L a function of cosThGJ, Phi
      if(_OnlyEven==kFALSE){
	if(_UseI0){
	_Sum +=       Form("+ SUM(L[1|%d],M[0|%d<L+1]){H_0_L_M*K_L*ReY_L_M(%s,%s,Y_L_M)}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	}
	if(_UseI12){
	  _Sum +=       Form("+ SUM(L[0|%d],M[0|%d<L+1]){H_1_L_M*K_L*ReY_L_M(%s,%s,Y_L_M)*COS2PHI}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	  _Sum +=        Form("+ SUM(L[1|%d],M[1|%d<L+1]){H_2_L_M*K_L*ImY_L_M(%s,%s,Y_L_M)*SIN2PHI}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	}
	if(_UseI3){
	  _Sum+=        Form("+ SUM(L[1|%d],M[1|%d<L+1]){H_3_L_M*K_L*ImY_L_M(%s,%s,Y_L_M)*PCIRC}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());

	}
      }
      else{//increment L by 2
	if(_UseI0){
	  _Sum +=       Form("+ SUM(L[2|%d:2],M[0|%d<L+1]){H_0_L_M*K_L*ReY_L_M(%s,%s,Y_L_M)}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	}
	if(_UseI12){
	  _Sum +=       Form("+ SUM(L[0|%d:2],M[0|%d<L+1]){H_1_L_M*K_L*ReY_L_M(%s,%s,Y_L_M)*COS2PHI}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	  _Sum +=        Form("+ SUM(L[2|%d:2],M[1|%d<L+1]){H_2_L_M*K_L*ImY_L_M(%s,%s,Y_L_M)*SIN2PHI}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	}
	if(_UseI3){
	  _Sum+=        Form("+ SUM(L[2|%d],M[1|%d<L+1]){H_3_L_M*K_L*ImY_L_M(%s,%s,Y_L_M)*PCIRC}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());

	}

      }

      P2S0AmpLoader::LoadPartialWaves(*_Setup,_Lmax,_Mmax,_Nref,_OnlyEven,_NegativeM); //generate formula vars for Hs in terms of partial waves
      _IsAmplitudes=kTRUE;

      //remove parameter circual polarisation from formula
      //to prevent integral resumming
      if(_PolCirc.length()==0&&_UseI3){
	auto helcirc = Form("PCIRC*parPcirc");
	auto 	tsum = TString(_Sum.data());
	tsum.ReplaceAll("PCIRC",helcirc);
	_Sum = tsum;
      }
       //remove parameter linear polarisation from formula
      //to prevent integral resumming
      if(_PolLin.length()==0&&_UseI12){
	auto plin = Form("COS2PHI*parPlin");
	auto tsum = TString(_Sum.data());
	tsum.ReplaceAll("COS2PHI",plin);

	plin = Form("SIN2PHI*parPlin");
	tsum = TString(tsum.Data());
	tsum.ReplaceAll("SIN2PHI",plin);
	 
	_Sum = tsum;
      }


      return _Sum;
    }
    
    std::string  PhotoTwoSpin0Amps::ConfigureMoments(){
       //Now expansion in spherical harmonic moments
      //H_Alpha_L_M are fit parameters => the moments
      //K_L = =TMath::Sqrt(2*L.+1.)/TMath::Sqrt(4*TMath::Pi())
      //Y_L_M are spherical harmonic functions depending on CosTh and Phi
      //Note the factor tau from JPAC paper is cared for in
      //PolarisedSphHarmonicMoments function
      _Sum = Form("H_0_0_0[2]*K_0*K_0"); //constant == 2*Y_0_0*K_0 = 2*K_0*K_0  as Y_0_0=1/sqrt(4pi)=K_0;
      //SUM(L[1|4]) => Sum over indice L between 1->4
      //H_0_L_M[0,-1,1] => create parameters e.g. H_0_1_1 initial value 0, range -1 to 1
      //M[0|4<L+1] sum over M up to 4 but <= L+1
      //COS2PHI => formula defined in  PolarisedSphHarmonicMoments
      //ReY_L_M(cosThGJ,Phi,Y_L_M)} => real part of Y^M_L a function of cosThGJ, Phi
      if(_OnlyEven==kFALSE){
	if(_UseI0){
	_Sum +=       Form("+ SUM(L[1|%d],M[0|%d<L+1]){H_0_L_M[0,-2,2]*K_L*ReY_L_M(%s,%s,Y_L_M)}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	}
	if(_UseI12){
	  _Sum +=       Form("+ SUM(L[0|%d],M[0|%d<L+1]){H_1_L_M[0,-2,2]*K_L*ReY_L_M(%s,%s,Y_L_M)*COS2PHI}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	  _Sum+=        Form("+ SUM(L[1|%d],M[1|%d<L+1]){H_2_L_M[0,-2,2]*K_L*ImY_L_M(%s,%s,Y_L_M)*SIN2PHI}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	}
	if(_UseI3){
	  _Sum+=        Form("+ SUM(L[1|%d],M[1|%d<L+1]){H_3_L_M[0,-2,2]*K_L*ImY_L_M(%s,%s,Y_L_M)*PCIRC}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());

	}
      }
      else{//increment L by 2
	if(_UseI0){
	  _Sum +=       Form("+ SUM(L[2|%d:2],M[0|%d<L+1]){H_0_L_M[0,-2,2]*K_L*ReY_L_M(%s,%s,Y_L_M)}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	}
	if(_UseI12){
	  _Sum +=       Form("+ SUM(L[0|%d:2],M[0|%d<L+1]){H_1_L_M[0,-2,2]*K_L*ReY_L_M(%s,%s,Y_L_M)*COS2PHI}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	  _Sum+=        Form("+ SUM(L[2|%d:2],M[1|%d<L+1]){H_2_L_M[0,-2,2]*K_L*ImY_L_M(%s,%s,Y_L_M)*SIN2PHI}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	}
	if(_UseI3){
	  _Sum+=        Form("+ SUM(L[2|%d],M[1|%d<L+1]){H_3_L_M[0,-2,2]*K_L*ImY_L_M(%s,%s,Y_L_M)*PCIRC}",2*_Lmax,2*_Mmax,_DecayAngleCosTh.data(),_DecayAnglePhi.data());
	}

      }

      //remove parameter circual polarisation from formula
      //to prevent integral resumming
      if(_PolCirc.length()==0&&_UseI3){
	auto helcirc = Form("PCIRC*parPcirc");
	auto tsum = TString(_Sum.data());
	tsum.ReplaceAll("PCIRC",helcirc);
	_Sum = tsum;
      }
      //remove parameter linear polarisation from formula
      //to prevent integral resumming
      if(_PolLin.length()==0&&_UseI12){
	auto plin = Form("COS2PHI*parPlin");
	auto tsum = TString(_Sum.data());
	tsum.ReplaceAll("COS2PHI",plin);

	plin = Form("SIN2PHI*parPlin");
	tsum = TString(tsum.Data());
	tsum.ReplaceAll("SIN2PHI",plin);
	 
	_Sum = tsum;
      }
     
      _IsAmplitudes=kFALSE;
      return _Sum;

    }
  
    void PhotoTwoSpin0Amps::LoadModelPDF(Long64_t Nevents){
      //Nevents in case this is a toy generator
      //  ComponentsPdfParser  parser = PolarisedSphHarmonicMoments("AngularDist","trucosThGJ","truphiGJ","truphiCM","trucosThCM",_Lmax*2,0,_Lmax*2);
      ComponentsPdfParser  parser = PolarisedSphHarmonicMoments();
      //std::cout<<"PhotoTwoSpin0Amps::LoadModelPDF got a parser"<<std::endl;
      auto level = RooMsgService::instance().globalKillBelow();
      RooMsgService::instance().setGlobalKillBelow(RooFit::DEBUG) ;
      _Setup->ParserPDF(_Sum,parser);
      //dynamic_cast<RooHSEventsPDF*>(_Setup->WS().pdf("PWA"))->SetConstInt();
      // std::cout<<"PhotoTwoSpin0Amps::LoadModelPDF loaded parser"<<std::endl;
      _Setup->LoadSpeciesPDF(GetName(),Nevents); //100000 events
      RooMsgService::instance().setGlobalKillBelow(level) ;
 
       // std::cout<<"PhotoTwoSpin0Amps::LoadModelPDF loaded a pdf"<<std::endl;
     //the range of moments depends on L, so here reset the ranges
      auto& pars = _Setup->Parameters();

      // std::cout<<" PhotoTwoSpin0Amps::LoadModelPDF static_range_cast not available until 6.28 c++14"<<std::endl;exit(0);
      for(auto par:static_range_cast<RooRealVar *>(pars)){
      	TString parname = par->GetName();
      	if(parname.Contains("H_")){
      	  auto parts = parname.Tokenize('_');
      	  auto L = TString(parts->At(2)->GetName()).Atoi();
      	  par->setMin(-2/sqrt(2*L+1));par->setMax(2/sqrt(2*L+1));
      	}
      }
    }
    
    void PhotoTwoSpin0Amps::PrintModel(){
      const auto pars = _Setup->Parameters();
      cout<<"PhotoTwoSpin0Amps::PrintModel() "<<GetName()<<" : Parameters "<<endl;
      for(auto par:pars){
	auto rpar  = dynamic_cast<RooRealVar*>(par);
	if(rpar!=nullptr)cout<<rpar->GetName()<<" = "<<rpar->getVal()<<" +- "<<rpar->getError()<<endl;
      }
      //if(doingMoments==kFALSE){
      const auto forms = _Setup->Formulas();
      forms.Print();
      for(auto par:forms){
	auto rpar  = dynamic_cast<RooFormulaVar*>(par);
	if(rpar==nullptr) continue;
	cout<<rpar->GetName()<<" = "<<rpar->getVal()<<endl;
      }
      cout<<"PhotoTwoSpin0Amps::PrintModel() : Formulas "<<endl;
      for(auto par:forms){
	auto rpar  = dynamic_cast<RooFormulaVar*>(par);
	if(rpar==nullptr) continue;
	TString sform(rpar->formula().GetTitle());
	cout<<rpar->GetName()<<" = "<<sform<<endl;
	}
      cout<<"PhotoTwoSpin0Amps::PrintModel() : Summation "<<endl;
      cout<<GetSummation()<<endl;

    }
    
    
    ComponentsPdfParser PhotoTwoSpin0Amps::PolarisedSphHarmonicMoments(){
      TString name = GetName();
      TString cth = _DecayAngleCosTh;
      TString phi = _DecayAnglePhi;
      TString phiPol = _PolPhi;
      TString Pol = _PolLin;
      TString BeamHel = _BeamHelicity;
      TString PolCirc = _PolCirc;
      TString PolLin = _PolLin;
      if(PolLin.Length()==0)_constLinPol=kTRUE;
      
      Int_t Lmax = _Lmax*2; //L for moments is 2xL for partial waves
      Int_t Mmax = _Mmax*2;
      
      ComponentsPdfParser mp(name.Data());
   //From eqn a15a,b  https://arxiv.org/pdf/1906.04841.pdf 
      mp.AddFunctionTemplate("RooHSSphHarmonic","Y_*_0(*,*,*,*,1)");//tau=1 for M=0
      mp.AddFunctionTemplate("RooHSSphHarmonic","Y_*_*(*,*,*,*,2)");//tau=2
      mp.AddFunctionTemplate("RooHSSphHarmonicRe","ReY_*_*(*,*,Y_*_*)");
      mp.AddFunctionTemplate("RooHSSphHarmonicIm","ImY_*_*(*,*,Y_*_*)");

      
      if(_UseI3){

	if(_constLinPol==kFALSE){
	  if(PolCirc.Length()==0){
	    //make circ pol a parameter
	    mp.SetVars(Form("%s,%s,%s,%s,%s",cth.Data(),phi.Data(),phiPol.Data(),PolLin.Data(),BeamHel.Data()));
	    //mp.AddParameter("parPcirc[0,0,1]");
	  }
	  else{
	    mp.SetVars(Form("%s,%s,%s,%s,%s,%s",cth.Data(),phi.Data(),phiPol.Data(),PolLin.Data(),BeamHel.Data(),PolCirc.Data()));
	  }
	}
	else{
	  if(PolCirc.Length()==0){
	    mp.SetVars(Form("%s,%s,%s,%s",cth.Data(),phi.Data(),phiPol.Data(),BeamHel.Data()));
	  }
	  else{
	    mp.SetVars(Form("%s,%s,%s,%s,%s",cth.Data(),phi.Data(),phiPol.Data(),BeamHel.Data(),PolCirc.Data()));
	  }
	}
      }
      else{
	if(_constLinPol==kFALSE){
	  mp.SetVars(Form("%s,%s,%s,%s",cth.Data(),phi.Data(),phiPol.Data(),PolLin.Data()));
	}
	else{
	  mp.SetVars(Form("%s,%s,%s",cth.Data(),phi.Data(),phiPol.Data()));
	}
     }
      
        //Create   RooHSSphHarmonic functions, this only needs done here because RooHSSphHarmonic are not created by the parser otherwise...
      //Note extra factor 2 when M!=0 (tau in paper)
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1!0]){Y_L_M(%s,%s,L,M,2)}+SUM(L[0|%d]){Y_L_0(%s,%s,L,0,1)}",Lmax,0,Mmax,cth.Data(),phi.Data(),Lmax,cth.Data(),phi.Data())));


      if(_constLinPol==kFALSE){ //linpol observable
	mp.AddFormula(Form("COS2PHI=@%s[]*cos(2*@%s[])",Pol.Data(),phiPol.Data())); //Note +ve (due to -ve in Eqn A15a)
	mp.AddFormula(Form("SIN2PHI=-@%s[]*sin(2*@%s[])",Pol.Data(),phiPol.Data())); //Note -ve
      }
      else if(PolLin.Length()==0){//lin pol parameter
	mp.AddFormula(Form("COS2PHI=cos(2*@%s[])",phiPol.Data())); //Note +ve
	mp.AddFormula(Form("SIN2PHI=-sin(2*@%s[])",phiPol.Data())); //Note -ve
      }
      else{//lin pol constant
	mp.AddFormula(Form("COS2PHI=%s*cos(2*@%s[])",Pol.Data(),phiPol.Data())); //Note +ve
	mp.AddFormula(Form("SIN2PHI=-%s*sin(2*@%s[])",Pol.Data(),phiPol.Data())); //Note -ve
      }
      
      for(Int_t iL=0;iL<=Lmax;iL++)
	mp.AddConstant(Form("K_%d[%lf]",iL,TMath::Sqrt(2*iL+1.)/TMath::Sqrt(4*TMath::Pi())));


      //Check if we want circular polarisation as parameter or observable
      if(_UseI3){
	//case helicity is continous variable
	if(_HelicityIsCat==kFALSE){
	  if(PolCirc.Length()==0){
	    mp.AddFormula(Form("PCIRC=(@%s[]/TMath::Abs(@%s[]))",BeamHel.Data(),BeamHel.Data()));
	    mp.AddParameter("parPcirc[0.5,0,1]");
	  }
	  else{  
	    mp.AddFormula(Form("PCIRC=(@%s[]/TMath::Abs(@%s[]))*@%s[]",BeamHel.Data(),BeamHel.Data(),PolCirc.Data()));
	  }
	}
	//case helciity is category 
	if(_HelicityIsCat==kTRUE){
	  if(PolCirc.Length()==0){
	    mp.AddFormula(Form("PCIRC=@%s[]",BeamHel.Data()));
	    mp.AddParameter("parPcirc[0.5,0,1]");

	  }
	  else{
	  mp.AddFormula(Form("PCIRC=@%s[]*@%s[]",BeamHel.Data(),PolCirc.Data()));
	  }
	}
      }

     //Check if we want linear polarisation as parameter or observable
      if(_UseI12){
	//case helicity is continous variable
	if(PolLin.Length()==0){
	  mp.AddParameter("parPlin[0.5,0,1]");
	}
     }

      
      return mp;
      
    }


    
  }
}

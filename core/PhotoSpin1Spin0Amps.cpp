#include "PhotoSpin1Spin0Amps.h"
#include "PhotoSpin1Spin0AmpLoader.h"

namespace HS{
  namespace FIT{


    std::string PhotoSpin1Spin0Amps::ConfigurePWAs(){
      _Sum = "H_0_0_0_0_0*K_0*K_0";
      //   _Sum = "H_0_0_0_0_0*K_0*K_0";//constraint == H_0_0_0*Y_0_0*K_0   as K_0=Y_0_0=1/sqrt(4pi);
      //SUM(L[1|4]) => Sum over indice L between 1->4
      //H_0_L_M[0,-1,1] => create parameters e.g. H_0_1_1 initial value 0, range -1 to 1
      //M[0|4<L+1] sum over M up to 4 but <= L+1
      //COS2PHI => formula defined in  PolarisedSphHarmonicMoments
      //ReY_L_M(cosThGJ,Phi,Y_L_M)} => real part of Y^M_L a function of cosThGJ, Phi


 //To remove H0_0000 from sum, 2 cases
      //l=0, L=1,2,....
      //l=2, L=0,1,2,....
      	_Sum +=       Form("+ SUM(L[1|%d],l[0|0],m[0|0],M[0|%d<L+1]){H_0_l_m_L_M*K_L*K_l*ReDab_l_m_L_M}",2*_Jmax,2*_Mmax);

      	_Sum +=       Form("+ SUM(L[0|%d],l[2|2],m[0|2<L+1<l+1],M[0|%d<L+1]){H_0_l_m_L_M*K_L*K_l*ReDab_l_m_L_M}",2*_Jmax,2*_Mmax);


	_Sum +=       Form("+ SUM(L[0|%d],l[0|2:2],m[0|2<L+1<l+1],M[0|%d<L+1]){H_1_l_m_L_M*K_L*K_l*ReDab_l_m_L_M*COS2PHI}",2*_Jmax,2*_Mmax);

	_Sum +=       Form("+ SUM(L[0|%d],l[0|2:2],m[1|2<L+1<l+1],M[0|0<L+1]){H_2_l_m_L_M*K_L*K_l*ImDab_l_m_L_M*SIN2PHI}",2*_Jmax);

	_Sum +=       Form("+ SUM(L[0|%d],l[0|2:2],m[0|2<L+1<l+1],M[1|%d<L+1]){H_2_l_m_L_M*K_L*K_l*ImDab_l_m_L_M*SIN2PHI}",2*_Jmax,2*_Mmax);
        
      

    

	PS1S0AmpLoader::LoadPartialWaves(*_Setup,_Jmax,_Lmax,_Mmax,_mmax,_Nref,_OnlyEven,_NegativeM); //generate formula vars for Hs in terms of partial waves
      _IsAmplitudes=kTRUE;
      return _Sum;
    }
    
    std::string  PhotoSpin1Spin0Amps::ConfigureMoments(){
       //Now expansion in spherical harmonic moments
      //H_Alpha_L_M are fit parameters => the moments
      //K_L = =TMath::Sqrt(2*L.+1.)/TMath::Sqrt(4*TMath::Pi())
      //Y_L_M are spherical harmonic functions depending on CosTh and Phi
      //Note the factor tau from JPAC paper is cared for in
      //PolarisedSphHarmonicMoments function
      _Sum = "H_0_0_0_0_0[2]*K_0*K_0";  //constant  = 2*K_0*K_0 and D0_00 = 1
      //SUM(L[1|4]) => Sum over indice L between 1->4
      //H_0_L_M[0,-1,1] => create parameters e.g. H_0_1_1 initial value 0, range -1 to 1
      //M[0|4<L+1] sum over M up to 4 but <= L+1
      //COS2PHI => formula defined in  PolarisedSphHarmonicMoments
      //ReY_L_M(cosThGJ,Phi,Y_L_M)} => real part of Y^M_L a function of cosThGJ, Phi
      
      //To remove H0_0000 from sum, 2 cases
      //l=0, L=1,2,....
      //l=2, L=0,1,2,....
	_Sum +=       Form("+ SUM(L[1|%d],l[0|0],m[0|0],M[0|%d<L+1]){H_0_l_m_L_M[0,-2,2]*K_L*K_l*ReDab_l_m_L_M}",2*_Jmax,2*_Mmax);

	_Sum +=       Form("+ SUM(L[0|%d],l[2|2],m[0|2<L+1<l+1],M[0|%d<L+1]){H_0_l_m_L_M[0,-2,2]*K_L*K_l*ReDab_l_m_L_M}",2*_Jmax,2*_Mmax);

	_Sum +=       Form("+ SUM(L[0|%d],l[0|2:2],m[0|2<L+1<l+1],M[0|%d<L+1]){H_1_l_m_L_M[0,-2,2]*K_L*K_l*ReDab_l_m_L_M*COS2PHI}",2*_Jmax,2*_Mmax);

	_Sum+=        Form("+ SUM(L[0|%d],l[0|2:2],m[1|2<L+1<l+1],M[0|0<L+1]){H_2_l_m_L_M[0,-2,2]*K_L*K_l*ImDab_l_m_L_M*SIN2PHI}",2*_Jmax);

		_Sum+=        Form("+ SUM(L[0|%d],l[0|2:2],m[0|2<L+1<l+1],M[1|%d<L+1]){H_2_l_m_L_M[0,-2,2]*K_L*K_l*ImDab_l_m_L_M*SIN2PHI}",2*_Lmax,2*_Mmax);
    
  
      _IsAmplitudes=kFALSE;
      return _Sum;

    }
    
    void PhotoSpin1Spin0Amps::LoadModelPDF(Long64_t Nevents){
      //Nevents in case this is a toy generator
      //  ComponentsPdfParser  parser = PolarisedSphHarmonicMoments("AngularDist","trucosThGJ","truphiGJ","truphiCM","trucosThCM",_Lmax*2,0,_Lmax*2);
      ComponentsPdfParser  parser = PolWignerDFunctionProductMoments();
      //std::cout<<"PhotoTwoSpin0Amps::LoadModelPDF got a parser"<<std::endl;
      auto level = RooMsgService::instance().globalKillBelow();
      RooMsgService::instance().setGlobalKillBelow(RooFit::DEBUG) ;
      _Setup->ParserPDF(_Sum,parser);
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
    
    void PhotoSpin1Spin0Amps::PrintModel(){
      const auto pars = _Setup->Parameters();
      cout<<"PhotoSpin1Spin0Amps::PrintModel() "<<GetName()<<" : Parameters "<<endl;
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
      cout<<"PhotoSpin1Spin0Amps::PrintModel() : Formulas "<<endl;
      for(auto par:forms){
	auto rpar  = dynamic_cast<RooFormulaVar*>(par);
	if(rpar==nullptr) continue;
	TString sform(rpar->GetTitle());
	cout<<rpar->GetName()<<" = "<<sform<<endl;
	}
      cout<<"PhotoSpin1Spin0Amps::PrintModel() : Summation "<<endl;
      cout<<GetSummation()<<endl;

    }
    
    
    ComponentsPdfParser PhotoSpin1Spin0Amps::PolWignerDFunctionProductMoments(){
      TString name = GetName();
      TString cthGJ = _DecayAngleCosThGJ;
      TString phiGJ = _DecayAnglePhiGJ;
      TString cthHF = _DecayAngleCosThHF;
      TString phiHF = _DecayAnglePhiHF;
      TString phiPol = _PolPhi;
      TString Pol = _Polarisation;
      Int_t Lmax = _Lmax*2; //L for moments is 2xL for partial waves
      Int_t Mmax = _Mmax*2;
      Int_t lmax = _lmax*2;
      Int_t mmax = _mmax*2;
      Int_t Jmax = _Jmax*2;
      
      ComponentsPdfParser mp(name.Data());
     if(_constPol==kFALSE){
       mp.SetVars(Form("%s,%s,%s,%s,%s,%s",cthGJ.Data(),phiGJ.Data(),cthHF.Data(),phiHF.Data(),phiPol.Data(),Pol.Data()));
     }
     else{
       mp.SetVars(Form("%s,%s,%s,%s,%s",cthGJ.Data(),phiGJ.Data(),cthHF.Data(),phiHF.Data(),phiPol.Data()));
     }
      //Adapted From eqn a15a,b  https://arxiv.org/pdf/1906.04841.pdf 

      mp.AddFunctionTemplate("RooHSDWigner","Da_*_*_*(*,*,*,*,*)");
      mp.AddFunctionTemplate("RooHSDWigner","Db_*_*_*(*,*,*,*,*)");
      mp.AddFunctionTemplate("RooHSDWignerProduct","Dab_*_*_*_*(*,*)");
      mp.AddFunctionTemplate("RooHSDWignerProductRe","ReDab_*_*_*_*(*,*,*)");
      mp.AddFunctionTemplate("RooHSDWignerProductIm","ImDab_*_*_*_*(*,*,*)");

      //Create   RooHSSphHarmonic functions, this only needs done here because RooHSSphHarmonic are not created by the parser otherwise...
      //Note extra factor 2 when M!=0 (tau in paper)
         mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1]){Db_L_M_0(%s,%s,L,M,0)}",lmax,0,mmax,cthHF.Data(),phiHF.Data())));
         mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1],N[%d|%d>-L-1<L+1]){Da_L_M_N(%s,%s,L,M,N)}",Jmax,0,Mmax,0,mmax,cthGJ.Data(),phiHF.Data())));
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d<L+1],l[0|%d],m[%d|%d<L+1<l+1]){Dab_l_m_L_M(Da_L_M_m,Db_l_m_0)}",Jmax,0,Mmax,lmax,0,mmax)));

      if(_constPol==kFALSE){
	mp.AddFormula(Form("COS2PHI=@%s[]*cos(2*@%s[])",Pol.Data(),phiPol.Data())); //Note +ve
	mp.AddFormula(Form("SIN2PHI=-@%s[]*sin(2*@%s[])",Pol.Data(),phiPol.Data())); //Note -ve
      }
      else{
	mp.AddFormula(Form("COS2PHI=%s*cos(2*@%s[])",Pol.Data(),phiPol.Data())); //Note +ve
	mp.AddFormula(Form("SIN2PHI=%s*sin(2*@%s[])",Pol.Data(),phiPol.Data())); //Note -ve
      }
      for(Int_t iJ=0;iJ<=Jmax;iJ++)
	mp.AddConstant(Form("K_%d[%lf]",iJ,TMath::Sqrt(2*iJ+1.)/TMath::Sqrt(4*TMath::Pi())));
      for(Int_t il=0;il<=lmax;il++)
	mp.AddConstant(Form("K_%d[%lf]",il,TMath::Sqrt(2*il+1.)/TMath::Sqrt(4*TMath::Pi())));
     
      return mp;
      
    }


    
  }
}

#pragma once

#include "ComponentsPdfParser.h"
#include "TString.h"
//#include ""
//#include ""

namespace HS{
  namespace FIT{

    
    //e.g.PolarisedSphHarmonicMoments("Moments","CosTh","Phi","PolPhi","Pol",3,2);
    ComponentsPdfParser SphHarmonicMoments(TString name,TString cth,TString phi,Int_t Lmax,Int_t Mmin,Int_t Mmax){
      ComponentsPdfParser mp(name.Data());
      mp.SetVars(Form("%s,%s",cth.Data(),phi.Data()));
      
      //From eqn a15a,b  https://arxiv.org/pdf/1906.04841.pdf 
      mp.AddFunctionTemplate("RooHSSphHarmonic","Y_*_0(*,*,*,*,1)");//tau=1 for M=0
      mp.AddFunctionTemplate("RooHSSphHarmonic","Y_*_*(*,*,*,*,2)");//tau=2
      mp.AddFunctionTemplate("RooHSSphHarmonicRe","ReY_*_*(*,*,Y_*_*)");
      //mp.AddFunctionTemplate("RooHSSphHarmonicIm","ImY_*_*(*,*,Y_*_*)");

      //Create   RooHSSphHarmonic functions, this only needs done here because RooHSSphHarmonic are not created by the parser otherwise...
      //mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[1|%d],M[%d|%d>-L-1<L+1]){Y_L_M(%s,%s,L,M,2)}+SUM(L[0|%d]){Y_L_0(%s,%s,L,0,1)}",Lmax,Mmin,Mmax,cth.Data(),phi.Data(),0,cth.Data(),phi.Data())));
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1!0]){Y_L_M(%s,%s,L,M,2)}+SUM(L[0|%d]){Y_L_0(%s,%s,L,0,1)}",Lmax,Mmin,Mmax,cth.Data(),phi.Data(),Lmax,cth.Data(),phi.Data())));

      for(Int_t iL=0;iL<=Lmax;iL++)
	mp.AddConstant(Form("K_%d[%lf]",iL,TMath::Sqrt(2*iL+1.)/TMath::Sqrt(4*TMath::Pi())));
      

       return mp;
      
    }

   //e.g.PolarisedSphHarmonicMoments("Moments","CosTh","Phi","PolPhi","Pol",3,2);
   ComponentsPdfParser PolarisedSphHarmonicMoments(TString name,TString cth,TString phi,TString phiPol,TString Pol,Int_t Lmax,Int_t Mmin,Int_t Mmax){
      ComponentsPdfParser mp(name.Data());
      mp.SetVars(Form("%s,%s,%s,%s",cth.Data(),phi.Data(),phiPol.Data(),Pol.Data()));
      
      //From eqn a15a,b  https://arxiv.org/pdf/1906.04841.pdf 
      mp.AddFunctionTemplate("RooHSSphHarmonic","Y_*_0(*,*,*,*,1)");//tau=1 for M=0
      mp.AddFunctionTemplate("RooHSSphHarmonic","Y_*_*(*,*,*,*,2)");//tau=2
      mp.AddFunctionTemplate("RooHSSphHarmonicRe","ReY_*_*(*,*,Y_*_*)");
      mp.AddFunctionTemplate("RooHSSphHarmonicIm","ImY_*_*(*,*,Y_*_*)");

      //Create   RooHSSphHarmonic functions, this only needs done here because RooHSSphHarmonic are not created by the parser otherwise...
      //Note extra factor 2 when M!=0 (tau in paper)
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1!0]){Y_L_M(%s,%s,L,M,2)}+SUM(L[0|%d]){Y_L_0(%s,%s,L,0,1)}",Lmax,Mmin,Mmax,cth.Data(),phi.Data(),Lmax,cth.Data(),phi.Data())));


      
      mp.AddFormula(Form("COS2PHI=@%s[]*cos(2*@%s[])",Pol.Data(),phiPol.Data()));
      mp.AddFormula(Form("SIN2PHI=-@%s[]*sin(2*@%s[])",Pol.Data(),phiPol.Data()));
       //mp.AddFormula(Form("SIN2PHI=@%s[]*sin(2*@%s[])",Pol.Data(),phiPol.Data()));

      for(Int_t iL=0;iL<=Lmax;iL++)
	mp.AddConstant(Form("K_%d[%lf]",iL,TMath::Sqrt(2*iL+1.)/TMath::Sqrt(4*TMath::Pi())));
      
      //////////////////
      //string sum;
      //sum+="H_0_0_0[1]"; //constant == 1
      //sum+="+SUM(L[1|4],M[0|4<L+1]){H_0_L_M[0,-1,1]*Y_L_M_Re(CosTheta,Phi,Y_L_M)}";
      //sum+="+SUM(L[1|4],M[0|4<L+1]){H_1_L_M[0,-1,1]*Y_L_M_Re(CosTheta,Phi,Y_L_M)*COS2PHI}";
      //sum+="+SUM(L[1|4],M[1|4<L+1]){H_2_L_M[0,-1,1]*Y_L_M_Im(CosTheta,Phi,Y_L_M)*SIN2PHI}"; //NB M=0 =>H2 =0 eqn D14
      
      return mp;
      
    }

     //e.g.ComponentsPdfParser  parser=HS::FIT::WignerDFunctionMoments("Moments","Theta","CosTh","Phi",4,0,4,0,4);
    //Note if cth="" then Theta assumed to be actual data not formula,
    //if cth has a string then this is assumed to be the tree data and Theta caclulated via formula

    ComponentsPdfParser WignerDFunctionMoments(TString name,TString theta,TString cth,TString phi,Int_t Lmax,Int_t Mmin,Int_t Mmax,Int_t Nmin,Int_t Nmax){
      
      ComponentsPdfParser mp(name.Data());
      if(cth.Length()!=0) //using formula for Theta calcualted from cosTh
	mp.SetVars(Form("%s,%s",cth.Data(),phi.Data()));
      else
	mp.SetVars(Form("%s,%s",theta.Data(),phi.Data()));
      
      mp.AddFunctionTemplate("RooHSDWigner","D_*_*_0(*,*,*,*,*)");//tau=1 for M=0
      mp.AddFunctionTemplate("RooHSDWigner","D_*_*_*(*,*,*,*,*)");//tau=1 as Sum M>-L +L
      mp.AddFunctionTemplate("RooHSDWignerRe","ReD_*_*_*(*,*,D_*_*_*)");
      //mp.AddFunctionTemplate("RooHSWignerDIm","ImD_*_*(*,*,D_*_*_*)");


      //TString theta= Form("ThetaFrom%s",cth.Data());
      //Create   RooHSWignerD functions, this only needs done here because RooHSWignerD are not created by the parser otherwise...
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1],N[%d|%d>-L-1<L+1]){D_L_M_N(%s,%s,L,M,N)}+SUM(L[0|%d]){D_L_0_0(%s,%s,L,0,0)}",Lmax,Mmin,Mmax,Nmin,Nmax,theta.Data(),phi.Data(),Lmax,theta.Data(),phi.Data())));

      //     mp.AddFormula(Form("%s=TMath::ACos(@%s[])",theta.Data(),cth.Data()));


      for(Int_t iL=0;iL<=Lmax;iL++)
	mp.AddConstant(Form("K_%d[%lf]",iL,(2*iL+1.)/(4*TMath::Pi())));
      
      
      return mp;
      
    }
    /////////////////////////////////////////////////////////////////////  
    ComponentsPdfParser WignerDFunctionProductMoments(TString name,TString theta1,TString cth1,TString phi1,TString theta2,TString cth2,TString phi2,Int_t Lmax1,Int_t Mmin1,Int_t Mmax1,Int_t Nmin1,Int_t Nmax1,Int_t Lmax2){

     Int_t Mmin2=Nmin1;
     Int_t Mmax2=Nmax1;


     ComponentsPdfParser mp(name.Data());

     if(cth1.Length()!=0)
	mp.SetVars(Form("%s,%s,%s,%s",cth1.Data(),phi1.Data(),cth2.Data(),phi2.Data()));
      else
	mp.SetVars(Form("%s,%s,%s,%s",theta1.Data(),phi1.Data(),theta2.Data(),phi2.Data()));
	
      
      mp.AddFunctionTemplate("RooHSDWigner","Da_*_*_*(*,*,*,*,*)");
      mp.AddFunctionTemplate("RooHSDWigner","Db_*_*_*(*,*,*,*,*)");
      mp.AddFunctionTemplate("RooHSDWignerProduct","Dab_*_*_*_*(*,*)");

      
  
      //create the required RooHSWignerD objects
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1]){Db_L_M_0(%s,%s,L,M,0)}",Lmax2,Mmin2,Mmax2,theta2.Data(),phi2.Data())));
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1],N[%d|%d>-L-1<L+1]){Da_L_M_N(%s,%s,L,M,N)}",Lmax1,Mmin1,Mmax1,Nmin1,Nmax1,theta1.Data(),phi1.Data())));


      //Create   RooHSWignerDProduct functions, this only needs done here because they are not created by the parser otherwise...
      auto pdfstring=mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d<L+1],l[0|%d],m[%d|%d0<L+1<l+1]){Dab_l_m_L_M(Da_L_M_m,Db_l_m_0)}",Lmax1,Mmin1,Mmax1,Lmax2,Nmin1,Nmax1)));
      
      auto LMAX=Lmax1>Lmax2?Lmax1:Lmax2;
      for(Int_t iL=0;iL<=LMAX;iL++)
	mp.AddConstant(Form("K_%d[%lf]",iL,(2*iL+1.)/(4*TMath::Pi())));
      
      return mp;
      
    }
  /////////////////////////////////////////////////////////////////////  
    ComponentsPdfParser WignerDFunctionProductMomentsSym(TString name,TString theta1,TString cth1,TString phi1,TString theta2,TString cth2,TString phi2,TString theta1sym,TString cth1sym,TString phi1sym,TString theta2sym,TString cth2sym,TString phi2sym,Int_t Lmax1,Int_t Mmin1,Int_t Mmax1,Int_t Nmin1,Int_t Nmax1,Int_t Lmax2){

     Int_t Mmin2=Nmin1;
     Int_t Mmax2=Nmax1;


     ComponentsPdfParser mp(name.Data());
 
     if(cth1.Length()!=0)
	mp.SetVars(Form("%s,%s,%s,%s,%s,%s,%s,%s",cth1.Data(),phi1.Data(),cth2.Data(),phi2.Data(),cth1sym.Data(),phi1sym.Data(),cth2sym.Data(),phi2sym.Data()));
      else
	mp.SetVars(Form("%s,%s,%s,%s,%s,%s,%s,%s",theta1.Data(),phi1.Data(),theta2.Data(),phi2.Data(),theta1sym.Data(),phi1sym.Data(),theta2sym.Data(),phi2sym.Data()));
	
      
      mp.AddFunctionTemplate("RooHSDWigner","Da_*_*_*(*,*,*,*,*)");
      mp.AddFunctionTemplate("RooHSDWigner","Db_*_*_*(*,*,*,*,*)");
      mp.AddFunctionTemplate("RooHSDWignerProduct","Dab_*_*_*_*(*,*)");      
  
      //create the required RooHSWignerD objects
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1]){Db_L_M_0(%s,%s,L,M,0)}",Lmax2,Mmin2,Mmax2,theta2.Data(),phi2.Data())));
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1],N[%d|%d>-L-1<L+1]){Da_L_M_N(%s,%s,L,M,N)}",Lmax1,Mmin1,Mmax1,Nmin1,Nmax1,theta1.Data(),phi1.Data())));

      //symmetrised functions
      mp.AddFunctionTemplate("RooHSDWigner","symDa_*_*_*(*,*,*,*,*)");
      mp.AddFunctionTemplate("RooHSDWigner","symDb_*_*_*(*,*,*,*,*)");
      mp.AddFunctionTemplate("RooHSDWignerProduct","symDab_*_*_*_*(*,*)");

      //create the required symmetrised RooHSWignerD objects
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1]){symDb_L_M_0(%s,%s,L,M,0)}",Lmax2,Mmin2,Mmax2,theta2sym.Data(),phi2sym.Data())));
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1],N[%d|%d>-L-1<L+1]){symDa_L_M_N(%s,%s,L,M,N)}",Lmax1,Mmin1,Mmax1,Nmin1,Nmax1,theta1sym.Data(),phi1sym.Data())));



      //Create   RooHSWignerDProduct functions, this only needs done here because they are not created by the parser otherwise...
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d<L+1],l[0|%d],m[%d|%d0<L+1<l+1]){Dab_l_m_L_M(Da_L_M_m,Db_l_m_0)}",Lmax1,Mmin1,Mmax1,Lmax2,Nmin1,Nmax1)));
      //symmeterised versions
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d<L+1],l[0|%d],m[%d|%d0<L+1<l+1]){symDab_l_m_L_M(symDa_L_M_m,symDb_l_m_0)}",Lmax1,Mmin1,Mmax1,Lmax2,Nmin1,Nmax1)));

      
      auto LMAX=Lmax1>Lmax2?Lmax1:Lmax2;
      for(Int_t iL=0;iL<=LMAX;iL++)
	mp.AddConstant(Form("K_%d[%lf]",iL,(2*iL+1.)/(4*TMath::Pi())));
      
      return mp;
      
    }

  }
}

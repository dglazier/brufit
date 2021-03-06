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

   
       return std::move(mp);
      
    }

   ComponentsPdfParser PolarisedSphHarmonicMoments(TString name,TString cth,TString phi,TString phiPol,TString Pol,Int_t Lmax,Int_t Mmin,Int_t Mmax){
      ComponentsPdfParser mp(name.Data());
      mp.SetVars(Form("%s,%s,%s,%s",cth.Data(),phi.Data(),phiPol.Data(),Pol.Data()));
      
      //From eqn a15a,b  https://arxiv.org/pdf/1906.04841.pdf 
      mp.AddFunctionTemplate("RooHSSphHarmonic","Y_*_0(*,*,*,*,1)");//tau=1 for M=0
      mp.AddFunctionTemplate("RooHSSphHarmonic","Y_*_*(*,*,*,*,2)");//tau=2
      mp.AddFunctionTemplate("RooHSSphHarmonicRe","ReY_*_*(*,*,Y_*_*)");
      mp.AddFunctionTemplate("RooHSSphHarmonicIm","ImY_*_*(*,*,Y_*_*)");

      //Create   RooHSSphHarmonic functions, this only needs done here because RooHSSphHarmonic are not created by the parser otherwise...
      mp.ConstructPDF(mp.ReplaceSummations(Form("SUM(L[0|%d],M[%d|%d>-L-1<L+1!0]){Y_L_M(%s,%s,L,M,2)}+SUM(L[0|%d]){Y_L_0(%s,%s,L,0,1)}",Lmax,Mmin,Mmax,cth.Data(),phi.Data(),Lmax,cth.Data(),phi.Data())));


      
      mp.AddFormula(Form("COS2PHI=@%s[]*cos(2*@%s[])",Pol.Data(),phiPol.Data()));
      mp.AddFormula(Form("SIN2PHI=-@%s[]*sin(2*@%s[])",Pol.Data(),phiPol.Data()));

      //////////////////
      //string sum;
      //sum+="H_0_0_0[1]"; //constant == 1
      //sum+="+SUM(L[1-4],M[0-L]){H_0_L_M[0,-1,1]*Y_L_M_Re(CosTheta,Phi,Y_L_M)}";
      //sum+="+SUM(L[0-4],M[0-L]){H_1_L_M[0,-1,1]*Y_L_M_Re(CosTheta,Phi,Y_L_M)*COS2PHI}";
      //sum+="+SUM(L[1-4],M[1-L]){H_2_L_M[0,-1,1]*Y_L_M_Im(CosTheta,Phi,Y_L_M)*SIN2PHI}"; //NB M=0 =>H2 =0 eqn D14
      
      // mp.ConstructPDF(sum);
      return std::move(mp);
      
    }

    
  }
}

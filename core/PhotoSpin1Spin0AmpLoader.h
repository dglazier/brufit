#pragma once

#include <TMath.h>
#include <TError.h>
#include <Math/SpecFuncMathMore.h>
#include "Setup.h"

namespace PS1S0AmpLoader{
  using namespace std;
  using namespace HS::FIT;
  
  //Declare functions
  //ClebschGordan(ilpr,L,il,impr,M,im);
  ///////////////////////////////////////////////////////////////////////////
  inline Double_t ClebschGordan(Int_t l1, Int_t l2, Int_t l3, Int_t m1, Int_t m2, Int_t m3){
    using TMath::Power;
    using ROOT::Math::wigner_3j;
    using TMath::Sqrt;
  
    return (Power(-1.0,l1-l2+m3)*Sqrt(2.0*l3+1.0)*wigner_3j(2*l1,2*l2,2*l3,2*m1,2*m2,-2*m3));
     
  }


  
 ///////////////////////////////////////////////////////////////////////////
inline TString BruTerm(int reflsign,double factor,int j,int l,int m,int jpr,int lpr,int mpr,int alpha,int negm=1){
  if(negm==0&&(m<0||mpr<0)) return "XXX";
  char refl='a';
  if(reflsign==-1)refl='b';

  int reflfactor = 1;
  if(alpha==1||alpha==2) reflfactor = reflsign;



  if(j==jpr&&l==lpr&&m==mpr)return Form("%d*%lf*(@%c_%d_%d_%d[]*@%c_%d_%d_%d[])",reflfactor,factor,refl,j,l,m,refl,jpr,lpr,mpr);
    
  TString func="cos"; //for real part
  // if(alpha==2) func="sin"; //for imaginery part
  
  return Form("%lf*(@%c_%d_%d_%d[]*@%c_%d_%d_%d[]*%s(@%cphi_%d_%d_%d[]-@%cphi_%d_%d_%d[]))",reflfactor*factor,refl,j,l,m,refl,jpr,lpr,mpr,func.Data(),refl,j,l,m,refl,jpr,lpr,mpr);
}


 ///////////////////////////////////////////////////////////////////////////
  TString Simplify(TString moment){
    //Look for like terms and sum together (some may cancel)
   
    auto expr = TString(moment(moment.First("=")+1,moment.Length()));
    //  cout<<expr<<endl;
    auto terms = (expr.Tokenize('+'));

    auto MatchTerm = [](TString& aname,TString other){
						      //Check if other is matched to aname
						      auto isNeg = false;
						      auto isOtherNeg = false;
						      auto name = aname;
						      bool isChanging=false;
						      if(name[0]=='-' ){
							isNeg=true;
							name = TString(name(1,name.Length()));
						      }
						      if(other.Contains(name)){
							if(other[0]=='-' ){
							  isOtherNeg=true;
							  other = TString(other(1,other.Length()));
							}

							if(isOtherNeg==isNeg){//same sign
							  name.Prepend("2*"); //multiply by 2
							  if(isNeg==true) name.Prepend("-");
							}
							else{//different signs, term cancel!
							  name="";
							}
							isChanging=true;
							aname = name;
						      }
						      return isChanging;
    };
    auto MatchTermUnOrder = [](TString& aname,TString other){
							     //Check if other is matched to aname
							     auto isNeg = false;
							     auto isOtherNeg = false;
							     auto name = aname;
							     bool isChanging=false;
							     if(name[0]=='-' ){
							       isNeg=true;
							       name = TString(name(1,name.Length()));
							     }
							     
							     auto firstAt = name.First('@');
							     auto firstClose = name.First(']');
							     auto deltaBra = name.First('[')-firstAt+2;
							     auto firstStar=name.First('*');
							     //auto secondAt =
							     TString afirst = name(firstAt,deltaBra);
							     TString asecond = name(firstClose+2,deltaBra);
							     TString numerical = name(0,firstStar);
							     if(afirst==asecond) //already got by MatchTerm
							       return isChanging;

							     auto reverse = asecond+"*"+afirst;
							    
							     if(other.Contains(reverse)&&other.Contains(numerical)){
							       if(other[0]=='-' ){
								 isOtherNeg=true;
								 other = TString(other(1,other.Length()));
							       }

							       if(isOtherNeg==isNeg){//same sign
								 name.Prepend("2*"); //multiply by 2
								 if(isNeg==true) name.Prepend("-");
							       }
							       else{//different signs, term cancel!
								 name="";
							       }
							       isChanging=true;
							       aname = name;
							     }
							     return isChanging;
    };

    TString newExpr="";
    //Loop over all the terms in the H sum
    for(int i=0;i<terms->GetEntries();++i){
      auto name = TString(terms->At(i)->GetName());
      if(name=="") continue;
     
      for(int j=i+1;j<terms->GetEntries();++j){
	//Loop over all the other terms in the H sum
 	auto othername = TString(terms->At(j)->GetName());
	// see if term i matches term j
	auto changed = MatchTerm(name,othername);
	if(changed==false) changed =  MatchTermUnOrder(name,othername);

	if(changed){
	  //there was a match so remove j as i now includes it
	  dynamic_cast<TObjString*>(terms->At(j))->SetString("");
	  break;
	}
      }
     
      if(name!=""){//if i has not cancelled j, add i back into string
	           //(may also include j part now)
	newExpr+=name;
	newExpr+='+';

      }
    }	
    //remove the last +
    newExpr.Remove(TString::kTrailing,'+');

    if(newExpr!="") //Add back HXYZ=
      newExpr.Prepend(moment(0,moment.First("=")+1));
    else{
      newExpr.Prepend(moment(0,moment.First("=")+1));
     }
    return newExpr;
  }

  //////////////////////////////////////////////////////////////////////
 inline TString CGMatrixReflectivity(TString Moment,Int_t J,Int_t M, Int_t S, Int_t Lambda,Int_t _Jmax, Int_t _Lmax, Int_t _Mmax, Int_t _mmax,Int_t alpha,Bool_t negm=1){
  
  //Sum_l,l',m'm' CLM CL0 rho_lm.rho_l'm'*
  //where CLM and CL0 are clebsch Gordan coeficients
  
  //sdm is hermitian so rho_(lm)(lm')=rho_(lm')(lm)*

  Int_t col =0;

  Int_t is = 1;
  Int_t ispr = 1;
  Int_t nwaves=1; //S0
  

   for(Int_t i=1;i<=_Jmax;i++)
    {
      Int_t nwaves_m=0;
      Int_t nwaves_l=3;
      
      Int_t m=i;
      if(i>_Mmax) m = _Mmax; //restrict m terms
      nwaves_m=m +m*negm+ 1;
      
      nwaves+=nwaves_l*nwaves_m;
    }

  Int_t nsdmes=nwaves*(nwaves+1)/2;
  Int_t nmoments = (2*_Lmax+1)*(_Lmax+1);
  
  cout<<"NWAVES "<<nwaves<<"  NREAL "<<nwaves*2-1<<"  NSDMES "<<nsdmes<<" NMOMENTS "<<nmoments<<endl;

  //  exit(0);

  vector<Double_t> Are(nsdmes);
  vector<Double_t> Aim(nsdmes);
  
  int sdmeRow=0;
  int sdmeCol=0;

  TString sumTotal(Form("%s_%d_%d_%d_%d_%d=",Moment.Data(),alpha,S,Lambda,J,M));
  
  for(Int_t ij=0;ij<=_Jmax;ij++){
    for(Int_t il=std::max(0,ij-1);il<=std::min(_Lmax,ij+1);il++){
      for(Int_t im=-ij;im<=ij;im++){
	if(TMath::Abs(im)>_Mmax) continue;
	for(Int_t ilambda=-is;ilambda<=is;ilambda++){
	  if(TMath::Abs(ilambda)>_mmax) continue;
	  for(Int_t ijpr=0;ijpr<=_Jmax;ijpr++){
	    for(Int_t ilpr=std::max(0,ijpr-1);ilpr<=std::min(_Lmax,ijpr+1);ilpr++){
	      for(Int_t impr=-ijpr;impr<=ijpr;impr++){
		if(TMath::Abs(impr)>_Mmax) continue;
		for(Int_t ilambdapr=-ispr;ilambdapr<=ispr;ilambdapr++){
		  if(TMath::Abs(ilambdapr)>_mmax) continue; 
		  //  if((ilambda+ilambdapr)!=Lambda) continue;
		  //if((im+impr)!=M) continue;
	  

		  auto CJM = ClebschGordan(ij,J,ijpr,im,M,impr);
		  auto CJLambda = ClebschGordan(ij,J,ijpr,ilambda,Lambda,ilambdapr);
		  auto CS0 = ClebschGordan(is,S,ispr,0,0,0);
		  auto CSLambda = ClebschGordan(is,S,ispr,ilambda,Lambda,ilambdapr);
		  auto Cllambdapr = ClebschGordan(ilpr,ispr,ijpr,0,ilambdapr,ilambdapr);
		  auto Cllambda = ClebschGordan(il,is,ij,0,ilambda,ilambda);
		  auto factor = TMath::Sqrt((2.*ilpr+1)*(2.*il+1))* TMath::Sqrt((2.*is+1)/(2.*ispr+1))/(2*ijpr+1);
		  
		  auto conjCJM = ClebschGordan(ijpr,J,ij,impr,M,im);
		  auto conjCJLambda = ClebschGordan(ijpr,J,ij,ilambdapr,Lambda,ilambda);
		  auto conjCS0 = ClebschGordan(ispr,S,is,0,0,0);
		  auto conjCSLambda = ClebschGordan(ispr,S,is,ilambdapr,Lambda,ilambda);
		  auto conjCllambda = ClebschGordan(ilpr,ispr,ijpr,0,ilambdapr,ilambdapr);
		  auto conjCllambdapr = ClebschGordan(il,is,ij,0,ilambda,ilambda);
		  auto conjfactor = TMath::Sqrt((2.*ilpr+1)*(2.*il+1))* TMath::Sqrt((2.*ispr+1)/(2.*is+1))/(2*ij+1);
		  
		  int mmprimesign = 1;
		  if(TMath::Abs((im-impr)%2) == 1) mmprimesign=-1;
		  int msign = 1;
		  if(TMath::Abs((im)%2) == 1) msign=-1;
		  int mprimesign = 1;
		  if(TMath::Abs((impr)%2) == 1) mprimesign=-1;
		  int naturality = 1;
		  if(TMath::Abs((ij-il)%2) == 1) naturality =-1;
		  int naturalityprime = 1;
		  if(TMath::Abs((ijpr-ilpr)%2) == 1) naturalityprime=-1;
		  int nnprime =1;
		  if(TMath::Abs((ij+il+ijpr+ilpr)%2)==1) nnprime=-1;
		  //cout<<iL<<" "<<im<<" "<<ilpr<<" "<<impr<<" factor "<<factor<<endl;
		  auto ccfactor=CJM*CJLambda*CS0*CSLambda*Cllambdapr*Cllambda*factor; //+ conjCJM*conjCJLambda*conjCS0*conjCSLambda*conjCllambdapr*conjCllambda*conjfactor;//CCFa1.a2*+conjCCFa2*a1*
	
		  //factor*=CM*C0;
		  if(ccfactor!=0){
		    // if(sdmeRow<=sdmeCol){
		    // if(ij==ijpr&&im==impr&&il==ilpr&&is==ispr)ccfactor/=2;
		      //cout<<il<<" "<<im<<" "<<ilpr<<" "<<impr<<" factor "<<factor<<" row "<<sdmeRow<< " col "<<sdmeCol<<endl;
    
		      TString sum =""; 
		      int refl=1;

		      // std::cout<<"j: "<<ij<<"  l: "<<il<<"  m: "<<im<<"  lambda: "<<ilambda<<"  "<<"jpr: "<<ijpr<<"  lpr: "<<ilpr<<"  mpr: "<<impr<<"  lambdapr: "<<ilambdapr<<"  "<<CJM<<"  "<<CJLambda<<"  "<<CS0<<"  "<<CSLambda<<"  "<<Cllambdapr<<"  "<<Cllambda<<std::endl;
		      
		       if(alpha==0){
			 sum+= BruTerm(refl,ccfactor,ij,il,im,ijpr,ilpr,impr,alpha,negm);
			   sum+= " + ";
			   sum+= BruTerm(refl,nnprime*mmprimesign*ccfactor,ij,il,-im,ijpr,ilpr,-impr,alpha,negm);
			   
			   int refl=-1;
			   sum+= " + ";
			   sum+= BruTerm(refl,ccfactor,ij,il,im,ijpr,ilpr,impr,alpha,negm);
			   sum+= " + " ;
			   sum+= BruTerm(refl,nnprime*mmprimesign*ccfactor,ij,il,-im,ijpr,ilpr,-impr,alpha,negm);
			 }
			 else if(alpha==1){
			   sum+= BruTerm(refl,naturality*msign*ccfactor,ij,il,-im,ijpr,ilpr,impr,alpha,negm);
			   sum+= " + ";
			   sum+= BruTerm(refl,naturalityprime*mprimesign*ccfactor,ij,il,im,ijpr,ilpr,-impr,alpha,negm);
			   
			   int refl=-1;
			   sum+= " + ";
			   sum+= BruTerm(refl,naturality*msign*ccfactor,ij,il,-im,ijpr,ilpr,impr,alpha,negm);
			   sum+= " + ";
			   sum+= BruTerm(refl,naturalityprime*mprimesign*ccfactor,ij,il,im,ijpr,ilpr,-impr,alpha,negm);
			 }
			 else if(alpha==2){
			   sum+= BruTerm(refl,naturality*msign*ccfactor,ij,il,-im,ijpr,ilpr,impr,alpha,negm);
			   sum+= " + ";
			   sum+= BruTerm(refl,-1*naturalityprime*mprimesign*ccfactor,ij,il,im,ijpr,ilpr,-impr,alpha,negm);
			   int refl=-1;
			   sum+= " + ";
			   sum+= BruTerm(refl,naturality*msign*ccfactor,ij,il,-im,ijpr,ilpr,impr,alpha,negm);
			   sum+= " + ";
			   sum+= BruTerm(refl,-1*naturalityprime*mprimesign*ccfactor,ij,il,im,ijpr,ilpr,-impr,alpha,negm);
			 }
		      
		      //if not including -ve m remove XXX terms
		      sum.ReplaceAll("+ -1*XXX","");
		      sum.ReplaceAll("-1*XXX","");
		      sum.ReplaceAll("+ 1*XXX","");
		      sum.ReplaceAll("1*XXX","");
		      sum.ReplaceAll("+ XXX","");
		      sum.ReplaceAll("XXX +","");
		      sum.ReplaceAll("XXX","");
		      
		      while(sum.EndsWith("+")==kTRUE)
			sum.Chop();
		      //	if(sum.EndsWith("+")==kFALSE){
		      sum+="+"; //in case no -ve will get extra +
		      //}
		      //  else{
		      // }
		      sumTotal+=sum;
		      col++;
		    }
		  }
		  if(negm==1 ){
		    sdmeRow++;
		    if(sdmeRow==nwaves){sdmeRow=0;sdmeCol++;}
		  }
		  else {
		    if(im>=0&&impr>=0)sdmeRow++;
		    if(sdmeRow==nwaves){sdmeRow=0;sdmeCol++;}
		  }
		  
		  //	}
	      }
	    }
	  }
	}
      }
      //l and m
    }
    
  }
  cout<<"calculated "<<col <<" elements"<<endl;
  //Tidy up string
  sumTotal.ReplaceAll(" ","");
  while(sumTotal.EndsWith("+")==kTRUE)
    sumTotal.Chop();

  
  sumTotal.Capacity(sumTotal.Length());

  TString copySum;
  for(int ic=0;ic<sumTotal.Length();ic++)
    copySum+=sumTotal[ic];
    
  cout<<"chekc size "<<sumTotal.Sizeof()<<" "<<sumTotal.Length()<<endl;
  cout<<sumTotal<<endl;

  
    //auto simplified = sumTotal;
    auto simplified = Simplify(sumTotal);
    //    cout<<"PhotoLoader moment = "<<simplified<<endl;
    return simplified;

}


  

  ////////////////////////////////////////////////////////////////////////////////////////
    void LoadPartialWaves(Setup& setup,Int_t _Jmax, Int_t _Lmax,Int_t _Mmax,Int_t _mmax,Int_t nRefl=2,Bool_t onlyEven=true,Bool_t negm=1){
    /// Generate the H equation in terms of amplitudes
    /// a_l_m is used for +ve reflectivty amps and b_l_m -ve
    /// The final amplitude is constrained to = sqrt(1 - others_squared )
    /// This keeps all amplitudes in range 0-1

    gErrorIgnoreLevel = kFatal;
    auto level = RooMsgService::instance().globalKillBelow();
    RooMsgService::instance().setGlobalKillBelow(RooFit::INFO) ;

  //First calculate number of waves for each reflectivity
    Int_t nwaves=1;
    for(Int_t k=1;k<=_Jmax;k++)
      {
	Int_t nwaves_m=0;
	Int_t nwaves_l=3;
	Int_t m = k;
	if(k>_Mmax) m = _Mmax; //restrict m terms
	nwaves_m=m +m*negm+ 1;
	cout<<"nwaves "<<k<<" "<<nwaves_l<<" "<<nwaves_m<<endl;
	nwaves+=nwaves_l*nwaves_m;
      }


    ///////////////////////////////////////
    ///Create each wave (LoadParameter)
    //Now load each wave
   
    Int_t counter=0;
    cout<<" nwaves "<<nwaves<<endl;
    TString forNormalise="TMath::Sqrt(1.0";
    
    for(Int_t ij=0;ij<=_Jmax;ij++){
      for(Int_t il=std::max(0,ij-1);il<=std::min(_Lmax,ij+1);il++){
	for(Int_t im=-ij;im<=ij;im++){
	  if(TMath::Abs(im)>_Mmax) continue;
	  if(negm==0&&im<0) continue;
	  if(ij==0&&il==0) continue;
	  cout<<counter << "j "<<ij<<" l "<<il<<" m "<<im<<endl;
	  // cout<<counter<< " "<<nwaves<<endl;
	  counter++;
	  if(counter!=nwaves){//not final parameter for normalisation
	    setup.LoadParameter(Form("a_%d_%d_%d[0,0,1]",ij,il,im));
	    forNormalise+=Form("-@a_%d_%d_%d[]*@a_%d_%d_%d[]",ij,il,im,ij,il,im);
	  }
	  setup.LoadParameter(Form("aphi_%d_%d_%d[0,%lf,%lf]",ij,il,im,-TMath::Pi(),TMath::Pi()));
	  if(nRefl==2){
	    setup.LoadParameter(Form("b_%d_%d_%d[0,0,1]",ij,il,im));
	    forNormalise+=Form("-@b_%d_%d_%d[]*@b_%d_%d_%d[]",ij,il,im,ij,il,im);
	    setup.LoadParameter(Form("bphi_%d_%d_%d[0,%lf,%lf]",ij,il,im,-TMath::Pi(),TMath::Pi()));
	  }

	  if(counter==nwaves){//final parameter for normalisation
	    forNormalise+=")";
	    // setup.LoadParameter(Form("a_%d_%d[0,0,1]",il,im));
	    setup.LoadFormula(Form("a_%d_%d_%d=%s",ij,il,im,forNormalise.Data()));
	    //and for connvenience i.e. do not need to know last wave name
	    setup.LoadFormula(Form("normalise=%s",forNormalise.Data()));
	    
	    // std::cout<<" NORMALISATION "<<Form("a_%d_%d_%d=%s",ij,il,im,forNormalise.Data())<<std::endl;exit(0);
	  }
	}
      }
    }
  
    //////////////////////////////////////////////////////////////////
    /////Now create Moment formula as funciotn of partial waves
    //setup.LoadFormula(CGMatrixReflectivity("H",0,0,lmax,mmax,0,1));
   // Create moment as a function of partial waves
   for(int iJ=0; iJ<=2*_Jmax;iJ++)
     {
       for(int iM=0; iM<=iJ; iM++)
	 {
	   for(int iS=0; iS<=1; iS++)
	     {
	       //if iS==0: S=0, iLambda=0
	       //if iS==1: S=2, iLambda = 0,1,2
	       if(iS==0)
		 {
		   int  iLambda =0;
		   TString moment = CGMatrixReflectivity("H",iJ, iM, iS, iLambda, _Jmax, _Lmax, _Mmax, _mmax, 0, 1);
		   if( moment.EndsWith("=")==kFALSE) setup.LoadFormula(moment);
		   else {
		     moment.ReplaceAll("=","[0]");
		     setup.LoadConstant(moment);
		   }
		   TString moment1 = CGMatrixReflectivity("H",iJ, iM, iS, iLambda, _Jmax, _Lmax, _Mmax, _mmax, 1, 1);
		   if( moment1.EndsWith("=")==kFALSE) setup.LoadFormula(moment1);
		   else {
		     moment1.ReplaceAll("=","[0]");
		     setup.LoadConstant(moment1);
		   }
		   TString moment2 = CGMatrixReflectivity("H",iJ, iM, iS, iLambda, _Jmax, _Lmax, _Mmax, _mmax, 2, 1);
		   if( moment2.EndsWith("=")==kFALSE) setup.LoadFormula(moment2);
		   else {
		     moment2.ReplaceAll("=","[0]");
		     setup.LoadConstant(moment2);
		   }
		 }
	       if(iS==1)
		 {
		   for(int iLambda=0;iLambda<=2;iLambda++)
		     {
		       TString moment = CGMatrixReflectivity("H",iJ, iM, iS*2, iLambda, _Jmax, _Lmax, _Mmax, _mmax, 0, 1);
		       if (moment.EndsWith("=")==kFALSE) setup.LoadFormula(moment);
		       else {
			 moment.ReplaceAll("=","[0]");
			 setup.LoadConstant(moment);
		       }
		      
		       TString moment1 = CGMatrixReflectivity("H",iJ, iM, iS*2, iLambda, _Jmax, _Lmax, _Mmax, _mmax, 1, 1);
		       if (moment1.EndsWith("=")==kFALSE) setup.LoadFormula(moment1);
		       else {
			 moment1.ReplaceAll("=","[0]");
			 setup.LoadConstant(moment1);
		       }
		       TString moment2 = CGMatrixReflectivity("H",iJ, iM, iS*2, iLambda, _Jmax, _Lmax, _Mmax, _mmax, 2, 1);
		       if (moment2.EndsWith("=")==kFALSE) setup.LoadFormula(moment2);
		       else {
			 moment2.ReplaceAll("=","[0]");
			 setup.LoadConstant(moment2);
		       }
		     }
		 }
	     }
	 }
     }
   
    
   
   std::cout<<"Number of waves was "<<counter<<" compared to "<<nwaves <<endl;
   
   
   gErrorIgnoreLevel = kInfo ;
   RooMsgService::instance().setGlobalKillBelow(level) ;
  }
 
}


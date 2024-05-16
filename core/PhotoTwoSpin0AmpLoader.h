#pragma once

#include <TMath.h>
#include <TError.h>
#include <TF1.h>
#include <Math/SpecFuncMathMore.h>
#include "Setup.h"

namespace P2S0AmpLoader{
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
  inline TString BruTermCircle(int reflsign,double factor,int l,int m,int lpr,int mpr,int alpha,int negm){
   
    if(negm==0&&(m<0||mpr<0)) return "XXX";
    char refl='a';
    if(reflsign==-1)refl='b';

    int reflfactor = 1;
    if(alpha==1||alpha==2) reflfactor = reflsign; //-ve in Eqn D8b&c canceled by -ve in A9b 

    if(factor == 1) factor=1.0;
    if(factor == -1) factor=-1.0;
  
    if(l==lpr&&m==mpr)return Form("%d*%0.16f*(@%c_%d_%d[]*@%c_%d_%d[])",reflfactor,factor,refl,l,m,refl,lpr,mpr);
      
    TString func="cos"; //for real part
    if(alpha==3) func="sin"; //for imaginery part

    // cout<<"FACTOR "<<factor<<endl;
     return Form("%0.16f*(TMath::Abs(@%c_%d_%d[])*TMath::Abs(@%c_%d_%d[])*%s(@%cphi_%d_%d[]-@%cphi_%d_%d[]))",reflfactor*factor,refl,l,m,refl,lpr,mpr,func.Data(),refl,l,m,refl,lpr,mpr);
   }

TString Simplify(TString moment){

  auto expr = TString(moment(moment.First("=")+1,moment.Length()));
  auto terms = (expr.Tokenize('+'));

  auto SwitchWaves = [](const TString& pws){
    auto pw1 = pws(pws.First('_'),pws.First('[') - pws.First('_'));
    auto pwtemp = TString(pws(pws.First(']'),pws.Length()-pws.First(']')));
    auto pw2 = pwtemp(pwtemp.First('_'),pwtemp.First('[') - pwtemp.First('_'));
    auto ptemp = pws;
    ptemp.ReplaceAll(pw1,"PW1PART");
    ptemp.ReplaceAll(pw2,"PW2PART");
    ptemp.ReplaceAll("PW1PART",pw2);
    ptemp.ReplaceAll("PW2PART",pw1);
    return ptemp;
  };
  auto NumericalPart = [](const TString& pname){
    auto pnum = TString(pname(0,pname.First('(')-1));
    TF1 fnum("fnum",pnum);//use TF1 to convert arithmetic to single double	
    Double_t anum = fnum.Eval(0);
    return anum;
  };
  auto MatchTermWithConj = [&SwitchWaves,&NumericalPart](TString& aname,TString other){
		     auto isNeg = false;
		     auto isOtherNeg = false;
		     auto name = aname;
		     //store numerical factors
		     auto anum = NumericalPart(name);
		     auto onum = NumericalPart(other);
		     //cout<<"number " <<Form("%0.16f",anum)<<" other num "<<Form("%0.16f",onum)<<" dif "<<Form("%0.16f",anum-onum)<<" sum "<<Form("%0.16f",anum+onum)<<endl;
		     
		     //remove numerical factor
		     auto pwname = TString(name(name.First('('),name.Last(')')-name.First('(')+1));
		     //get switched pws tp get conjugate
		     auto conjname=SwitchWaves(pwname);
 
		     bool isChanging=false;
		     if(other.Contains(pwname)||other.Contains(conjname)){
		       auto conjfactor = 1;
		       if(pwname.Contains("sin")){
			 //check if same 
			 if(other.Contains(pwname)) conjfactor=1;
			 //or conjugate
			 else  conjfactor = -1; //sin changes sign when phi=-ve phi
		       }
		       //sum numerical parts
		       auto nsum = anum + onum*conjfactor;
		       if(TMath::Abs(nsum)<1E-15){
			 name=""; //exactly cancelled
		       }
		       else{//sum terms
			 name = Form("%0.16f * %s",nsum,pwname.Data());
		       }
		      
		       isChanging=true;
		       aname = name;
		     }

		     return isChanging;
		   };

  TString newExpr="";
  for(int i=0;i<terms->GetEntries();++i){
    auto name = TString(terms->At(i)->GetName());
    if(name=="") continue;
    //cout<<"Testing "<<name<<endl;
    for(int j=i+1;j<terms->GetEntries();++j){
      auto othername = TString(terms->At(j)->GetName());
      if(othername.Length()==0) continue;
      auto changed = MatchTermWithConj(name,othername);
      //cout<<" with "<<othername<<" "<<changed<<endl;
      if(changed){
	//cout<<"new name is "<<name<<" from "<<terms->At(i)->GetName()<<" and "<<othername<<endl;
	dynamic_cast<TObjString*>(terms->At(j))->SetString("");
	break;
      }
    }
    //    if(name == terms->At(i)->GetName())cout<<"Testing "<<name<<endl;
    if(name!=""){
      newExpr+=name;
      if(i!=terms->GetEntries()-1){
	newExpr+='+';
      }
    }
  }
 
    if(newExpr.EndsWith("+"))  newExpr = newExpr.Strip(TString::kTrailing,'+');
    newExpr.ReplaceAll(" ","");//remove whitespace so can iterate
    newExpr.Prepend(moment(0,moment.First("=")+1));
  
  
  return newExpr;
}
TString SimplifyAll(TString moment){
  TString target = moment;
  TString changed = "";
  auto count =0;
  //iterate until all terms combined together
  while (1){
    changed  = Simplify(target);
    //cout<<"target "<<target <<endl;
    //cout<<"changed "<<changed<<endl<<endl;
    if(changed==target) break;
    target = changed;
    count++;
  }
  //  cout<<"SimplifyAll interations "<<count<<endl;
  return target;
}
  /*///////////////////////////////////////////////////////////////////////////
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
  */
  //////////////////////////////////////////////////////////////////////
  inline TString CGMatrixReflectivity(TString Moment,Int_t L,Int_t M,
				      Int_t lmax, Int_t mmax,Int_t alpha,
				      Bool_t useNegRef,Bool_t onlyEven,Bool_t negm){
  
    //Sum_l,l',m'm' CLM CL0 rho_lm.rho_l'm'*
    //where CLM and CL0 are clebsch Gordan coeficients
  
    //sdm is hermitian so rho_(lm)(lm')=rho_(lm')(lm)*

    Int_t col =0;

    Int_t nwaves=1; //S0
  
    for(Int_t i=1;i<=lmax;i++){
      Int_t m=i;
      if(i>mmax) m = mmax; //restrict m terms
      nwaves+=m +m*negm+ 1;
    }
    Int_t nsdmes=nwaves*(nwaves+1)/2;
    Int_t nmoments = (2*lmax+1)*(lmax+1);
  
  
    vector<Double_t> Are(nsdmes);
    vector<Double_t> Aim(nsdmes);
  
    int sdmeRow=0;
    int sdmeCol=0;

    TString sumTotal(Form("%s_%d_%d_%d=",Moment.Data(),alpha,L,M));
 
    for(Int_t il=0;il<=lmax;il++){
      for(Int_t im=-il;im<=il;im++){
	if(TMath::Abs(im)>mmax) continue;
	for(Int_t ilpr=0;ilpr<=lmax;ilpr++){
	  for(Int_t impr=-ilpr;impr<=ilpr;impr++){
	    if(TMath::Abs(impr)>mmax) continue;
	    if(onlyEven == true){
	      if(il%2==1) continue;
	      if(ilpr%2==1) continue;
	    }
	    auto CM = ClebschGordan(ilpr,L,il,impr,M,im);
	    auto C0 = ClebschGordan(ilpr,L,il,0,0,0);
	    auto factor = TMath::Sqrt((2.*ilpr+1)/(2.*il+1));
	  
	  
	    auto conjCM = ClebschGordan(il,L,ilpr,im,M,impr);
	    auto conjC0 = ClebschGordan(il,L,ilpr,0,0,0);
	    auto conjfactor = TMath::Sqrt((2.*il+1)/(2.*ilpr+1));

	    int mmprimesign = 1;
	    if(TMath::Abs((im-impr)%2) == 1) mmprimesign=-1;
	    int msign = 1;
	    if(TMath::Abs((im)%2) == 1) msign=-1;
	    int mprimesign = 1;
	    if(TMath::Abs((impr)%2) == 1) mprimesign=-1;
	    auto ccfactor=CM*C0*factor; //2* as repeat -m

	 
	    if(ccfactor!=0){ //sum terms with non-zero cj coefficents
	 
	      TString sum =""; 
	      int refl=1;

	      if(alpha==0){
		//from eqn D8a
		sum+= BruTermCircle(refl,ccfactor,il,im,ilpr,impr,alpha,negm);
		
		sum+= " + ";
		sum+= BruTermCircle(refl,mmprimesign*ccfactor,il,-im,ilpr,-impr,alpha,negm);

		if(useNegRef==kTRUE){
		
		  refl=-1;
		  sum+= " + ";
		  sum+= BruTermCircle(refl,ccfactor,il,im,ilpr,impr,alpha,negm);
		  sum+= " + " ;
		  sum+= BruTermCircle(refl,mmprimesign*ccfactor,il,-im,ilpr,-impr,alpha,negm);
		}
	      }
	      else if(alpha==1){
		sum+= BruTermCircle(refl,msign*ccfactor,il,-im,ilpr,impr,alpha,negm);
		sum+= " + ";
		sum+= BruTermCircle(refl,mprimesign*ccfactor,il,im,ilpr,-impr,alpha,negm);
		
		if(useNegRef==kTRUE){
		  refl=-1;
		  sum+= " + ";
		  sum+= BruTermCircle(refl,msign*ccfactor,il,-im,ilpr,impr,alpha,negm);
		  sum+= " + ";
		  sum+= BruTermCircle(refl,mprimesign*ccfactor,il,im,ilpr,-impr,alpha,negm);
		}		
	      }
	      else if(alpha==2){
		sum+= BruTermCircle(refl,msign*ccfactor,il,-im,ilpr,impr,alpha,negm);
		sum+= " + ";
		//extra -1 in equation D8c
		sum+= BruTermCircle(refl,-1*mprimesign*ccfactor,il,im,ilpr,-impr,alpha,negm);

		
		if(useNegRef==kTRUE){
		  refl=-1;
		  sum+= " + ";
		  sum+= BruTermCircle(refl,msign*ccfactor,il,-im,ilpr,impr,alpha,negm);
		  sum+= " + ";
		  sum+= BruTermCircle(refl,-1*mprimesign*ccfactor,il,im,ilpr,-impr,alpha,negm);
		}
		
	      }
	      
	      if(alpha==3){
		//from eqn D8a
		//note mmprimesign cancels -ve in ccfactor for opposite SDME entries
		ccfactor*=-1;//-ve in Eqn A9d
		sum+= BruTermCircle(refl,ccfactor,il,im,ilpr,impr,alpha,negm);
		sum+= " +- ";//include + for easy splitting in simplify
		sum+= BruTermCircle(refl,mmprimesign*ccfactor,il,-im,ilpr,-impr,alpha,negm);

		if(useNegRef==kTRUE){
		
		  refl=-1;
		  sum+= " + ";
		  sum+= BruTermCircle(refl,ccfactor,il,im,ilpr,impr,alpha,negm);
		  sum+= " +- " ;//include + for easy splitting in simplify
		  sum+= BruTermCircle(refl,mmprimesign*ccfactor,il,-im,ilpr,-impr,alpha,negm);
		}
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
	      
	      sum+="+"; //in case no -ve will get extra +
	      //std::cout<<" sum "<<sum<<std::endl;
	      sumTotal+=sum;
	      col++;
	      }
	    }
	    /* if(negm==1 ){ */
	    /*   sdmeRow++; */
	    /*   if(sdmeRow==nwaves){sdmeRow=0;sdmeCol++;} */
	    /* } */
	    /* else { */
	    /*   if(im>=0&&impr>=0)sdmeRow++; */
	    /*   if(sdmeRow==nwaves){sdmeRow=0;sdmeCol++;} */
	    /* } */

	  
	}
	
      }
      
    }
     
    //Tidy up string
    sumTotal.ReplaceAll(" ","");
    while(sumTotal.EndsWith("+")==kTRUE)
      sumTotal.Chop();

  
    sumTotal.Capacity(sumTotal.Length());

    TString copySum;
    for(int ic=0;ic<sumTotal.Length();ic++)
      copySum+=sumTotal[ic];
    
    //auto simplified = sumTotal;
    //auto simplified = Simplify(sumTotal);
    auto simplified = SimplifyAll(sumTotal);
    //cout<<"PhotoLoader moment = "<<simplified<<endl;
    return simplified;

  }

  ////////////////////////////////////////////////////////////////////////////////////////
  void LoadPartialWaves(Setup& setup,Int_t lmax, Int_t mmax,Int_t nRefl=2,Bool_t onlyEven=true,Bool_t negm=1){
    /// Generate the H equation in terms of amplitudes
    /// a_l_m is used for +ve reflectivty amps and b_l_m -ve
    /// The final amplitude is constrained to = sqrt(1 - others_squared )
    /// This keeps all amplitudes in range 0-1

    gErrorIgnoreLevel = kFatal;
    auto level = RooMsgService::instance().globalKillBelow();
    RooMsgService::instance().setGlobalKillBelow(RooFit::INFO) ;

    //First calculate number of waves for each reflectivity
    Int_t nwaves=1;
    for(Int_t i=1;i<=lmax;i++){
      Int_t m=i;
      if(i>mmax) m = mmax; //restrict m terms
      nwaves+=m +m*negm+ 1;
    }


    ///////////////////////////////////////
    ///Create each wave (LoadParameter)
    //Now load each wave
    Int_t counter=0;
    TString forNormalise="TMath::Sqrt(1.0";
   
    // for(Int_t il=lmax;il<=lmax;il++){
    // for(Int_t il=0;il<=lmax;il++){
    for(Int_t il=0;il<=lmax;il++){
      for(Int_t im=-il;im<=il;im++){
	if(TMath::Abs(im)>mmax) continue;
	if(negm==0&&im<0) continue;
	counter++;
	if(onlyEven == true) if(il%2==1) continue;
	if(counter!=nwaves){//not final parameter for normalisation
	  setup.LoadParameter(Form("a_%d_%d[0,-1,1]",il,im));
	  forNormalise+=Form("-@a_%d_%d[]*@a_%d_%d[]",il,im,il,im);
	}

	setup.LoadParameter(Form("aphi_%d_%d[0,%lf,%lf]",il,im,-10*TMath::Pi(),10*TMath::Pi()));

	if(nRefl==2){
	  setup.LoadParameter(Form("b_%d_%d[0,-1,1]",il,im));
	  forNormalise+=Form("-@b_%d_%d[]*@b_%d_%d[]",il,im,il,im);
	  setup.LoadParameter(Form("bphi_%d_%d[0,%lf,%lf]",il,im,-10*TMath::Pi(),10*TMath::Pi()));
	}
      
	if(counter==nwaves){//final parameter for normalisation
	  forNormalise+=")";
	  //	setup.LoadParameter(Form("a_%d_%d[0,0,1]",il,im));
	  setup.LoadFormula(Form("a_%d_%d=%s",il,im,forNormalise.Data()));
	  //and for connvenience i.e. do not need to know last wave name
	  setup.LoadFormula(Form("normalise=%s",forNormalise.Data()));
	
	}
      }
    }
    //////////////////////////////////////////////////////////////////
    /////Now create Moment formula as funciotn of partial waves
    //setup.LoadFormula(CGMatrixReflectivity("H",0,0,lmax,mmax,0,1));
    ///for(int iL=2*lmax;iL<=2*lmax;iL++)
    for(int iL=0;iL<=2*lmax;iL++)
      for(int iM=0;iM<=iL;iM++)
	{
	  auto momf = CGMatrixReflectivity("H",iL,iM,lmax,mmax,0,(nRefl==2),onlyEven,negm);
	  if(momf[momf.Length()-1]!='=')setup.LoadFormula(momf);
	  else {
	    momf.ReplaceAll("=","[0]");
	    setup.LoadConstant(momf);
	  }
	  momf=CGMatrixReflectivity("H",iL,iM,lmax,mmax,1,(nRefl==2),onlyEven,negm);
	  if(momf[momf.Length()-1]!='=')setup.LoadFormula(momf);
	  else {
	    momf.ReplaceAll("=","[0]");
	    setup.LoadConstant(momf);
	  }
	 
	  if(iM!=0){
	    momf=CGMatrixReflectivity("H",iL,iM,lmax,mmax,2,(nRefl==2),onlyEven,negm);
	    if(momf[momf.Length()-1]!='=')setup.LoadFormula(momf);
	 
	    else {
	      momf.ReplaceAll("=","[0]");
	      setup.LoadConstant(momf);
	    }
	  }
	  if(iM!=0){
	    momf=CGMatrixReflectivity("H",iL,iM,lmax,mmax,3,(nRefl==2),onlyEven,negm);
	    if(momf[momf.Length()-1]!='=')setup.LoadFormula(momf);
	    
	    else {
	      momf.ReplaceAll("=","[0]");
	      setup.LoadConstant(momf);
	    }
	  }
	}
    gErrorIgnoreLevel = kInfo ;
    RooMsgService::instance().setGlobalKillBelow(level) ;
   }

}

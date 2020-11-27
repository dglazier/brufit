#include "PlotResults.h"
#include "RooHSEventsPDF.h"
#include "RooMcmc.h"
#include "FitManager.h"
#include <RooPlot.h>
#include <RooMsgService.h>
#include <TCanvas.h> 
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TDirectory.h>

namespace HS{
  namespace FIT{
    using namespace RooFit;

    AutocorrPlot::AutocorrPlot(Setup *setup, RooMcmc* mcmc)
    {
      fCanvases->SetName(TString("Autocorrelation plot"));

      auto tree = mcmc->GetTree();
  
      auto& pars = setup->ParsAndYields();
      std::cout<<"The parameters are: " << pars<<std::endl;
      Int_t Npars = pars.size();
      std::cout<<"The number of parameters is:  "<<Npars<<std::endl;      const   Int_t Nentries =200;// MCMCTree->GetEntries();

      auto canName = "Autocorrelation Plot";
      auto canvas = new TCanvas(canName, canName);
      fCanvases->Add(canvas);
      
      auto mg = new TMultiGraph();
      auto leg = new TLegend(0.7, 0.5, 0.9, 0.9);

      Int_t param_iter = 0;
      vector<Double_t> param(Npars);
      vector<Double_t> chain_mean(Npars);
   

      for (RooAbsArg* ipar : pars)
	{//Loop over parameters
	  
	  Double_t entryTree[Nentries];
	  Double_t lag[Nentries]; 
	  Double_t autocorr[Nentries];	
	
	  tree->SetBranchAddress(ipar->GetName(), &param[param_iter]);
	  //	ResultTree->SetBranchAddress(ipar->GetName(), &chain_mean[param_iter]);        
	  //std::cout<<ipar->GetName()<<std::endl;
	  
	  Double_t chain_total=0;
	  for (int ientry = 0; ientry<Nentries; ientry++)
	    {//Loop over entries of the tree
	    	     	     	    
	      tree->GetEntry(ientry);
	      entryTree[ientry] = param[param_iter];
	      lag[ientry] = ientry;
	      //Calculate the average of the chain
	      chain_total+=param[param_iter];
	      chain_mean[param_iter] = chain_total/Nentries;
	      
	    }//Close loop over entries
	  	
	  Double_t parD = 0;	
	  for (int ilagEntryD= 0; ilagEntryD<Nentries; ilagEntryD++)
	    {
	      parD += (entryTree[ilagEntryD]-chain_mean[param_iter])*(entryTree[ilagEntryD]-chain_mean[param_iter]);	    
	    }
	
	  for (int ilag=0; ilag<Nentries; ilag++)
	    {//compute the autocorr value as a function of lag
	      Double_t parN=0;
	      Double_t par= 0;
	      for(int ilagEntryN =0; ilagEntryN<(Nentries-ilag); ilagEntryN++)
		{
		  parN += (entryTree[ilagEntryN]-chain_mean[param_iter])*(entryTree[ilagEntryN+ilag]-chain_mean[param_iter]);           
		}
	    
	      par = parN/parD;
	      autocorr[ilag] = par;
	      if(ilag==100){std::cout<<chain_mean[param_iter]<<"  "<<entryTree[ilag]<<"   "<<autocorr[ilag]<<std::endl;}

	    }        

	  TGraph *gr1 = new TGraph(Nentries, lag, autocorr);
	  gr1->SetTitle(ipar->GetName());
	  gr1->SetLineWidth(2);
	  gr1->SetLineColor(param_iter%9+1);
	  gr1->SetMarkerColor(param_iter%9+1);
	  mg->Add(gr1);
	  leg->AddEntry(gr1,ipar->GetName());	 
	  
	  param_iter++;	
	}//Close loop over params
      
      mg->Draw("AL"); 
      mg->SetTitle("Autocorrelation plot");
      mg->GetXaxis()->SetTitle("lag");
      mg->GetYaxis()->SetTitle("autocorrelation");
      leg->Draw("same");

      canvas->Modified();
      canvas->Update();
      canvas->Draw("");		 
    }
  }//FIT
}//HS

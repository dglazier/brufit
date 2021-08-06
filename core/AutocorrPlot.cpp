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

    AutocorrPlot::AutocorrPlot(Setup *setup, RooMcmc* mcmc, TList* canvases)
    {
    
      auto tree =mcmc->GetTree()->CopyTree("",""); //make a copy
      RemoveNegativeInNames(tree);
      Int_t burnIn  = mcmc->GetNumBurnInSteps();
  
      auto& pars = setup->ParsAndYields();
      std::cout<<"The parameters are: " << pars<<std::endl;
      Int_t Npars = pars.size();
      std::cout<<"The number of parameters is:  "<<Npars<<std::endl;    
      const   Int_t Nentries = tree->GetEntries()-burnIn;

      auto canName = "Autocorrelation Plot";
      auto canvas = new TCanvas(canName, canName);
      canvases->Add(canvas);
      
      auto mg = new TMultiGraph();
      auto leg = new TLegend(0.7, 0.5, 0.9, 0.9);

      Int_t param_iter = 0;
      vector<Double_t> param(Npars);
      vector<Double_t> chain_mean(Npars);

      for (RooAbsArg* ipar : pars)
	{//Loop over parameters
	  if(ipar->isConstant()==kTRUE) continue;

	  Double_t entryTree[Nentries];
	  Double_t lag[Nentries]; 
	  Double_t autocorr[Nentries];	
	
	  tree->SetBranchAddress(CheckForNegatives(ipar->GetName()), &param[param_iter]);
	  Double_t time = 0;
	  Double_t timePar = 0;
	  Double_t chain_total=0;
 
	  for (int ientry = 0; ientry<Nentries; ientry++)
	    {//Loop over entries of the tree
	    	     	     	    
	      tree->GetEntry(ientry+burnIn);
	      entryTree[ientry] = param[param_iter];
	      lag[ientry] = ientry;
	      //Calculate the average of the chain
	      chain_total+=param[param_iter];
	      chain_mean[param_iter] = chain_total/Nentries;
	      
	    }//Close loop over entries
	  	
	  Double_t parD = 0;	
	  for (int ilagEntryD= 0; ilagEntryD<Nentries; ilagEntryD++)
	    { //This is really a loop over entries
	      parD += (entryTree[ilagEntryD]-chain_mean[param_iter])*(entryTree[ilagEntryD]-chain_mean[param_iter]);	    
	    }
	
	  for (int ilag=0; ilag<Nentries; ilag++)
	    {//compute the autocorr value as a function of lag. This is really a loop over lag
	      Double_t parN=0;
	      Double_t par= 0;
	      for(int ilagEntryN =0; ilagEntryN<(Nentries-ilag); ilagEntryN++)
		{//This is a loop over entries of the tree at fixed tau
		  parN += (entryTree[ilagEntryN]-chain_mean[param_iter])*(entryTree[ilagEntryN+ilag]-chain_mean[param_iter]);           
		}
	    
	      par = parN/parD;
	      autocorr[ilag] = par;	      
	      
	    }//End of loop ilag

	  Int_t timepar = std::min(5000, (Nentries-burnIn));


	  for(Int_t itime=0; itime<timepar; itime++)
	    {
	      timePar += autocorr[itime];
	    }
	  time = 1+2*timePar;        
	  std::cout<<"The autocorr time for " <<ipar->GetName()<< " is: "<<time<<std::endl;
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


      delete tree;
      
      canvas->Modified();
      canvas->Update();
      canvas->Draw("");		 
    }
  }//FIT
}//HS

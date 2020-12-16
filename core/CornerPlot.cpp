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

    CornerPlot::CornerPlot(Setup *setup, RooMcmc *mcmc)
    {
      fCanvases->SetName(TString("Corner Plot"));
      auto vars=setup->FitVars();

      auto tree = mcmc->GetTree();

      auto& pars = setup->ParsAndYields();
      std::cout<<"The Parameters are: "<<std::endl;
      Int_t Npars = pars.size();
      std::cout<<"The number of parameters is: "<<Npars<<std::endl;

      vector<Double_t> param(Npars);
      vector<Double_t> chain_mean(Npars);
      Int_t param_iter=0;

      for (RooAbsArg* ipar : pars)
	{//Loop over parameters once to set values from tree
	  tree->SetBranchAddress(ipar->GetName(), &param[param_iter++]);
	}//close loop over setting params

      for(auto var : vars)
	{
      auto canName = "Corner Plot";
      auto canvas = new TCanvas(canName, canName);
      
      fCanvases->Add(canvas);
      canvas->Divide(Npars, Npars);

      /*
      tree->UseCurrentStyle();
      gStyle->SetOptStat(0);
      gStyle->SetHistFillColor(1); 
      gStyle->SetTitleSize(0.15, "t");
      gStyle->SetTitleY(1.02);*/
      Int_t counter=0;//Counter for new line on canvas      

      for (RooAbsArg* ipar : pars)
	{//Loop over parameters twice to draw corner plot
	 //fix the first parameter 
         //loop over a second parameter to draw the corresponding hists
	  Int_t int_counter =0; //Counter for across canvas

	  for (RooAbsArg* ipar2 : pars)
	    {//second loop 
	      
	      if(ipar==ipar2)
		{
		  canvas->cd(Npars*counter+1+int_counter);
		  tree->Draw(ipar->GetName());
		  counter++;
		  
       		  break;
		}

	      else if(ipar!=ipar2)
		{
		  canvas->cd(Npars*counter+1+int_counter);
		  
		  TString DrawPar = ipar->GetName();
		  TString DrawPar2 = ipar2->GetName();
		  TString Draw2D = DrawPar + ":" + DrawPar2;
		  tree->Draw(Draw2D,"", "col");
		  int_counter++;

		}
	      
	    }//End Second loop ipar2
	}//End First loop ipar
      canvas->Modified();
      canvas->Update();
      canvas->Draw();
      
      
	}//Loop over vars

    }//CornerPlot
  }//FIT
}//HS


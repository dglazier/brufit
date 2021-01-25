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
#include <TStyle.h>
#include <TLine.h>

namespace HS{
  namespace FIT{
    using namespace RooFit;

    CornerFullPlot::CornerFullPlot(Setup *setup, RooMcmc *mcmc, TList* canvases)
    {
      auto myStyle = TStyle{*gStyle};
      myStyle.SetName("MCMCStyle");
      
      myStyle.SetOptStat(0);
      myStyle.SetTitleY(1.0);
      myStyle.SetTitleTextColor(2);
      myStyle.SetTitleSize(0.2, "t");
      myStyle.SetLabelSize(0.125, "xy");
      myStyle.SetNdivisions(4, "xy");

      auto defStyle=gStyle;
      gStyle=&myStyle;
  
      
 
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

      auto canName = "Corner Full Plot";
      auto canvas = new TCanvas(canName, canName);
      
      canvases->Add(canvas);
      canvas->Divide(Npars, Npars);

      tree->UseCurrentStyle();

      for(auto var : vars)
	{
	  
	  Int_t counter=0;//Counter for new line on canvas      
	  
	  for (RooAbsArg* ipar : pars)
	    {//Loop over parameters twice to draw corner plot
	      //fix the first parameter 
	      //loop over a second parameter to draw the corresponding hists
	      Int_t int_counter =1; //Counter for across canvas

	      for (RooAbsArg* ipar2 : pars)
		{//second loop 
	      
		  if(ipar==ipar2)
		    {
		      auto can = canvas->cd(Npars*counter+int_counter); 
		      can->SetBorderSize(0);
		      can->SetTopMargin(0.0);
		      can->SetBottomMargin(0.1);
		      can->SetLeftMargin(0.19);
		      can->SetRightMargin(0.01);

		      TString DrawParInd1 = ipar->GetName();
		      TString histname= TString("hmcmc_")+ipar->GetName();
		      DrawParInd1+=">>";
		      DrawParInd1+=histname;

		      tree->Draw(DrawParInd1);
		      auto htemp=dynamic_cast<TH1*>(gDirectory->FindObject(histname));
		      Double_t mean = htemp->GetMean();
		      gDirectory->ls();

		      TLine *line = new TLine(mean,0,mean, htemp->GetMaximum());
		      line->SetLineColor(kRed);
		      line->Draw(); 

		      counter++;
		  
		      break;
		    }

		  else if(ipar!=ipar2)
		    {
		      auto can = canvas->cd(Npars*counter+int_counter);
		      can->SetBorderSize(0);
		      can->SetTopMargin(0);
		      can->SetBottomMargin(0.1);
		      can->SetLeftMargin(0.19);
		      can->SetRightMargin(0.01);
		      
		      TString DrawPar = ipar->GetName();
		      TString DrawPar2 = ipar2->GetName();

		      Double_t maxY =  tree->GetMaximum(DrawPar);		 
		      Double_t minY =  tree->GetMinimum(DrawPar);
		      Double_t maxX =  tree->GetMaximum(DrawPar2);		 
		      Double_t minX =  tree->GetMinimum(DrawPar2);

		      TString histname="hist_";
		      histname+=ipar->GetName();
		      histname+=ipar2->GetName();

		      auto hist =  TH2F{histname,histname, 50,minX, maxX,50, minY,maxY };
		      TString Draw2D = DrawPar + ":" + DrawPar2;
		      TString Draw2D1 = DrawPar + ":" + DrawPar2 + ">>"+histname;
		      
		      tree->Draw(Draw2D1);
		      tree->Draw(Draw2D, "", "col");

		      Double_t meanX = hist.GetMean(1);
		      Double_t meanY = hist.GetMean(2);
	          
		      TLine* lineV = new TLine(meanX,minY, meanX, maxY);
		      TLine* lineH = new TLine(minX, meanY, maxX, meanY);
		      lineV->SetLineColor(kRed);
		      lineH->SetLineColor(kRed);
		      lineV->Draw();
		      lineH->Draw();
		  
		     
		      /* In case of needing the mode..
		      //Start Mode Line Two
		      TH1F* histY = new TH1F("histY", "histY", 100,-100000,100000);
		      TString DrawParInd1Y = ipar->GetName();
		      TString DrawParIndY = DrawParInd1Y + ">>histY";
		      tree->Draw(DrawParIndY);
		      Double_t meanY = histY->GetMean();
		      Double_t sigmaY = histY->GetRMS();
		      TH1F *hist2Y = new TH1F("hist2Y", "hist2Y", 100, meanY-3*sigmaY, meanY+3*sigmaY);
		      TString DrawParFinY = DrawParInd1Y + ">>hist2Y";
		      tree->Draw(DrawParFinY);
		      Int_t binmaxY = hist2Y->GetMaximumBin();
		      Double_t modeY = hist2Y->GetBinCenter(binmaxY);

		      TH1F* histX = new TH1F("histX", "histX", 100,-100000,100000);
		      TString DrawParInd1X = ipar2->GetName();
		      TString DrawParIndX = DrawParInd1X + ">>histX";
		      tree->Draw(DrawParIndX);
		      Double_t meanX = histX->GetMean();
		      Double_t sigmaX = histX->GetRMS();
		      TH1F *hist2X = new TH1F("hist2X", "hist2X", 100, meanX-3*sigmaX, meanX+3*sigmaX);
		      TString DrawParFinX = DrawParInd1X + ">>hist2X";
		      tree->Draw(DrawParFinX);
		      Int_t binmaxX = hist2X->GetMaximumBin();
		      Double_t modeX = hist2X->GetBinCenter(binmaxX);

		  
		      TString DrawPar = ipar->GetName();
		      TString DrawPar2 = ipar2->GetName();
		      TString title = DrawPar + ":" + DrawPar2;
		      TString Draw2D = DrawPar + ":" + DrawPar2;
		      tree->Draw(Draw2D, "", "col");
        
		      TLine *lineV = new TLine(modeX, meanY-3*sigmaY, modeX, meanY+3*sigmaY);
		      TLine *lineH = new TLine(meanX-3*sigmaX, modeY, meanX+3*sigmaX, modeY);
		      lineV->SetLineColor(kRed);
		      lineH->SetLineColor(kRed);
		      tree->Draw(Draw2D,"", "col");
		      lineV->Draw();
		      lineH->Draw(); 
		  
		      //End Mode Line Two
		      */
		      int_counter++;

		    }
	      
		}//End Second loop ipar2
	    }//End First loop ipar
	  canvas->Modified();
	  canvas->Update();
	  canvas->Draw();

     
      
	}//Loop over vars

      //change style back
      gStyle=defStyle;

    }//CornerFullPlot
  }//FIT
}//HS


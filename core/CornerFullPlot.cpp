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

      auto tree =mcmc->GetTree()->CopyTree("",""); //make a copy
      RemoveNegativeInNames(tree);
      
      auto& pars = setup->ParsAndYields();
      //  std::cout<<"The Parameters are: "<<std::endl;
      Int_t Npars = pars.size();
      //std::cout<<"The number of parameters is: "<<Npars<<std::endl;
      auto canName = "Corner Full Plot";
      auto canvas = new TCanvas(canName, canName);
      
      canvases->Add(canvas);
      canvas->Divide(Npars, Npars);

  	  
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
		  can->SetBottomMargin(0.);
		  can->SetLeftMargin(0.);
		  can->SetRightMargin(0.0);

		  if(ipar->isConstant()==kFALSE){
		    TString DrawParInd1 = CheckForNegatives(ipar->GetName());
		    TString histname= TString("cornerf_")+ipar->GetName();
		    DrawParInd1+=">>";
		    DrawParInd1+=histname;

		    tree->Draw(DrawParInd1);
		    auto htemp=dynamic_cast<TH1*>(gDirectory->FindObject(histname));
		    if(htemp!=nullptr){
		      Double_t mean = htemp->GetMean();

		      TLine *line = new TLine(mean,0,mean, htemp->GetMaximum());
		      line->SetLineColor(kRed);
		      line->Draw(); 
		    }
		  }
		  counter++;
		  
		  break;
		}

	      else if(ipar!=ipar2)
		{
		  auto can = canvas->cd(Npars*counter+int_counter);
		  can->SetBorderSize(0);
		  can->SetTopMargin(0);
		  can->SetBottomMargin(0.);
		  can->SetLeftMargin(0.);
		  can->SetRightMargin(0.0);
		      
		  if(ipar->isConstant()==kFALSE&&ipar2->isConstant()==kFALSE){
		    TString DrawPar = CheckForNegatives(ipar->GetName());
		    TString DrawPar2 = CheckForNegatives(ipar2->GetName());

		    Double_t maxY =  tree->GetMaximum(DrawPar);		 
		    Double_t minY =  tree->GetMinimum(DrawPar);
		    Double_t maxX =  tree->GetMaximum(DrawPar2);		 
		    Double_t minX =  tree->GetMinimum(DrawPar2);

		    TString histname="cornerf2_";
		    histname+=ipar->GetName();
		    histname+=ipar2->GetName();

		    auto hist =  new TH2F{histname,histname, 50,minX, maxX,50, minY,maxY };
		    // TString Draw2D = DrawPar + ":" + DrawPar2;
		    TString Draw2D1 = DrawPar + ":" + DrawPar2 + ">>"+histname;
		      
		    tree->Draw(Draw2D1,"","col");
		    //tree->Draw(Draw2D, "", "col");

		    Double_t meanX = hist->GetMean(1);
		    Double_t meanY = hist->GetMean(2);
	          
		    TLine* lineV = new TLine(meanX,minY, meanX, maxY);
		    TLine* lineH = new TLine(minX, meanY, maxX, meanY);
		    lineV->SetLineColor(kRed);
		    lineH->SetLineColor(kRed);
		    lineV->Draw();
		    lineH->Draw();
		  
		  }
		  int_counter++;
		  can->Modified();
		  can->Update();


		}
	      
	    }//End Second loop ipar2
	}//End First loop ipar
      canvas->Modified();
      canvas->Update();
      canvas->Draw();

      delete tree;
 
      //change style back
      gStyle=defStyle;

    }//CornerFullPlot



  }//FIT
}//HS


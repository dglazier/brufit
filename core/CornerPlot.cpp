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

    CornerPlot::CornerPlot(Setup *setup, RooMcmc *mcmc, TList* canvases)
    {
      auto myStyle = TStyle{*gStyle};
      myStyle.SetName("MCMCStyle");
      
      myStyle.SetFillColor(0);
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

      auto burnIn=mcmc->GetNumBurnInSteps();

      auto& pars = setup->ParsAndYields();
      // std::cout<<"The Parameters are: "<<std::endl;
      Int_t Npars = pars.size();
      //std::cout<<"The number of parameters is: "<<Npars<<std::endl;
      
      //Don't need canvas for each var
      auto canName = "Corner Plot";
      //create canvas, going to be given to canvases for ownership
      auto canvas = new TCanvas(canName, canName); 
      canvas->Divide(Npars, Npars);
    
      canvases->Add(canvas);
       
   
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
		      
		  if(ipar->isConstant()==kFALSE && setup->IsParSetConst(ipar->GetName())==kFALSE){
		    TString DrawParInd1 = CheckForNegatives(ipar->GetName());
		    TString DrawParInd = DrawParInd1 + ">>";
		    TString histname= "corner_";
		    histname+=ipar->GetName();
		    DrawParInd+=histname;
		      
		    tree->Draw(DrawParInd,"","",tree->GetEntries()-burnIn,burnIn);
		    auto hist=dynamic_cast<TH1*>(gDirectory->FindObject(histname));
		    if(hist!=nullptr){

		      Double_t mean = hist->GetMean();
		      Double_t sigma = hist->GetRMS();
		      auto hist2 = TH1F{histname+"range", ipar->GetName(), 100, mean-3*sigma, mean+3*sigma};
		   
		      TString DrawParFin = DrawParInd1 + ">>"+histname+"range";
		      tree->Draw(DrawParFin,"","",tree->GetEntries()-burnIn,burnIn); 
		      hist2.DrawCopy();

		      //need to let ROOT delete line
		      auto* line = new TLine{mean,0,mean, hist->GetMaximum()};
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
		  can->SetBottomMargin(0.1);
		  can->SetLeftMargin(0.19);
		  can->SetRightMargin(0.01);
		  
		  if(ipar->isConstant()==kFALSE && ipar2->isConstant()==kFALSE
		     && setup->IsParSetConst(ipar->GetName())==kFALSE && setup->IsParSetConst(ipar2->GetName())==kFALSE){
		    
		    TString DrawPar = CheckForNegatives(ipar->GetName());
		    TString DrawPar2 = CheckForNegatives(ipar2->GetName());
		    TString title = DrawPar + ":" + DrawPar2;
		    TString Draw2D = DrawPar + ":" + DrawPar2 + ">>";
		    TString histname="corner2_";
		    histname+=ipar->GetName();
		    histname+="_";
		    histname+=ipar2->GetName();
		    Draw2D +=histname;
	
		    tree->Draw(Draw2D, "", "goff",tree->GetEntries()-burnIn,burnIn);
		    auto hist=dynamic_cast<TH2*>(gDirectory->FindObject(histname));
		    Double_t meanX = hist->GetMean(1);
		    Double_t meanY = hist->GetMean(2);
		    Double_t rmsX = hist->GetRMS(1);
		    Double_t rmsY = hist->GetRMS(2);
		    delete hist;
		
		    auto hist2 = new TH2F{histname,histname, 50, meanX-3*rmsX, meanX+3*rmsX, 50, meanY-3*rmsY, meanY+3*rmsY};
		    hist2->SetTitle(title);

		    TString Draw2DFin = DrawPar + ":" + DrawPar2 + ">>"+histname;
		    tree->Draw(Draw2DFin,"", "col",tree->GetEntries()-burnIn,burnIn);
		  
		    //need to let ROOT delete line
		    auto lineV = new TLine{meanX,meanY-3*rmsY, meanX, meanY+3*rmsY};
		    auto lineH = new TLine{meanX-3*rmsX, meanY, meanX+3*rmsX, meanY};
		    lineV->SetLineColor(kRed);
		    lineH->SetLineColor(kRed);
		    lineV->Draw();
		    lineH->Draw();
		  
		  }
		  can->Modified();
		  can->Update();

		  int_counter++;


		}
	      
	    }//End Second loop ipar2
	}//End First loop ipar
      canvas->Modified();
      canvas->Update();
      canvas->Draw();
    
      delete tree;
 
      
      //change style back
      gStyle=defStyle;


    }//CornerPlot
  }//FIT
}//HS


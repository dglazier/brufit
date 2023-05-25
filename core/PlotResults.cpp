#include "PlotResults.h"
#include "RooHSEventsPDF.h"
#include "RooMcmc.h"
#include <RooMsgService.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TDirectory.h>

namespace HS{
  namespace FIT{
    
    PlotResults::PlotResults(const Setup *setup,const RooDataSet* data,const TString& tag,const TString& opt):fPlotOptions{opt}{

      using namespace RooFit;
      //stop plots being deleted when current directory is deleted
      RooPlot::AddDirectory(false);

      Bool_t IsBatch=gROOT->IsBatch();
      if(fPlotOptions.Contains("goff")==kTRUE) gROOT->SetBatch(kTRUE);
      
      //cout<<"PlotResults::PlotResults "<<fCanvases.get()<<" "<<setup<<" "<<endl;
      //fCanvases->SetOwner();
      fCanvases->SetName(TString("RFPlots")+setup->GetName());

  	
      auto vars=setup->FitVars();
      auto *model=setup->Model();
      // cout<<"model "<<model<<endl;
      //For RooHSEventsPDF flag plotting so can calc partial integrals
      //for 1 observable case;
      RooHSEventsPDF::SetIsPlotting(kTRUE);

     
      for( auto *var : vars){
	
	auto canName=tag+"_"+ var->GetName();
	auto canvas=new TCanvas(canName,canName);
	fCanvases->Add(canvas);
	
	canvas->Divide(2,1);
	canvas->cd(1);

	RooPlot* frame = var->frame();
       	data->plotOn(frame, DataError(RooAbsData::SumW2) ) ; 

	const auto& pdfs = setup->constPDFs();

	if(pdfs.getSize()>1) //only draw components if >1
	  for(Int_t ic=0;ic<pdfs.getSize();ic++)
	    model->plotOn(frame,Components(pdfs[ic]),LineStyle(kDashed),LineColor(ic%8+1),Precision(1E-2));

	
      	model->plotOn(frame,LineColor(kRed)) ;

	//	model->getParameters(data)->Print("v");
	//	std::cout<<"PlotOn  "<<model->expectedEvents(model->getObservables(data))<<std::endl;
       	//model->plotOn(frame,LineColor(kRed),Precision(4E-2)) ;
	
	model->paramOn(frame,
		       Layout(0.1, 0.2, 0.9),
		       Format("NE",AutoPrecision(1))); //show fit parameters
	
	
	frame->SetTitle(TString("Fit components for ")+var->GetName());

	frame->Draw();
	fRooPlots.push_back(rooplot_uptr{frame});//keep it live
	
	canvas->SetTheta(frame->chiSquare()); //Save Chi2 as theta variable in TCanvas. Retrieve with TCanvas::GetTheta().
	auto level = RooMsgService::instance().globalKillBelow();
	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
 
	///////////////////////////////////////////
	//Residual distributions
	auto halfCanvas=canvas->cd(2);
	halfCanvas->Divide(1,2);
	halfCanvas->cd(1);
	
	auto hres=roohist_uptr(frame->residHist());
	hres->Draw();
     	fRooHists.push_back(std::move(hres));//keep it live

	//Pull distributions
	halfCanvas->cd(2);
	
	auto hpull=roohist_uptr(frame->pullHist());
	hpull->Draw();
	fRooHists.push_back(std::move(hpull));//keep it live
	//////////////////////////////////////////////

	
	RooMsgService::instance().setGlobalKillBelow(level);

	canvas->Modified();
	canvas->Update();
	if(fPlotOptions.Contains("goff")==kFALSE)canvas->Draw("");


     }

      //Turn off plotting in RooHSEventsPDF
      RooHSEventsPDF::SetIsPlotting(kFALSE);
      if(fPlotOptions.Contains("goff")==kTRUE) gROOT->SetBatch(IsBatch);

     }

    void PlotResults::Write(){
      fCanvases->Write();
    }

    TString PlotResults::CheckForNegatives(TString name){
      if(name.Contains('-')){
	name.ReplaceAll("-","neg");
      }
      return name;
    }
    
    void PlotResults::RemoveNegativeInNames(TTree* tree){
      
      auto branches=tree->GetListOfBranches();
      for(auto br:*branches){
	TString name=br->GetName();
	if(name.Contains('-')){
	  name.ReplaceAll("-","neg");
	  dynamic_cast<TBranch*>(br)->SetName(name);
	}
      }
      
    }
  }//namespace FIT

}//namespace HS

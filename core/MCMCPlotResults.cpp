#include "PlotResults.h"
#include "RooHSEventsPDF.h"
#include "RooMcmc.h"
#include "FitManager.h"
#include <RooPlot.h>
#include <RooMsgService.h>
#include <TCanvas.h> 
#include <TDirectory.h>

namespace HS{
  namespace FIT{
    using namespace RooFit;

    ////////Can use this to get the name of the mcmc file
    // TString fileName=setup->GetOutDir()+setup->GetName()+"/Results"+setup->GetTitle()+"HSRooMcmcSeq"+".root";


//I think this fileName will need to be constructed in the  void FitManager::PlotDataModel() function and passed to the MCMCPlotResults :
    MCMCPlotResults::MCMCPlotResults(Setup *setup, const RooDataSet* data, const TString& tag, RooMcmc* mcmc)
    {

    cout<<"MCMCPlotResults  "<<fCanvases.get()<<" "<<setup<<" "<<endl;
    fCanvases->SetName(TString("RFPlots")+setup->GetName());

    auto vars=setup->FitVars();
    auto model=setup->Model();
    cout<<"model "<<model<<endl;

    RooHSEventsPDF_IsPlotting=kTRUE;

    //Get the tree from the mcmc
    auto tree = mcmc->GetTree();
    auto burnIn=mcmc->GetNumBurnInSteps();

    for(auto var : vars)
      {
	auto canName = tag+"_"+var->GetName();
	auto canvas = new TCanvas(canName, canName);
	fCanvases->Add(canvas);

	RooPlot* frame = var->frame();
	data->plotOn(frame, DataError(RooAbsData::SumW2));

	//needs a loop over the parameter values from the mcmc
	//and to draw the model for every nth entry of the mcmc

	//Start here
	auto& pars = setup->ParsAndYields();
	std::cout<<"The parameters are: " << pars<<std::endl;      
	
	//loop over mcmc tree samples
	Int_t Nentries = tree->GetEntries();
	Int_t NthDraw = (Nentries-burnIn)/25;
	Int_t mod = 0; //mod<NthDraw!
	Int_t Npars = pars.size();
	std::cout<<"The number of parameters is:  "<<Npars<<std::endl;
	vector<Double_t> params(Npars);
	Int_t param_index = 0;
 

	for (int ientry = burnIn; ientry<Nentries; ientry++)
	  {//Loop over entries of the tree
	    if(ientry%NthDraw==mod)
	      {//Draw a selection of the entries
		
		Int_t param_index = 0;
		for(RooAbsArg* ipar : pars)
		  {//Loop over parameters
		    
		    tree->SetBranchAddress(ipar->GetName(), &params[param_index]);
		    tree->GetEntry(ientry);
		    //  std::cout<<param_index<<"  "<<ipar->GetName()<<"  "<<params[param_index]<<std::endl;
		    string string1 = ipar->GetName();	     
		    string string2 = "_str";
		    string ipar_str = string1 + string2;
		    
		    if(ipar_str.find("Yld") != std::string::npos)
		      {//If yield, Set Yields
			setup->SetYldVal(ipar->GetName(), params[param_index]);
			//	std::cout<< ipar->GetName()<<"  " <<params[param_index]<<std::endl;
		      }//Close if yields
		    else
		      {//If par, Set pars
			setup->SetParVal(ipar->GetName(), params[param_index]);
			//std::cout<<ipar->GetName()<<"  "<<params[param_index]<<std::endl;
		      }//Close if pars
		    
		    if(param_index==(Npars-1))
		      {//if the parameters have all been set
			//plot the model
			setup->TotalPDF();
			auto model = setup->Model();
			model->plotOn(frame,LineColor(kRed), LineWidth(1)) ;
			const auto& pdfs = setup->PDFs();
			for (Int_t ic = 0; ic<pdfs.getSize(); ic++)
			  {
			    model->plotOn(frame, Components(pdfs[ic]),LineWidth(1), LineColor(ic%8+3));
			  }
			frame->Draw();
		      }
		    param_index++;	 
		  }//Close loop over params
				

	      }//Close if selection of entries    
	  }//Close loop over entries
	
	
	//End here
	canvas->Modified();
	canvas->Update();
	canvas->Draw("");
      }//loop over vars
    }//MCMCPlotResults

  }//FIT
}//H

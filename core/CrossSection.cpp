#include "CrossSection.h"
#include "RooHSEventsPDF.h"
#include "RooHSEventsHistPDF.h"
#include "RooComponentsPDF.h"
#include <TH1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TDirectory.h>
#include <RooRandom.h>
#include <RooStats/RooStatsUtils.h>

#include <memory>

#include <utility>

namespace HS{
namespace FIT{
	
	void CrossSection::Run(){
		
		LoadFitResult();
		
		CreateCurrSetup();
		
		cout << " CrossSection::Run() " << Bins().GetBins().GetNAxis() << " binning variable(s) provided for calculation." << endl;
		if(Bins().GetBins().GetNAxis()>2){
			cout << " CrossSection::Run() More than two different types of bins. Cross section calculation only supports two bins (beam energy and angle)." << endl;
			return;
		}
		if(Bins().GetBins().GetNAxis()>1 && fBeamEnergyBinName==""){
			cout << "More than two bin variables and none of them was set as beam energy. Cross section cannot be calculated." << endl;
			return;
		}
		
		FillEventsPDFs();
		
		CalcFlux();
		CalcYield();
		CalcAcceptanceCorrection();
		CalcCrossSection();
		
	}
	
	void CrossSection::SaveResults(){
		
		TString fileName=Form("%s%s/ResultsCrossSection.root",fCurrSetup->GetOutDir().Data(),GetCurrName().Data());
		cout << "Save to " << fileName << endl;
		auto outfile=std::make_unique<TFile> (fileName,"recreate");
		SetName("cs");
		Write();
	}
	
	void CrossSection::LoadFitResult(){
		if(fResultFileName==TString())
			return;
		
		TString resultFile=fResultDir+Bins().BinName(GetDataBin(GetFiti()))+"/"+fResultFileName;
		std::unique_ptr<TFile> fitFile{TFile::Open(resultFile)};
		std::unique_ptr<RooDataSet> result{dynamic_cast<RooDataSet*>( fitFile->Get(Minimiser::FinalParName()))};
		
		//Set the values of the paramteres to those in the given result
		if(result.get()){
			auto newPars = SetUp().ParsAndYields();
			auto* resAll = result->get(); //get all result info
			auto* resPars=resAll->selectCommon(newPars); //just select pars and yieds
			newPars.assignFast(*resPars); //set values to results
			cout<<"ToyManager::LoadFitResult setting values from fit results "<<resultFile<<" : "<<endl;
			newPars.Print("v");
		}
	}
	
	void CrossSection::SetBeamEnergyBinLimits(TString bin){
		fBeamEnergyBinName = bin;
		vector<Double_t> limits;
		TAxis a = (TAxis)Bins().GetBins().GetAxis(Bins().GetBins().GetAxisi(bin));
		Int_t nbins = a.GetNbins();
		for(Int_t i=1;i<=(nbins+1);i++)
			limits.push_back(a.GetBinLowEdge(i));
		SetBeamEnergyBinLimits(limits);
	}
	
	void CrossSection::CalcFlux(){
		auto fluxfile=std::make_unique<TFile> (fFluxfile,"read");
		TH1F* hFlux = (TH1F*) fluxfile->Get(fFluxhistname)->Clone("flux");
		Int_t nbins = fBeamEnergyBinLimits.size()-1;
		cout << "CrossSection::CalcFlux for " << nbins << " bins." << endl;
		if(nbins==1){ //only one beam energy bin, we can integrate directly
			fBeamEnergyValue = fBeamEnergyBinLimits[0]+(fBeamEnergyBinLimits[1]-fBeamEnergyBinLimits[0])/2.;
			Double_t integral = hFlux->Integral(hFlux->FindBin(fBeamEnergyBinLimits[0]),hFlux->FindBin(fBeamEnergyBinLimits[1])-1);
			cout << "CrossSection::CalcFlux Integrated flux from " << fBeamEnergyBinLimits[0] << " to " << fBeamEnergyBinLimits[1] << " = " << integral << endl;
			fFlux = integral;
		}
		if(nbins>1){ //need to find correct bin first
			TAxis a = (TAxis)Bins().GetBins().GetAxis(Bins().GetBins().GetAxisi(fBeamEnergyBinName));
			//find correct bin limits, feels clumsy, is there a better way to do this?
			Double_t binvalue = -1.;
			TString currname = GetCurrName();
			TObjArray* token = currname.Tokenize("_");
			for(auto i:*token){
				TString namebuffer = ((TObjString*)i)->String();
				if(namebuffer.Contains(fBeamEnergyBinName)){ // pick out part of name that contains fBeamEnergyBinName
					namebuffer.ReplaceAll(fBeamEnergyBinName,""); //remove fBeamEnergyBinName from bin, only central value is left
					binvalue = namebuffer.Atof();
					break;
				}
			}
			Int_t binnumber = a.FindBin(binvalue);
			Double_t lowedge = a.GetBinLowEdge(binnumber);
			Double_t upedge = a.GetBinUpEdge(binnumber);
			
			Double_t integral = hFlux->Integral(hFlux->FindBin(lowedge),hFlux->FindBin(upedge)-1);
			cout << "CrossSection::CalcFlux Integrated flux from " << lowedge << " to " << upedge << " = " << integral << endl;
			fFlux = integral;
			fBeamEnergyValue = binvalue;
		}
	}
	
	void CrossSection::CalcYield(){
		UInt_t idata=GetDataBin(GetFiti());
		Double_t sumofweightsData=Data().Get(idata)->sumEntries();
		fYield = sumofweightsData;
	}
	
	void CrossSection::CalcAcceptanceCorrection(){
		auto pdfs=fCurrSetup->PDFs();
		if(pdfs.getSize()>1)
			cout<< "FitManager::CalcAcceptanceCorrection() Found more than one pdf!!! Last is used as acceptance!!!" << endl;
		for(Int_t ip=0;ip<pdfs.getSize();ip++){
			auto pdf=dynamic_cast<RooHSEventsPDF*>( &pdfs[ip]);
			pdf->Print();
			if(pdf){
				if(Bins().FileNames(pdf->GetName()).size()==0)
					continue;
// 				
				Double_t integralAccepted=pdf->unnormalisedIntegral(1,"");
				Double_t integralGenerated=pdf->unnormalisedIntegral(2,"");
				
				if(integralGenerated)
					cout << "FitManager::CalcAcceptanceCorrection() accepted=" << integralAccepted << " generated=" << integralGenerated << " ratio=" << integralAccepted/integralGenerated << endl;
				else
					cout << "FitManager::CalcAcceptanceCorrection() accepted=" << integralAccepted << " generated=" << integralGenerated << " Can't calculate acceptance!!!" << endl;
				
				if(integralGenerated)
					fAcceptance = (integralAccepted/integralGenerated);
				
			}
		}
	}
	
	void CrossSection::CalcCrossSection(){ //TODO include errors
		if(fAcceptance==0)
			cout << "CrossSection::CalcCrossSection() Acceptance is 0!! Cannot normalise cross section!!" << endl;
		else
			fCrossSection = fYield/fAcceptance;
		
		if(fFlux==0)
			cout << "CrossSection::CalcCrossSection() Flux is 0!! Cannot normalise cross section!!" << endl;
		else
			fCrossSection/=fFlux;
		
		if(fTargetThickness==0)
			cout << "CrossSection::CalcCrossSection() Target thickness is 0!! Cannot normalise cross section!!" << endl;
		else
			fCrossSection/=fTargetThickness;
		
		if(fBranchingRatio==0)
			cout << "CrossSection::CalcCrossSection() Branching ratio is 0!! Cannot normalise cross section!!" << endl;
		else
			fCrossSection/=fBranchingRatio;
		
		Double_t binwidth = 0.;
		Int_t naxis = Bins().GetBins().GetNAxis();
		for(Int_t i=0;i<naxis;i++){
			TAxis a = (TAxis)Bins().GetBins().GetAxis(i);
			TString axisname = a.GetName();
			//need to make sure to get the correct binning variable, allow only two bins in Run(),
			//if it is not fBeamEnergyBinName it must be correct
			if(axisname==fBeamEnergyBinName)
				continue;
			//find correct bin limits, feels clumsy, is there a better way to do this?
			Double_t binvalue = -1.;
			TString currname = GetCurrName();
			TObjArray* token = currname.Tokenize("_");
			for(auto i:*token){
				TString namebuffer = ((TObjString*)i)->String();
				if(namebuffer.Contains(axisname)){ // pick out part of name that contains axisname
					namebuffer.ReplaceAll(axisname,""); //remove axisname from bin, only central value is left
					binvalue = namebuffer.Atof();
					break;
				}
			}
			Int_t binnumber = a.FindBin(binvalue);
			Double_t lowedge = a.GetBinLowEdge(binnumber);
			Double_t upedge = a.GetBinUpEdge(binnumber);
			binwidth = upedge-lowedge;
			fBinValue = binvalue;
		}
		if(binwidth==0)
			cout << "CrossSection::CalcCrossSection() Bin width is 0!! Cannot normalise cross section!!" << endl;
		else
			fCrossSection/=binwidth;
		
	}
	
	void CrossSection::DrawResults(){
		cout << "CrossSection::DrawResults()" << endl;
		CrossSection* a;
		Int_t nbins = Bins().GetSize();
		Double_t csbuffer[nbins];
		Double_t binningbuffer[nbins];
		Double_t ebinning[nbins];
		std::set<Double_t> ebinningset;
		for(Int_t i=0;i<nbins;i++){ //loop over all bins and fill buffer
			TString fileName=Form("%s%s/ResultsCrossSection.root",SetUp().GetOutDir().Data(),Bins().BinName(i).Data());
			auto file=std::make_unique<TFile> (fileName,"read");
			a = (CrossSection*)file->Get("cs");
			cout << a->GetBeamEnergyValue() << " " << a->GetBinValue() << " " << a->GetCrossSection() << endl;
			
			csbuffer[i] = a->GetCrossSection();
			binningbuffer[i] = a->GetBinValue();
			ebinning[i] = a->GetBeamEnergyValue();
			ebinningset.insert(a->GetBeamEnergyValue());
		}
		
		// now sort all results into correct Ebin arrays for easy plotting
		Int_t ebins = ebinningset.size();
		Double_t cs[ebins][nbins/ebins]; //in total nbins of which ebins energy, therefore nbins/ebins is the angular binning
		Double_t binning[ebins][nbins/ebins];
		Int_t ecounter=0;
		for(auto i:ebinningset){
			Int_t counter=0;
			for(Int_t j=0;j<nbins;j++){
				if(i==ebinning[j]){
					cs[ecounter][counter] = csbuffer[j];
					binning[ecounter][counter] = binningbuffer[j];
					counter++;
				}
			}
			ecounter++;
		}
			
		TGraphErrors* gResults[ebins];
		TCanvas* cResults = new TCanvas("results","results");
		cResults->DivideSquare(ebins);
		for(Int_t e=0; e<ebins;e++){
			cResults->cd(e+1);
			gResults[e] = new TGraphErrors(nbins/ebins,binning[e],cs[e]);
			gResults[e]->SetNameTitle(TString::Format("Results%d",e),TString::Format("Results%d",e));
			gResults[e]->Draw("AP");
			gResults[e]->SetMarkerStyle(20);
		}
		
	}
	
}
}
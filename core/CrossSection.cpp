#include "CrossSection.h"
#include "RooHSEventsPDF.h"
#include "RooHSEventsHistPDF.h"
#include "RooComponentsPDF.h"
#include <TH1.h>
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
		
		FillEventsPDFs();
		
		cout << " CrossSection::Run() " << Bins().GetBins().GetNAxis() << " bin(s) provided for calculation." << endl;
		if(Bins().GetBins().GetNAxis()>2){
			cout << " CrossSection::Run() More than two different types of bins. Cross section calculation only supports two bins (beam energy and angle)." << endl;
			return;
		}
		
		CalcYield();
		CalcAcceptanceCorrection();
		CalcCrossSection();
		
	}
	
	void CrossSection::SaveResults(){
		
		TString fileName=Form("%s%s/CrossSection.root",fCurrSetup->GetOutDir().Data(),GetCurrName().Data());
		cout << "Save to " << fileName << endl;
		auto outfile=std::make_unique<TFile> (fileName,"recreate");
		SetName("cs");
		Write();
	}
	
	void CrossSection::LoadFitResult(){
		if(fResultFileName==TString())
			return;
		
		TString resultFile=fResultOutDir+Bins().BinName(GetDataBin(GetFiti()))+"/"+fResultFileName;
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
		vector<Double_t> limits;
		TAxis a = (TAxis)Bins().GetBins().GetAxis(Bins().GetBins().GetAxisi(bin));
		Int_t nbins = a.GetNbins();
		for(Int_t i=1;i<=(nbins+1);i++)
			limits.push_back(a.GetBinLowEdge(i));
		SetBeamEnergyBinLimits(limits);
	}
	
	void CrossSection::LoadFlux(TString filename, TString histname){
		fFlux.clear(); // clear current flux before setting new
		auto fluxfile=std::make_unique<TFile> (filename,"read");
		TH1F* hFlux = (TH1F*) fluxfile->Get(histname)->Clone("flux");
		Int_t nbins = fBeamEnergyBinLimits.size()-1;
		cout << "CrossSection::LoadFlux for " << nbins << " bins." << endl;
		for(Int_t i=0; i<nbins;i++){
			Double_t integral = hFlux->Integral(hFlux->FindBin(fBeamEnergyBinLimits[i]),hFlux->FindBin(fBeamEnergyBinLimits[i+1])-1);
			cout << "CrossSection::LoadFlux Integrated flux from " << fBeamEnergyBinLimits[i] << " to " << fBeamEnergyBinLimits[i+1] << " = " << integral << endl;
			fFlux.push_back(integral);
		}
	}
	
	void CrossSection::CalcYield(){
		UInt_t idata=GetDataBin(GetFiti());
		Double_t sumofweightsData=Data().Get(idata)->sumEntries();
		fYield = sumofweightsData;
	}
	
	void CrossSection::CalcAcceptanceCorrection(){
		UInt_t idata=GetDataBin(GetFiti());
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
	
	void CrossSection::CalcCrossSection(){
		cout << "CrossSection::CalcCrossSection() not yet implemented" << endl;
		
		if(fAcceptance==0)
			cout << "CrossSection::CalcCrossSection() Acceptance is 0!! Cannot calculate cross section!!" << endl;
		else
			fCrossSection = fYield/fAcceptance;
		
// 		if(fFlux==0)
// 			cout << "CrossSection::CalcCrossSection() Flux is 0!! Cannot normalise cross section!!" << endl;
// 		else
// 			fCrossSection/=fFlux;
		
		if(fTargetThickness==0)
			cout << "CrossSection::CalcCrossSection() Target thickness is 0!! Cannot normalise cross section!!" << endl;
		else
			fCrossSection/=fTargetThickness;
	}
	
	void CrossSection::DrawResults(){
		cout << "CrossSection::DrawResults() not yet implemented" << endl;
	}
	
}
}
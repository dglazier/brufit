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
		
		TString fileName=Form("%s%s/ResultsCrossSection.root",fCurrSetup->GetOutDir().Data(),GetCurrName().Data());
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
		fBeamEnergyBinName = bin;
		vector<Double_t> limits;
		TAxis a = (TAxis)Bins().GetBins().GetAxis(Bins().GetBins().GetAxisi(bin));
		Int_t nbins = a.GetNbins();
		for(Int_t i=1;i<=(nbins+1);i++)
			limits.push_back(a.GetBinLowEdge(i));
		SetBeamEnergyBinLimits(limits);
	}
	
	void CrossSection::LoadFlux(TString filename, TString histname){
		auto fluxfile=std::make_unique<TFile> (filename,"read");
		TH1F* hFlux = (TH1F*) fluxfile->Get(histname)->Clone("flux");
		Int_t nbins = fBeamEnergyBinLimits.size()-1;
		cout << "CrossSection::LoadFlux for " << nbins << " bins." << endl;
		if(nbins==1){ //only one beam energy bin, we can integrate directly
			Double_t integral = hFlux->Integral(hFlux->FindBin(fBeamEnergyBinLimits[0]),hFlux->FindBin(fBeamEnergyBinLimits[1])-1);
			cout << "CrossSection::LoadFlux Integrated flux from " << fBeamEnergyBinLimits[0] << " to " << fBeamEnergyBinLimits[1] << " = " << integral << endl;
			fFlux = integral;
		}
		if(nbins>1){ //need to find correct bin first
			
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
			if(axisname==fBeamEnergyBinName) //need to make sure to get the correct binning variable, allow only two bins in Run(), if it is not fBeamEnergyBinName it must be correct
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
		}
		if(binwidth==0)
			cout << "CrossSection::CalcCrossSection() Bin width is 0!! Cannot normalise cross section!!" << endl;
		else
			fCrossSection/=binwidth;
	}
	
	void CrossSection::DrawResults(){
		cout << "CrossSection::DrawResults() not yet implemented" << endl;
	}
	
}
}
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
		
		LoadResult();
		
		CreateCurrSetup();
		
		FillEventsPDFs();
		
		CalcAcceptanceCorrection();
		
		SaveResults();
		
	}
	
	void CrossSection::SaveResults(){
		
		TString fileName="test.root";//Form("%s%s/%s.root",fCurrSetup->GetOutDir().Data(),GetCurrName().Data(),GetCurrTitle().Data());
		auto outfile=std::make_unique<TFile> (fileName,"recreate");
		fAcceptanceTree->Write();
	}
	
	void CrossSection::LoadResult(){
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
			cout<<"ToyManager::LoadResult setting values from fit results "<<resultFile<<" : "<<endl;
			newPars.Print("v");
		}
	}
	
	void CrossSection::CalcAcceptanceCorrection(){
		UInt_t idata=GetDataBin(GetFiti());
		auto pdfs=fCurrSetup->PDFs();
		
		fAcceptanceTree = new TTree("acceptanceCorrection","acceptanceCorrection");
		
		cout<< "FitManager::CalcAcceptanceCorrection() Found " << pdfs.getSize() << "pdfs." << endl;
		for(Int_t ip=0;ip<pdfs.getSize();ip++){
			auto pdf=dynamic_cast<RooHSEventsPDF*>( &pdfs[ip]);
			pdf->Print();
			if(pdf){
				if(Bins().FileNames(pdf->GetName()).size()==0)
					continue;
// 				
				Double_t integralAccepted=pdf->unnormalisedIntegral(1,"");
				Double_t integralGenerated=pdf->unnormalisedIntegral(2,"");
				Double_t sumofweightsData=Data().Get(idata)->sumEntries();
				Double_t acceptanceCorrectedYield=0;
				if(integralGenerated)
					acceptanceCorrectedYield=sumofweightsData/(integralAccepted/integralGenerated);
				
				if(integralGenerated)
					cout << "FitManager::CalcAcceptanceCorrection() accepted=" << integralAccepted << " generated=" << integralGenerated << " ratio=" << integralAccepted/integralGenerated << endl;
				else
					cout << "FitManager::CalcAcceptanceCorrection() accepted=" << integralAccepted << " generated=" << integralGenerated << " Can't calculate acceptance!!!" << endl;
				cout << "FitManager::CalcAcceptanceCorrection() data yield=" << sumofweightsData << endl;
				cout << "FitManager::CalcAcceptanceCorrection() acceptance corrected data yield=" << acceptanceCorrectedYield << endl;
				
				fAcceptanceTree->Branch(pdf->GetName()+TString("_acc"),&integralAccepted);
				fAcceptanceTree->Branch(pdf->GetName()+TString("_gen"),&integralGenerated);
				fAcceptanceTree->Branch(pdf->GetName()+TString("_yld"),&sumofweightsData);
				fAcceptanceTree->Branch(pdf->GetName()+TString("_yldcorr"),&acceptanceCorrectedYield);
				fAcceptanceTree->Fill();
			}
		}
	}
	
}
}
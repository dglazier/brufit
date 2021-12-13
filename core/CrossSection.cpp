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

	Bool_t CrossSection::Run(){

		LoadFitResult();

		CreateCurrSetup();

		cout << " CrossSection::Run() " << Bins().GetBins().GetNAxis() << " binning variable(s) provided for calculation." << endl;
		if(Bins().GetBins().GetNAxis()>2){
			cout << " CrossSection::Run() More than two different types of bins. Cross section calculation only supports two bins (beam energy and angle)." << endl;
			return kTRUE;//note could be false, will leave as true so same behviour as before bool return implemented
		}
		if(Bins().GetBins().GetNAxis()>1 && fBeamEnergyBinName==""){
			cout << "More than two bin variables and none of them was set as beam energy. Cross section cannot be calculated." << endl;
			return kTRUE;
		}

		FillEventsPDFs();

		CalcFlux();
		CalcYield();
		CalcAcceptanceCorrection();
		CalcCrossSection();
		return kTRUE;
	}

	void CrossSection::SaveResults(){

		TString fileName=Form("%s%s/ResultsCrossSection.root",fCurrSetup->GetOutDir().Data(),GetCurrName().Data());
		cout << "Save to " << fileName << endl;
		auto outfile=std::unique_ptr<TFile> (new TFile{fileName,"recreate"});
		SetName("cs");
		if(fSampleAcceptance) fAcceptanceTree->Write();
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
			cout<<"CrossSection::LoadFitResult setting values from fit results "<<resultFile<<" : "<<endl;
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
		std::unique_ptr<TFile> fluxfile{TFile::Open(fFluxfile)};
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
		dset_uptr ds = Data().Get(idata);
		Double_t sumofweightsData=ds->sumEntries();
		//need to calculate sum of weights squared by hand. Use same code fragment as RooFit for sumEntries()
		Double_t sumofweights2Data(0), carry(0);
		Int_t numentries = ds->numEntries();
		for (Int_t i=0 ; i<numentries ; i++) {
			ds->get(i) ;
			Double_t y = ds->weight()*ds->weight() - carry;
			Double_t t = sumofweights2Data + y;
			carry = (t - sumofweights2Data) - y;
			sumofweights2Data = t;
		}
		fYield = sumofweightsData;
		fYield_err = TMath::Sqrt(sumofweights2Data);
		cout<< "CrossSection::CalcYield() Sum of weights = " << fYield << "+/-" << fYield_err << endl;
	}

	void CrossSection::CalcAcceptanceCorrection(){
		Double_t acceptance = 0;
		Double_t acceptance_err = 0;

		auto pdfs=fCurrSetup->PDFs();
		if(pdfs.getSize()>1)
			cout<< "CrossSection::CalcAcceptanceCorrection() Found more than one pdf!!! Last is used as acceptance!!!" << endl;
		for(Int_t ip=0;ip<pdfs.getSize();ip++){
			auto pdf=dynamic_cast<RooHSEventsPDF*>( &pdfs[ip]);
			pdf->Print();
			if(pdf){
				if(Bins().FileNames(pdf->GetName()).size()==0)
					continue;
				if(fSampleAcceptance){
					cout << "CrossSection::CalcAcceptanceCorrection() Sample acceptance from MCMC tree" << endl;
					TString resultFile=fResultDir+Bins().BinName(GetDataBin(GetFiti()))+"/"+fResultFileName;
					std::unique_ptr<TFile> fitFile{TFile::Open(resultFile)};

					//Get MCMC result tree and put in data set
					std::unique_ptr<TTree> resultTree{dynamic_cast<TTree*>( fitFile->Get("MCMCTree"))};//Set the values of the paramteres to those in the given result
					if(!resultTree.get()){
						cout<<"CrossSection::CalcAcceptanceCorrection Couldn't load fit result!" << endl;
						return;
					}
					auto newPars = fCurrSetup->ParsAndYields();
// 					newPars.Print("v");
					RooDataSet mcmcDS("mcmcDS","mcmcDS",resultTree.get(),newPars);
					mcmcDS.Print("v");

					Int_t numentries = mcmcDS.numEntries();
					//fAcceptanceTree->SetNameTitle("acc","acceptance"); //output tree
					//fAcceptanceTree->Branch("acc",&acceptance);
					fAcceptanceTree.reset(new TTree("acc","acceptance"));//changed to shared pointer for copy constructor
					fAcceptanceTree->Branch("acc",&acceptance);
					for(Int_t i=0; i<numentries; i++){
						auto* resAll = mcmcDS.get(i); //get all result info
						newPars.assignFast(*resAll); //set values to results
// 						newPars.Print("v");

						Double_t integralAccepted=pdf->unnormalisedIntegral(1,"");
						Double_t integralGenerated=pdf->unnormalisedIntegral(2,"");
						cout << "CrossSection::CalcAcceptanceCorrection() accepted=" << integralAccepted << " generated=" << integralGenerated << " ratio=" << integralAccepted/integralGenerated << endl;
						acceptance = (integralAccepted/integralGenerated);
						fAcceptanceTree->Fill();
					}
					TH1F hacc("h","h",1,0,1); //dummy for easy mean and stddev calculation
					fAcceptanceTree->Draw("acc>>h","","goff");
					acceptance = hacc.GetMean();
					acceptance_err = hacc.GetStdDev();
				}
				else{
					Double_t integralAccepted=pdf->unnormalisedIntegral(1,"");
					Double_t integralGenerated=pdf->unnormalisedIntegral(2,"");

					if(integralGenerated){
						cout << "CrossSection::CalcAcceptanceCorrection() accepted=" << integralAccepted << " generated=" << integralGenerated << " ratio=" << integralAccepted/integralGenerated << endl;
						acceptance = (integralAccepted/integralGenerated);
					}
					else
						cout << "CrossSection::CalcAcceptanceCorrection() accepted=" << integralAccepted << " generated=" << integralGenerated << " Can't calculate acceptance!!!" << endl;
				}
				fAcceptance = acceptance;
				fAcceptance_err = acceptance_err;
			}
		}
	}

	void CrossSection::CalcCrossSection(){ //TODO include errors
		if(fAcceptance==0)
			cout << "CrossSection::CalcCrossSection() Acceptance is 0!! Cannot normalise cross section!!" << endl;
		else{
			fCrossSection = fYield/fAcceptance;
			fCrossSection_err = fYield_err/fAcceptance;

		}

		if(fFlux==0)
			cout << "CrossSection::CalcCrossSection() Flux is 0!! Cannot normalise cross section!!" << endl;
		else{
			fCrossSection/=fFlux;
			fCrossSection_err/=fFlux;
		}

		if(fTargetThickness==0)
			cout << "CrossSection::CalcCrossSection() Target thickness is 0!! Cannot normalise cross section!!" << endl;
		else{
			fCrossSection/=fTargetThickness;
			fCrossSection_err/=fTargetThickness;
		}

		if(fBranchingRatio==0)
			cout << "CrossSection::CalcCrossSection() Branching ratio is 0!! Cannot normalise cross section!!" << endl;
		else{
			fCrossSection/=fBranchingRatio;
			fCrossSection_err/=fBranchingRatio;
		}

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
			for(auto j:*token){
				TString namebuffer = ((TObjString*)j)->String();
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
		else{
			fCrossSection/=binwidth;
			fCrossSection_err/=binwidth;
		}

	}

	void CrossSection::DrawResults(TString outputfile){
		cout << "CrossSection::DrawResults()" << endl;
		CrossSection* a;
		Int_t nbins = Bins().GetSize();
		Double_t csbuffer[nbins];
		Double_t cserrbuffer[nbins];
		Double_t yieldbuffer[nbins];
		Double_t yielderrbuffer[nbins];
		Double_t acccorryieldbuffer[nbins];
		Double_t acccorryielderrbuffer[nbins];
		Double_t fluxnormyieldbuffer[nbins];
		Double_t fluxnormielderrbuffer[nbins];
		Double_t binningbuffer[nbins];
		Double_t ebinning[nbins];
		std::set<Double_t> ebinningset;
		for(Int_t i=0;i<nbins;i++){ //loop over all bins and fill buffer
			TString fileName=Form("%s%s/ResultsCrossSection.root",SetUp().GetOutDir().Data(),Bins().BinName(i).Data());
			std::unique_ptr<TFile> file{TFile::Open(fileName)};
			if(file==nullptr){ //check file exists
				cout << fileName << " couldn't be opened. Does it exist?" << endl;
				csbuffer[i] = 0;
				cserrbuffer[i] = 0;
				binningbuffer[i] = 0;
				ebinning[i] = 0;
				continue;
			}
			a = (CrossSection*)file->Get("cs");
			cout << a->GetBeamEnergyValue() << " " << a->GetBinValue() << " " << a->GetCrossSection() << "+/-" << a->GetCrossSection_err() << " " << a->GetAcceptance() << "+/-" << a->GetAcceptance_err() << endl;

			csbuffer[i] = a->GetCrossSection();
			cserrbuffer[i] = a->GetCrossSection_err();

			Double_t thickness = a->GetTargetThickness();
			Double_t br = a->GetBranchingRatio();
			yieldbuffer[i] = a->GetYield();
			yielderrbuffer[i] = a->GetYield_err();
			acccorryieldbuffer[i] = yieldbuffer[i]/a->GetAcceptance();
			acccorryielderrbuffer[i] = yielderrbuffer[i]/a->GetAcceptance();
			fluxnormyieldbuffer[i] = yieldbuffer[i]/a->GetFlux();
			fluxnormielderrbuffer[i] = yielderrbuffer[i]/a->GetFlux();

			binningbuffer[i] = a->GetBinValue();
			ebinning[i] = a->GetBeamEnergyValue();
			ebinningset.insert(a->GetBeamEnergyValue());
		}

		// now sort all results into correct Ebin arrays for easy plotting
		Int_t ebins = ebinningset.size();
		Double_t cs[ebins][nbins/ebins]; //in total nbins of which ebins energy, therefore nbins/ebins is the angular binning
		Double_t cs_err[ebins][nbins/ebins]; //in total nbins of which ebins energy, therefore nbins/ebins is the angular binning
		Double_t yield[ebins][nbins/ebins];
		Double_t yield_err[ebins][nbins/ebins];
		Double_t acccorryield[ebins][nbins/ebins];
		Double_t acccorryield_err[ebins][nbins/ebins];
		Double_t fluxnormyield[ebins][nbins/ebins];
		Double_t fluxnormield_err[ebins][nbins/ebins];
		Double_t binning[ebins][nbins/ebins];
		Double_t energybins[ebins];
		Int_t ecounter=0;
		for(auto i:ebinningset){
			Int_t counter=0;
			for(Int_t j=0;j<nbins;j++){
				if(i==ebinning[j]){
					cs[ecounter][counter] = csbuffer[j];
					cs_err[ecounter][counter] = cserrbuffer[j];
					yield[ecounter][counter] = yieldbuffer[j];
					yield_err[ecounter][counter] = yielderrbuffer[j];
					acccorryield[ecounter][counter] = acccorryieldbuffer[j];
					acccorryield_err[ecounter][counter] = acccorryielderrbuffer[j];
					fluxnormyield[ecounter][counter] = fluxnormyieldbuffer[j];
					fluxnormield_err[ecounter][counter] = fluxnormielderrbuffer[j];
					binning[ecounter][counter] = binningbuffer[j];
					energybins[ecounter]=i;
					counter++;
				}
			}
			ecounter++;
		}

		TString axisname;
		Int_t naxis = Bins().GetBins().GetNAxis();
		for(Int_t i=0;i<naxis;i++){
			TAxis axis = (TAxis)Bins().GetBins().GetAxis(i);
			//need to make sure to get the correct binning variable, allow only two bins in Run(),
			//if it is not fBeamEnergyBinName it must be correct
			if(axis.GetName()==fBeamEnergyBinName)
				continue;
			axisname = axis.GetName();
		}

		TGraphErrors* gResults[ebins];
		TCanvas* cResults = new TCanvas("results","results");
		cResults->DivideSquare(ebins);
		for(Int_t e=0; e<ebins;e++){
			cResults->cd(e+1);
			gResults[e] = new TGraphErrors(nbins/ebins,binning[e],cs[e],0,cs_err[e]);
			gResults[e]->SetNameTitle(TString::Format("Results%d",e),fBeamEnergyBinName+TString::Format("=%f",energybins[e]));
			gResults[e]->Draw("AP");
			gResults[e]->SetMarkerStyle(20);
			gResults[e]->GetXaxis()->SetTitle(axisname);
		}

		TGraphErrors* gYields[ebins];
		TCanvas* cYields = new TCanvas("yields","yields");
		cYields->DivideSquare(ebins);
		for(Int_t e=0; e<ebins;e++){
			cYields->cd(e+1);
			gYields[e] = new TGraphErrors(nbins/ebins,binning[e],yield[e],0,yield_err[e]);
			gYields[e]->SetNameTitle(TString::Format("Yields%d",e),fBeamEnergyBinName+TString::Format("=%f",energybins[e]));
			gYields[e]->Draw("AP");
			gYields[e]->SetMarkerStyle(20);
			gYields[e]->GetXaxis()->SetTitle(axisname);
		}

		TGraphErrors* gAccCorrYields[ebins];
		TCanvas* cAccCorrYields = new TCanvas("acccorryields","acccorryields");
		cAccCorrYields->DivideSquare(ebins);
		for(Int_t e=0; e<ebins;e++){
			cAccCorrYields->cd(e+1);
			gAccCorrYields[e] = new TGraphErrors(nbins/ebins,binning[e],acccorryield[e],0,acccorryield_err[e]);
			gAccCorrYields[e]->SetNameTitle(TString::Format("AccCorrYields%d",e),fBeamEnergyBinName+TString::Format("=%f",energybins[e]));
			gAccCorrYields[e]->Draw("AP");
			gAccCorrYields[e]->SetMarkerStyle(20);
			gAccCorrYields[e]->GetXaxis()->SetTitle(axisname);
		}

		TGraphErrors* gFluxNormYields[ebins];
		TCanvas* cFluxNormYields = new TCanvas("fluxnormyields","fluxnormyields");
		cFluxNormYields->DivideSquare(ebins);
		for(Int_t e=0; e<ebins;e++){
			cFluxNormYields->cd(e+1);
			gFluxNormYields[e] = new TGraphErrors(nbins/ebins,binning[e],fluxnormyield[e],0,fluxnormield_err[e]);
			gFluxNormYields[e]->SetNameTitle(TString::Format("FluxNormYields%d",e),fBeamEnergyBinName+TString::Format("=%f",energybins[e]));
			gFluxNormYields[e]->Draw("AP");
			gFluxNormYields[e]->SetMarkerStyle(20);
			gFluxNormYields[e]->GetXaxis()->SetTitle(axisname);
		}

		if(outputfile != ""){
			cout << "Save to " << outputfile << endl;
			auto outfile=std::unique_ptr<TFile> (new TFile{outputfile,"recreate"});
			if(outfile){
				cResults->Write();
				cYields->Write();
				cAccCorrYields->Write();
				cFluxNormYields->Write();
				for(Int_t e=0; e<ebins;e++){
					gResults[e]->Write();
					gYields[e]->Write();
					gAccCorrYields[e]->Write();
					gFluxNormYields[e]->Write();
				}
			}
		}

	}

}
}

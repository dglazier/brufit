////////////////////////////////////////////////////////////////
///
///Class:               CrossSection
///Description:
///

#pragma once

#include <utility>


#include "FitManager.h"
#include "Setup.h"
#include "Weights.h"

namespace HS{
namespace FIT{

	using tree_uptr =std::unique_ptr<TTree>;
	using tree_shptr =std::shared_ptr<TTree>;
	using strings_t = std::vector<TString>;

	class CrossSection  : public FitManager{

	public:
		CrossSection()=default;
		CrossSection(const CrossSection&)=default;
		CrossSection(const FitManager& fm,TString outDir="",TString resultFile=""):FitManager(fm),fResultDir(std::move(std::move(outDir))),fResultFileName(std::move(std::move(resultFile))){};
		CrossSection(CrossSection&&)=default;
		~CrossSection() override =default;
		CrossSection& operator=(const CrossSection& other) = default;
		CrossSection& operator=(CrossSection&& other) = delete;

		Bool_t Run() override;
		void SaveResults() override;

		void LoadFitResult();

		void SetBeamEnergyBinLimits(std::vector<Double_t> limits){fBeamEnergyBinLimits = limits;};
		void SetBeamEnergyBinLimits(TString bin);
		void SetFlux(TString filename, TString histname){fFluxfile=filename; fFluxhistname=histname;};
		void SetTargetThickness(Double_t n){fTargetThickness=n;};
		void SetBranchingRatio(Double_t n){fBranchingRatio=n;};
		void SampleAcceptance(Bool_t b = kTRUE){fSampleAcceptance=b;};


		void SetResultDir(TString name){fResultDir=std::move(name);};
		void SetResultFileName(TString name){fResultFileName=std::move(name);};

		void CalcFlux();
		void CalcYield();
		void CalcAcceptanceCorrection();
		void CalcCrossSection();
		void DrawResults(TString outputfile = "");

		Double_t GetFlux(){return fFlux;};
		Double_t GetTargetThickness(){return fTargetThickness;};
		Double_t GetBranchingRatio(){return fBranchingRatio;};
		Double_t GetAcceptance(){return fAcceptance;};
		Double_t GetAcceptance_err(){return fAcceptance_err;};
		Double_t GetYield(){return fYield;};
		Double_t GetYield_err(){return fYield_err;};
		Double_t GetCrossSection(){return fCrossSection;};
		Double_t GetCrossSection_err(){return fCrossSection_err;};
		Double_t GetBinValue(){return fBinValue;};
		Double_t GetBeamEnergyValue(){return fBeamEnergyValue;};
		Double_t GetBeamEnergyNBins(){return fBeamEnergyBinLimits.size()-1;};

	protected:

	private:

		TString fResultDir;
		TString fResultFileName;

		std::vector<Double_t> fBeamEnergyBinLimits = {0};
		TString fBeamEnergyBinName = "";
		TString fFluxfile = "";
		TString fFluxhistname = "";
		Double_t fFlux = 0.;
		Double_t fTargetThickness = 1.; //inverse barn
		Double_t fBranchingRatio = 1;
		Double_t fAcceptance = 0.;
		Double_t fAcceptance_err = 0.;
		Double_t fYield = 0.;
		Double_t fYield_err = 0.;
		Double_t fCrossSection = 0.;
		Double_t fCrossSection_err = 0.;
		Double_t fBinValue = 0.;
		Double_t fBeamEnergyValue = 0.;

		Bool_t fSampleAcceptance = kFALSE;
		tree_shptr fAcceptanceTree;

	ClassDefOverride(HS::FIT::CrossSection,1);
	}; //class CrossSection

}//namespace FIT
}//namespace HS

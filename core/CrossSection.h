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
	using strings_t = std::vector<TString>;

	class CrossSection  : public FitManager{
	
	public:
		CrossSection()=default;
		CrossSection(const CrossSection&)=default;
		CrossSection(const FitManager& fm,TString outDir="",TString resultFile=""):FitManager(fm),fResultOutDir(std::move(std::move(outDir))),fResultFileName(std::move(std::move(resultFile))){};
		CrossSection(CrossSection&&)=default;
		~CrossSection() override =default;
		CrossSection& operator=(const CrossSection& other) = default;
		CrossSection& operator=(CrossSection&& other) = default;
		
		void Run() override;
		void SaveResults() override;
		
		void LoadFitResult();
		
		void SetBeamEnergyBinLimits(std::vector<Double_t> limits){fBeamEnergyBinLimits = limits;};
		void SetBeamEnergyBinLimits(TString bin);
		void LoadFlux(TString filename, TString histname);
		void SetTargetThickness(Double_t n){fTargetThickness=n;};
		void SetBranchingRatio(Double_t n){fBranchingRatio=n;};
		
		
		void SetResultOutdir(TString name){fResultOutDir=std::move(name);};
		void SetResultFileName(TString name){fResultFileName=std::move(name);};
		
		void CalcYield();
		void CalcAcceptanceCorrection();
		void CalcCrossSection();
		void DrawResults();
		
		Double_t GetFlux(){return fFlux;};
		Double_t GetTargetThickness(){return fTargetThickness;};
		Double_t GetBranchingRatio(){return fBranchingRatio;};
		Double_t GetAcceptance(){return fAcceptance;};
		Double_t GetYield(){return fYield;};
		Double_t GetCrossSection(){return fCrossSection;};
		
		
	protected:
	
	private:

		TString fResultOutDir;
		TString fResultFileName;
		
		std::vector<Double_t> fBeamEnergyBinLimits = {0};
		TString fBeamEnergyBinName = "";
		Double_t fFlux = 0.;
		Double_t fTargetThickness = 1.; //inverse barn
		Double_t fBranchingRatio = 1;
		Double_t fAcceptance = 0.;
		Double_t fYield = 0.;
		Double_t fCrossSection = 0.;

	ClassDefOverride(HS::FIT::CrossSection,1);
	}; //class CrossSection
	
}//namespace FIT
}//namespace HS


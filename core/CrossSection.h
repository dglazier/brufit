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
		
		void LoadResult();
		
		void SetBeamEnergyBinLimits(std::vector<Double_t> limits){fBeamEnergyBinLimits = limits;};
		void SetBeamEnergyBinLimits(TString bin);
		void LoadFlux(TString filename, TString histname);
		void SetTargetThickness(Double_t n){fTargetThickness=n;};
		
		void SetResultOutdir(TString name){fResultOutDir=std::move(name);}
		void SetResultFileName(TString name){fResultFileName=std::move(name);}
		
		void CalcAcceptanceCorrection();
		
		void CalcCrossSection();
		void DrawResults();
		
	protected:
	
	private:

		TString fResultOutDir;
		TString fResultFileName;
		
		TTree* fAcceptanceTree=nullptr;
		
		std::vector<Double_t> fBeamEnergyBinLimits = {0};
		std::vector<Double_t> fFlux = {0};
		Double_t fTargetThickness = 1.; //inverse barn

	ClassDefOverride(HS::FIT::CrossSection,1);
	}; //class CrossSection
	
}//namespace FIT
}//namespace HS


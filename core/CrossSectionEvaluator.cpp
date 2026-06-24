#include "CrossSectionEvaluator.h"
#include "BruEventsPDF.h"
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TDirectory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMultiVarGaussian.h>
#include <ROOT/TProcessExecutor.hxx>
#include <ROOT/TSeq.hxx>
#include <iostream>
#include <algorithm>

namespace HS {
namespace FIT {

    Bool_t CrossSectionEvaluator::LoadDefaultBins() {
        if (Bins().GetBins().GetNAxis() > 0) return kTRUE; 

        TString binFile = fResultDir + "/DataBinsConfig.root";
        std::unique_ptr<TFile> f(TFile::Open(binFile));
        if (!f || f->IsZombie()) {
            std::cerr << " -> Error: Could not open " << binFile << std::endl;
            return kFALSE;
        }

        for (auto key : *f->GetListOfKeys()) {
            TObject* obj = f->Get(key->GetName());
            if (auto b = dynamic_cast<HS::FIT::Bins*>(obj)) {
                Bins().GetBins() = *b; 
                return kTRUE;
            }
        }
        return kFALSE;
    }

    void CrossSectionEvaluator::InitTree() {
        if (fOutTree) { delete fOutTree; fOutTree = nullptr; }
        
        fOutTree = new TTree("CrossSection", "Cross Section Results");
        fOutTree->Branch("GlobalBin", &fOutGlobalBin, "GlobalBin/I");
        fOutTree->Branch("CrossSection", &fOutCrossSection, "CrossSection/D");
        fOutTree->Branch("CrossSection_err", &fOutCrossSectionErr, "CrossSection_err/D");
        fOutTree->Branch("Yield", &fOutYield, "Yield/D");
        fOutTree->Branch("Yield_err", &fOutYieldErr, "Yield_err/D");
        fOutTree->Branch("Acceptance", &fOutAcceptance, "Acceptance/D");
        fOutTree->Branch("Acceptance_err", &fOutAcceptanceErr, "Acceptance_err/D");

        Int_t naxis = Bins().GetBins().GetNAxis();
        fBinCenters.resize(naxis, 0.0);
        fBinWidths.resize(naxis, 0.0);

        for (Int_t i = 0; i < naxis; i++) {
            TString axisName = Bins().GetBins().GetAxis(i).GetName();
            fOutTree->Branch(axisName + "_center", &fBinCenters[i], axisName + "_center/D");
            fOutTree->Branch(axisName + "_width", &fBinWidths[i], axisName + "_width/D");
        }
    }

    // Run() evaluates strictly the SINGLE CURRENT BIN assigned by Here::Go
    Bool_t CrossSectionEvaluator::Run() {
        if (!fCSFormula) return kFALSE;
        if (!LoadDefaultBins()) return kFALSE;

        Int_t globalBinIndex = GetDataBin(GetFiti());
        
        std::cout << "\n==================================================" << std::endl;
        std::cout << "--- Processing Bin " << globalBinIndex << " [" << GetCurrName() << "] ---" << std::endl;

        // 1. Create transient physics Setup for this bin
        CreateCurrSetup();
        
        // 2. Load True Parameters into active Setup
        if (!LoadFitResult(globalBinIndex)) {
            std::cerr << " -> Failed to load fit result." << std::endl;
            return kFALSE;
        }

        // 3. Connect Data for THIS bin natively
        fCurrDataSet = std::move(Data().Get(GetFiti()));
        FillEventsPDFs();

        InitTree();

        CSData binData;
        binData.kinematics = Bins().GetBins().GetBinDimensions(globalBinIndex);

        CalcYield(globalBinIndex, binData);
        CalcAcceptance(globalBinIndex, binData);

        auto result = fCSFormula(binData);

        fOutGlobalBin = globalBinIndex;
        fOutYield = binData.yield;
        fOutYieldErr = binData.yield_err;
        fOutAcceptance = binData.acceptance;
        fOutAcceptanceErr = binData.acceptance_err;
        fOutCrossSection = result.first;
        fOutCrossSectionErr = result.second;

        Int_t naxis = Bins().GetBins().GetNAxis();
        for (Int_t i = 0; i < naxis; i++) {
            TString axisName = Bins().GetBins().GetAxis(i).GetName();
            fBinCenters[i] = binData.kinematics.at(axisName).center;
            fBinWidths[i] = binData.kinematics.at(axisName).width;
        }

        std::cout << " -> Bin Summary:" << std::endl;
        std::cout << "    * Yield:      " << fOutYield << " +/- " << fOutYieldErr << std::endl;
        std::cout << "    * Acceptance: " << fOutAcceptance << " +/- " << fOutAcceptanceErr << std::endl;
        std::cout << "    * CrossSec:   " << fOutCrossSection << " +/- " << fOutCrossSectionErr << std::endl;
        std::cout << "==================================================\n" << std::endl;

        fOutTree->Fill();
        return kTRUE;
    }

    Bool_t CrossSectionEvaluator::LoadFitResult(Int_t globalBinIndex) {
        if (fResultFileName == TString()) return kFALSE;

        TString resultFile = fResultDir + Bins().BinName(globalBinIndex) + "/" + fResultFileName;
        std::unique_ptr<TFile> fitFile(TFile::Open(resultFile));
        if (!fitFile || fitFile->IsZombie()) return kFALSE;

        std::unique_ptr<RooDataSet> result{dynamic_cast<RooDataSet*>(fitFile->Get(Minimiser::FinalParName()))};
        if (result.get()) {
            auto newPars = fCurrSetup->ParsAndYields(); 
            auto* resAll = result->get(); 
            auto* resPars = resAll->selectCommon(newPars); 
            newPars.assignFast(*resPars); 
            delete resPars;
        } else return kFALSE;
        
        return kTRUE;
    }

    void CrossSectionEvaluator::CalcYield(Int_t globalBinIndex, CSData& binData) {
        if (!fCurrDataSet) return;

        if (!fCurrDataSet->isWeighted()) {
            binData.yield = fCurrDataSet->numEntries();
            binData.yield_err = TMath::Sqrt(binData.yield);
        } else {
            Double_t sumofweightsData = fCurrDataSet->sumEntries();
            Double_t sumofweights2Data(0), carry(0);
            
            Int_t numentries = fCurrDataSet->numEntries();
            for (Int_t i = 0 ; i < numentries ; i++) {
                fCurrDataSet->get(i);
                Double_t w = fCurrDataSet->weight();
                Double_t y = (w * w) - carry;
                Double_t t = sumofweights2Data + y;
                carry = (t - sumofweights2Data) - y;
                sumofweights2Data = t;
            }
            binData.yield = sumofweightsData;
            binData.yield_err = TMath::Sqrt(sumofweights2Data);
        }
    }

    void CrossSectionEvaluator::CalcAcceptance(Int_t globalBinIndex, CSData& binData) {
        auto pdfs = fCurrSetup->PDFs();
        bru::BruEventsPDF* pdf = nullptr;

        for (Int_t ip = 0; ip < pdfs.getSize(); ip++) {
            if ((pdf = dynamic_cast<bru::BruEventsPDF*>(&pdfs[ip]))) break; 
        }
        if (!pdf) return;

        // --- 1. Compute Central Values ---
        Double_t integralAccepted = pdf->unnormalisedIntegral(1, ""); 
        Double_t integralGenerated = 0;

        if (fAccMode == AcceptanceMode::kFullPDF) {
            integralGenerated = pdf->unnormalisedIntegral(2, "") * fGenScale;
        } else if (fAccMode == AcceptanceMode::kPhaseSpaceScale) {
            Double_t genYield = static_cast<Double_t>(pdf->GetNMCGenEntries());
            integralGenerated = genYield * fGenScale;
        }

        if (integralGenerated > 0) {
            binData.acceptance = integralAccepted / integralGenerated;
            std::cout << " -> Computed Central Acceptance: " << binData.acceptance 
                      << " (Acc=" << integralAccepted << ", GenScaled=" << integralGenerated << ")" << std::endl;
        } else {
            std::cerr << " -> Generated integral is zero! Acceptance will be 0." << std::endl;
            return;
        }

        // --- 2. Compute MCMC/Minuit Parameter Shape Variance ---
        Double_t mcmc_variance = 0.0;
        
        if (fSampleAcceptance) {
            auto newPars = fCurrSetup->ParsAndYields();
            
            std::vector<std::vector<Double_t>> paramSamples;
            std::vector<RooRealVar*> paramVars;
            std::vector<TString> paramNames;
            Int_t nSamplesToProcess = 0;

            // ISOLATION SCOPE
            {
                TString resultFile = fResultDir + Bins().BinName(globalBinIndex) + "/" + fResultFileName;
                std::unique_ptr<TFile> fitFile(TFile::Open(resultFile));
                
                if (fitFile && !fitFile->IsZombie()) {
                    TTree* mcmcTree = dynamic_cast<TTree*>(fitFile->Get("MCMCTree"));
                    RooFitResult* minuitResult = dynamic_cast<RooFitResult*>(fitFile->Get("MinuitResult"));
                    std::unique_ptr<RooDataSet> sampleDS;

                    if (mcmcTree) {
                        std::cout << "    -> Extracting parameter samples from MCMCTree." << std::endl;
                        sampleDS.reset(new RooDataSet("mcmcDS", "mcmcDS", newPars, RooFit::Import(*mcmcTree)));
                    } else if (minuitResult && minuitResult->floatParsFinal().getSize() > 0) {
                        std::cout << "    -> Extracting " << fNSamplesMinuit << " parameter samples from Minuit Covariance." << std::endl;
                        RooArgList floatPars = minuitResult->floatParsFinal();
                        RooMultiVarGaussian mvg("mvg", "mvg", floatPars, *minuitResult);
                        RooArgSet genVars(floatPars); 
                        sampleDS.reset(static_cast<RooDataSet*>(mvg.generate(genVars, fNSamplesMinuit)));
                    }

                    if (sampleDS && sampleDS->numEntries() > 0) {
                        Int_t nTotal = sampleDS->numEntries();
                        nSamplesToProcess = std::min(nTotal, fNSamplesMinuit);
                        Int_t step = std::max(1, nTotal / nSamplesToProcess);
                        
                        std::cout << "    -> Spacing extraction for " << nSamplesToProcess << " samples (step size: " << step << ")." << std::endl;

                        const RooArgSet* firstRow = sampleDS->get(0);
                        for (auto* arg : *firstRow) {
                            if (arg) {
                                paramNames.push_back(arg->GetName());
                                paramVars.push_back(dynamic_cast<RooRealVar*>(newPars.find(arg->GetName())));
                            }
                        }

                        paramSamples.resize(nSamplesToProcess, std::vector<Double_t>(paramNames.size()));

                        Int_t extractedCount = 0;
                        for (Int_t i = 0; i < nTotal && extractedCount < nSamplesToProcess; i += step) {
                            const RooArgSet* row = sampleDS->get(i);
                            for (size_t p = 0; p < paramNames.size(); p++) {
                                paramSamples[extractedCount][p] = row->getRealValue(paramNames[p]);
                            }
                            extractedCount++;
                        }
                        nSamplesToProcess = extractedCount; 
                    }
                } 
            } 
            
            // Execute ROOT::TProcessExecutor Multiprocessing Pool
            if (nSamplesToProcess > 0) {
                std::cout << "    -> Parallelizing " << nSamplesToProcess << " acceptance integrals using ROOT::TProcessExecutor..." << std::endl;
                
                RooArgSet* originalPars = (RooArgSet*)newPars.snapshot();
                ROOT::TProcessExecutor pool(fNThreads);
                
                auto sample_func = [&](Int_t i) -> Double_t {
                    
                    for (size_t p = 0; p < paramNames.size(); p++) {
                        if (paramVars[p]) paramVars[p]->setVal(paramSamples[i][p]);
                    }
                    pdf->getVal(); // Flush RooFit caches
                    
                    Double_t tmpAcc = pdf->unnormalisedIntegral(1, "");
                    Double_t tmpGen = (fAccMode == AcceptanceMode::kFullPDF) 
                                      ? pdf->unnormalisedIntegral(2, "") * fGenScale
                                      : integralGenerated; 
                    
                    if (tmpGen > 0) return tmpAcc / tmpGen;
                    return 0.0;
                };
                
                auto accSamples = pool.Map(sample_func, ROOT::TSeqI(0, nSamplesToProcess));
                
                Double_t mean = 0;
                for (Double_t v : accSamples) mean += v;
                mean /= accSamples.size();
                for (Double_t v : accSamples) mcmc_variance += (v - mean) * (v - mean);
                mcmc_variance /= (accSamples.size() > 1 ? accSamples.size() - 1 : 1);

                std::cout << "    -> MCMC Parameter Acceptance Variance: " << mcmc_variance << std::endl;

                // Restore central parameters
                newPars.assignFast(*originalPars);
                delete originalPars;
                pdf->getVal(); 
            } 
        } 
        
        // --- 3. Compute Base MC Statistical Variance ---
        Double_t stat_variance = 0.0;
        Double_t rawGenCount = static_cast<Double_t>(pdf->GetNMCGenEntries());
        
        if (rawGenCount > 0) {
            Double_t rawAcceptance = integralAccepted / (integralGenerated / fGenScale); 
            Double_t raw_err = 0.0;
            
            if (rawAcceptance > 0 && rawAcceptance < 1.0) {
                raw_err = TMath::Sqrt(rawAcceptance * (1.0 - rawAcceptance) / rawGenCount);
            } else {
                raw_err = rawAcceptance / TMath::Sqrt(rawGenCount);
            }
            
            stat_variance = (raw_err / fGenScale) * (raw_err / fGenScale);
            std::cout << "    -> Base MC Statistical Variance: " << stat_variance << std::endl;
        }

        binData.acceptance_err = TMath::Sqrt(mcmc_variance + stat_variance);
    }

    void CrossSectionEvaluator::SaveResults() {
        // Correctly resolving the private field error by utilizing the Setup() public accessor
        TString fileName = Form("%s%s/ResultsCrossSection.root", SetUp().GetOutDir().Data(), GetCurrName().Data());
        std::cout << "--- Saving bin result to " << fileName << " ---" << std::endl;
        
        auto outfile = std::unique_ptr<TFile>(new TFile{fileName, "recreate"});
        if (outfile && fOutTree) {
            fOutTree->SetDirectory(outfile.get());
            fOutTree->Write();
        }
        fOutTree = nullptr; // File takes ownership cleanly
    }

} // namespace FIT
} // namespace HS

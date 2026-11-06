#include "BruMcmc.h"
#include "BruMappedRhat.h" 
#include "BruComponentsPDF.h"
#include "BruMetropolisHastings.h"

#include <TROOT.h>
#include <TIterator.h>
#include <TDirectory.h>
#include <TLeaf.h>
#include <TTreeIndex.h>
#include <RooStats/SequentialProposal.h>
#include <RooStats/ProposalHelper.h>
#include <TRobustEstimator.h>
#include <TPrincipal.h>
#include <TMatrixD.h>
#include <TH1D.h>
#include <TMatrixDSym.h>
#include <RooGlobalFunc.h>
#include <RooProduct.h>
#include <RooRealVar.h>
#include <RooAddition.h>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace HS{
  namespace FIT{

    using namespace RooFit;
    using namespace RooStats;
    
    void BruMcmc::InitModel() {
        std::cout << "BruMcmc::InitModel" << std::endl;
        fPdf = fSetup->Model();
        
        fPOI.removeAll();
        fPOI.add(fSetup->Parameters());
        fPOI.add(fSetup->Yields());
        
        fNuisParams.removeAll();
        fConditionalObs.removeAll();
        fGlobalObs.removeAll();
        fPriorPdf = nullptr; 

        if (fSetup->ApplyMCVariance()) {
            std::cout << "BruMcmc: Globally activating exact MC Variance caching." << std::endl;
            RooArgSet* comps = fPdf->getComponents();
            for (auto* obj : *comps) {
                if (obj->InheritsFrom("bru::BruEventsPDF")) {
                    static_cast<bru::BruEventsPDF*>(obj)->TrackMCVariance(kTRUE);
                }
            }
            delete comps;
        }
    }

    void BruMcmc::Run(Setup &setup,RooAbsData &fitdata){
      fSetup=&setup;
      fData=&fitdata;
      std::cout<<"BruMcmc::Run"<<std::endl;
      
      InitModel();
      SetupBasicUsage();

      MakeChain();
      CalculateChainNLLStats();
    }
    
    RooAbsReal* BruMcmc::BuildNLL(RooAbsData* data,
        std::unique_ptr<RooAbsReal>& baseNll,
        std::unique_ptr<RooRealVar>& alphaVar,
        std::unique_ptr<RooProduct>& correctedNll,
        std::unique_ptr<RooAbsReal>& constraintNll,
        std::unique_ptr<RooAddition>& totalNll) 
    {
        auto foptions = fSetup->FitOptions();
        TObject* opt=nullptr;
        if((opt=foptions.find("Save"))!=nullptr) foptions.Remove(opt);
        if((opt=foptions.find("SumW2Error"))!=nullptr) foptions.Remove(opt);
        
        auto cmd1 = RooFit::ConditionalObservables(fConditionalObs);
        foptions.Add(dynamic_cast<RooCmdArg*>(&cmd1));
        
        auto cmd2 = RooFit::Constrain(RooArgSet()); 
        foptions.Add(dynamic_cast<RooCmdArg*>(&cmd2));

        baseNll.reset(fPdf->createNLL(*data, foptions));
        baseNll->constOptimizeTestStatistic(RooAbsArg::Activate, false);

        Double_t combinedScale = 1.0;
        fLastBetaFactor = 1.0;
        
        if (!_isBatchMode && fSetup->ApplyMCVariance()) {
            if (fCachedBetaVal > 0.0) {
                // Hard-staged routing: reuse frozen Beta factor during Phase 3 tuning & Phase 4 sampling
                fLastBetaFactor = fCachedBetaVal;
                combinedScale *= fCachedBetaVal;
            } else {
                Double_t max_sigma_rel = 0.0;
                TString domName = "None";

                baseNll->getVal(); 

                RooArgSet* allComps = fPdf->getComponents();
                for (auto* obj : *allComps) {
                    if (obj->InheritsFrom("bru::BruEventsPDF")) {
                        auto* evPdf = static_cast<bru::BruEventsPDF*>(obj);
                        
                        Double_t sig = evPdf->GetRelativeVariance();
                        if (sig > max_sigma_rel) { 
                            max_sigma_rel = sig; 
                            domName = evPdf->GetName(); 
                        }
                    }
                }
                delete allComps;
                
                if (max_sigma_rel > 0.0) {
                    Double_t N_data = data->sumEntries(); 
                    Double_t betaVal = 1.0 / (1.0 + N_data * (max_sigma_rel * max_sigma_rel));
                    fLastBetaFactor = betaVal;
                    combinedScale *= betaVal;
                    
                    // Lock Beta once Phase 2 mapping has finished and Phase 3/4 begins
                    if (!_isCovarianceMode) {
                        fCachedBetaVal = betaVal;
                    }
                    
                    std::cout << "\n=========================================" << std::endl;
                    std::cout << " [BruMcmc Likelihood Scaling] " << std::endl;
                    std::cout << " -> Dominant PDF : " << domName << " (Rel Err: " << max_sigma_rel * 100.0 << " %)" << std::endl;
                    std::cout << " -> Beta Factor  : " << betaVal << (fCachedBetaVal > 0.0 ? " [LOCKED]" : "") << std::endl;
                    std::cout << "=========================================\n" << std::endl;
                }
            }
        }
        
        if(data->isNonPoissonWeighted() && fCorrectForWeights){
            Double_t SumW = SumWeights();
            Double_t SumW2 = SumWeights2();
            Double_t alphaVal = SumW / SumW2;
            combinedScale *= alphaVal;

            if (!_isBatchMode) {
                std::cout << " -> Weights Alpha: " << alphaVal << std::endl;
                std::cout << " -> Final Scale  : " << combinedScale << "\n" << std::endl;
            }
        }

        RooAbsReal* activeDataNll = baseNll.get();
        
        if (combinedScale != 1.0) {
            TString NllName = baseNll->GetName();
            NllName.ReplaceAll("-", "m");
            NllName.ReplaceAll("+", "p");
            baseNll->SetName(NllName);
        
            alphaVar.reset(new RooRealVar("alpha_weight", "alpha_weight", combinedScale));
            alphaVar->setConstant(kTRUE);

            correctedNll.reset(new RooProduct("scaled_data_nll", Form("%lf * %s", combinedScale, baseNll->GetName()), RooArgList(*alphaVar, *baseNll)));
            activeDataNll = correctedNll.get();
        }

        if (fPriorPdf) {
            RooArgSet emptySet;
            constraintNll.reset(fPriorPdf->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(emptySet)));
            totalNll.reset(new RooAddition("total_nll", "Total NLL with Unscaled Constraints", RooArgList(*activeDataNll, *constraintNll)));
            return totalNll.get();
        }

        return activeDataNll;
    }
    
    void BruMcmc::BuildBatchedNLLs(int numBatches) {
      std::cout << "BruMcmc: Deep searching model for BruComponentsPDFs to pre-calculate phase space..." << std::endl;
      
      std::unique_ptr<RooArgSet> comps(fPdf->getComponents());
      for (auto* arg : *comps) {
          auto* bruPdf = dynamic_cast<bru::BruComponentsPDF*>(arg);
          if (bruPdf) {
              std::cout << " -> Forcing initial integration for component: " << bruPdf->GetName() << std::endl;
              bruPdf->analyticalIntegral(1, ""); 
          }
      }
      
      ClearBatches(); 
      int totalEvents = fData->numEntries();
      int batchSize = totalEvents / numBatches;

      std::cout << "BruMcmc: Slicing data into " << numBatches << " batches..." << std::endl;

      for (int b = 0; b < numBatches; ++b) {
          int startIndex = b * batchSize;
          int endIndex = (b == numBatches - 1) ? totalEvents : (b + 1) * batchSize;

          auto subData = std::unique_ptr<RooAbsData>(fData->reduce(RooFit::EventRange(startIndex, endIndex)));
          auto cache = std::unique_ptr<NLLCache>(new NLLCache());
          
          cache->finalNll = BuildNLL(subData.get(), cache->baseNll, cache->alphaVar, cache->correctedNll, cache->constraintNll, cache->totalNll);
          cache->finalNll->constOptimizeTestStatistic(RooAbsArg::Activate, false);
          
          fBatchedNLLPointers.push_back(cache->finalNll);
          fBatchedData.push_back(std::move(subData));
          fBatchedNLLCache.push_back(std::move(cache));
      }
    } 
    
    void BruMcmc::ClearBatches() {
      fBatchedNLLPointers.clear();
      fBatchedNLLCache.clear();
      fBatchedData.clear();
      fBatchSwapFreq = 0;
      std::cout << "BruMcmc: Batched memory cleared." << std::endl;
    }
    
    bool BruMcmc::RunHastings(RooAbsReal* nll) {
        fParams.reset(nll->getParameters(*fData)); 
        RemoveConstantParameters(fParams.get());

        BruMetropolisHastings mh;
        mh.SetFunction(*nll);
        mh.SetParameters(*(fParams.get()));
        if (fChainParams.getSize() > 0) mh.SetChainParameters(fChainParams);
        mh.SetProposalFunction(*fPropFunc);
        mh.SetNumIters(fNumIters);

        if (fBatchSwapFreq > 0 && !fBatchedNLLPointers.empty()) {
            mh.SetBatchedFunctions(fBatchedNLLPointers, fBatchSwapFreq);
        }

        fChain.reset(mh.ConstructChain()); 
        fChainAcceptance = mh.GetAcceptance();
        fLastEvalErrorCount = mh.GetEvalErrorCount();
	
        if(fChain == nullptr){
            if(fTreeMCMC){ delete fTreeMCMC; fTreeMCMC=nullptr; }
            return false;
        }
        return true;
    }

    void BruMcmc::SaveChainToTree() {
        if(fChainData.get()){ fChainData.reset(); }
        if(fTreeMCMC){ delete fTreeMCMC; fTreeMCMC=nullptr; }
     
        const RooDataSet* internalData = fChain->GetAsConstDataSet();
        if(!internalData) return;

        auto saveDir = gDirectory;
        if(fOutFile) fOutFile->cd();
        
        fTreeMCMC = RooStats::GetAsTTree("MCMCTree","MCMCTree", *internalData);
        
        if(fChain->Size() > fNumBurnInSteps){
            fChainData.reset(dynamic_cast<RooDataSet*>(internalData->reduce(RooFit::EventRange(fNumBurnInSteps, fChain->Size()), RooFit::Name("mcmcChain"))) );
        } else {
            fChainData.reset(dynamic_cast<RooDataSet*>(internalData->Clone("mcmcChain")));
        }
        saveDir->cd();
    }

    void BruMcmc::CalculateChainNLLStats() {
        fMeanNLL  = 1e12;
        fSigmaNLL = 0.0;
        fMinNLL   = 1e12;

        if (!fChainData || fChainData->numEntries() == 0) return;

        const Double_t invBeta = (fLastBetaFactor > 1e-9) ? (1.0 / fLastBetaFactor) : 1.0;
        Double_t sumW  = 0.0;
        Double_t sumWN = 0.0;
        Double_t minVal = 1e12;

        for (int i = 0; i < fChainData->numEntries(); ++i) {
            const RooArgSet* row = fChainData->get(i);
            if (!row) continue;
            Double_t w = fChainData->weight();
            if (w <= 0.0) w = 1.0;

            Double_t rawNll = row->getRealValue("nll_MarkovChain_local_", 1e12);
            if (rawNll >= 1e11 && fChain && i + fNumBurnInSteps < fChain->Size()) {
                rawNll = fChain->NLL(i + fNumBurnInSteps);
            }
            if (rawNll >= 1e11) continue;

            Double_t physNll = rawNll * invBeta;
            if (physNll < minVal) minVal = physNll;
            sumW  += w;
            sumWN += w * physNll;
        }

        if (sumW <= 0.0) return;

        fMeanNLL = sumWN / sumW;
        fMinNLL  = minVal;

        Double_t sumVar = 0.0;
        for (int i = 0; i < fChainData->numEntries(); ++i) {
            const RooArgSet* row = fChainData->get(i);
            if (!row) continue;
            Double_t w = fChainData->weight();
            if (w <= 0.0) w = 1.0;

            Double_t rawNll = row->getRealValue("nll_MarkovChain_local_", 1e12);
            if (rawNll >= 1e11 && fChain && i + fNumBurnInSteps < fChain->Size()) {
                rawNll = fChain->NLL(i + fNumBurnInSteps);
            }
            if (rawNll >= 1e11) continue;

            Double_t physNll = rawNll * invBeta;
            Double_t diff = physNll - fMeanNLL;
            sumVar += w * diff * diff;
        }

        fSigmaNLL = std::sqrt(sumVar / sumW);
    }

    // ======================================================================
    // Stationarity Gate: Locate where Phase 2 NLL enters the stationary basin
    // ======================================================================
    size_t BruMcmc::FindStationaryBurnIn(TTree* tree, Int_t nPars) const {
        if (!tree) return 0;
        Long64_t nEntries = tree->GetEntries();
        if (nEntries < 400) return static_cast<size_t>(nEntries / 2);

        TBranch* nllBranch = tree->GetBranch("nll_MarkovChain_local_");
        if (!nllBranch) return static_cast<size_t>(nEntries / 2);

        Double_t nllVal = 0.0;
        tree->ResetBranchAddresses();
        tree->SetBranchAddress("nll_MarkovChain_local_", &nllVal);

        std::vector<Double_t> nllSeries(nEntries, 0.0);
        for (Long64_t i = 0; i < nEntries; ++i) {
            tree->GetEntry(i);
            nllSeries[i] = nllVal;
        }
        tree->ResetBranchAddresses();

        // Characterize the stationary basin using the final 25% of the chain
        Long64_t tailStart = static_cast<Long64_t>(nEntries * 0.75);
        Long64_t tailCount = nEntries - tailStart;
        if (tailCount <= 10) return static_cast<size_t>(nEntries / 2);

        Double_t sumTail = 0.0;
        for (Long64_t i = tailStart; i < nEntries; ++i) sumTail += nllSeries[i];
        Double_t meanTail = sumTail / tailCount;

        Double_t varTail = 0.0;
        for (Long64_t i = tailStart; i < nEntries; ++i) {
            Double_t d = nllSeries[i] - meanTail;
            varTail += d * d;
        }
        Double_t sigmaTail = std::sqrt(varTail / (tailCount - 1));

        // Allow 3.5 sigma thermal fluctuation above the stationary tail mean
        Double_t thermalFloor = 0.5 * std::sqrt(static_cast<Double_t>(std::max(1, nPars)));
        Double_t threshold = meanTail + std::max(3.5 * sigmaTail, thermalFloor);

        // Scan forward in blocks from 25% to 70% of the chain to find where it settled
        Long64_t minBurn = static_cast<Long64_t>(nEntries * 0.25);
        Long64_t maxBurn = static_cast<Long64_t>(nEntries * 0.70);
        Long64_t blockSize = std::max<Long64_t>(200, nEntries / 100);

        Long64_t chosenBurn = maxBurn;
        for (Long64_t idx = minBurn; idx <= maxBurn; idx += blockSize) {
            Long64_t bEnd = std::min(nEntries, idx + blockSize);
            Double_t bMean = 0.0;
            for (Long64_t j = idx; j < bEnd; ++j) bMean += nllSeries[j];
            bMean /= (bEnd - idx);

            if (bMean <= threshold) {
                chosenBurn = idx;
                break;
            }
        }

        std::cout << "BruMcmc: Stationarity gate selected burn-in = " << chosenBurn << " / " << nEntries
                  << " (" << Form("%.1f%%", 100.0 * chosenBurn / nEntries) << " discarded, retaining "
                  << (nEntries - chosenBurn) << " stationary steps | Tail NLL = "
                  << Form("%.2f +/- %.2f", meanTail, sigmaTail) << ")" << std::endl;

        return static_cast<size_t>(chosenBurn);
    }

    Bool_t BruMcmc::MakeChain() {
      fSuccess = kFALSE; 
      
      if (!fData || !fPdf) return kFALSE;
      if (fPOI.getSize() == 0) return kFALSE;

      std::unique_ptr<RooAbsReal> baseNll;
      std::unique_ptr<RooRealVar> alphaVar;
      std::unique_ptr<RooProduct> correctedNll;
      std::unique_ptr<RooAbsReal> constraintNll;
      std::unique_ptr<RooAddition> totalNll;
        
      RooAbsReal* finalNll = nullptr;

      if (fBatchSwapFreq > 0 && !fBatchedNLLPointers.empty()) {
          finalNll = fBatchedNLLPointers[0]; 
      } else {
          finalNll = BuildNLL(fData, baseNll, alphaVar, correctedNll, constraintNll, totalNll);
      }

      if (!finalNll) return kFALSE;

      std::unique_ptr<BruSequentialProposal> defaultPropFunc;
      if (fPropFunc == nullptr) {
          defaultPropFunc.reset(new BruSequentialProposal(fNorm));
          fPropFunc = defaultPropFunc.get();
      }

      bool success = RunHastings(finalNll);

      if (defaultPropFunc) fPropFunc = nullptr; 

      if (!success) {
          CleanMakeChain();
          return kFALSE; 
      }

      SaveChainToTree();
     
      if (baseNll) {
          baseNll->constOptimizeTestStatistic(RooAbsArg::DeActivate, false);
      }

      CleanMakeChain();
       
      fSuccess = kTRUE; 
      return kTRUE;
    }
   
    // ======================================================================
    // Helper 1: Extract Means (Pass 1)
    // ======================================================================
    std::vector<double> BruMcmc::ExtractChainMeans(TTree* tree, size_t burnin, Int_t Nentries, Int_t Npars, 
                                                   std::vector<Double_t>& params, const std::vector<bool>& isCyclic, 
                                                   const std::vector<double>& minVal, const std::vector<double>& maxVal) {
        std::vector<double> sumSin(Npars, 0.0), sumCos(Npars, 0.0);
        std::vector<double> means(Npars, 0.0);

        for (int ientry = burnin; ientry < Nentries + static_cast<Int_t>(burnin); ientry++) {
            tree->GetEntry(ientry);
            for (int p = 0; p < Npars; p++) {
                if (isCyclic[p]) {
                    double len = maxVal[p] - minVal[p];
                    double angle = (params[p] - minVal[p]) / len * 2.0 * TMath::Pi() - TMath::Pi();
                    sumSin[p] += TMath::Sin(angle);
                    sumCos[p] += TMath::Cos(angle);
                } else {
                    means[p] += params[p];
                }
            }
        }

        std::vector<double> circMean(Npars, 0.0);
        for (int p = 0; p < Npars; p++) {
            if (isCyclic[p]) {
                double meanAngle = TMath::ATan2(sumSin[p], sumCos[p]);
                circMean[p] = (meanAngle + TMath::Pi()) / (2.0 * TMath::Pi()) * (maxVal[p] - minVal[p]) + minVal[p];
            } else {
                circMean[p] = means[p] / Nentries; 
            }
        }
        return circMean;
    }

    // ======================================================================
    // Helper 2: Calculate Window-Detrended Empirical Covariance (Pass 2)
    // ======================================================================
    TMatrixDSym BruMcmc::CalculateEmpiricalCovariance(TTree* tree, size_t burnin, Int_t Nentries, Int_t Npars, 
                                                      std::vector<Double_t>& params, const std::vector<double>& means, 
                                                      const std::vector<bool>& isCyclic, const std::vector<double>& minVal, 
                                                      const std::vector<double>& maxVal) {
        TMatrixDSym covMatSym(Npars);
        covMatSym.Zero();

        // Use 4 large sub-windows to detrend residual macroscopic drift while preserving
        // the natural oscillation period of slow polarized moments (H_1, H_2).
        const Int_t numWindows = (Nentries >= 8000) ? 4 : 1;
        const Int_t winSize = Nentries / numWindows;

        for (Int_t w = 0; w < numWindows; ++w) {
            size_t wStart = burnin + static_cast<size_t>(w * winSize);
            Int_t wEntries = (w == numWindows - 1) ? (Nentries - w * winSize) : winSize;
            if (wEntries <= 1) continue;

            std::vector<double> localMeans = (numWindows > 1)
                ? ExtractChainMeans(tree, wStart, wEntries, Npars, params, isCyclic, minVal, maxVal)
                : means;

            for (size_t ientry = wStart; ientry < wStart + static_cast<size_t>(wEntries); ++ientry) {
                tree->GetEntry(ientry);
                std::vector<double> delta(Npars, 0.0);

                for (int p = 0; p < Npars; p++) {
                    if (isCyclic[p]) {
                        double len = maxVal[p] - minVal[p];
                        delta[p] = std::remainder(params[p] - localMeans[p], len); 
                    } else {
                        delta[p] = params[p] - localMeans[p];
                    }
                }

                for(int i = 0; i < Npars; i++) {
                    for(int j = i; j < Npars; j++) {
                        covMatSym(i, j) += delta[i] * delta[j];
                    }
                }
            }
        }

        Double_t denom = std::max(1, Nentries - numWindows);
        for(int i = 0; i < Npars; i++) {
            for(int j = i; j < Npars; j++) {
                covMatSym(i, j) /= denom;
                covMatSym(j, i) = covMatSym(i, j); 
            }
        }
        return covMatSym;
    }

    // ======================================================================
    // Helper 3: Decouple Yields (Off-diagonal yield covariances are strictly 0)
    // ======================================================================
    void BruMcmc::DecoupleYields(TMatrixDSym& covMat, const std::vector<bool>& isYield) {
        int Npars = covMat.GetNrows();
        for(int i = 0; i < Npars; i++) {
            for(int j = i + 1; j < Npars; j++) {
                if (isYield[i] || isYield[j]) {
                    covMat(i, j) = 0.0;
                    covMat(j, i) = 0.0;
                }
            }
        }
    }

    // ======================================================================
    // Helper 4: Minimum Relative Variance Floor & Off-Diagonal Shrinkage
    //           (Preserves exact diagonal variances & leaves yields untouched)
    // ======================================================================
    void BruMcmc::ApplyLinearShrinkage(TMatrixDSym& covMat, const std::vector<bool>& isYield, Double_t alpha, Double_t minRelCov) {
        int Npars = covMat.GetNrows();
        std::vector<Double_t> physVars;
        physVars.reserve(Npars);

        for (int i = 0; i < Npars; i++) {
            if (!isYield[i] && covMat(i, i) > 0.0) {
                physVars.push_back(covMat(i, i));
            }
        }
        if (physVars.empty()) return;

        std::sort(physVars.begin(), physVars.end());
        Double_t medianPhysVar = physVars[physVars.size() / 2];
        Double_t epsilon = std::max(medianPhysVar * minRelCov, 1e-8);

        // 1. Apply floor ONLY to pinned/collapsed physics parameters (< epsilon);
        //    never clamp large physical variances (such as low-L polarized moments H_1_0_0).
        for (int i = 0; i < Npars; i++) {
            if (isYield[i]) continue;
            if (covMat(i, i) < epsilon) {
                covMat(i, i) = epsilon;
            }
        }

        if (alpha <= 0.0) return;

        // 2. Diagonal blending: shrinks noisy off-diagonal elements by (1 - alpha)
        //    while preserving 100% of every diagonal variance.
        TMatrixDSym diagMat(Npars);
        diagMat.Zero();
        for (int i = 0; i < Npars; i++) {
            if (isYield[i]) {
                diagMat(i, i) = covMat(i, i); 
            } else {
                diagMat(i, i) = std::max(covMat(i, i), epsilon); 
            }
        }
        
        if (alpha >= 1.0) {
            covMat = diagMat; 
        } else {
            covMat *= (1.0 - alpha);
            diagMat *= alpha;
            covMat += diagMat;
        }
    }

    // ======================================================================
    // Refactored Master Function
    // ======================================================================
    TMatrixDSym BruMcmc::MakeMcmcCovarianceMatrix(TTree* tree, size_t burnin, Bool_t decoupleYields, Double_t shrinkageAlpha, Double_t minRelCov) {
        auto pars = fSetup->NonConstParsAndYields();
        Int_t Npars = pars.size();
        
        Int_t Nentries = 0;
        if (tree && tree->GetEntries() > static_cast<Long64_t>(burnin)) {
            Nentries = tree->GetEntries() - burnin;
        }
        
        if (Nentries <= 0) {
            std::cerr << "BruMcmc: FATAL ERROR - Insufficient accepted steps (" << Nentries << ") to build matrix!" << std::endl;
            std::cerr << "BruMcmc: Returning Identity Matrix to prevent Cholesky collapse." << std::endl;
            TMatrixDSym identityMat(Npars);
            for (int i=0; i<Npars; i++) identityMat(i,i) = 1.0;
            return identityMat;
        }

        std::vector<Double_t> params(Npars);

        int pindex = 0;
        tree->ResetBranchAddresses();
        for(RooAbsArg* ipar : pars) {
            if(ipar->isConstant()) continue;
            if(tree->SetBranchAddress(ipar->GetName(), &params[pindex]) == 0) pindex++;
        }
        Npars = pindex; 

        std::vector<bool> isCyclic(Npars, false);
        std::vector<bool> isYield(Npars, false); 
        std::vector<double> maxVal(Npars, 0.0), minVal(Npars, 0.0);
        
        for (int i = 0; i < Npars; ++i) {
            auto var = dynamic_cast<RooRealVar*>(pars[i]);
            if (fCyclicPars.contains(*var)) {
                isCyclic[i] = true;
                maxVal[i] = var->getMax();
                minVal[i] = var->getMin();
            }
            if (fSetup->Yields().contains(*var)) {
                isYield[i] = true; 
            }
        }

        // If shrinkageAlpha < 0.0, apply mild dimension-adaptive off-diagonal shrinkage
        // (e.g. 2.7% at d=40, 5.1% at d=77) to suppress sample noise in N*(N-1)/2 elements.
        Double_t effAlpha = shrinkageAlpha;
        if (effAlpha < 0.0) {
            effAlpha = std::clamp(static_cast<Double_t>(Npars) / 1500.0, 0.02, 0.08);
        }

        std::cout << "BruMcmc: Calculating Empirical Covariance for " << Nentries
                  << " accepted steps (Shrinkage alpha = " << effAlpha << ")..." << std::endl;

        std::vector<double> means = ExtractChainMeans(tree, burnin, Nentries, Npars, params, isCyclic, minVal, maxVal);
        TMatrixDSym covMatSym = CalculateEmpiricalCovariance(tree, burnin, Nentries, Npars, params, means, isCyclic, minVal, maxVal);
        
        if (decoupleYields) {
            DecoupleYields(covMatSym, isYield);
        }
        
        if (effAlpha > 0.0 || minRelCov > 0.0) {
            ApplyLinearShrinkage(covMatSym, isYield, effAlpha, minRelCov);
        }

        if (decoupleYields) {
            DecoupleYields(covMatSym, isYield);
        }
        
        covMatSym.Print();
        tree->ResetBranchAddresses();
        return covMatSym;
    }
  
    void BruMcmc::AddEntryBranch(){
      Long64_t entry = 0;
      Double_t weight = 1.0; 
      
      if (fTreeMCMC) {
          auto entryBranch = fTreeMCMC->Branch("entry", &entry, "entry/L");
          auto weightBranch = fTreeMCMC->Branch("weight", &weight, "weight/D");
          
          const RooDataSet* internalData = (fChain && fChain->Size() == fTreeMCMC->GetEntries())
                                           ? fChain->GetAsConstDataSet()
                                           : fChainData.get();
          
          for(entry = 0; entry < fTreeMCMC->GetEntries(); entry++) {
              if (internalData && entry < internalData->numEntries()) {
                  internalData->get(entry); 
                  weight = internalData->weight(); 
              }
              entryBranch->Fill();
              weightBranch->Fill();
          }
      }
    } 
    
    void BruMcmc::Result(){
      AddEntryBranch();
      AddFormulaToMCMCTree();
      CalculateChainNLLStats();

      auto activeList = fSetup->NonConstParsAndYields();
      for(auto* arg : activeList){
        auto* targetPar = dynamic_cast<RooRealVar*>(arg);
        if (!targetPar) continue;
        
        TString pName = targetPar->GetName();
        
        Double_t sumW   = 0.0;
        Double_t sumWX  = 0.0;
        
        for (int entry = 0; entry < fChainData->numEntries(); ++entry) {
            const RooArgSet* row = fChainData->get(entry); 
            Double_t weight = fChainData->weight();
            Double_t val    = row->getRealValue(pName);    
            
            sumW   += weight;
            sumWX  += weight * val;
        }
        
        Double_t mean = sumW > 0 ? (sumWX / sumW) : 0.0;
        Double_t sumVariance = 0.0;

        for (int entry = 0; entry < fChainData->numEntries(); ++entry) {
            const RooArgSet* row = fChainData->get(entry);
            Double_t weight = fChainData->weight();
            Double_t val    = row->getRealValue(pName);
            
            sumVariance += weight * (val - mean) * (val - mean);
        }
        
        Double_t sigma = sumW > 0 ? std::sqrt(sumVariance / sumW) : 0.0;

        std::cout << pName << " " << mean << " +- " << sigma << std::endl;
        
        targetPar->setVal(mean);
        targetPar->setError(sigma);
      }
    } 

    void BruMcmc::AddFormulaToMCMCTree(){
      if (!fTreeMCMC) return;
      fTreeMCMC->ResetBranchAddresses();

      std::cout<<"BruMcmc::AddFormulaToMCMCTree()"<<std::endl;
      auto formulas=fSetup->ParameterFormulas(); 
      if(!formulas.getSize()) return;

      _formVals.assign(formulas.getSize(), 0.0);
      _formBranches.assign(formulas.getSize(), nullptr);

      Int_t iform=0;

      auto parLeaves=fTreeMCMC->GetListOfLeaves();
      
      for(auto* formu_abs:formulas){
        auto* formu=dynamic_cast<RooFormulaVar*>(formu_abs);
        TString formuName=formu->GetName();
        _formVals[iform]=0;
        _formBranches[iform]=nullptr;
        _formBranches[iform]=fTreeMCMC->Branch(formuName,&_formVals[iform],formuName+"/D");
        iform++;
      }

      Long64_t Nmcmc=fTreeMCMC->GetEntries();
      Int_t Nleaf=parLeaves->GetEntries();
      RooArgSet allPars(fSetup->ParsAndYields());
 
      for(Int_t entry=0;entry<Nmcmc;entry++){
        fTreeMCMC->GetEntry(entry);
 
        for(Int_t ibr=0;ibr<Nleaf;ibr++){
          auto *leaf=dynamic_cast<TLeaf*>(parLeaves->At(ibr));	
          auto* brVar=dynamic_cast<RooRealVar*>(allPars.find(leaf->GetName()));
          if(brVar!=nullptr) brVar->setVal(leaf->GetValue());
        }
        iform=0;
        for(auto* formu_abs:formulas){
          auto* formu=dynamic_cast<RooFormulaVar*>(formu_abs);
          _formVals[iform]=formu->getValV();
          _formBranches[iform]->Fill();
          iform++;
        }
      }  
      std::cout<<"BruMcmc::AddFormulaToMCMCTree() done"<<std::endl;
    }
    
    Double_t  BruMcmc::SumWeights(){
      Double_t sumw(0), carry(0);
      Int_t i ;
      for (i=0 ; i<fData->numEntries() ; i++) {
        fData->get(i) ;
 
        Double_t y = fData->weight() - carry;
        Double_t t = sumw + y;
        carry = (t - sumw) - y;
        sumw = t;
      }
      return sumw;
    }
    
    Double_t  BruMcmc::SumWeights2(){
      Double_t sumw(0), carry(0);
      Int_t i ;
      for (i=0 ; i<fData->numEntries() ; i++) {
        fData->get(i) ;
 
        Double_t y = fData->weight()*fData->weight() - carry;
        Double_t t = sumw + y;
        carry = (t - sumw) - y;
        sumw = t;
      }
      return sumw;
    }

    void BruMcmc::SetupBasicUsage()
    {
      fPropFunc = nullptr;
      if (fTreeMCMC) {
          delete fTreeMCMC;
          fTreeMCMC = nullptr;
      }
      TString fileName=fSetup->GetOutDir()+fSetup->GetName()+"/Results"+fSetup->GetTitle()+GetName()+GetTag()+".root";
      fOutFile.reset(TFile::Open(fileName,"recreate"));
    }
    
    file_uptr BruMcmc::SaveInfo(){
      auto saveDir= gDirectory;
      if (fOutFile) fOutFile->cd();
      
      if (fTreeMCMC && fChainData) {
          fTreeMCMC->SetDirectory(fOutFile.get());
          std::cout<<"BruMcmc::SaveInfo() Result"<<std::endl;
          Result(); 
          fTreeMCMC->Write();
          std::cout<<"BruMcmc::SaveInfo() written MCMC"<<std::endl;
          delete fTreeMCMC; fTreeMCMC=nullptr;
      } else {
          std::cout << "BruMcmc::SaveInfo() WARNING: No MCMC Tree to save." << std::endl;
      }

      RooArgSet saveArgs(fSetup->Parameters());
      saveArgs.add(fSetup->Yields());
      
      RooRealVar Nllval("NLL", "NLL", NLL());
      if (fChain || fChainData) {
          saveArgs.add(Nllval);
      }
     
      RooDataSet saveDS(FinalParName(),TString(GetName())+"Results",saveArgs);
      saveDS.add(saveArgs);
      saveDS.Write();
      TTree* treeDS=RooStats::GetAsTTree(ResultTreeName(),ResultTreeName(),saveDS);
      if (treeDS) {
          treeDS->Write();
          delete treeDS; treeDS=nullptr;
      }

      std::cout<<"BruMcmc::SaveInfo() Done to "<< (fOutFile ? fOutFile->GetName() : "null") <<std::endl;
      if (saveDir) saveDir->cd();
      return std::move(fOutFile);
    }

    void BruMcmc::SaveStepInfo(){
      std::cout<<"BruMcmc::SaveStepInfo() "<<fOutFile.get()<<" "<<fTreeMCMC<<" "<< (fOutFile ? fOutFile->GetName() : "null") <<std::endl;
      auto saveDir= gDirectory;
      if (fOutFile) fOutFile->cd();
      
      if (fTreeMCMC) {
          fTreeMCMC->SetDirectory(fOutFile.get());
          AddEntryBranch();
          AddFormulaToMCMCTree();
          fTreeMCMC->Write();
          delete fTreeMCMC; fTreeMCMC=nullptr;
      }
      if (saveDir) saveDir->cd();
    }

    void BruMcmc::SetParVals(RooArgSet* toThesePars){
      for( auto &pory: fSetup->ParsAndYields()){
        if( dynamic_cast<RooRealVar*>(pory)){
          dynamic_cast<RooRealVar*>(pory)->setVal(dynamic_cast<RooRealVar*>(toThesePars->find(pory->GetName()))->getVal());
        }
      }
    }
 
    void BruMcmcSeq::Run(Setup &setup,RooAbsData &fitdata){
        fSetup=&setup;
        fData=&fitdata;
        SetData(fitdata);
        InitModel();
        SetupBasicUsage();
         
        RooStats::SequentialProposal sp(fNorm);
        SetProposalFunction(sp);
        MakeChain();
        CalculateChainNLLStats();
    }

    void BruMcmcSeqHelper::Run(Setup &setup,RooAbsData &fitdata){
         fSetup=&setup;
         fData=&fitdata;
         SetData(fitdata);
         InitModel();
         SetupBasicUsage();

         SetProposalFunction(_proposal);
         MakeChain();
         CalculateChainNLLStats();
    }

    // =======================================================
    // Per-Parameter Sample Count Scaling for BruMcmcCovariance
    // Formula: N_phase_i = Nuser_i * N_base_i * N_parameters
    // Base defaults: N_base = {500, 5000, 250}
    // =======================================================
    void BruMcmcCovariance::ConfigureSampleCounts(Int_t nActivePars) {
        if (_fUserScaleVec.empty()) return;

        Int_t effectivePars = std::max(1, nActivePars);
        _fNumItersVec.resize(_fUserScaleVec.size());

        for (size_t i = 0; i < _fUserScaleVec.size(); ++i) {
            Int_t basePerPar = (i < _baseSamplesPerPar.size()) ? _baseSamplesPerPar[i] : 250;
            Double_t userMult = _fUserScaleVec[i];

            if (_scaleSamplesByPars) {
                // Legacy safety guard: if a script passes raw sample counts (> 500, e.g. {20000, 200000, 10000}),
                // normalize by the 40-parameter (Lmax=4) baseline so it scales proportionally with nActivePars.
                if (userMult > 500.0) {
                    userMult = userMult / (static_cast<Double_t>(basePerPar) * 40.0);
                }
                Long64_t scaled = std::llround(userMult * static_cast<Double_t>(basePerPar) * sqrt(effectivePars));
                _fNumItersVec[i] = static_cast<Int_t>(std::clamp<Long64_t>(scaled, 500LL, 5000000LL));
            } else {
                if (userMult > 500.0) {
                    _fNumItersVec[i] = static_cast<Int_t>(std::llround(userMult));
                } else {
                    Long64_t scaled = std::llround(userMult * static_cast<Double_t>(basePerPar) * 40.0);
                    _fNumItersVec[i] = static_cast<Int_t>(std::clamp<Long64_t>(scaled, 500LL, 5000000LL));
                }
            }
        }

        _iNIter = 0;
        if (_fNumItersVec.size() >= 3) {
            std::cout << "BruMcmcCovariance: Configured sample counts for " << effectivePars
                      << " active parameters -> Phase 1 (1D): " << _fNumItersVec[0]
                      << ", Phase 2 (ND): " << _fNumItersVec[1]
                      << ", Phase 4 (Cov): " << _fNumItersVec[2] << std::endl;
        }
    }

    // =======================================================
    // Modular Phase Implementations for BruMcmcCovariance
    // =======================================================
    void BruMcmcCovariance::RestorePhase2Params() {
        if (!fPhase2EndParams) return;
        for (auto* arg : fSetup->ParsAndYields()) {
            auto* masterVar = dynamic_cast<RooRealVar*>(arg);
            auto* snapVar   = dynamic_cast<RooRealVar*>(fPhase2EndParams->find(arg->GetName()));
            if (masterVar && snapVar && !masterVar->isConstant()) {
                masterVar->setVal(snapVar->getVal());
            }
        }
    }

    Bool_t BruMcmcCovariance::ExecutePhase1_BurnIn(Int_t maxRetries, Int_t numBatches) {
        std::cout << "\n*** Starting Phase 1: Burn-in ***" << std::endl;
        
        // Ensure Phase 1 starts cleanly in 1D sequential mode with a 1D-appropriate acceptance window
        _propSeq.ClearParStepScales();
        _propSeq.SetIsSequential(kTRUE);
        _propSeq.SetScale(fNorm);
        _propSeq.SetAcceptanceRange(fMinAcc, fMaxAcc);
        _propSeq.SetTargetAccept(fTargetAcc);

        _isBatchMode = (numBatches > 1) ? kTRUE : kFALSE;
        if (_isBatchMode) {
            BuildBatchedNLLs(numBatches);
            SetStochasticSwapping(500); 
            
            for (auto* y : static_range_cast<RooRealVar*>(fSetup->Yields())) {
                if (y && !y->isConstant()) {
                    y->setVal(y->getVal() / numBatches);
                }
            }
        }

        ChangeNIter();
        SetTag("1DStep");
        SetupBasicUsage();
        SetProposalFunction(_propSeq);
        
        Bool_t made = MakeChain();
        Int_t retries = 0;
        
        while((made == kFALSE || (fChain && fChain->Size() < fNumIters * 0.8)) && retries < maxRetries) {
            std::cout << "\n*** BruMcmcCovariance: 1DStep Aborted Early! Retrying (" << retries + 1 << "/" << maxRetries << ") with adapted scale ***\n" << std::endl;
            made = MakeChain();
            retries++;
        }

        if (_isBatchMode) {
            ClearBatches(); 
            SetStochasticSwapping(0); 
            for (auto* y : static_range_cast<RooRealVar*>(fSetup->Yields())) {
                if (y && !y->isConstant()) {
                    y->setVal(y->getVal() * numBatches);
                }
            }
        }
        _isBatchMode = kFALSE;

        // Extract 1D parameter step widths from the converged second half of Phase 1
        // so Phase 2 (NDStep) steps each parameter according to its actual sensitivity
        if (fTreeMCMC && fTreeMCMC->GetEntries() > 100) {
            auto allPars = fSetup->NonConstParsAndYields();
            Int_t maxPars = allPars.size();
            size_t p1Burn = fTreeMCMC->GetEntries() / 2;
            Int_t p1Entries = fTreeMCMC->GetEntries() - static_cast<Long64_t>(p1Burn);

            std::vector<Double_t> pVals(maxPars, 0.0);
            std::vector<RooRealVar*> activeVars;
            activeVars.reserve(maxPars);
            RooArgList activeParList;

            Int_t pindex = 0;
            fTreeMCMC->ResetBranchAddresses();
            for (RooAbsArg* ipar : allPars) {
                if (ipar->isConstant()) continue;
                auto* var = dynamic_cast<RooRealVar*>(ipar);
                if (var && fTreeMCMC->SetBranchAddress(var->GetName(), &pVals[pindex]) == 0) {
                    activeVars.push_back(var);
                    activeParList.add(*var);
                    pindex++;
                }
            }
            Int_t nPars = pindex;

            if (nPars > 0 && p1Entries > 10) {
                std::vector<bool> isCyclic(nPars, false);
                std::vector<double> minVal(nPars, 0.0), maxVal(nPars, 0.0);

                for (int i = 0; i < nPars; ++i) {
                    if (fCyclicPars.contains(*activeVars[i])) {
                        isCyclic[i] = true;
                        minVal[i] = activeVars[i]->getMin();
                        maxVal[i] = activeVars[i]->getMax();
                    }
                }

                std::vector<double> p1Means = ExtractChainMeans(fTreeMCMC, p1Burn, p1Entries, nPars, pVals, isCyclic, minVal, maxVal);
                TMatrixDSym p1Cov = CalculateEmpiricalCovariance(fTreeMCMC, p1Burn, p1Entries, nPars, pVals, p1Means, isCyclic, minVal, maxVal);
                fTreeMCMC->ResetBranchAddresses();

                std::vector<Double_t> p1Sigmas(nPars, 0.0);
                for (int i = 0; i < nPars; ++i) {
                    Double_t len = activeVars[i]->getMax() - activeVars[i]->getMin();
                    Double_t sig = (p1Cov(i, i) > 0.0) ? std::sqrt(p1Cov(i, i)) : 0.05 * len;
                    // Compensate for the 4x batching scale on yields during Phase 1
                    if (fSetup->Yields().contains(*activeVars[i]) && numBatches > 1) {
                        sig *= std::sqrt(static_cast<Double_t>(numBatches));
                    }
                    p1Sigmas[i] = sig;
                }
                _propSeq.SetParStepScales(activeParList, p1Sigmas);
                _propSeq.SetScale(2.382 / std::sqrt(static_cast<Double_t>(std::max(1, nPars))));
            } else {
                fTreeMCMC->ResetBranchAddresses();
            }
        }

        return made;
    }

    Bool_t BruMcmcCovariance::ExecutePhase2_Mapping(Int_t maxRetries) {
        _isCovarianceMode = kTRUE;
        ChangeNIter();
        SaveStepInfo();
        SetTag("NDStep");
        SetupBasicUsage();
        SetProposalFunction(_propSeq);
        _propSeq.SetIsSequential(kFALSE);
        
        Bool_t made = MakeChain();
        Int_t retries = 0;
        
        while((made == kFALSE || (fChain && fChain->Size() < fNumIters * 0.8)) && retries < maxRetries) {
            std::cout << "\n*** BruMcmcCovariance: NDStep Aborted Early! Retrying (" << retries + 1 << "/" << maxRetries << ") with adapted scale ***\n" << std::endl;
            made = MakeChain();
            retries++;
        }
        _isCovarianceMode = kFALSE;
        return made;
    }

    // =======================================================
    // HELPER: Build Agnostic Diagnostic Groups
    // =======================================================
    std::vector<RooArgList> BruMcmcCovariance::BuildDiagnosticGroups(const RooArgSet& activePars) {
        std::vector<RooArgList> paramGroups;
        RooArgList yields = fSetup->Yields();
        RooArgList physicsPars;
        
        for (auto* p : activePars) {
            if (!yields.contains(*p)) physicsPars.add(*p);
        }

        if (yields.getSize() > 0) paramGroups.push_back(yields);

        int chunkSize = 5;
        RooArgList currentChunk;
        for (int i = 0; i < physicsPars.getSize(); ++i) {
            currentChunk.add(*physicsPars.at(i));
            if (currentChunk.getSize() == chunkSize || i == physicsPars.getSize() - 1) {
                paramGroups.push_back(currentChunk);
                currentChunk.removeAll();
            }
        }
        return paramGroups;
    }

    // =======================================================
    // PHASE 3A: Bracketed Acceptance Pre-Tuning
    // =======================================================
    Bool_t BruMcmcCovariance::ExecuteTuning_Acceptance(const RooArgSet& activePars, Int_t maxRetries) {
        std::cout << "\n*** Phase 3A: Adaptive Acceptance Tuning ***" << std::endl;
        _isTuningMode = kTRUE;
        
        Int_t originalIters = fNumIters;
        Int_t tuningIters = 500;
        SetNumIters(tuningIters);

        Double_t currentScale = _propCov.StepSizeFactor();
        Int_t currentGibbs = _propCov.GetGibbsBlockSize();

        // Steer into [1.5%, 8.0%] (targeting 3.5%) before Phase 3B starts
        Double_t effMinAcc = (_tuneMode == McmcTuneMode::kMappedRhat || fMaxAcc > 0.5) ? 0.01 : fMinAcc;
        Double_t effMaxAcc = (_tuneMode == McmcTuneMode::kMappedRhat || fMaxAcc > 0.5) ? 0.1 : fMaxAcc;
        Double_t targetBandCenter = (_tuneMode == McmcTuneMode::kMappedRhat || fMaxAcc > 0.5) ? 0.035 : fTargetAcc;
        const Double_t minValidAcc = std::max(0.001, fAbsoluteMinAcc);

        for (int step = 0; step < maxRetries; ++step) {
            RestorePhase2Params();
            _propCov.SetScale(currentScale);
            _propCov.SetGibbsBlockSize(currentGibbs);
            SetProposalFunction(_propCov);

            Bool_t made = MakeChain();
            Int_t actualSteps = fChain ? fChain->Size() : 0;
            Double_t acc = fChainAcceptance;
            Bool_t bailedOut = (!made || actualSteps < tuningIters * 0.8);

            std::cout << " -> Step " << step + 1 << "/" << maxRetries 
                      << " (Scale: " << currentScale << ", Gibbs: " << currentGibbs 
                      << ") | Acceptance: " << acc * 100.0 << "% (Target Band: " 
                      << effMinAcc * 100.0 << "% - " << effMaxAcc * 100.0 << "%)" << std::endl;

            if (!bailedOut && acc >= effMinAcc && acc <= effMaxAcc) {
                std::cout << "--> [SUCCESS] Optimal acceptance band achieved." << std::endl;
                RestorePhase2Params();
                SetNumIters(originalIters);
                _isTuningMode = kFALSE;
                return kTRUE;
            }

            if (bailedOut || acc < minValidAcc) {
                currentScale *= 0.5;
                if (currentScale < 0.02) {
                    if (currentGibbs == 0) currentGibbs = 15;
                    else if (currentGibbs > 5) currentGibbs = std::max(5, currentGibbs / 2);
                    _propCov.SetGibbsBlockSize(currentGibbs);
                    currentScale = 0.25; 
                    std::cout << "    [FALLBACK] Global steps blocked by boundaries. Engaging Gibbs Block Size: " << currentGibbs << std::endl;
                }
            } else {
                Double_t ratio = acc / targetBandCenter;
                Double_t multiplier = 1.0;
                if (ratio > 1.0) {
                    multiplier = std::min(2.2, std::pow(ratio, 0.7));
                } else {
                    multiplier = std::max(0.45, std::pow(ratio, 0.6));
                }
                currentScale *= multiplier;
            }
        }

        std::cout << "--> [WARNING] Acceptance tuning exhausted retries. Proceeding with best-effort config." << std::endl;
        RestorePhase2Params();
        SetNumIters(originalIters);
        _isTuningMode = kFALSE;
        return kFALSE; 
    }

    // =======================================================
    // PHASE 3B: Grouped Diagnostic Verification (R-hat + ESS)
    // =======================================================
    Bool_t BruMcmcCovariance::ExecuteTuning_MappedRhat(const RooArgSet& activePars, Int_t maxRetries) {
        std::cout << "\n*** Phase 3B: Grouped Diagnostic Verification (R-hat + ESS) ***" << std::endl;
        _isTuningMode = kTRUE;

        Int_t originalIters = fNumIters;
        Int_t nActivePars = activePars.getSize();

        // Scale trial window with dimension: 1000 for d=40, 1925 for d=77, 2500 for d=126
        Int_t tuneWindow = 5000;//std::clamp(25 * nActivePars, _rhatWindowSize, 2500);
        // Realistic RWM ESS target in a short trial window (~15.0) so Phase 3B does not over-inflate
        Double_t targetESS = std::max(15.0, static_cast<Double_t>(tuneWindow) / std::max(15.0, 2.5 * nActivePars));
        // Allow acceptance down to fAbsoluteMinAcc (0.1%) in high dimensions
        const Double_t minValidAcc = std::max(0.001, fAbsoluteMinAcc);

        SetNumIters(tuneWindow);

        std::vector<RooArgList> paramGroups = BuildDiagnosticGroups(activePars);
        BruMappedRhat diagHelper;

        Int_t initialGibbs = _propCov.GetGibbsBlockSize();
        Double_t initialScale = _propCov.StepSizeFactor();

        Double_t bestScale = initialScale;
        Int_t bestGibbs = initialGibbs;
        Double_t bestScore = 1e9;
        Double_t bestGlobalMeanRhat = 9999.0;

        auto runTuningStage = [&](Int_t stageGibbs, Double_t startScale, Int_t stageRetries) -> Bool_t {
            Int_t currentGibbs = stageGibbs;
            Double_t overrideScale = startScale;
            Double_t lastGoodScale = 0.0;
            Double_t minBailoutScale = 1e9;

            for (int retry = 0; retry < stageRetries; ++retry) {
                RestorePhase2Params();
                _propCov.SetGibbsBlockSize(currentGibbs);
                _propCov.SetScale(overrideScale);
                SetProposalFunction(_propCov);
                SetNumIters(tuneWindow);

                Bool_t made = MakeChain();
                Int_t actualSteps = fChain ? fChain->Size() : 0;
                Double_t acc = fChainAcceptance;

                double errFrac = static_cast<double>(fLastEvalErrorCount) / (actualSteps > 0 ? actualSteps : 1);
                std::cout << "    [DIAGNOSTIC] there were " << errFrac << " fraction of likelihood error evaluations" << std::endl;

                Bool_t bailedOut = (!made || actualSteps < tuneWindow * 0.85 || fTreeMCMC == nullptr);
                if (bailedOut) {
                    std::cout << "    [DIAGNOSTIC] Chain aborted early (" << actualSteps << "/" << tuneWindow << " steps)." << std::endl;
                }

                Double_t maxRhat = 0.0;
                Double_t sumRhat = 0.0;
                Double_t minESS = 1e9;
                Double_t sumESS = 0.0;
                Int_t nPhysGroups = 0;

                if (!bailedOut && fTreeMCMC != nullptr) {
                    for (size_t i = 0; i < paramGroups.size(); ++i) {
                        std::pair<Double_t, Double_t> res = diagHelper.CalculateDiagnostics(fTreeMCMC, paramGroups[i]);
                        if (i > 0 || paramGroups.size() == 1) {
                            if (res.first > maxRhat) maxRhat = res.first;
                            if (res.second < minESS) minESS = res.second;
                            sumRhat += res.first;
                            sumESS += res.second;
                            nPhysGroups++;
                        }
                    }
                } else {
                    maxRhat = BruMappedRhat::kConvergenceFailure;
                    minESS = 0.0;
                }

                Double_t meanRhat = (nPhysGroups > 0) ? (sumRhat / nPhysGroups) : BruMappedRhat::kConvergenceFailure;
                Double_t meanESS  = (nPhysGroups > 0) ? (sumESS / nPhysGroups) : 0.0;

                std::cout << " -> Step " << retry + 1 << "/" << stageRetries
                          << " (Scale: " << overrideScale << ", Gibbs: " << currentGibbs << ", Win: " << tuneWindow << ")\n"
                          << "    | Acc: " << acc * 100.0 << "%"
                          << " | Mean R-hat: " << (meanRhat > 100 ? "FAIL" : Form("%.4f", meanRhat))
                          << " | Worst R-hat: " << (maxRhat > 100 ? "FAIL" : Form("%.4f", maxRhat))
                          << " | Mean ESS: " << Form("%.1f", meanESS)
                          << " | Worst ESS: " << Form("%.1f", minESS)
                          << " (Target: >" << Form("%.1f", targetESS) << ")" << std::endl;

                if (!bailedOut && acc >= minValidAcc && meanRhat < 50.0) {
                    Double_t score = meanRhat + 0.35 * (maxRhat - meanRhat) - 0.002 * meanESS;
                    if (score < bestScore) {
                        bestScore = score;
                        bestScale = overrideScale;
                        bestGibbs = currentGibbs;
                    }
                    if (currentGibbs == 0 && meanRhat < bestGlobalMeanRhat) {
                        bestGlobalMeanRhat = meanRhat;
                    }
                }

                if (!bailedOut && acc >= minValidAcc &&
                    (((maxRhat <= _rhatTarget) && (meanESS >= (targetESS * 0.75))) )) {
                    std::cout << "--> [SUCCESS] All parameter blocks achieved stationarity and mixing!" << std::endl;
                    bestScale = overrideScale;
                    bestGibbs = currentGibbs;
                    return kTRUE;
                }

                // Bracketed geometric search: only bisect downward on bailout or acc < minValidAcc (0.1%)
                if (bailedOut || acc < minValidAcc) {
                    minBailoutScale = std::min(minBailoutScale, overrideScale);
                    if (lastGoodScale > 0.0) {
                        overrideScale = std::sqrt(lastGoodScale * minBailoutScale);
                        std::cout << "    [DIAGNOSTIC] Bailout/low-acc at ceiling. Bisecting scale to " << overrideScale
                                  << " in bracket [" << lastGoodScale << ", " << minBailoutScale << "]." << std::endl;
                    } else {
                        overrideScale *= 0.6;
                        std::cout << "    [DIAGNOSTIC] Bailout at initial scale. Deflating scale to " << overrideScale << std::endl;
                    }
                } else {
                    lastGoodScale = std::max(lastGoodScale, overrideScale);
                    if (minBailoutScale < 1e8) {
                        Double_t nextScale = std::sqrt(lastGoodScale * minBailoutScale);
                        if (std::abs(nextScale - overrideScale) / overrideScale < 0.05) {
                            std::cout << "    [DIAGNOSTIC] Bracket converged around scale = " << bestScale << "." << std::endl;
                            break;
                        }
                        overrideScale = nextScale;
                        std::cout << "    [DIAGNOSTIC] Bisecting upward toward bailout ceiling: scale = " << overrideScale << std::endl;
                    } else {
                        overrideScale *= 1.35;
                        std::cout << "    [DIAGNOSTIC] Inflating step scale to " << overrideScale << " for broader mixing." << std::endl;
                    }
                }
            }
            return kFALSE;
        };

        // Stage 1: Global Covariance Jumps (or user-specified initialGibbs)
        if (runTuningStage(initialGibbs, initialScale, maxRetries)) {
            _propCov.SetGibbsBlockSize(bestGibbs);
            _propCov.SetScale(bestScale);
            RestorePhase2Params();
            SetNumIters(originalIters);
            _isTuningMode = kFALSE;
            return kTRUE;
        }

        // Stage 2: Conditional Gibbs Fallback ONLY if Global jumps completely failed (Mean R-hat > 1.50)
        if (initialGibbs == 0 && bestGlobalMeanRhat > 1.50) {
            Int_t fallbackGibbs = std::clamp(nActivePars / 4, 5, 25);
            std::cout << "\n===========================================================\n"
                      << " [FALLBACK] Global jumps failed to mix (Best Mean R-hat = " << bestGlobalMeanRhat << ").\n"
                      << "            Engaging Gibbs Block Size = " << fallbackGibbs << ".\n"
                      << "===========================================================\n" << std::endl;
            if (runTuningStage(fallbackGibbs, 0.35, std::max(3, maxRetries / 2))) {
                _propCov.SetGibbsBlockSize(bestGibbs);
                _propCov.SetScale(bestScale);
                RestorePhase2Params();
                SetNumIters(originalIters);
                _isTuningMode = kFALSE;
                return kTRUE;
            }
        }

        std::cout << "\n===========================================================\n"
                  << " [INFO] Diagnostic tuning complete. Locking best observed\n"
                  << "        config (Scale = " << bestScale << ", Gibbs = " << bestGibbs << ").\n"
                  << "===========================================================\n" << std::endl;

        _propCov.SetGibbsBlockSize(bestGibbs);
        _propCov.SetScale(bestScale);
        RestorePhase2Params();
        SetNumIters(originalIters);
        _isTuningMode = kFALSE;
        return (bestGlobalMeanRhat <= 1.35);
    }

    // =======================================================
    // PHASE 4: Official Run
    // =======================================================
    Bool_t BruMcmcCovariance::ExecutePhase4_Official() {
        RestorePhase2Params();
        if (_fNumItersVec.size() >= 3) {
            SetNumIters(_fNumItersVec[2]);
        }
        SetProposalFunction(_propCov);
        std::cout << "\n*** Launching Official Covariance Chain (Scale: "
                  << _propCov.StepSizeFactor() << ", Gibbs: " << _propCov.GetGibbsBlockSize() << ") ***\n" << std::endl;
        _isResultMode = kTRUE;
        _propCov.TuneCovarianceStep(kFALSE);
        Bool_t made = MakeChain();
        if(!made) std::cerr << "BruMcmcCovariance: Official covariance chain failed." << std::endl;
        _isResultMode = kFALSE;
        return made;
    }

    void BruMcmcCovariance::WriteDiagnosticsAndCovariance(const TMatrixDSym& finalCovMat,
                                                          const std::vector<Double_t>& rhatVals,
                                                          const std::vector<Double_t>& essVals) {
        if (!fOutFile) return;
        auto saveDir = gDirectory;
        fOutFile->cd();

        finalCovMat.Write("PosteriorCovariance");

        Double_t finalAcc = fChainAcceptance;
        Double_t meanNll  = fMeanNLL;
        Double_t sigmaNll = fSigmaNLL;
        Double_t minNll   = fMinNLL;
        Double_t maxRhat  = fMaxRhat;
        Double_t minEss   = fMinESS;
        Double_t meanEss  = fMeanESS;

        auto allPars = fSetup->NonConstParsAndYields();
        std::vector<RooArgList> paramGroups = BuildDiagnosticGroups(allPars);
        std::vector<Double_t> localRhat = rhatVals;
        std::vector<Double_t> localEss  = essVals;
        localRhat.resize(paramGroups.size(), 0.0);
        localEss.resize(paramGroups.size(), 0.0);

        TTree* diagTree = new TTree("MCDiagnostics", "MCMC Convergence Diagnostics");
        diagTree->Branch("Acceptance", &finalAcc, "Acceptance/D");
        diagTree->Branch("MeanNLL", &meanNll, "MeanNLL/D");
        diagTree->Branch("SigmaNLL", &sigmaNll, "SigmaNLL/D");
        diagTree->Branch("MinNLL", &minNll, "MinNLL/D");
        diagTree->Branch("MaxRhat", &maxRhat, "MaxRhat/D");
        diagTree->Branch("MinESS", &minEss, "MinESS/D");
        diagTree->Branch("MeanESS", &meanEss, "MeanESS/D");

        for (size_t i = 0; i < paramGroups.size(); ++i) {
            TString groupName = (i == 0 && fSetup->Yields().getSize() > 0) ? "Yields" : Form("Physics_Block_%zu", i);
            diagTree->Branch(groupName + "_Rhat", &localRhat[i], groupName + "_Rhat/D");
            diagTree->Branch(groupName + "_ESS", &localEss[i], groupName + "_ESS/D");
        }
        diagTree->Fill();
        diagTree->Write();
        delete diagTree;

        if (saveDir) saveDir->cd();
        std::cout << "--> Diagnostics successfully written to 'MCDiagnostics' tree.\n" << std::endl;
    }

    // ==============================================================================================
    // MAIN RUN: MULTI-TIERED ESCALATION & BAILOUT ARCHITECTURE
    // ==============================================================================================
    void BruMcmcCovariance::Run(Setup &setup, RooAbsData &fitdata) {
        fData = &fitdata;
        fSetup = &setup;
        InitModel();

        auto allPars = fSetup->NonConstParsAndYields();
        ConfigureSampleCounts(allPars.size());

        // Reset per-chain cached state
        fCachedBetaVal = -1.0;
        fLastBetaFactor = 1.0;
        fMeanNLL = 1e12;
        fSigmaNLL = 0.0;
        fMinNLL = 1e12;
        fMaxRhat = BruMappedRhat::kConvergenceFailure;
        fMinESS = 0.0;
        fMeanESS = 0.0;
        fWasRescued = kFALSE;
        fPhase2EndParams.reset();
        
        _propSeq.SetCyclicParameters(fCyclicPars);
        _propCov.SetCyclicParameters(fCyclicPars);
        _propSeq.SetScale(fNorm);
        _propCov.SetScale(fNorm);
        _propCov.SetYields(fSetup->Yields()); 
        _propCov.TuneCovarianceStep(kFALSE);
	
        _propSeq.SetGibbsBlockSize(_NGibbs);
        _propCov.SetGibbsBlockSize(_NGibbs);

        if (_doSeq) ExecutePhase1_BurnIn(10, 4);
        
        std::unique_ptr<RooDataSet> phase2DataBackup;
        
        if (_doND) {
            std::cout << "\n*** Phase 2: Mapping Phase (Gibbs Size: " << _NGibbs << ") ***" << std::endl;
            
            Double_t originalTarget = fTargetAcc;
            Double_t originalMin = fMinAcc;
            Double_t originalMax = fMaxAcc;
            
            SetDesiredAcceptance(_coldTargetAcc - 0.05, _coldTargetAcc + 0.05, _coldTargetAcc);
            
            _propSeq.SetAcceptanceRange(fMinAcc, fMaxAcc);
            _propSeq.SetTargetAccept(fTargetAcc);
            
            ExecutePhase2_Mapping(10);
            
            SetDesiredAcceptance(originalMin, originalMax, originalTarget);
            _propSeq.SetAcceptanceRange(fMinAcc, fMaxAcc);
            _propSeq.SetTargetAccept(fTargetAcc);
        }
        
        if (fTreeMCMC != nullptr && _doCov) {
            
            // Snapshot converged Phase 2 parameters so Phase 3 trial chains never corrupt Phase 4's start
            RooArgSet allParsAndYields(fSetup->ParsAndYields());
            fPhase2EndParams.reset(dynamic_cast<RooArgSet*>(allParsAndYields.snapshot()));

            if (fChainData) {
                phase2DataBackup.reset(dynamic_cast<RooDataSet*>(fChainData->Clone("phase2_backup")));
            }

            ChangeNIter();
            
            // Locate stationary burn-in point from Phase 2 NLL trajectory and apply adaptive shrinkage
            size_t actualBurnIn = FindStationaryBurnIn(fTreeMCMC, allPars.size());
            std::unique_ptr<TMatrixDSym> covMat(new TMatrixDSym(
                MakeMcmcCovarianceMatrix(fTreeMCMC, actualBurnIn, kTRUE, -1.0, 1e-4)));

            _propCov.SetCovariance(*covMat, allPars);
            
            SaveStepInfo();
            fTreeMCMC = nullptr; 
            
            SetTag("");
            SetupBasicUsage();
            
            _propCov.SetGibbsBlockSize(0);
            _propCov.TuneCovarianceStep(kFALSE);
            _propCov.SetScale(0.25);

            SetProposalFunction(_propCov);

            Int_t maxGlobalRetries = 1; 
            Int_t globalRetryCount = 0;
            Bool_t globalConvergence = kFALSE;

            while (!globalConvergence && globalRetryCount < maxGlobalRetries) {
                
                if (globalRetryCount > 0) {
                    std::cout << "\n=======================================================" << std::endl;
                    std::cout << " GLOBAL POSTERIOR REFINEMENT: Iteration " << globalRetryCount + 1 << " / " << maxGlobalRetries << std::endl;
                    std::cout << "=======================================================\n" << std::endl;
                }

                if (fTuneCovStep || _tuneCovStep) {
                    ExecuteTuning_Acceptance(allPars, 8);
                }
                if (_tuneMode == McmcTuneMode::kMappedRhat) {
                    ExecuteTuning_MappedRhat(allPars, _rhatMaxRetries);
                } 
		
                Bool_t phase4Success = ExecutePhase4_Official();
                
                // --- THE ULTIMATE CRITICAL RESCUE ---
                if (!phase4Success || (fTreeMCMC && fTreeMCMC->GetEntries() < fNumBurnInSteps + 50)) {
                    std::cout << "\n===========================================================" << std::endl;
                    std::cout << " [CRITICAL RESCUE] Official Covariance chain catastrophically failed." << std::endl;
                    std::cout << " Abandoning Covariance Matrix. Falling back to Phase 2 Mapping events." << std::endl;
                    std::cout << "===========================================================\n" << std::endl;
                    
                    if (fTreeMCMC) { delete fTreeMCMC; fTreeMCMC = nullptr; }
                    
                    if (phase2DataBackup && phase2DataBackup->numEntries() > 50) {
                        Int_t totalEvents = phase2DataBackup->numEntries();
                        Int_t fallbackBurnIn = totalEvents / 2;
                        std::cout << " -> Recovering " << totalEvents - fallbackBurnIn << " events from mapping phase (50% burn-in assumed)." << std::endl;

                        fChainData.reset(dynamic_cast<RooDataSet*>(phase2DataBackup->reduce(
                            RooFit::EventRange(fallbackBurnIn, totalEvents), RooFit::Name("mcmcChain"))));

                        fTreeMCMC = RooStats::GetAsTTree("MCMCTree", "MCMCTree", *fChainData);
                        fChainAcceptance = _coldTargetAcc;
                        
                        fNumBurnInSteps = 0;
                        fWasRescued = kTRUE;
                        globalConvergence = kTRUE;
                    } else {
                        std::cout << " -> Phase 2 mapping data was also severely truncated. Cannot safely rescue." << std::endl;
                    }
                }
                
                if (fTreeMCMC != nullptr) {
                    CalculateChainNLLStats();

                    std::cout << "\n*** Calculating Sub-Set Diagnostics (R-hat & ESS) ***" << std::endl;
                    BruMappedRhat diagHelper;
                    std::vector<RooArgList> paramGroups = BuildDiagnosticGroups(allPars);
                    
                    Double_t maxRhat = 0.0;
                    Double_t minESS = 1e9;
                    Double_t sumESS = 0.0;
                    Int_t nEssGroups = 0;
                    
                    std::vector<Double_t> rhatVals(paramGroups.size());
                    std::vector<Double_t> essVals(paramGroups.size());

                    for (size_t i = 0; i < paramGroups.size(); ++i) {
                        std::pair<Double_t, Double_t> res = diagHelper.CalculateDiagnostics(fTreeMCMC, paramGroups[i]);
                        rhatVals[i] = res.first;
                        essVals[i] = res.second;
                        
                        if (res.first > maxRhat) maxRhat = res.first;
                        if (res.second < minESS) minESS = res.second;
                        sumESS += res.second;
                        nEssGroups++;
                        
                        TString groupName = (i == 0 && fSetup->Yields().getSize() > 0) ? "Yields" : Form("Physics_Block_%zu", i);
                        std::cout << " -> " << groupName << " [" << paramGroups[i].getSize() << " pars]" 
                                  << " | R-hat: " << Form("%.4f", rhatVals[i]) 
                                  << " | ESS: " << Form("%.1f", essVals[i]) << std::endl;
                    }

                    fMaxRhat = maxRhat;
                    fMinESS  = (nEssGroups > 0) ? minESS : 0.0;
                    fMeanESS = (nEssGroups > 0) ? (sumESS / nEssGroups) : 0.0;

                    if (globalConvergence) {
                        std::cout << "--> [SUCCESS] Posterior rescued and finalized." << std::endl;
                    } else if (maxRhat != BruMappedRhat::kConvergenceFailure && maxRhat <= _rhatTarget) {
                        std::cout << "--> [SUCCESS] Official chain converged perfectly! Escaping refinement loop." << std::endl;
                        globalConvergence = kTRUE;
                    } else {
                        if (globalRetryCount < maxGlobalRetries - 1) {
                            std::cout << "--> [WARNING] Chain failed global convergence (Worst R-hat: " << maxRhat << ")." << std::endl;
                            
                            std::cout << "--> Rebuilding covariance matrix from recent sampling to improve proposal..." << std::endl;
                            TMatrixDSym rescuedMat = MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kTRUE, -1.0, 1e-4);
                            _propCov.SetCovariance(rescuedMat, allPars);
                            
                            if (fTreeMCMC) { delete fTreeMCMC; fTreeMCMC = nullptr; }
                        } else {
                            std::cout << "--> [WARNING] Max global retries exhausted. Saving best-effort posterior." << std::endl;
                        }
                    }

                    if (globalConvergence || globalRetryCount == maxGlobalRetries - 1) {
                        TMatrixDSym finalCovMat = MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kTRUE, 0.0, 0.0);
                        fLastPosteriorCov.ResizeTo(finalCovMat.GetNrows(), finalCovMat.GetNcols());
                        fLastPosteriorCov = finalCovMat;
                        fLastRhatVals = rhatVals;
                        fLastEssVals  = essVals;

                        WriteDiagnosticsAndCovariance(finalCovMat, rhatVals, essVals);
                    }
                } 
                
                if (globalConvergence) break; 
                globalRetryCount++;
            }
        }
    }
   
  } //namespace FIT
} //namespace HS

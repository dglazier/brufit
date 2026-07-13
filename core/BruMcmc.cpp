#include "BruMcmc.h"
#include "BruMappedRhat.h" // <--- Include the diagnostic helper
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

        // --> NEW: Globally activate MC Variance tracking right at initialization!
        // This ensures the exact interference cross-matrix is perfectly populated 
        // during the Phase 1 DoFirstIntegrations() sweep.
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
    }
    
    RooAbsReal* BruMcmc::BuildNLL(RooAbsData* data,
        std::unique_ptr<RooAbsReal>& baseNll,
        std::unique_ptr<RooRealVar>& alphaVar,
        std::unique_ptr<RooProduct>& correctedNll,
        std::unique_ptr<RooAbsReal>& constraintNll,
        std::unique_ptr<RooAddition>& totalNll) 
    {
        // 1. DATA NLL OPTIONS (Strictly isolated from Priors)
        auto foptions = fSetup->FitOptions();
        TObject* opt=nullptr;
        if((opt=foptions.find("Save"))!=nullptr) foptions.Remove(opt);
        if((opt=foptions.find("SumW2Error"))!=nullptr) foptions.Remove(opt);
        
        auto cmd1 = RooFit::ConditionalObservables(fConditionalObs);
        foptions.Add(dynamic_cast<RooCmdArg*>(&cmd1));
        
        // CRITICAL: Force RooFit NOT to internalize constraints in the data NLL
        auto cmd2 = RooFit::Constrain(RooArgSet()); 
        foptions.Add(dynamic_cast<RooCmdArg*>(&cmd2));

        // Build the pure physics NLL
        baseNll.reset(fPdf->createNLL(*data, foptions));
        baseNll->constOptimizeTestStatistic(RooAbsArg::Activate, false);

        Double_t combinedScale = 1.0;
        
        // --- 2. MC VARIANCE SCALING (BETA) DEEP SEARCH ---
        if (!_isBatchMode && fSetup->ApplyMCVariance()) {
            Double_t max_sigma_rel = 0.0;
            TString domName = "None";

            // Force evaluation. Because TrackMCVariance was activated in InitModel, 
            // the cross-matrix was flawlessly cached during the initial burn-in phase.
            baseNll->getVal(); 

            // Extract the variance 
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
            
            // Replaced arbitrary magic numbers with rigorous > 0 logic
            if (max_sigma_rel > 0.0) {
                Double_t N_data = data->sumEntries(); 
                Double_t betaVal = 1.0 / (1.0 + N_data * (max_sigma_rel * max_sigma_rel));
                combinedScale *= betaVal;
                
                std::cout << "\n=========================================" << std::endl;
                std::cout << " [BruMcmc Likelihood Scaling] " << std::endl;
                std::cout << " -> Dominant PDF : " << domName << " (Rel Err: " << max_sigma_rel * 100.0 << " %)" << std::endl;
                std::cout << " -> Beta Factor  : " << betaVal << std::endl;
                std::cout << "=========================================\n" << std::endl;
            }
        }
        
        // --- 3. ORIGINAL WEIGHT SCALING (ALPHA) ---
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

        // --- 4. SCALE ONLY THE DATA NLL ---
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

        // --- 5. RE-ADD THE CONSTRAINTS (UNSCALED) ---
        if (fPriorPdf) {
            RooArgSet emptySet;
            // Build the NLL penalty for the prior. CloneData(kFALSE) stops redundant memory usage.
            constraintNll.reset(fPriorPdf->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(emptySet)));
            
            // Final NLL = (Scaled Data NLL) + (Unscaled Constraint NLL)
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
      
      ClearBatches(); // Ensure we are clean
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
   
    TMatrixDSym BruMcmc::MakeMcmcCovarianceMatrix(TTree* tree, size_t burnin, Bool_t decoupleYields) {
      auto pars = fSetup->NonConstParsAndYields();
      Int_t Npars = pars.size();
      Int_t Nentries = tree->GetEntries() - burnin;
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
      std::vector<double> sumSin(Npars, 0.0), sumCos(Npars, 0.0);
      std::vector<double> means(Npars, 0.0);
        
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
        
      for (int ientry = burnin; ientry < Nentries + burnin; ientry++) {
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
              means[p] /= Nentries; 
          }
      }

      std::cout << "BruMcmc: Calculating Empirical Covariance for " << Nentries << " accepted steps..." << std::endl;

      TMatrixDSym covMatSym(Npars);
      for (int i = 0; i < Npars; i++) {
          for (int j = 0; j < Npars; j++) covMatSym(i, j) = 0.0;
      }

      for (int ientry = burnin; ientry < Nentries + burnin; ientry++) {
          tree->GetEntry(ientry);
          std::vector<double> delta(Npars, 0.0);

          for (int p = 0; p < Npars; p++) {
              if (isCyclic[p]) {
                  double len = maxVal[p] - minVal[p];
                  delta[p] = std::remainder(params[p] - circMean[p], len); 
              } else {
                  delta[p] = params[p] - means[p];
              }
          }

          for(int i = 0; i < Npars; i++) {
              for(int j = i; j < Npars; j++) {
                  covMatSym(i, j) += delta[i] * delta[j];
              }
          }
      }

      for(int i = 0; i < Npars; i++) {
          for(int j = i; j < Npars; j++) {
              covMatSym(i, j) /= (Nentries - 1);
              
              if (decoupleYields && i != j && (isYield[i] || isYield[j])) {
                  covMatSym(i, j) = 0.0;
              }
              
              covMatSym(j, i) = covMatSym(i, j); // Mirror
          }
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
          
          const RooDataSet* internalData = fChain ? fChain->GetAsConstDataSet() : nullptr;
          
          for(entry = 0; entry < fTreeMCMC->GetEntries(); entry++) {
              if (internalData) {
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
  
      for(auto* arg : *fParams){
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
      fTreeMCMC->ResetBranchAddresses();

      std::cout<<"BruMcmc::AddFormulaToMCMCTree()"<<std::endl;
      auto formulas=fSetup->ParameterFormulas(); 
      if(!formulas.getSize()) return;

      _formVals.reserve(formulas.getSize());
      _formBranches.reserve(formulas.getSize());

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
 
      for(Int_t entry=0;entry<Nmcmc;entry++){
        
        fTreeMCMC->GetEntry(entry);
 
        for(Int_t ibr=0;ibr<Nleaf;ibr++){
          auto *leaf=dynamic_cast<TLeaf*>(parLeaves->At(ibr));	
          auto* brVar=dynamic_cast<RooRealVar*>(fParams->find(leaf->GetName()));
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
      
      RooRealVar Nllval("NLL", "NLL", fChain ? NLL() : 0.0);
      if (fChain) {
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
    }

    void BruMcmcSeqHelper::Run(Setup &setup,RooAbsData &fitdata){
         fSetup=&setup;
         fData=&fitdata;
         SetData(fitdata);
         InitModel();
         SetupBasicUsage();

         SetProposalFunction(_proposal);
         MakeChain();
    }

    // =======================================================
    // Modular Phase Implementations for BruMcmcCovariance
    // =======================================================
    Bool_t BruMcmcCovariance::ExecutePhase1_BurnIn(Int_t maxRetries, Int_t numBatches) {
        std::cout << "\n*** Starting Phase 1: Burn-in ***" << std::endl;
        
        _isBatchMode = (numBatches > 1) ? kTRUE : kFALSE;
        if (_isBatchMode) {
            BuildBatchedNLLs(numBatches);
            SetStochasticSwapping(500); 
            
            // Scale down yields for batching
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
        while(made == kFALSE && retries < maxRetries) {
            std::cout << "\n*** BruMcmcCovariance: 1DStep Failed! Retrying (" << retries + 1 << "/" << maxRetries << ") ***\n" << std::endl;
            _propSeq.SetScale(fNorm); 
            made = MakeChain();
            retries++;
        }

        if (_isBatchMode) {
            ClearBatches(); 
            SetStochasticSwapping(0); 
            // Scale yields back up
            for (auto* y : static_range_cast<RooRealVar*>(fSetup->Yields())) {
                if (y && !y->isConstant()) {
                    y->setVal(y->getVal() * numBatches);
                }
            }
        }
        _isBatchMode = kFALSE;
        return made;
    }

    Bool_t BruMcmcCovariance::ExecutePhase2_Mapping(Int_t maxRetries) {
        std::cout << "\n*** Starting Phase 2: Covariance Mapping (Full Data) ***" << std::endl;
        _isCovarianceMode = kTRUE;
        ChangeNIter();
        SaveStepInfo();
        SetTag("NDStep");
        SetupBasicUsage();
        SetProposalFunction(_propSeq);
        _propSeq.SetIsSequential(kFALSE);
        
        Bool_t made = MakeChain();
        Int_t retries = 0;
        while(made == kFALSE && retries < maxRetries) {
            std::cout << "\n*** BruMcmcCovariance: NDStep Failed! Retrying (" << retries + 1 << "/" << maxRetries << ") ***\n" << std::endl;
            _propSeq.SetScale(fNorm); 
            made = MakeChain();
            retries++;
        }
        _isCovarianceMode = kFALSE;
        return made;
    }
// =======================================================
    // PHASE 3: Diagnostic-Driven Tuning (R-hat + ESS)
    // =======================================================
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
    // PHASE 3A: Pure Acceptance Tuning
    // =======================================================
    Bool_t BruMcmcCovariance::ExecuteTuning_Acceptance(const RooArgSet& activePars, Int_t maxRetries) {
        std::cout << "\n*** Phase 3A: Initial Acceptance Tuning ***" << std::endl;
        _isTuningMode = kTRUE;
        Int_t officialIters = fNumIters;
        SetNumIters(500); 
        
        Bool_t tuned = kFALSE;
        Int_t retries = 0;

        while(!tuned && retries < maxRetries) {
            MakeChain(); 
            
            std::cout << " -> Step " << retries + 1 << "/" << maxRetries 
                      << " | Acceptance: " << Form("%.2f%%", fChainAcceptance * 100.0) << std::endl;

            if (fChainAcceptance > fMinAcc && fChainAcceptance < fMaxAcc) {
                std::cout << "--> [SUCCESS] Optimal acceptance scale achieved." << std::endl;
                tuned = kTRUE;
            } else {
                Double_t safeAcc = fChainAcceptance > 0 ? fChainAcceptance : 0.01;
                Double_t scaleFactor = safeAcc / fTargetAcc;
                
                if (scaleFactor < 0.1) scaleFactor = 0.1;
                if (scaleFactor > 4.0) scaleFactor = 4.0;
                
                _propCov.ApplyNewScale(_propCov.StepSizeFactor() * scaleFactor);
                
                if (fTreeMCMC) { delete fTreeMCMC; fTreeMCMC = nullptr; }
                retries++;
            }
        }
        
        SetNumIters(officialIters);
        _isTuningMode = kFALSE;
        if (fTreeMCMC) { delete fTreeMCMC; fTreeMCMC = nullptr; }
        return tuned;
    }

    // =======================================================
    // PHASE 3B: Grouped Diagnostic Tuning (R-hat + ESS)
    // =======================================================
    Bool_t BruMcmcCovariance::ExecuteTuning_MappedRhat(const RooArgSet& activePars, Int_t maxRetries) {
        std::cout << "\n*** Phase 3B: Grouped Diagnostic Verification (R-hat + ESS) ***" << std::endl;
        _isTuningMode = kTRUE;
        
        Int_t officialIters = fNumIters;
        Int_t tuneWindow = 1000; 
        Bool_t tuned = kFALSE;
        Int_t retries = 0;
        
        BruMappedRhat diagHelper;
        std::vector<RooArgList> groups = BuildDiagnosticGroups(activePars);

        while (!tuned && retries < maxRetries) {
            SetNumIters(tuneWindow);
            MakeChain(); 
            
            Double_t maxRhat = 0.0;
            Double_t minESS = 1e9;
            
            // Evaluate all groups to find the worst performing sub-space
            for (const auto& group : groups) {
                std::pair<Double_t, Double_t> diag = diagHelper.CalculateDiagnostics(fTreeMCMC, group, tuneWindow);
                if (diag.first == BruMappedRhat::kConvergenceFailure) { maxRhat = BruMappedRhat::kConvergenceFailure; break; }
                if (diag.first > maxRhat) maxRhat = diag.first;
                if (diag.second < minESS) minESS = diag.second;
            }

            Double_t acc = fChainAcceptance;
	    Double_t targetESS = tuneWindow * _essFraction;

            std::cout << " -> Step " << retries + 1 << "/" << maxRetries 
                      << "\n    | Acc: " << Form("%.1f%%", acc * 100.0) 
                      << " | Worst R-hat: " << (maxRhat == BruMappedRhat::kConvergenceFailure ? "FAIL" : Form("%.4f", maxRhat))
                      << " | Worst ESS: " << Form("%.1f", minESS) << " (Target: >" << targetESS << ")" << std::endl;

            if (maxRhat != BruMappedRhat::kConvergenceFailure && maxRhat <= _rhatTarget && minESS >= targetESS) {
                std::cout << "--> [SUCCESS] Topology is symmetric and ESS is healthy across all groups." << std::endl;
                tuned = kTRUE;
                
            } else {
                Double_t scaleModifier = 1.0;
                if (acc < 0.02) {
                    std::cout << "    [DIAGNOSTIC] Boundary trap (Acc < 5%). Shrinking scale." << std::endl;
                    scaleModifier = 0.5;
                } else if (acc > 10.40) {
                    std::cout << "    [DIAGNOSTIC] Crawling (Acc > 40%). Expanding scale." << std::endl;
                    scaleModifier = 1.5;
                } else {
                    std::cout << "    [DIAGNOSTIC] Convergence failed. Valley is wide. Expanding scale to force larger jumps." << std::endl;
                    scaleModifier = 1.5;
		    // tuneWindow *= 1.5; // Allow more time to cross the valley
                }
                
                _propCov.ApplyNewScale(_propCov.StepSizeFactor() * scaleModifier);
                if (fTreeMCMC) { delete fTreeMCMC; fTreeMCMC = nullptr; }
                retries++;
            }
        }

        if (!tuned) std::cout << "--> [WARNING] Diagnostic tuning exhausted. Proceeding with best-effort scale." << std::endl;

        SetNumIters(officialIters);
        _isTuningMode = kFALSE;
        if (fTreeMCMC) { delete fTreeMCMC; fTreeMCMC = nullptr; } 
        return tuned;
    }

    // =======================================================
    // PHASE 4: Official Run
    // =======================================================
    Bool_t BruMcmcCovariance::ExecutePhase4_Official() {
        std::cout << "\n*** Launching Official Covariance Chain ***\n" << std::endl;
        _isResultMode = kTRUE;
        Bool_t made = MakeChain();
        if(!made) std::cerr << "BruMcmcCovariance: Official covariance chain failed." << std::endl;
        _isResultMode = kFALSE;
        return made;
    }

    // =======================================================
    // MAIN RUN (With Global Weighted-Average Refinement)
    // =======================================================
    void BruMcmcCovariance::Run(Setup &setup, RooAbsData &fitdata) {
        fData = &fitdata;
        fSetup = &setup;
        InitModel();
        
        _propSeq.SetCyclicParameters(fCyclicPars);
        _propCov.SetCyclicParameters(fCyclicPars);
        _propSeq.SetScale(fNorm);
        _propCov.SetScale(fNorm);

        if (_doSeq) ExecutePhase1_BurnIn(10, 4);
        if (_doND) ExecutePhase2_Mapping(10);
        
  if (fTreeMCMC != nullptr && _doCov) {
            ChangeNIter();
            
            // 1. Safely extract the sample weight BEFORE the tree is saved and destroyed
            Double_t totalSamplesWeight = fTreeMCMC->GetEntries(); 

            // 2. Build initial matrix
            std::unique_ptr<TMatrixDSym> covMat(new TMatrixDSym(MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kTRUE)));
            for(int i = 0; i < covMat->GetNrows(); i++) {
                for(int j = 0; j < covMat->GetNcols(); j++) {
                    if (i != j) (*covMat)(i, j) *= 0.85; 
                }
            }

            auto allPars = fSetup->NonConstParsAndYields();
            _propCov.SetCovariance(*covMat, allPars);
            
            // 3. Save the Phase 2 tree to disk. 
            // This closes the ROOT file and destroys the tree in memory!
            SaveStepInfo();
            fTreeMCMC = nullptr; // <--- CRITICAL C++ SAFETY: Nullify the dangling pointer
            
            SetTag("");
            SetupBasicUsage();
            SetProposalFunction(_propCov);

            Int_t maxGlobalRetries = 3; 
            Int_t globalRetryCount = 0;
            Bool_t globalConvergence = kFALSE;

              // --- THE GLOBAL REFINEMENT LOOP ---
            while (!globalConvergence && globalRetryCount < maxGlobalRetries) {
                
                if (globalRetryCount > 0) {
                    std::cout << "\n=======================================================" << std::endl;
                    std::cout << " GLOBAL POSTERIOR REFINEMENT: Iteration " << globalRetryCount + 1 << " / " << maxGlobalRetries << std::endl;
                    std::cout << "=======================================================\n" << std::endl;
                }

                if (_tuneCovStep) {
                    ExecuteTuning_Acceptance(allPars, 10);
                    ExecuteTuning_MappedRhat(allPars, 10);
                }
                
                if (ExecutePhase4_Official() && fTreeMCMC != nullptr) {
                    TMatrixDSym finalCovMat = MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kFALSE);
                    
                    std::cout << "\n*** Calculating Sub-Set Diagnostics (R-hat & ESS) ***" << std::endl;
                    BruMappedRhat diagHelper;
                    std::vector<RooArgList> paramGroups = BuildDiagnosticGroups(allPars);
                    
                    Double_t finalAcc = fChainAcceptance;
                    Double_t maxRhat = 0.0;
                    Double_t minESS = 1e9;
                    
                    std::vector<Double_t> rhatVals(paramGroups.size());
                    std::vector<Double_t> essVals(paramGroups.size());

                    for (size_t i = 0; i < paramGroups.size(); ++i) {
                        std::pair<Double_t, Double_t> res = diagHelper.CalculateDiagnostics(fTreeMCMC, paramGroups[i]);
                        rhatVals[i] = res.first;
                        essVals[i] = res.second;
                        
                        if (res.first > maxRhat) maxRhat = res.first;
                        if (res.second < minESS) minESS = res.second;
                        
                        TString groupName = (i == 0 && fSetup->Yields().getSize() > 0) ? "Yields" : Form("Physics_Block_%zu", i);
                        std::cout << " -> " << groupName << " [" << paramGroups[i].getSize() << " pars]" 
                                  << " | R-hat: " << Form("%.4f", rhatVals[i]) 
                                  << " | ESS: " << Form("%.1f", essVals[i]) << std::endl;
                    }

                    // --- GLOBAL AUDIT & HAARIO MATRIX UPDATE ---
                    if (maxRhat != BruMappedRhat::kConvergenceFailure && maxRhat <= _rhatTarget) {
                        std::cout << "--> [SUCCESS] Official chain converged perfectly! Escaping refinement loop." << std::endl;
                        globalConvergence = kTRUE;
                    } else {
                        if (globalRetryCount < maxGlobalRetries - 1) {
                            std::cout << "--> [WARNING] Chain failed global convergence (Worst R-hat: " << maxRhat << ")." << std::endl;
                            std::cout << "--> Applying Haario Weighted-Average to refine Covariance Matrix shape..." << std::endl;
                            
                            // Haario Update: C_new = (W_old/W_tot)*C_old + (W_new/W_tot)*C_empirical
                            Double_t currentSamples = fTreeMCMC->GetEntries();
                            Double_t newTotal = totalSamplesWeight + currentSamples;
                            
                            TMatrixDSym currentBaseMat = _propCov.GetBaseMatrix();
                            TMatrixDSym empiricalMat = MakeMcmcCovarianceMatrix(fTreeMCMC, 50, kTRUE);
                            
                            currentBaseMat *= (totalSamplesWeight / newTotal);
                            empiricalMat *= (currentSamples / newTotal);
                            currentBaseMat += empiricalMat; // Safely blend the topologies
                            
                            _propCov.SetCovariance(currentBaseMat, allPars);
                            totalSamplesWeight = newTotal;
                            
                            if (fTreeMCMC) { delete fTreeMCMC; fTreeMCMC = nullptr; }
                        } else {
                            std::cout << "--> [WARNING] Max global retries exhausted. Saving best-effort posterior." << std::endl;
                        }
                    }

                    // --- SAVE AND EXIT ---
                    if (globalConvergence || globalRetryCount == maxGlobalRetries - 1) {
                        if (fOutFile) {
                            fOutFile->cd();
                            finalCovMat.Write("PosteriorCovariance");
                            
                            TTree* diagTree = new TTree("MCDiagnostics", "MCMC Convergence Diagnostics");
                            diagTree->Branch("Acceptance", &finalAcc, "Acceptance/D");
                            for (size_t i = 0; i < paramGroups.size(); ++i) {
                                TString groupName = (i == 0 && fSetup->Yields().getSize() > 0) ? "Yields" : Form("Physics_Block_%zu", i);
                                diagTree->Branch(groupName + "_Rhat", &rhatVals[i], groupName + "_Rhat/D");
                                diagTree->Branch(groupName + "_ESS", &essVals[i], groupName + "_ESS/D");
                            }
                            diagTree->Fill(); 
                            diagTree->Write(); 
                            delete diagTree; 
                            
                            std::cout << "--> Diagnostics successfully written to 'MCDiagnostics' tree.\n" << std::endl;
                        }
                    }
                } else {
                    std::cerr << "BruMcmcCovariance: Official covariance chain failed critically." << std::endl;
                    break; 
                }
                globalRetryCount++;
            }
        }
    }
   
   
  } //namespace FIT
} //namespace HS

// #include "BruMcmc.h"
// #include "BruComponentsPDF.h"
// #include "BruMetropolisHastings.h"
// #include "BruMappedRhat.h"

// #include <TROOT.h>
// #include <TIterator.h>
// #include <TDirectory.h>
// #include <TLeaf.h>
// #include <TTreeIndex.h>
// #include <RooStats/SequentialProposal.h>
// #include <RooStats/ProposalHelper.h>
// #include <TRobustEstimator.h>
// #include <TPrincipal.h>
// #include <TMatrixD.h>
// #include <TH1D.h>
// #include <TMatrixDSym.h>
// #include <RooGlobalFunc.h>
// #include <RooProduct.h>
// #include <RooRealVar.h>
// #include <RooAddition.h>
// #include <cmath>

// namespace HS{
//   namespace FIT{

//     using namespace RooFit;
//     using namespace RooStats;
    
//     void BruMcmc::InitModel() {
//         std::cout << "BruMcmc::InitModel" << std::endl;
//         fPdf = fSetup->Model();
        
//         fPOI.removeAll();
//         fPOI.add(fSetup->Parameters());
//         fPOI.add(fSetup->Yields());
        
//         fNuisParams.removeAll();
//         fConditionalObs.removeAll();
//         fGlobalObs.removeAll();
//         fPriorPdf = nullptr; 

//         // --> NEW: Globally activate MC Variance tracking right at initialization!
//         // This ensures the exact interference cross-matrix is perfectly populated 
//         // during the Phase 1 DoFirstIntegrations() sweep.
//         if (fSetup->ApplyMCVariance()) {
//             std::cout << "BruMcmc: Globally activating exact MC Variance caching." << std::endl;
//             RooArgSet* comps = fPdf->getComponents();
//             for (auto* obj : *comps) {
//                 if (obj->InheritsFrom("bru::BruEventsPDF")) {
//                     static_cast<bru::BruEventsPDF*>(obj)->TrackMCVariance(kTRUE);
//                 }
//             }
//             delete comps;
//         }
//     }

//     void BruMcmc::Run(Setup &setup,RooAbsData &fitdata){
//       fSetup=&setup;
//       fData=&fitdata;
//       std::cout<<"BruMcmc::Run"<<std::endl;
      
//       InitModel();
//       SetupBasicUsage();

//       MakeChain();
//     }
    
//     RooAbsReal* BruMcmc::BuildNLL(RooAbsData* data,
//         std::unique_ptr<RooAbsReal>& baseNll,
//         std::unique_ptr<RooRealVar>& alphaVar,
//         std::unique_ptr<RooProduct>& correctedNll,
//         std::unique_ptr<RooAbsReal>& constraintNll,
//         std::unique_ptr<RooAddition>& totalNll) 
//     {
//         // 1. DATA NLL OPTIONS (Strictly isolated from Priors)
//         auto foptions = fSetup->FitOptions();
//         TObject* opt=nullptr;
//         if((opt=foptions.find("Save"))!=nullptr) foptions.Remove(opt);
//         if((opt=foptions.find("SumW2Error"))!=nullptr) foptions.Remove(opt);
        
//         auto cmd1 = RooFit::ConditionalObservables(fConditionalObs);
//         foptions.Add(dynamic_cast<RooCmdArg*>(&cmd1));
        
//         // CRITICAL: Force RooFit NOT to internalize constraints in the data NLL
//         auto cmd2 = RooFit::Constrain(RooArgSet()); 
//         foptions.Add(dynamic_cast<RooCmdArg*>(&cmd2));

//         // Build the pure physics NLL
//         baseNll.reset(fPdf->createNLL(*data, foptions));
//         baseNll->constOptimizeTestStatistic(RooAbsArg::Activate, false);

//         Double_t combinedScale = 1.0;
        
//         // --- 2. MC VARIANCE SCALING (BETA) DEEP SEARCH ---
//         if (!_isBatchMode && fSetup->ApplyMCVariance()) {
//             Double_t max_sigma_rel = 0.0;
//             TString domName = "None";

//             // Force evaluation. Because TrackMCVariance was activated in InitModel, 
//             // the cross-matrix was flawlessly cached during the initial burn-in phase.
//             baseNll->getVal(); 

//             // Extract the variance 
//             RooArgSet* allComps = fPdf->getComponents();
//             for (auto* obj : *allComps) {
//                 if (obj->InheritsFrom("bru::BruEventsPDF")) {
//                     auto* evPdf = static_cast<bru::BruEventsPDF*>(obj);
                    
//                     Double_t sig = evPdf->GetRelativeVariance();
//                     if (sig > max_sigma_rel) { 
//                         max_sigma_rel = sig; 
//                         domName = evPdf->GetName(); 
//                     }
//                 }
//             }
//             delete allComps;
            
//             // Replaced arbitrary magic numbers with rigorous > 0 logic
//             if (max_sigma_rel > 0.0) {
//                 Double_t N_data = data->sumEntries(); 
//                 Double_t betaVal = 1.0 / (1.0 + N_data * (max_sigma_rel * max_sigma_rel));
//                 combinedScale *= betaVal;
                
//                 std::cout << "\n=========================================" << std::endl;
//                 std::cout << " [BruMcmc Likelihood Scaling] " << std::endl;
//                 std::cout << " -> Dominant PDF : " << domName << " (Rel Err: " << max_sigma_rel * 100.0 << " %)" << std::endl;
//                 std::cout << " -> Beta Factor  : " << betaVal << std::endl;
//                 std::cout << "=========================================\n" << std::endl;
//             }
//         }
        
//         // --- 3. ORIGINAL WEIGHT SCALING (ALPHA) ---
//         if(data->isNonPoissonWeighted() && fCorrectForWeights){
//             Double_t SumW = SumWeights();
//             Double_t SumW2 = SumWeights2();
//             Double_t alphaVal = SumW / SumW2;
//             combinedScale *= alphaVal;

//             if (!_isBatchMode) {
//                 std::cout << " -> Weights Alpha: " << alphaVal << std::endl;
//                 std::cout << " -> Final Scale  : " << combinedScale << "\n" << std::endl;
//             }
//         }

//         // --- 4. SCALE ONLY THE DATA NLL ---
//         RooAbsReal* activeDataNll = baseNll.get();
        
//         if (combinedScale != 1.0) {
//             TString NllName = baseNll->GetName();
//             NllName.ReplaceAll("-", "m");
//             NllName.ReplaceAll("+", "p");
//             baseNll->SetName(NllName);
        
//             alphaVar.reset(new RooRealVar("alpha_weight", "alpha_weight", combinedScale));
//             alphaVar->setConstant(kTRUE);

//             correctedNll.reset(new RooProduct("scaled_data_nll", Form("%lf * %s", combinedScale, baseNll->GetName()), RooArgList(*alphaVar, *baseNll)));
//             activeDataNll = correctedNll.get();
//         }

//         // --- 5. RE-ADD THE CONSTRAINTS (UNSCALED) ---
//         if (fPriorPdf) {
//             RooArgSet emptySet;
//             // Build the NLL penalty for the prior. CloneData(kFALSE) stops redundant memory usage.
//             constraintNll.reset(fPriorPdf->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(emptySet)));
            
//             // Final NLL = (Scaled Data NLL) + (Unscaled Constraint NLL)
//             totalNll.reset(new RooAddition("total_nll", "Total NLL with Unscaled Constraints", RooArgList(*activeDataNll, *constraintNll)));
            
//             return totalNll.get();
//         }

//         return activeDataNll;
//     }
    
//     void BruMcmc::BuildBatchedNLLs(int numBatches) {
//       // =======================================================
//       // THE FIX: Deep Search the PDF Tree
//       // Force the Master PDF to calculate and cache the integrals ONCE.
//       // We must iterate over all components because the BruComponentsPDF 
//       // is usually hidden inside a RooSimultaneous or RooProdPdf!
//       // =======================================================
//       std::cout << "BruMcmc: Deep searching model for BruComponentsPDFs to pre-calculate phase space..." << std::endl;
      
//       std::unique_ptr<RooArgSet> comps(fPdf->getComponents());
//       for (auto* arg : *comps) {
//           auto* bruPdf = dynamic_cast<bru::BruComponentsPDF*>(arg);
//           if (bruPdf) {
//               std::cout << " -> Forcing initial integration for component: " << bruPdf->GetName() << std::endl;
              
//               // This triggers CheckChange (which returns true because _Last is 0)
//               // and runs DoFirstIntegrations(), permanently caching the math!
//               bruPdf->analyticalIntegral(1, ""); 
//           }
//       }
      
//       ClearBatches(); // Ensure we are clean
//       int totalEvents = fData->numEntries();
//       int batchSize = totalEvents / numBatches;

//       std::cout << "BruMcmc: Slicing data into " << numBatches << " batches..." << std::endl;

//       for (int b = 0; b < numBatches; ++b) {
//           int startIndex = b * batchSize;
//           int endIndex = (b == numBatches - 1) ? totalEvents : (b + 1) * batchSize;

//           // Inside the batch loop...
//           auto subData = std::unique_ptr<RooAbsData>(fData->reduce(RooFit::EventRange(startIndex, endIndex)));
//           auto cache = std::unique_ptr<NLLCache>(new NLLCache());
          
//           // Elegantly pass the slice into your custom BuildNLL
//           cache->finalNll = BuildNLL(subData.get(), cache->baseNll, cache->alphaVar, cache->correctedNll, cache->constraintNll, cache->totalNll);
          
//           // Turn off expensive profiling for burn-in batches
//           cache->finalNll->constOptimizeTestStatistic(RooAbsArg::Activate, false);
          
//           fBatchedNLLPointers.push_back(cache->finalNll);
//           fBatchedData.push_back(std::move(subData));
//           fBatchedNLLCache.push_back(std::move(cache));
//       }
//     } 
    
//     void BruMcmc::ClearBatches() {
//       fBatchedNLLPointers.clear();
//       fBatchedNLLCache.clear();
//       fBatchedData.clear();
//       fBatchSwapFreq = 0;
//       std::cout << "BruMcmc: Batched memory cleared." << std::endl;
//     }
    
//     // ---------------------------------------------------------
//     // MODULAR METHOD 2: Run the stepping loop
//     // ---------------------------------------------------------
//     bool BruMcmc::RunHastings(RooAbsReal* nll) {
//         fParams.reset(nll->getParameters(*fData)); 
//         RemoveConstantParameters(fParams.get());

//         BruMetropolisHastings mh;
//         mh.SetFunction(*nll);
//         mh.SetParameters(*(fParams.get()));
//         if (fChainParams.getSize() > 0) mh.SetChainParameters(fChainParams);
//         mh.SetProposalFunction(*fPropFunc);
//         mh.SetNumIters(fNumIters);

//         // ==========================================================
//         // --- STOCHASTIC GRADIENT BATCHING ---
//         // Pass the pre-compiled NLL slices down to the MH engine
//         // so it can rotate the physics landscape during the burn-in!
//         // ==========================================================
//         if (fBatchSwapFreq > 0 && !fBatchedNLLPointers.empty()) {
//             mh.SetBatchedFunctions(fBatchedNLLPointers, fBatchSwapFreq);
//         }

//         fChain.reset(mh.ConstructChain()); 
//         fChainAcceptance = mh.GetAcceptance();

//         if(fChain == nullptr){
//             if(fTreeMCMC){ delete fTreeMCMC; fTreeMCMC=nullptr; }
//             return false;
//         }
//         return true;
//     }

//     // ---------------------------------------------------------
//     // MODULAR METHOD 3: Format and Extract Data
//     // ---------------------------------------------------------
//     void BruMcmc::SaveChainToTree() {
//         if(fChainData.get()){ fChainData.reset(); }
//         if(fTreeMCMC){ delete fTreeMCMC; fTreeMCMC=nullptr; }
     
//         const RooDataSet* internalData = fChain->GetAsConstDataSet();
//         if(!internalData) return;

//         auto saveDir = gDirectory;
//         if(fOutFile) fOutFile->cd();
        
//         fTreeMCMC = RooStats::GetAsTTree("MCMCTree","MCMCTree", *internalData);
        
//         if(fChain->Size() > fNumBurnInSteps){
// 	        fChainData.reset(dynamic_cast<RooDataSet*>(internalData->reduce(RooFit::EventRange(fNumBurnInSteps, fChain->Size()), RooFit::Name("mcmcChain"))) );
//         } else {
// 	        fChainData.reset(dynamic_cast<RooDataSet*>(internalData->Clone("mcmcChain")));
//         }
//         saveDir->cd();
//     }

//     // ---------------------------------------------------------
//     // MAIN ENTRY POINT 
//     // ---------------------------------------------------------
//     Bool_t BruMcmc::MakeChain() {

//       fSuccess = kFALSE; // Assume failure until it completes successfully
      
//       if (!fData || !fPdf) return kFALSE;
//       if (fPOI.getSize() == 0) return kFALSE;

//       std::unique_ptr<RooAbsReal> baseNll;
//       std::unique_ptr<RooRealVar> alphaVar;
//       std::unique_ptr<RooProduct> correctedNll;
//       std::unique_ptr<RooAbsReal> constraintNll;
//       std::unique_ptr<RooAddition> totalNll;
        
//       RooAbsReal* finalNll = nullptr;

//       // ==========================================================
//       // --- NLL ROUTING ---
//       // If stochastic swapping is active, use the pre-compiled batch cache
//       // Otherwise, build a fresh NLL using the current fData pointer
//       // ==========================================================
//       if (fBatchSwapFreq > 0 && !fBatchedNLLPointers.empty()) {
//           finalNll = fBatchedNLLPointers[0]; 
//       } else {
// 	      finalNll = BuildNLL(fData, baseNll, alphaVar, correctedNll, constraintNll, totalNll);
//       }

//       if (!finalNll) return kFALSE;

//       std::unique_ptr<BruSequentialProposal> defaultPropFunc;
//       if (fPropFunc == nullptr) {
//           defaultPropFunc.reset(new BruSequentialProposal(fNorm));
//           fPropFunc = defaultPropFunc.get();
//       }

//       // Execute the Hastings stepping engine
//       bool success = RunHastings(finalNll);

//       if (defaultPropFunc) fPropFunc = nullptr; 

//       if (!success) {
//           CleanMakeChain();
//           return kFALSE; 
//       }

//       SaveChainToTree();
     
//       // Only deactivate the test statistic optimization if we built it dynamically here
//       if (baseNll) {
//           baseNll->constOptimizeTestStatistic(RooAbsArg::DeActivate, false);
//       }

//       CleanMakeChain();
       
//       fSuccess = kTRUE; // We made it to the end! Mark as successful.
//       return kTRUE;
//     }
   
//     ////////////////////////////////////////////////////////
//     TMatrixDSym BruMcmc::MakeMcmcCovarianceMatrix(TTree* tree, size_t burnin, Bool_t decoupleYields) {
//       auto pars = fSetup->NonConstParsAndYields();
//       Int_t Npars = pars.size();
//       Int_t Nentries = tree->GetEntries() - burnin;
//       std::vector<Double_t> params(Npars);

//       int pindex = 0;
//       tree->ResetBranchAddresses();
//       for(RooAbsArg* ipar : pars) {
//           if(ipar->isConstant()) continue;
//           if(tree->SetBranchAddress(ipar->GetName(), &params[pindex]) == 0) pindex++;
//       }
//       Npars = pindex; 

//       std::vector<bool> isCyclic(Npars, false);
//       std::vector<bool> isYield(Npars, false); // Track yields to decouple them later
//       std::vector<double> maxVal(Npars, 0.0), minVal(Npars, 0.0);
//       std::vector<double> sumSin(Npars, 0.0), sumCos(Npars, 0.0);
//       std::vector<double> means(Npars, 0.0);
        
//       for (int i = 0; i < Npars; ++i) {
//           auto var = dynamic_cast<RooRealVar*>(pars[i]);
          
//           if (fCyclicPars.contains(*var)) {
//               isCyclic[i] = true;
//               maxVal[i] = var->getMax();
//               minVal[i] = var->getMin();
//           }
//           if (fSetup->Yields().contains(*var)) {
//               isYield[i] = true; // Flag as a yield
//           }
//       }
        
//       // --- PASS 1: Calculate Means (Circular & Arithmetic) ---
//       for (int ientry = burnin; ientry < Nentries + burnin; ientry++) {
//           tree->GetEntry(ientry);
//           for (int p = 0; p < Npars; p++) {
//               if (isCyclic[p]) {
//                   double len = maxVal[p] - minVal[p];
//                   double angle = (params[p] - minVal[p]) / len * 2.0 * TMath::Pi() - TMath::Pi();
//                   sumSin[p] += TMath::Sin(angle);
//                   sumCos[p] += TMath::Cos(angle);
//               } else {
//                   means[p] += params[p];
//               }
//           }
//       }
        
//       std::vector<double> circMean(Npars, 0.0);
//       for (int p = 0; p < Npars; p++) {
//           if (isCyclic[p]) {
//               double meanAngle = TMath::ATan2(sumSin[p], sumCos[p]);
//               circMean[p] = (meanAngle + TMath::Pi()) / (2.0 * TMath::Pi()) * (maxVal[p] - minVal[p]) + minVal[p];
//           } else {
//               means[p] /= Nentries; // Finalize arithmetic mean
//           }
//       }

//       std::cout << "BruMcmc: Calculating Empirical Covariance for " << Nentries << " accepted steps..." << std::endl;

//       // --- PASS 2: Calculate Covariance Matrix ---
//       TMatrixDSym covMatSym(Npars);
//       for (int i = 0; i < Npars; i++) {
//           for (int j = 0; j < Npars; j++) covMatSym(i, j) = 0.0;
//       }

//       for (int ientry = burnin; ientry < Nentries + burnin; ientry++) {
//           tree->GetEntry(ientry);
//           std::vector<double> delta(Npars, 0.0);

//           for (int p = 0; p < Npars; p++) {
//               if (isCyclic[p]) {
//                   double len = maxVal[p] - minVal[p];
//                   delta[p] = std::remainder(params[p] - circMean[p], len); 
//               } else {
//                   delta[p] = params[p] - means[p];
//               }
//           }

//           // Build upper triangle
//           for(int i = 0; i < Npars; i++) {
//               for(int j = i; j < Npars; j++) {
//                   covMatSym(i, j) += delta[i] * delta[j];
//               }
//           }
//       }

//       // Finalize Matrix: Divide by (N-1) and mirror to lower triangle
//       for(int i = 0; i < Npars; i++) {
//           for(int j = i; j < Npars; j++) {
//               covMatSym(i, j) /= (Nentries - 1);
              
//               // DECOUPLE YIELD DRAG: Only zero out correlations if decoupleYields is kTRUE
//               if (decoupleYields && i != j && (isYield[i] || isYield[j])) {
//                   covMatSym(i, j) = 0.0;
//               }
              
//               covMatSym(j, i) = covMatSym(i, j); // Mirror
//           }
//       }
     
//       covMatSym.Print();
//       tree->ResetBranchAddresses();
//       return covMatSym;
//     } 
//     ////////////////////////////////////////////////////////
  
//     /////////////////////////////////////////////////////////
//     void BruMcmc::AddEntryBranch(){
      
//       Long64_t entry = 0;
//       Double_t weight = 1.0; // Default to 1 just in case
      
//       if (fTreeMCMC) {
//           auto entryBranch = fTreeMCMC->Branch("entry", &entry, "entry/L");
//           auto weightBranch = fTreeMCMC->Branch("weight", &weight, "weight/D");
          
//           // Grab the raw dataset that has the weights preserved
//           const RooDataSet* internalData = fChain ? fChain->GetAsConstDataSet() : nullptr;
          
//           for(entry = 0; entry < fTreeMCMC->GetEntries(); entry++) {
//               if (internalData) {
//                   internalData->get(entry); // Load the row
//                   weight = internalData->weight(); // Extract the hidden weight
//               }
              
//               entryBranch->Fill();
//               weightBranch->Fill();
//           }
//       }
//     } 
    
//     void BruMcmc::Result(){
//       AddEntryBranch();
//       AddFormulaToMCMCTree();
  
//       // Use a modern C++ range-based loop over the RooArgSet
//       for(auto* arg : *fParams){

//         auto* targetPar = dynamic_cast<RooRealVar*>(arg);
//         if (!targetPar) continue;
        
//         TString pName = targetPar->GetName();
        
//         // ==========================================================
//         // --- STABLE TWO-PASS ALGORITHM FOR UNCERTAINTIES ---
//         // ==========================================================
//         Double_t sumW   = 0.0;
//         Double_t sumWX  = 0.0;
        
//         // PASS 1: Calculate the Mean safely
//         for (int entry = 0; entry < fChainData->numEntries(); ++entry) {
//             const RooArgSet* row = fChainData->get(entry); // Get the actual row
//             Double_t weight = fChainData->weight();
//             Double_t val    = row->getRealValue(pName);    // Safely extract by name
            
//             sumW   += weight;
//             sumWX  += weight * val;
//         }
        
//         Double_t mean = sumW > 0 ? (sumWX / sumW) : 0.0;
//         Double_t sumVariance = 0.0;

//         // PASS 2: Calculate Variance using (x - mu)^2 to prevent cancellation
//         for (int entry = 0; entry < fChainData->numEntries(); ++entry) {
//             const RooArgSet* row = fChainData->get(entry);
//             Double_t weight = fChainData->weight();
//             Double_t val    = row->getRealValue(pName);
            
//             sumVariance += weight * (val - mean) * (val - mean);
//         }
        
//         Double_t sigma = sumW > 0 ? std::sqrt(sumVariance / sumW) : 0.0;
//         // ==========================================================

//         std::cout << pName << " " << mean << " +- " << sigma << std::endl;
        
//         targetPar->setVal(mean);
//         targetPar->setError(sigma);
//       }
//     } 

//     void BruMcmc::AddFormulaToMCMCTree(){
//       fTreeMCMC->ResetBranchAddresses();

//       std::cout<<"BruMcmc::AddFormulaToMCMCTree()"<<std::endl;
//       auto formulas=fSetup->ParameterFormulas(); 
//       if(!formulas.getSize()) return;

//       _formVals.reserve(formulas.getSize());
//       _formBranches.reserve(formulas.getSize());

//       Int_t iform=0;

//       auto parLeaves=fTreeMCMC->GetListOfLeaves();
      
//       for(auto* formu_abs:formulas){
//         auto* formu=dynamic_cast<RooFormulaVar*>(formu_abs);
//         TString formuName=formu->GetName();
//         _formVals[iform]=0;
//         _formBranches[iform]=nullptr;
//         _formBranches[iform]=fTreeMCMC->Branch(formuName,&_formVals[iform],formuName+"/D");
//         iform++;
//       }

//       Long64_t Nmcmc=fTreeMCMC->GetEntries();
//       Int_t Nleaf=parLeaves->GetEntries();
 
//       for(Int_t entry=0;entry<Nmcmc;entry++){
        
//         fTreeMCMC->GetEntry(entry);
 
//         for(Int_t ibr=0;ibr<Nleaf;ibr++){
//           auto *leaf=dynamic_cast<TLeaf*>(parLeaves->At(ibr));	
//           auto* brVar=dynamic_cast<RooRealVar*>(fParams->find(leaf->GetName()));
//           if(brVar!=nullptr) brVar->setVal(leaf->GetValue());
            
//         }
//         iform=0;
//         for(auto* formu_abs:formulas){
//           auto* formu=dynamic_cast<RooFormulaVar*>(formu_abs);
          
//           _formVals[iform]=formu->getValV();
//           _formBranches[iform]->Fill();
//           iform++;

//         }
//       }  
//     std::cout<<"BruMcmc::AddFormulaToMCMCTree() done"<<std::endl;
//      }
//     ///////////////////////////////////////////////
//     Double_t  BruMcmc::SumWeights(){
//       Double_t sumw(0), carry(0);
//       Int_t i ;
//       for (i=0 ; i<fData->numEntries() ; i++) {
//         fData->get(i) ;
 
//         Double_t y = fData->weight() - carry;
//         Double_t t = sumw + y;
//         carry = (t - sumw) - y;
//         sumw = t;
//       }
//       return sumw;
//     }
//     ///////////////////////////////////////////////////////
//     Double_t  BruMcmc::SumWeights2(){
//       Double_t sumw(0), carry(0);
//       Int_t i ;
//       for (i=0 ; i<fData->numEntries() ; i++) {
//         fData->get(i) ;
 
//         Double_t y = fData->weight()*fData->weight() - carry;
//         Double_t t = sumw + y;
//         carry = (t - sumw) - y;
//         sumw = t;
//       }
//       return sumw;
//     }

//     /////////////////////////////////////////////////////
//     void BruMcmc::SetupBasicUsage()
//     {
//       fPropFunc = nullptr;
     
//       TString fileName=fSetup->GetOutDir()+fSetup->GetName()+"/Results"+fSetup->GetTitle()+GetName()+GetTag()+".root";

//       fOutFile.reset(TFile::Open(fileName,"recreate"));
       
//      }
//     ///////////////////////////////////////////////////////////////
//     file_uptr BruMcmc::SaveInfo(){
//       auto saveDir= gDirectory;
//       if (fOutFile) fOutFile->cd();
      
//       if (fTreeMCMC && fChainData) {
//           fTreeMCMC->SetDirectory(fOutFile.get());
//           std::cout<<"BruMcmc::SaveInfo() Result"<<std::endl;
//           Result(); 
//           fTreeMCMC->Write();
//           std::cout<<"BruMcmc::SaveInfo() written MCMC"<<std::endl;
//           delete fTreeMCMC; fTreeMCMC=nullptr;
//       } else {
//           std::cout << "BruMcmc::SaveInfo() WARNING: No MCMC Tree to save." << std::endl;
//       }

//       RooArgSet saveArgs(fSetup->Parameters());
//       saveArgs.add(fSetup->Yields());
      
//       RooRealVar Nllval("NLL", "NLL", fChain ? NLL() : 0.0);
//       if (fChain) {
//           saveArgs.add(Nllval);
//       }
     
//       RooDataSet saveDS(FinalParName(),TString(GetName())+"Results",saveArgs);
//       saveDS.add(saveArgs);
//       saveDS.Write();
//       TTree* treeDS=RooStats::GetAsTTree(ResultTreeName(),ResultTreeName(),saveDS);
//       if (treeDS) {
//           treeDS->Write();
//           delete treeDS; treeDS=nullptr;
//       }

//       std::cout<<"BruMcmc::SaveInfo() Done to "<< (fOutFile ? fOutFile->GetName() : "null") <<std::endl;
//       if (saveDir) saveDir->cd();
//       return std::move(fOutFile);
//     }

//     ///////////////////////////////////////////////////////////////
//     void BruMcmc::SaveStepInfo(){
      
//       std::cout<<"BruMcmc::SaveStepInfo() "<<fOutFile.get()<<" "<<fTreeMCMC<<" "<< (fOutFile ? fOutFile->GetName() : "null") <<std::endl;
//       auto saveDir= gDirectory;
//       if (fOutFile) fOutFile->cd();
      
//       if (fTreeMCMC) {
//           fTreeMCMC->SetDirectory(fOutFile.get());
//           AddEntryBranch();
//           AddFormulaToMCMCTree();
//           fTreeMCMC->Write();
//           delete fTreeMCMC; fTreeMCMC=nullptr;
//       }
//     }

   
//      //////////////////////////////////////////////////////////////

//     void BruMcmc::SetParVals(RooArgSet* toThesePars){
//       for( auto &pory: fSetup->ParsAndYields()){
//         if( dynamic_cast<RooRealVar*>(pory)){
//           dynamic_cast<RooRealVar*>(pory)->setVal(dynamic_cast<RooRealVar*>(toThesePars->find(pory->GetName()))->getVal());
//         }
//       }
//     }
 
//    void BruMcmcSeq::Run(Setup &setup,RooAbsData &fitdata){
//      fSetup=&setup;
//     fData=&fitdata;
//     SetData(fitdata);
//     InitModel();
//     SetupBasicUsage();
     
//     RooStats::SequentialProposal sp(fNorm);
//     SetProposalFunction(sp);
//     MakeChain();
    
//    }

//    void BruMcmcSeqHelper::Run(Setup &setup,RooAbsData &fitdata){

//      fSetup=&setup;
//      fData=&fitdata;
//      SetData(fitdata);
//      InitModel();
//      SetupBasicUsage();

//      SetProposalFunction(_proposal);
//      MakeChain();

//    }

//     // =======================================================
//     // Modular Phase Implementations
//     // =======================================================
//     Bool_t BruMcmcCovariance::ExecutePhase1_BurnIn(Int_t maxRetries, Int_t numBatches) {
//         std::cout << "\n*** Starting Phase 1: Burn-in ***" << std::endl;
        
//         _isBatchMode = (numBatches > 1) ? kTRUE : kFALSE;
//         if (_isBatchMode) {
//             BuildBatchedNLLs(numBatches);
//             SetStochasticSwapping(500); 
            
//             // Scale down yields for batching
//             for (auto* y : static_range_cast<RooRealVar*>(fSetup->Yields())) {
//                 if (y && !y->isConstant()) {
//                     y->setVal(y->getVal() / numBatches);
//                 }
//             }
//         }

//         ChangeNIter();
//         SetTag("1DStep");
//         SetupBasicUsage();
//         SetProposalFunction(_propSeq);
        
//         Bool_t made = MakeChain();
//         Int_t retries = 0;
//         while(made == kFALSE && retries < maxRetries) {
//             std::cout << "\n*** BruMcmcCovariance: 1DStep Failed! Retrying (" << retries + 1 << "/" << maxRetries << ") ***\n" << std::endl;
//             _propSeq.SetScale(fNorm); 
//             made = MakeChain();
//             retries++;
//         }

//         if (_isBatchMode) {
//             ClearBatches(); 
//             SetStochasticSwapping(0); 
//             // Scale yields back up
//             for (auto* y : static_range_cast<RooRealVar*>(fSetup->Yields())) {
//                 if (y && !y->isConstant()) {
//                     y->setVal(y->getVal() * numBatches);
//                 }
//             }
//         }
//         _isBatchMode = kFALSE;
//         return made;
//     }

//     Bool_t BruMcmcCovariance::ExecutePhase2_Mapping(Int_t maxRetries) {
//         std::cout << "\n*** Starting Phase 2: Covariance Mapping (Full Data) ***" << std::endl;
//         _isCovarianceMode = kTRUE;
//         ChangeNIter();
//         SaveStepInfo();
//         SetTag("NDStep");
//         SetupBasicUsage();
//         SetProposalFunction(_propSeq);
//         _propSeq.SetIsSequential(kFALSE);
        
//         Bool_t made = MakeChain();
//         Int_t retries = 0;
//         while(made == kFALSE && retries < maxRetries) {
//             std::cout << "\n*** BruMcmcCovariance: NDStep Failed! Retrying (" << retries + 1 << "/" << maxRetries << ") ***\n" << std::endl;
//             _propSeq.SetScale(fNorm); 
//             made = MakeChain();
//             retries++;
//         }
//         _isCovarianceMode = kFALSE;
//         return made;
//     }
//     Bool_t BruMcmcCovariance::ExecutePhase3_Tuning(const RooArgSet& activePars, Int_t maxRetries) {
//       if (_tuneMode == McmcTuneMode::kMappedRhat) {
// 	return ExecuteTuning_MappedRhat(activePars);
//       } else {
// 	return ExecuteTuning_Acceptance(activePars, maxRetries);
//       }
//     }
//     // =======================================================
//     // PHASE 3A: Legacy Acceptance Tuning
//     // =======================================================
//     Bool_t BruMcmcCovariance::ExecuteTuning_Acceptance(const RooArgSet& activePars, Int_t maxRetries) {
//         std::cout << "\n*** Starting Acceptance-Based Covariance Tuning (250 steps) ***" << std::endl;
//         _isTuningMode = kTRUE;
//         Int_t officialIters = fNumIters;
//         SetNumIters(250); 
        
//         Bool_t tuned = kFALSE;
//         Int_t retries = 0;
//         TMatrixDSym tuneCovMat = _propCov.GetCovariance(); // Create deep copy

//         while(!tuned && retries < maxRetries) {
//             MakeChain(); 
            
//             if (fChainAcceptance > fMinAcc && fChainAcceptance < fMaxAcc) {
//                 std::cout << "--> Tuning Successful! Acceptance: " << fChainAcceptance << std::endl;
//                 tuned = kTRUE;
//             } else {
//                 Double_t currentScale = _propCov.StepSizeFactor();
//                 Double_t acc = fChainAcceptance > 0 ? fChainAcceptance : 0.01;
//                 currentScale *= (acc) / (fTargetAcc);
                
//                 std::cout << "--> Tuning Failed. Adjusting scale and retrying (" 
//                           << retries + 1 << "/" << maxRetries << ") with new scale " << currentScale << std::endl;

//                 _propCov.SetScale(currentScale); 
//                 _propCov.SetCovariance(tuneCovMat, activePars);
//                 retries++;
//             }
//         }
        
//         SetNumIters(officialIters);
//         _isTuningMode = kFALSE;
//         return tuned;
//     }

//     // =======================================================
//     // PHASE 3B: Dual-Metric Mapped R-hat Validation Tuning
//     // =======================================================
//     Bool_t BruMcmcCovariance::ExecuteTuning_MappedRhat(const RooArgSet& activePars) {
//         std::cout << "\n*** Starting Mapped R-hat Validation & Tuning Phase ***" << std::endl;
//         _isTuningMode = kTRUE;
        
//         Int_t officialIters = fNumIters;
//         Int_t currentWindowSize = _rhatWindowSize; 
//         Bool_t tuned = kFALSE;
//         Int_t retries = 0;
        
//         // Cache the baseline covariance matrix to apply scaling cleanly on retries
//         TMatrixDSym baseCovMat = _propCov.GetCovariance(); // Deep copy
//         BruMappedRhat diagnosticHelper; 

//         while (!tuned && retries < _rhatMaxRetries) {
//             std::cout << " -> Running validation window: " << currentWindowSize << " steps..." << std::endl;
            
//             SetNumIters(currentWindowSize);
            
//             // MakeChain builds the RAM tree internally. 
//             // We do NOT call SaveStepInfo() which would write it to disk and destroy the pointer.
//             MakeChain(); 
            
//             // Pass the RAM tree directly to the SIMD diagnostic engine
//             Double_t rHatValue = diagnosticHelper.CalculateSplitR(fTreeMCMC, activePars);
            
//             std::cout << " -> Result: Mapped R-hat = " 
//                       << (rHatValue == BruMappedRhat::kConvergenceFailure ? "DEAD CHAIN (W <= 0)" : Form("%.4f", rHatValue)) 
//                       << " | Acceptance = " << fChainAcceptance << std::endl;

//             // Use the named constant for explicit logic
//             if (rHatValue != BruMappedRhat::kConvergenceFailure && rHatValue <= _rhatTarget) {
//                 std::cout << "--> Mapped R-hat Tuning SUCCESSFUL! Locking covariance scale." << std::endl;
//                 tuned = kTRUE;
//             } else {
//                 Double_t scaleFactor = 1.0;
                
//                 // Dual-Metric Logic: Combine R-hat failure with Acceptance Compass
//                 if (fChainAcceptance < fMinAcc) {
//                     scaleFactor = 0.5; // Sampler is crashing into boundaries. Shrink aggressively.
//                     std::cout << "    [Action: Overstepping detected. Shrinking proposal scale.]" << std::endl;
//                 } else if (fChainAcceptance > fMaxAcc) {
//                     scaleFactor = 1.5; // Sampler is barely moving. Grow aggressively.
//                     std::cout << "    [Action: Diffusive crawl detected. Expanding proposal scale.]" << std::endl;
//                 } else {
//                     scaleFactor = 0.85; // Good acceptance but bad mixing -> Topological trap.
//                     currentWindowSize *= 2; // Expand the window to gather more robust statistics.
//                     std::cout << "    [Action: Local trap detected. Perturbing scale and expanding validation window.]" << std::endl;
//                 }

//                 // Apply scaling and reset the matrix
//                 Double_t newTotalScale = _propCov.StepSizeFactor() * scaleFactor;
// 		_propCov.ApplyNewScale(newTotalScale);
                 
//                 // Clean up the temporary tree to prevent memory leaks during retries
//                 if (fTreeMCMC) { 
//                     delete fTreeMCMC; 
//                     fTreeMCMC = nullptr; 
//                 }
                
//                 retries++;
//             }
//         }

//         if (!tuned) {
//              std::cout << "--> WARNING: Mapped R-hat Tuning exhausted " << _rhatMaxRetries 
//                        << " retries without converging. Launching official chain anyway..." << std::endl;
//         }

//         SetNumIters(officialIters);
//         _isTuningMode = kFALSE;
//         return tuned;
//     }
//     // Bool_t BruMcmcCovariance::ExecutePhase3_Tuning(const RooArgSet& activePars, Int_t maxRetries) {
//     //     std::cout << "\n*** Starting Covariance Tuning Phase (250 steps) ***" << std::endl;
//     //     _isTuningMode = kTRUE;
//     //     Int_t officialIters = fNumIters;
//     //     SetNumIters(250); 
        
//     //     Bool_t tuned = kFALSE;
//     //     Int_t retries = 0;
//     //     auto tuneCovMat = _propCov.GetCovariance();

//     //     while(!tuned && retries < maxRetries) {
//     //         MakeChain(); 
            
//     //         if (fChainAcceptance > fMinAcc && fChainAcceptance < fMaxAcc) {
//     //             std::cout << "--> Tuning Successful! Acceptance: " << fChainAcceptance << std::endl;
//     //             tuned = kTRUE;
//     //         } else {
//     //             Double_t currentScale = _propCov.StepSizeFactor();
//     //             Double_t acc = fChainAcceptance > 0 ? fChainAcceptance : 0.01;
//     //             currentScale *= (acc) / (fTargetAcc);
                
//     //             std::cout << "--> Tuning Failed. Adjusting scale and retrying (" 
//     //                       << retries + 1 << "/" << maxRetries << ") with new scale " << currentScale << std::endl;

//     //             _propCov.SetScale(currentScale); 
//     //             _propCov.SetCovariance(tuneCovMat, activePars);
//     //             retries++;
//     //         }
//     //     }
        
//     //     SetNumIters(officialIters);
//     //     _isTuningMode = kFALSE;
//     //     return tuned;
//     // }

//     Bool_t BruMcmcCovariance::ExecutePhase4_Official() {
//         std::cout << "\n*** Launching Official Covariance Chain ***\n" << std::endl;
//         _isResultMode = kTRUE;
//         Bool_t made = MakeChain();
        
//         if(!made) {
//             std::cerr << "BruMcmcCovariance: Official covariance chain failed." << std::endl;
//         }
//         _isResultMode = kFALSE;
//         return made;
//     }

//     // Now, the main Run method is incredibly clean:
//     void BruMcmcCovariance::Run(Setup &setup, RooAbsData &fitdata) {
//         fData = &fitdata;
//         fSetup = &setup;
//         InitModel();
        
//         _propSeq.SetCyclicParameters(fCyclicPars);
//         _propCov.SetCyclicParameters(fCyclicPars);
//         _propSeq.SetScale(fNorm);
//         _propCov.SetScale(fNorm);

//         if (_doSeq) ExecutePhase1_BurnIn(10, 4);
//         if (_doND) ExecutePhase2_Mapping(10);
        
//         if (fTreeMCMC != nullptr && _doCov) {
//             ChangeNIter();
//             std::unique_ptr<TMatrixDSym> covMat(new TMatrixDSym(MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kTRUE)));
            
//             // Apply standard shrinkage
//             for(int i = 0; i < covMat->GetNrows(); i++) {
//                 for(int j = 0; j < covMat->GetNcols(); j++) {
//                     if (i != j) (*covMat)(i, j) *= 0.85; 
//                 }
//             }

//             _propCov.SetCovariance(*covMat, fSetup->NonConstParsAndYields());
//             SaveStepInfo();
//             SetTag("");
//             SetupBasicUsage();
//             SetProposalFunction(_propCov);

//             if (_tuneCovStep) ExecutePhase3_Tuning(fSetup->NonConstParsAndYields(), 10);
            
//             if (ExecutePhase4_Official() && fTreeMCMC != nullptr) {
//                 TMatrixDSym finalCovMat = MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kFALSE);
//                 if (fOutFile) {
//                     fOutFile->cd();
//                     finalCovMat.Write("PosteriorCovariance");
//                 }
//             }
//         }
//     }

//   }//namespace FIT
// }//namespace HS

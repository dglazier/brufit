#include "BruMcmc.h"
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
            
            if (max_sigma_rel > 1E-9) {
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
      // =======================================================
      // THE FIX: Deep Search the PDF Tree
      // Force the Master PDF to calculate and cache the integrals ONCE.
      // We must iterate over all components because the BruComponentsPDF 
      // is usually hidden inside a RooSimultaneous or RooProdPdf!
      // =======================================================
      std::cout << "BruMcmc: Deep searching model for BruComponentsPDFs to pre-calculate phase space..." << std::endl;
      
      std::unique_ptr<RooArgSet> comps(fPdf->getComponents());
      for (auto* arg : *comps) {
          auto* bruPdf = dynamic_cast<bru::BruComponentsPDF*>(arg);
          if (bruPdf) {
              std::cout << " -> Forcing initial integration for component: " << bruPdf->GetName() << std::endl;
              
              // This triggers CheckChange (which returns true because _Last is 0)
              // and runs DoFirstIntegrations(), permanently caching the math!
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

	  // Inside the batch loop...
	  auto subData = std::unique_ptr<RooAbsData>(fData->reduce(RooFit::EventRange(startIndex, endIndex)));
	  auto cache = std::unique_ptr<NLLCache>(new NLLCache());
	  
	  // Elegantly pass the slice into your custom BuildNLL
	  cache->finalNll = BuildNLL(subData.get(), cache->baseNll, cache->alphaVar, cache->correctedNll, cache->constraintNll, cache->totalNll);
	  
	  // Turn off expensive profiling for burn-in batches
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
    
    // ---------------------------------------------------------
    // MODULAR METHOD 2: Run the stepping loop
    // ---------------------------------------------------------
    bool BruMcmc::RunHastings(RooAbsReal* nll) {
        fParams.reset(nll->getParameters(*fData)); 
        RemoveConstantParameters(fParams.get());

        BruMetropolisHastings mh;
        mh.SetFunction(*nll);
        mh.SetParameters(*(fParams.get()));
        if (fChainParams.getSize() > 0) mh.SetChainParameters(fChainParams);
        mh.SetProposalFunction(*fPropFunc);
        mh.SetNumIters(fNumIters);

        // ==========================================================
        // --- STOCHASTIC GRADIENT BATCHING ---
        // Pass the pre-compiled NLL slices down to the MH engine
        // so it can rotate the physics landscape during the burn-in!
        // ==========================================================
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

    // ---------------------------------------------------------
    // MODULAR METHOD 3: Format and Extract Data
    // ---------------------------------------------------------
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

    // ---------------------------------------------------------
    // MAIN ENTRY POINT 
    // ---------------------------------------------------------
    Bool_t BruMcmc::MakeChain() {

      fSuccess = kFALSE; // Assume failure until it completes successfully
      
      if (!fData || !fPdf) return kFALSE;
      if (fPOI.getSize() == 0) return kFALSE;

      std::unique_ptr<RooAbsReal> baseNll;
      std::unique_ptr<RooRealVar> alphaVar;
      std::unique_ptr<RooProduct> correctedNll;
      std::unique_ptr<RooAbsReal> constraintNll;
      std::unique_ptr<RooAddition> totalNll;
        
      RooAbsReal* finalNll = nullptr;

      // ==========================================================
      // --- NLL ROUTING ---
      // If stochastic swapping is active, use the pre-compiled batch cache
      // Otherwise, build a fresh NLL using the current fData pointer
      // ==========================================================
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

      // Execute the Hastings stepping engine
      bool success = RunHastings(finalNll);

      if (defaultPropFunc) fPropFunc = nullptr; 

      if (!success) {
          CleanMakeChain();
          return kFALSE; 
      }

      SaveChainToTree();
     
      // Only deactivate the test statistic optimization if we built it dynamically here
      if (baseNll) {
          baseNll->constOptimizeTestStatistic(RooAbsArg::DeActivate, false);
      }

      CleanMakeChain();
       
      fSuccess = kTRUE; // We made it to the end! Mark as successful.
      return kTRUE;
    }
   
    ////////////////////////////////////////////////////////
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
      std::vector<bool> isYield(Npars, false); // Track yields to decouple them later
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
              isYield[i] = true; // Flag as a yield
          }
      }
        
      // --- PASS 1: Calculate Means (Circular & Arithmetic) ---
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
              means[p] /= Nentries; // Finalize arithmetic mean
          }
      }

      std::cout << "BruMcmc: Calculating Empirical Covariance for " << Nentries << " accepted steps..." << std::endl;

      // --- PASS 2: Calculate Covariance Matrix ---
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

          // Build upper triangle
          for(int i = 0; i < Npars; i++) {
              for(int j = i; j < Npars; j++) {
                  covMatSym(i, j) += delta[i] * delta[j];
              }
          }
      }

      // Finalize Matrix: Divide by (N-1) and mirror to lower triangle
      for(int i = 0; i < Npars; i++) {
          for(int j = i; j < Npars; j++) {
              covMatSym(i, j) /= (Nentries - 1);
              
              // DECOUPLE YIELD DRAG: Only zero out correlations if decoupleYields is kTRUE
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
    ////////////////////////////////////////////////////////
  
    /////////////////////////////////////////////////////////
    void BruMcmc::AddEntryBranch(){
      
      Long64_t entry = 0;
      Double_t weight = 1.0; // Default to 1 just in case
      
      if (fTreeMCMC) {
          auto entryBranch = fTreeMCMC->Branch("entry", &entry, "entry/L");
          auto weightBranch = fTreeMCMC->Branch("weight", &weight, "weight/D");
          
          // Grab the raw dataset that has the weights preserved
          const RooDataSet* internalData = fChain ? fChain->GetAsConstDataSet() : nullptr;
          
          for(entry = 0; entry < fTreeMCMC->GetEntries(); entry++) {
              if (internalData) {
                  internalData->get(entry); // Load the row
                  weight = internalData->weight(); // Extract the hidden weight
              }
              
              entryBranch->Fill();
              weightBranch->Fill();
          }
      }
    } 
    
    void BruMcmc::Result(){
      AddEntryBranch();
      AddFormulaToMCMCTree();
  
      // Use a modern C++ range-based loop over the RooArgSet
      for(auto* arg : *fParams){

        auto* targetPar = dynamic_cast<RooRealVar*>(arg);
        if (!targetPar) continue;
        
        TString pName = targetPar->GetName();
        
        // ==========================================================
        // --- STABLE TWO-PASS ALGORITHM FOR UNCERTAINTIES ---
        // ==========================================================
        Double_t sumW   = 0.0;
        Double_t sumWX  = 0.0;
        
        // PASS 1: Calculate the Mean safely
        for (int entry = 0; entry < fChainData->numEntries(); ++entry) {
            const RooArgSet* row = fChainData->get(entry); // Get the actual row
            Double_t weight = fChainData->weight();
            Double_t val    = row->getRealValue(pName);    // Safely extract by name
            
            sumW   += weight;
            sumWX  += weight * val;
        }
        
        Double_t mean = sumW > 0 ? (sumWX / sumW) : 0.0;
        Double_t sumVariance = 0.0;

        // PASS 2: Calculate Variance using (x - mu)^2 to prevent cancellation
        for (int entry = 0; entry < fChainData->numEntries(); ++entry) {
            const RooArgSet* row = fChainData->get(entry);
            Double_t weight = fChainData->weight();
            Double_t val    = row->getRealValue(pName);
            
            sumVariance += weight * (val - mean) * (val - mean);
        }
        
        Double_t sigma = sumW > 0 ? std::sqrt(sumVariance / sumW) : 0.0;
        // ==========================================================

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
    ///////////////////////////////////////////////
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
    ///////////////////////////////////////////////////////
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

    /////////////////////////////////////////////////////
    void BruMcmc::SetupBasicUsage()
    {
      fPropFunc = nullptr;
     
      TString fileName=fSetup->GetOutDir()+fSetup->GetName()+"/Results"+fSetup->GetTitle()+GetName()+GetTag()+".root";

      fOutFile.reset(TFile::Open(fileName,"recreate"));
       
     }
    ///////////////////////////////////////////////////////////////
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

    ///////////////////////////////////////////////////////////////
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

   
     //////////////////////////////////////////////////////////////

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
void BruMcmcCovariance::Run(Setup &setup, RooAbsData &fitdata) {

      fData = &fitdata;
      fSetup = &setup;
    
      InitModel();
    
      // Pass the cyclic list down to the proposal functions
      _propSeq.SetCyclicParameters(fCyclicPars);
      _propCov.SetCyclicParameters(fCyclicPars);

      // =======================================================
      // --- CRITICAL MULTI-CHAIN FIX (FRESH START) ---
      // Reset the step sizes so fresh random starts can take large leaps!
      _propSeq.SetScale(fNorm);
      _propCov.SetScale(fNorm);

      const Int_t maxRetries = 10; 
      const int numBatches = 4;

      // Scale newly randomized Yields DOWN to match the 1/4 batched data slice for Phase 1
      for (auto* y : static_range_cast<RooRealVar*>(fSetup->Yields())) {
        if (y && !y->isConstant()) {
          y->setVal(y->getVal() / numBatches);
          std::cout << "BruMcmc: Scaled Randomized Yield '" << y->GetName() 
                    << "' DOWN to " << y->getVal() << " for Phase 1 batched burn-in." << std::endl;
        }
      }
      // =======================================================

      // =======================================================
      // PHASE 1: Stochastic Burn-in (Batched Data)
      // =======================================================
      if (_doSeq == kTRUE) {
        _isBatchMode=kTRUE;
        BuildBatchedNLLs(numBatches);
        SetStochasticSwapping(500); 

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
      }
      _isBatchMode=kFALSE;

      // =======================================================
      // TRANSITION TO FULL DATA
      // Clean up batched memory and upscale the Yields BEFORE Phase 2b!
      // =======================================================
      ClearBatches(); 
      SetStochasticSwapping(0); 

      for (auto* y : static_range_cast<RooRealVar*>(fSetup->Yields())) {
        if (y && !y->isConstant()) {
          y->setVal(y->getVal() * numBatches);
          std::cout << "BruMcmc: Scaled Yield Coordinate '" << y->GetName() 
                    << "' UP to " << y->getVal() << " for Full Dataset." << std::endl;
        }
      }   

      // =======================================================
      // PHASE 2b: Stationary Covariance Mapping (FULL DATA)
      // =======================================================
      if (_doND == kTRUE) {
        std::cout << "\n*** Starting Phase 2b: Covariance Mapping (Full Data) ***" << std::endl;
        _isCovarianceMode=kTRUE;
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
      }
      _isCovarianceMode=kFALSE;
	
      // =======================================================
      // PHASE 3 & 4: Covariance Matrix Extraction, Tuning, & Official Run
      // =======================================================
      if (fTreeMCMC != nullptr && _doCov == kTRUE) {
         
        ChangeNIter();
        std::cout << "\n BruMcmcCovariance::Run Covariance Matrix Calculation" << std::endl;
         
        // 1. Build the matrix (kTRUE = decouple off-diagonal Yields for the proposal)
        std::unique_ptr<TMatrixDSym> covMat(new TMatrixDSym(MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kTRUE)));
         
        // 2. NO YIELD SCALING REQUIRED! The matrix was mapped on the Full Dataset.
        // We only apply a light 15% shrinkage to kill high-dimensional statistical noise.
        double shrinkage = 0.15; 
        for(int i = 0; i < covMat->GetNrows(); i++) {
          for(int j = 0; j < covMat->GetNcols(); j++) {
             if (i != j) {
                 (*covMat)(i, j) *= (1.0 - shrinkage);
             }
          }
        }

        _propCov.SetCovariance(*covMat, fSetup->NonConstParsAndYields());
         
        SaveStepInfo();
        SetTag("");
        SetupBasicUsage();
        SetProposalFunction(_propCov);

        // --- TUNING PHASE ---
        if (fTuneCovStep == kTRUE) {
          std::cout << "\n*** Starting Covariance Tuning Phase (250 steps) ***" << std::endl;
	  _isTuningMode=kTRUE;

          Int_t officialIters = fNumIters;
          SetNumIters(250); 
             
          Bool_t tuned = kFALSE;
          Int_t retries = 0;
             
          auto tuneCovMat = _propCov.GetCovariance();

          while(!tuned && retries < maxRetries) {
            MakeChain(); 
                 
            if (fChainAcceptance > fMinAcc && fChainAcceptance < fMaxAcc) {
              std::cout << "--> Tuning Successful! Acceptance: " << fChainAcceptance << std::endl;
              tuned = kTRUE;
            } else {
              Double_t currentScale = _propCov.StepSizeFactor();
              Double_t acc = fChainAcceptance > 0 ? fChainAcceptance : 0.01;
              currentScale *= (acc) / (fTargetAcc);
                     
              std::cout << "--> Tuning Failed should be between " << fMinAcc << "-" << fMaxAcc
                        << " (Acceptance " << fChainAcceptance << "). Adjusting scale and retrying (" 
                        << retries + 1 << "/" << maxRetries << ") will try new scale " << currentScale 
                        << " using correction " << (acc)/(fTargetAcc) << std::endl;

              _propCov.SetScale(currentScale); 
              _propCov.SetCovariance(tuneCovMat, fSetup->NonConstParsAndYields());
              retries++;
            }
          }
             
          SetNumIters(officialIters);
          _tuneCovStep = kFALSE; 
          std::cout << "\n*** Tuning Complete. Launching Official Covariance Chain ***\n" << std::endl;
        }
	_isTuningMode=kFALSE;
	_isResultMode=kTRUE;

        // --- OFFICIAL RUN ---
        Bool_t made = MakeChain();
         
        if(made == kFALSE) {
          std::cerr << " BruMcmcCovariance::Run : Official covariance chain failed." << std::endl;
        } else if (fTreeMCMC != nullptr) {
          
          std::cout << "\n*** Extracting Final Posterior Covariance Matrix ***" << std::endl;
          
          // Pass kFALSE to get the pure, fully correlated physical matrix
          TMatrixDSym finalCovMat = MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kFALSE);
          
          if (fOutFile) {
              fOutFile->cd();
              finalCovMat.Write("PosteriorCovariance");
              std::cout << "--> Matrix successfully written to file as 'PosteriorCovariance'" << std::endl;
          }
        }
      }
      _isResultMode=kFALSE;

    }
       
  }//namespace FIT
}//namespace HS

// #include "BruMcmc.h"
// #include "BruComponentsPDF.h"
// #include "BruMetropolisHastings.h"

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
// #include <RooProdPdf.h>
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
//         std::unique_ptr<RooAbsReal>& scaledDataNll,
//         std::unique_ptr<RooAbsReal>& constraintNll,
//         std::unique_ptr<RooAbsReal>& totalNll) 
// {
//     // 1. DATA NLL OPTIONS (Strictly isolated from Priors)
//     auto foptions = fSetup->FitOptions();
//     TObject* opt=nullptr;
//     if((opt=foptions.find("Save"))!=nullptr) foptions.Remove(opt);
//     if((opt=foptions.find("SumW2Error"))!=nullptr) foptions.Remove(opt);
    
//     auto cmd1 = RooFit::ConditionalObservables(fConditionalObs);
//     foptions.Add(dynamic_cast<RooCmdArg*>(&cmd1));
    
//     // CRITICAL: Force RooFit NOT to internalize constraints in the data NLL
//     auto cmd2 = RooFit::Constrain(RooArgSet()); 
//     foptions.Add(dynamic_cast<RooCmdArg*>(&cmd2));

//     // Build the pure physics NLL
//     baseNll.reset(fPdf->createNLL(*data, foptions));
//     baseNll->constOptimizeTestStatistic(RooAbsArg::Activate, false);

//     Double_t combinedScale = 1.0;
    
//     // --- 2. MC VARIANCE SCALING (BETA) DEEP SEARCH ---
//     if (!_isBatchMode && fSetup->ApplyMCVariance()) {
//         Double_t max_sigma_rel = 0.0;
//         TString domName = "None";

//         baseNll->getVal(); 

//         RooArgSet* allComps = fPdf->getComponents();
//         for (auto* obj : *allComps) {
//             if (obj->InheritsFrom("bru::BruEventsPDF")) {
//                 auto* evPdf = static_cast<bru::BruEventsPDF*>(obj);
//                 Double_t sig = evPdf->GetRelativeVariance();
//                 if (sig > max_sigma_rel) { 
//                     max_sigma_rel = sig; 
//                     domName = evPdf->GetName(); 
//                 }
//             }
//         }
//         delete allComps; 
        
//         if (max_sigma_rel > 1E-9) {
//             Double_t N_data = data->sumEntries(); 
//             Double_t betaVal = 1.0 / (1.0 + N_data * (max_sigma_rel * max_sigma_rel));
//             combinedScale *= betaVal;
            
//             std::cout << "\n=========================================" << std::endl;
//             std::cout << " [BruMcmc Likelihood Scaling] " << std::endl;
//             std::cout << " -> Dominant PDF : " << domName << " (Rel Err: " << max_sigma_rel * 100.0 << " %)" << std::endl;
//             std::cout << " -> Beta Factor  : " << betaVal << std::endl;
//             std::cout << "=========================================\n" << std::endl;
//         }
//     }
    
//     // --- 3. ORIGINAL WEIGHT SCALING (ALPHA) ---
//     if(data->isNonPoissonWeighted() && fCorrectForWeights){
//         Double_t SumW = SumWeights();
//         Double_t SumW2 = SumWeights2();
//         Double_t alphaVal = SumW / SumW2;
//         combinedScale *= alphaVal;

//         if (!_isBatchMode) {
//             std::cout << " -> Weights Alpha: " << alphaVal << std::endl;
//             std::cout << " -> Final Scale  : " << combinedScale << "\n" << std::endl;
//         }
//     }

//     // --- 4. SCALE ONLY THE DATA NLL ---
//     RooAbsReal* activeDataNll = baseNll.get();
    
//     if (combinedScale != 1.0) {
//         TString NllName = baseNll->GetName();
//         NllName.ReplaceAll("-", "m");
//         NllName.ReplaceAll("+", "p");
//         baseNll->SetName(NllName);
    
//         alphaVar.reset(new RooRealVar("alpha_weight", "alpha_weight", combinedScale));
//         alphaVar->setConstant(kTRUE);

//         scaledDataNll.reset(new RooProduct("scaled_data_nll", Form("%lf * %s", combinedScale, baseNll->GetName()), RooArgList(*alphaVar, *baseNll)));
//         activeDataNll = scaledDataNll.get();
//     }

//     // --- 5. RE-ADD THE CONSTRAINTS (UNSCALED) ---
//     if (fPriorPdf) {
//         RooArgSet emptySet;
//         // Build the NLL penalty for the prior. CloneData(kFALSE) stops redundant memory usage.
//         constraintNll.reset(fPriorPdf->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(emptySet)));
        
//         // Final NLL = (Scaled Data NLL) + (Unscaled Constraint NLL)
//         totalNll.reset(new RooAddition("total_nll", "Total NLL with Unscaled Constraints", RooArgList(*activeDataNll, *constraintNll)));
        
//         return totalNll.get();
//     }

//     return activeDataNll;
// }
// // RooAbsReal* BruMcmc::BuildNLL(RooAbsData* data,
// //         std::unique_ptr<RooAbsReal>& baseNll,
// //         std::unique_ptr<RooRealVar>& alphaVar,
// //         std::unique_ptr<RooAbsReal>& scaledDataNll,
// //         std::unique_ptr<RooAbsReal>& constraintNll,
// //         std::unique_ptr<RooAbsReal>& totalNll) 
// // {
// //     // 1. DATA NLL OPTIONS (Strictly isolated from Priors)
// //     auto foptions = fSetup->FitOptions();
// //     TObject* opt=nullptr;
// //     if((opt=foptions.find("Save"))!=nullptr) foptions.Remove(opt);
// //     if((opt=foptions.find("SumW2Error"))!=nullptr) foptions.Remove(opt);
    
// //     auto cmd1 = RooFit::ConditionalObservables(fConditionalObs);
// //     foptions.Add(dynamic_cast<RooCmdArg*>(&cmd1));
    
// //     // CRITICAL: Force RooFit NOT to internalize constraints in the data NLL
// //     auto cmd2 = RooFit::Constrain(RooArgSet()); 
// //     foptions.Add(dynamic_cast<RooCmdArg*>(&cmd2));

// //     // Build the pure physics NLL
// //     baseNll.reset(fPdf->createNLL(*data, foptions));
// //     baseNll->constOptimizeTestStatistic(RooAbsArg::Activate, false);

// //     Double_t combinedScale = 1.0;
    
// //     // --- 2. MC VARIANCE SCALING (BETA) DEEP SEARCH ---
// //     if (!_isBatchMode && fSetup->ApplyMCVariance()) {
// //         Double_t max_sigma_rel = 0.0;
// //         TString domName = "None";

// //         baseNll->getVal(); 

// //         RooArgSet* allComps = fPdf->getComponents();
// //         for (auto* obj : *allComps) {
// //             if (obj->InheritsFrom("bru::BruEventsPDF")) {
// //                 auto* evPdf = static_cast<bru::BruEventsPDF*>(obj);
// //                 Double_t sig = evPdf->GetRelativeVariance();
// //                 if (sig > max_sigma_rel) { 
// //                     max_sigma_rel = sig; 
// //                     domName = evPdf->GetName(); 
// //                 }
// //             }
// //         }
// //         delete allComps; 
        
// //         if (max_sigma_rel > 1E-9) {
// //             Double_t N_data = data->sumEntries(); 
// //             Double_t betaVal = 1.0 / (1.0 + N_data * (max_sigma_rel * max_sigma_rel));
// //             combinedScale *= betaVal;
            
// //             std::cout << "\n=========================================" << std::endl;
// //             std::cout << " [BruMcmc Likelihood Scaling] " << std::endl;
// //             std::cout << " -> Dominant PDF : " << domName << " (Rel Err: " << max_sigma_rel * 100.0 << " %)" << std::endl;
// //             std::cout << " -> Beta Factor  : " << betaVal << std::endl;
// //             std::cout << "=========================================\n" << std::endl;
// //         }
// //     }
    
// //     // --- 3. ORIGINAL WEIGHT SCALING (ALPHA) ---
// //     if(data->isNonPoissonWeighted() && fCorrectForWeights){
// //         Double_t SumW = SumWeights();
// //         Double_t SumW2 = SumWeights2();
// //         Double_t alphaVal = SumW / SumW2;
// //         combinedScale *= alphaVal;

// //         if (!_isBatchMode) {
// //             std::cout << " -> Weights Alpha: " << alphaVal << std::endl;
// //             std::cout << " -> Final Scale  : " << combinedScale << "\n" << std::endl;
// //         }
// //     }

// //     // --- 4. SCALE ONLY THE DATA NLL ---
// //     RooAbsReal* activeDataNll = baseNll.get();
    
// //     if (combinedScale != 1.0) {
// //         TString NllName = baseNll->GetName();
// //         NllName.ReplaceAll("-", "m");
// //         NllName.ReplaceAll("+", "p");
// //         baseNll->SetName(NllName);
    
// //         alphaVar.reset(new RooRealVar("alpha_weight", "alpha_weight", combinedScale));
// //         alphaVar->setConstant(kTRUE);

// //         scaledDataNll.reset(new RooProduct("scaled_data_nll", Form("%lf * %s", combinedScale, baseNll->GetName()), RooArgList(*alphaVar, *baseNll)));
// //         activeDataNll = scaledDataNll.get();
// //     }

// //     // --- 5. RE-ADD THE CONSTRAINTS (UNSCALED) ---
// //     if (fPriorPdf) {
// //         RooArgSet emptySet;
// //         // Build the NLL penalty for the prior. CloneData(kFALSE) stops redundant memory usage.
// //         constraintNll.reset(fPriorPdf->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(emptySet)));
        
// //         // Final NLL = (Scaled Data NLL) + (Unscaled Constraint NLL)
// //         totalNll.reset(new RooAddition("total_nll", "Total NLL with Unscaled Constraints", RooArgList(*activeDataNll, *constraintNll)));
        
// //         return totalNll.get();
// //     }

// //     return activeDataNll;
// // }
//     // RooAbsReal* BruMcmc::BuildNLL(RooAbsData* data,
//     //     std::unique_ptr<RooAbsPdf>& localProdPdf,
//     //     std::unique_ptr<RooAbsReal>& baseNll,
//     //     std::unique_ptr<RooRealVar>& alphaVar,
//     //     std::unique_ptr<RooProduct>& correctedNll) 
//     // {
//     //     RooAbsPdf * activePdf = fPdf;
//     //     if (fPriorPdf) {
//     //         TString prodName = TString("product_") + TString(fPdf->GetName()) + TString("_") + TString(fPriorPdf->GetName());
//     //         localProdPdf.reset(new RooProdPdf(prodName, prodName, RooArgList(*fPdf, *fPriorPdf)));
//     //         activePdf = localProdPdf.get();
//     //     }

//     //     std::unique_ptr<RooArgSet> constrainedParams(activePdf->getParameters(*data));

//     //     auto foptions = fSetup->FitOptions();
//     //     TObject* opt=nullptr;
//     //     if((opt=foptions.find("Save"))!=nullptr) foptions.Remove(opt);
//     //     if((opt=foptions.find("SumW2Error"))!=nullptr) foptions.Remove(opt);
       
//     //     auto cmd1 = RooFit::ConditionalObservables(fConditionalObs);
//     //     foptions.Add(dynamic_cast<RooCmdArg*>(&cmd1));
      
//     //     auto cmd2 = RooFit::Constrain(*constrainedParams);
//     //     foptions.Add(dynamic_cast<RooCmdArg*>(&cmd2));
 
//     //     baseNll.reset(activePdf->createNLL(*data, foptions));
//     //     baseNll->constOptimizeTestStatistic(RooAbsArg::Activate, false);
      
//     //     if(data->isNonPoissonWeighted() && fCorrectForWeights){
//     //         Double_t SumW = SumWeights();
//     //         Double_t SumW2 = SumWeights2();
//     //         Double_t alphaVal = SumW / SumW2;

//     //         TString NllName = baseNll->GetName();
//     //         NllName.ReplaceAll("-", "m");
//     //         NllName.ReplaceAll("+", "p");
//     //         baseNll->SetName(NllName);
        
//     //         alphaVar.reset(new RooRealVar("alpha_weight", "alpha_weight", alphaVal));
//     //         alphaVar->setConstant(kTRUE);

//     //         correctedNll.reset(new RooProduct("alphanll", Form("%lf * %s", alphaVal, baseNll->GetName()), RooArgList(*alphaVar, *baseNll)));
//     //         return correctedNll.get();
//     //     }

//     //     return baseNll.get();
//     // }
    
//  void BruMcmc::BuildBatchedNLLs(int numBatches) {
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

// 	  // Inside the batch loop...
// 	  auto subData = std::unique_ptr<RooAbsData>(fData->reduce(RooFit::EventRange(startIndex, endIndex)));
// 	  auto cache = std::unique_ptr<NLLCache>(new NLLCache());
	  
// 	  // Elegantly pass the slice into your custom BuildNLL
// 	  cache->finalNll = BuildNLL(subData.get(), cache->baseNll, cache->alphaVar, cache->scaledDataNll, cache->constraintNll, cache->totalNll);
	  
// 	  // Turn off expensive profiling for burn-in batches
// 	  cache->finalNll->constOptimizeTestStatistic(RooAbsArg::Activate, false);
	  
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
// bool BruMcmc::RunHastings(RooAbsReal* nll) {
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
//     // bool BruMcmc::RunHastings(RooAbsReal* nll) {
//     //     fParams.reset(nll->getParameters(*fData)); 
//     //     RemoveConstantParameters(fParams.get());

//     //     BruMetropolisHastings mh;
//     //     mh.SetFunction(*nll);
//     //     mh.SetParameters(*(fParams.get()));
//     //     if (fChainParams.getSize() > 0) mh.SetChainParameters(fChainParams);
//     //     mh.SetProposalFunction(*fPropFunc);
//     //     mh.SetNumIters(fNumIters);

//     //     fChain.reset(mh.ConstructChain()); 
//     // 	fChainAcceptance = mh.GetAcceptance();

//     //     if(fChain == nullptr){
//     //         if(fTreeMCMC){ delete fTreeMCMC; fTreeMCMC=nullptr; }
//     //         return false;
//     //     }
//     //     return true;
//     // }  

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

//       std::unique_ptr<RooAbsPdf> localProdPdf;
//       std::unique_ptr<RooAbsReal> baseNll;
//       std::unique_ptr<RooRealVar> alphaVar;
//       std::unique_ptr<RooProduct> correctedNll;
        
//       RooAbsReal* finalNll = nullptr;

//       // ==========================================================
//       // --- NLL ROUTING ---
//       // If stochastic swapping is active, use the pre-compiled batch cache
//       // Otherwise, build a fresh NLL using the current fData pointer
//       // ==========================================================
//       if (fBatchSwapFreq > 0 && !fBatchedNLLPointers.empty()) {
//           finalNll = fBatchedNLLPointers[0]; 
//       } else {
// 	//          finalNll = BuildNLL(fData, localProdPdf, baseNll, alphaVar, correctedNll);
// 	  finalNll = BuildNLL(fData, baseNll, alphaVar, scaledDataNll, constraintNll, totalNll);
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
//  TMatrixDSym BruMcmc::MakeMcmcCovarianceMatrix(TTree* tree, size_t burnin, Bool_t decoupleYields) {
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
//  void BruMcmc::AddEntryBranch(){
      
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
//  void BruMcmc::Result(){
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
//     // void BruMcmc::Result(){
//     //   AddEntryBranch();
//     //   RooArgList saveFloatFinalList(*fChainData->get()) ;

//     //   AddFormulaToMCMCTree();
  
//     //   for(Int_t i = 0; i < fParams->getSize(); i++){

//     //     auto* var = dynamic_cast<RooRealVar*>(saveFloatFinalList.at(i));
        
//     //     // ==========================================================
//     //     // --- Explicitly Calculate Weighted Mean and Sigma ---
//     //     // ==========================================================
//     //     Double_t sumW   = 0.0;
//     //     Double_t sumWX  = 0.0;
//     //     Double_t sumWX2 = 0.0;
        
//     //     for (int entry = 0; entry < fChainData->numEntries(); ++entry) {
//     //         fChainData->get(entry); // Loads the row into the dataset's internal buffer
//     //         Double_t weight = fChainData->weight();
//     //         Double_t val    = var->getVal();
            
//     //         sumW   += weight;
//     //         sumWX  += weight * val;
//     //         sumWX2 += weight * val * val;
//     //     }
        
//     //     Double_t mean = 0.0;
//     //     Double_t sigma = 0.0;
        
//     //     if (sumW > 0) {
//     //         mean = sumWX / sumW;
//     //         Double_t variance = (sumWX2 / sumW) - (mean * mean);
//     //         sigma = variance > 0 ? std::sqrt(variance) : 0.0;
//     //     }
//     //     // ==========================================================

//     //     std::cout << var->GetName() << " " << mean << " +- " << sigma << std::endl;
        
//     //     auto var2 = dynamic_cast<RooRealVar*>(fParams->find(var->GetName()));
//     //     if (var2) {
//     //         var2->setVal(mean);
//     //         var2->setError(sigma);
//     //     }
//     //   }
//     // }
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
// void BruMcmcCovariance::Run(Setup &setup, RooAbsData &fitdata) {

//       fData = &fitdata;
//       fSetup = &setup;
    
//       InitModel();
    
//       // Pass the cyclic list down to the proposal functions
//       _propSeq.SetCyclicParameters(fCyclicPars);
//       _propCov.SetCyclicParameters(fCyclicPars);

//       // =======================================================
//       // --- CRITICAL MULTI-CHAIN FIX (FRESH START) ---
//       // Reset the step sizes so fresh random starts can take large leaps!
//       _propSeq.SetScale(fNorm);
//       _propCov.SetScale(fNorm);

//       const Int_t maxRetries = 10; 
//       const int numBatches = 4;

//       // Scale newly randomized Yields DOWN to match the 1/4 batched data slice for Phase 1
//       for (auto* y : static_range_cast<RooRealVar*>(fSetup->Yields())) {
//         if (y && !y->isConstant()) {
//           y->setVal(y->getVal() / numBatches);
//           std::cout << "BruMcmc: Scaled Randomized Yield '" << y->GetName() 
//                     << "' DOWN to " << y->getVal() << " for Phase 1 batched burn-in." << std::endl;
//         }
//       }
//       // =======================================================

//       // =======================================================
//       // PHASE 1: Stochastic Burn-in (Batched Data)
//       // =======================================================
//       if (_doSeq == kTRUE) {
//         _isBatchMode=kTRUE;
//         BuildBatchedNLLs(numBatches);
//         SetStochasticSwapping(500); 

//         ChangeNIter();
//         SetTag("1DStep");
//         SetupBasicUsage();
//         SetProposalFunction(_propSeq);
        
//         Bool_t made = MakeChain();
//         Int_t retries = 0;
//         while(made == kFALSE && retries < maxRetries) {
//           std::cout << "\n*** BruMcmcCovariance: 1DStep Failed! Retrying (" << retries + 1 << "/" << maxRetries << ") ***\n" << std::endl;
//           _propSeq.SetScale(fNorm); 
//           made = MakeChain();
//           retries++;
//         }
//       }
//       _isBatchMode=kFALSE;

//       // =======================================================
//       // TRANSITION TO FULL DATA
//       // Clean up batched memory and upscale the Yields BEFORE Phase 2b!
//       // =======================================================
//       ClearBatches(); 
//       SetStochasticSwapping(0); 

//       for (auto* y : static_range_cast<RooRealVar*>(fSetup->Yields())) {
//         if (y && !y->isConstant()) {
//           y->setVal(y->getVal() * numBatches);
//           std::cout << "BruMcmc: Scaled Yield Coordinate '" << y->GetName() 
//                     << "' UP to " << y->getVal() << " for Full Dataset." << std::endl;
//         }
//       }   

//       // =======================================================
//       // PHASE 2b: Stationary Covariance Mapping (FULL DATA)
//       // =======================================================
//       if (_doND == kTRUE) {
//         std::cout << "\n*** Starting Phase 2b: Covariance Mapping (Full Data) ***" << std::endl;
//         _isCovarianceMode=kTRUE;
//         ChangeNIter();
//         SaveStepInfo();
//         SetTag("NDStep");
//         SetupBasicUsage();
//         SetProposalFunction(_propSeq);
//         _propSeq.SetIsSequential(kFALSE);
        
//         Bool_t made = MakeChain();
//         Int_t retries = 0;
//         while(made == kFALSE && retries < maxRetries) {
//           std::cout << "\n*** BruMcmcCovariance: NDStep Failed! Retrying (" << retries + 1 << "/" << maxRetries << ") ***\n" << std::endl;
//           _propSeq.SetScale(fNorm); 
//           made = MakeChain();
//           retries++;
//         }
//       }
//       _isCovarianceMode=kFALSE;
	
//       // =======================================================
//       // PHASE 3 & 4: Covariance Matrix Extraction, Tuning, & Official Run
//       // =======================================================
//       if (fTreeMCMC != nullptr && _doCov == kTRUE) {
         
//         ChangeNIter();
//         std::cout << "\n BruMcmcCovariance::Run Covariance Matrix Calculation" << std::endl;
         
//         // 1. Build the matrix (kTRUE = decouple off-diagonal Yields for the proposal)
//         std::unique_ptr<TMatrixDSym> covMat(new TMatrixDSym(MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kTRUE)));
         
//         // 2. NO YIELD SCALING REQUIRED! The matrix was mapped on the Full Dataset.
//         // We only apply a light 15% shrinkage to kill high-dimensional statistical noise.
//         double shrinkage = 0.15; 
//         for(int i = 0; i < covMat->GetNrows(); i++) {
//           for(int j = 0; j < covMat->GetNcols(); j++) {
//              if (i != j) {
//                  (*covMat)(i, j) *= (1.0 - shrinkage);
//              }
//           }
//         }

//         _propCov.SetCovariance(*covMat, fSetup->NonConstParsAndYields());
         
//         SaveStepInfo();
//         SetTag("");
//         SetupBasicUsage();
//         SetProposalFunction(_propCov);

//         // --- TUNING PHASE ---
//         if (fTuneCovStep == kTRUE) {
//           std::cout << "\n*** Starting Covariance Tuning Phase (250 steps) ***" << std::endl;
// 	  _isTuningMode=kTRUE;

//           Int_t officialIters = fNumIters;
//           SetNumIters(250); 
             
//           Bool_t tuned = kFALSE;
//           Int_t retries = 0;
             
//           auto tuneCovMat = _propCov.GetCovariance();

//           while(!tuned && retries < maxRetries) {
//             MakeChain(); 
                 
//             if (fChainAcceptance > fMinAcc && fChainAcceptance < fMaxAcc) {
//               std::cout << "--> Tuning Successful! Acceptance: " << fChainAcceptance << std::endl;
//               tuned = kTRUE;
//             } else {
//               Double_t currentScale = _propCov.StepSizeFactor();
//               Double_t acc = fChainAcceptance > 0 ? fChainAcceptance : 0.01;
//               currentScale *= (acc) / (fTargetAcc);
                     
//               std::cout << "--> Tuning Failed should be between " << fMinAcc << "-" << fMaxAcc
//                         << " (Acceptance " << fChainAcceptance << "). Adjusting scale and retrying (" 
//                         << retries + 1 << "/" << maxRetries << ") will try new scale " << currentScale 
//                         << " using correction " << (acc)/(fTargetAcc) << std::endl;

//               _propCov.SetScale(currentScale); 
//               _propCov.SetCovariance(tuneCovMat, fSetup->NonConstParsAndYields());
//               retries++;
//             }
//           }
             
//           SetNumIters(officialIters);
//           _tuneCovStep = kFALSE; 
//           std::cout << "\n*** Tuning Complete. Launching Official Covariance Chain ***\n" << std::endl;
//         }
// 	_isTuningMode=kFALSE;
// 	_isResultMode=kTRUE;

//         // --- OFFICIAL RUN ---
//         Bool_t made = MakeChain();
         
//         if(made == kFALSE) {
//           std::cerr << " BruMcmcCovariance::Run : Official covariance chain failed." << std::endl;
//         } else if (fTreeMCMC != nullptr) {
          
//           std::cout << "\n*** Extracting Final Posterior Covariance Matrix ***" << std::endl;
          
//           // Pass kFALSE to get the pure, fully correlated physical matrix
//           TMatrixDSym finalCovMat = MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kFALSE);
          
//           if (fOutFile) {
//               fOutFile->cd();
//               finalCovMat.Write("PosteriorCovariance");
//               std::cout << "--> Matrix successfully written to file as 'PosteriorCovariance'" << std::endl;
//           }
//         }
//       }
//       _isResultMode=kFALSE;

//     }
       
//   }//namespace FIT
// }//namespace HS

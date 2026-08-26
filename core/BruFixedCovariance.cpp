#include "BruFixedCovariance.h"
#include "BruMappedRhat.h" // NEW: Required for the diagnostic calculations
#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>
#include <RooRealVar.h>
#include <iostream>
#include <memory>
#include <algorithm>

namespace HS {
namespace FIT {

    // ========================================================================
    // BruCovarianceReader Implementation
    // ========================================================================
    Bool_t BruCovarianceReader::SaveThis(const TString& filePath, const TString& matrixName, 
                                         const TMatrixDSym& mat, const std::vector<std::string>& params) {
        std::unique_ptr<TFile> file(TFile::Open(filePath, "RECREATE"));
        if (!file || file->IsZombie()) {
            std::cerr << "BruCovarianceReader::SaveThis - ERROR: Cannot create file " << filePath << std::endl;
            return kFALSE;
        }

        // 1. Write the Matrix
        mat.Write(matrixName);

        // 2. Build and write the Parameter Name List
        TList nameList;
        nameList.SetName(matrixName + "_Names");
        for (const auto& name : params) {
            nameList.Add(new TObjString(name.c_str()));
        }
        
        nameList.Write(nameList.GetName(), TObject::kSingleKey);
        
        file->Close();
        std::cout << "BruCovarianceReader::SaveThis - SUCCESS: Matrix and parameter list saved to " << filePath << std::endl;
        return kTRUE;
    }

    Bool_t BruCovarianceReader::Load(const TString& filePath, const TString& matrixName, const TString& listName) {
        fSavedNames.clear();
        
        std::unique_ptr<TFile> file(TFile::Open(filePath, "READ"));
        if (!file || file->IsZombie()) {
            std::cerr << "BruCovarianceReader::Load - ERROR: Cannot open file " << filePath << std::endl;
            return kFALSE;
        }

        auto* rawMatrix = dynamic_cast<TMatrixDSym*>(file->Get(matrixName));
        if (!rawMatrix) {
            std::cerr << "BruCovarianceReader::Load - ERROR: Cannot find TMatrixDSym named " << matrixName << std::endl;
            return kFALSE;
        }
        fMatrix.ResizeTo(*rawMatrix);
        fMatrix = *rawMatrix; 

        // Only attempt to load the parameter list if a listName is provided.
        // Matrices saved directly via FitManager (Mode 3) do not have this list attached.
        if (listName != "") {
            auto* rawList = dynamic_cast<TList*>(file->Get(listName));
            if (!rawList) {
                std::cerr << "BruCovarianceReader::Load - ERROR: Cannot find TList of names " << listName << std::endl;
                return kFALSE;
            }

            for (int i = 0; i < rawList->GetSize(); ++i) {
                auto* objStr = dynamic_cast<TObjString*>(rawList->At(i));
                if (objStr) {
                    fSavedNames.push_back(objStr->GetString().Data());
                }
            }
        }
        
        return kTRUE;
    }

    Bool_t BruCovarianceReader::AlignAndValidate(const RooArgSet& currentPars) {
        Int_t nPars = currentPars.getSize();
        
        // 1. Check Dimensions
        if (fMatrix.GetNrows() != nPars || fMatrix.GetNcols() != nPars) {
            std::cerr << "\nBruCovarianceReader::AlignAndValidate - FATAL: Dimension mismatch." << std::endl;
            std::cerr << " -> Loaded Matrix : " << fMatrix.GetNrows() << "x" << fMatrix.GetNcols() << std::endl;
            std::cerr << " -> Model Expected: " << nPars << " non-constant parameters." << std::endl;
            return kFALSE;
        }

        // If no names were loaded (e.g. reading from a direct previous result), 
        // we assume the matrix ordering inherently matches the current model perfectly (1-to-1).
        if (fSavedNames.empty()) {
            std::cout << "BruCovarianceReader::AlignAndValidate - No parameter list provided. Assuming 1-to-1 direct mapping." << std::endl;
            return kTRUE;
        }

        // 2. Extract Target Names
        std::vector<std::string> targetNames;
        for (auto* arg : currentPars) {
            targetNames.push_back(arg->GetName());
        }

        // 3. Build Index Map [Target Index -> Source Matrix Index]
        std::vector<int> indexMap(nPars, -1);
        Bool_t mappingSuccessful = kTRUE;

        for (size_t t = 0; t < targetNames.size(); ++t) {
            auto it = std::find(fSavedNames.begin(), fSavedNames.end(), targetNames[t]);
            if (it != fSavedNames.end()) {
                indexMap[t] = std::distance(fSavedNames.begin(), it);
            } else {
                std::cerr << "BruCovarianceReader::AlignAndValidate - FATAL: Required model parameter '" 
                          << targetNames[t] << "' is MISSING from the loaded matrix." << std::endl;
                mappingSuccessful = kFALSE;
            }
        }

        if (!mappingSuccessful) {
            std::cerr << "\n--- PARAMETERS AVAILABLE IN LOADED MATRIX ---" << std::endl;
            for (const auto& name : fSavedNames) std::cerr << "  " << name << std::endl;
            std::cerr << "---------------------------------------------\n" << std::endl;
            return kFALSE;
        }

        // 4. Check if reordering is actually necessary
        bool needsReorder = false;
        for (size_t i = 0; i < indexMap.size(); ++i) {
            if (indexMap[i] != static_cast<int>(i)) { 
                needsReorder = true; 
                break; 
            }
        }

        // 5. Permute the Matrix Physically
        if (needsReorder) {
            std::cout << "BruCovarianceReader::AlignAndValidate - Matrix alignment required. Reordering rows/columns..." << std::endl;
            TMatrixDSym alignedMat(nPars);
            
            for (int i = 0; i < nPars; ++i) {
                for (int j = 0; j < nPars; ++j) {
                    // Map the target (i, j) to the original source coordinates
                    alignedMat(i, j) = fMatrix(indexMap[i], indexMap[j]);
                }
            }
            
            fMatrix = alignedMat;          // Replace with properly ordered matrix
            fSavedNames = targetNames;     // Sync internal name list
        }

        std::cout << "BruCovarianceReader::AlignAndValidate - SUCCESS: Matrix dimensions and parameter alignment complete." << std::endl;
        return kTRUE;
    }

    // ========================================================================
    // BruMcmcFixedCovariance Implementation
    // ========================================================================
    
    // Mode 1 Constructor: Static file mapping
    BruMcmcFixedCovariance::BruMcmcFixedCovariance(const TString& covFilePath, const TString& matrixName, 
                                                   std::vector<Int_t> Niters, Int_t Nburn, Float_t norm, 
                                                   float target, float accmin, float accmax)
        : BruMcmcCovariance(Niters, Nburn, norm, target, accmin, accmax),
          _loadMode(CovLoadMode::kStaticFile),
          fPathStr1(covFilePath),
          fMatrixName(matrixName),
          fListName(matrixName + "_Names") 
    {
        SetNameTitle("BruMcmcFixedCovariance", "BruMcmcFixedCovariance minimiser");
    }

    // Mode 2 & 3 Constructor: Dynamic file mapping per bin
    BruMcmcFixedCovariance::BruMcmcFixedCovariance(CovLoadMode mode, const TString& pathStr1, const TString& pathStr2, const TString& matrixName, 
                                                   std::vector<Int_t> Niters, Int_t Nburn, Float_t norm, 
                                                   float target, float accmin, float accmax)
        : BruMcmcCovariance(Niters, Nburn, norm, target, accmin, accmax),
          _loadMode(mode),
          fPathStr1(pathStr1),
          fPathStr2(pathStr2),
          fMatrixName(matrixName)
    {
        SetNameTitle("BruMcmcFixedCovariance", "BruMcmcFixedCovariance minimiser");
        
        // For Mode 3 (Previous Result), the matrix is saved natively by FitManager without a custom name list.
        if (_loadMode == CovLoadMode::kPreviousResult) {
            fListName = ""; 
        } else {
            fListName = matrixName + "_Names";
        }
    }

    void BruMcmcFixedCovariance::Run(Setup& setup, RooAbsData& fitdata) {
        fData = &fitdata;
        fSetup = &setup;
        
        InitModel();
        
        _propSeq.SetCyclicParameters(fCyclicPars);
        _propCov.SetCyclicParameters(fCyclicPars);
        _propSeq.SetScale(fNorm);
        _propCov.SetScale(fNorm);
        
        // CRITICAL: Ensure the yield mask is configured before any matrix tuning operations
        _propCov.SetYields(fSetup->Yields());

        auto activePars = fSetup->NonConstParsAndYields();

        // 1. Run Burn-in 
        if (!ExecutePhase1_BurnIn(10, 1)) {
            std::cerr << "BruMcmcFixedCovariance::Run - FATAL: Phase 1 Burn-in failed." << std::endl;
            return;
        }

        // Resolve the actual covariance file path based on the selected mode
        TString actualCovPath;
        if (_loadMode == CovLoadMode::kStaticFile) {
            actualCovPath = fPathStr1;
        } else if (_loadMode == CovLoadMode::kBinDirectory) {
            actualCovPath = fPathStr1 + "/" + setup.GetName() + "/" + fPathStr2;
        } else if (_loadMode == CovLoadMode::kPreviousResult) {
            actualCovPath = fPathStr1 + "/" + setup.GetName() + "/Results" + fPathStr2 + ".root";
        }

        // 2. Load and Align Matrix
        std::cout << "\n*** Phase 2: Loading & Aligning External Covariance Matrix ***" << std::endl;
        std::cout << "--> Looking for Matrix in: " << actualCovPath << std::endl;
        
        if (!_covReader.Load(actualCovPath, fMatrixName, fListName)) {
            std::cerr << "BruMcmcFixedCovariance::Run - FATAL: Could not load covariance matrix from " << actualCovPath << std::endl;
            return; 
        }
        
        if (!_covReader.AlignAndValidate(activePars)) return; 

        TMatrixDSym covMat = _covReader.GetMatrix();

        // 3. Apply Shrinkage
        std::cout << "--> Applying dimensional shrinkage factor of " << _shrinkage << " to off-diagonals." << std::endl;
        for(int i = 0; i < covMat.GetNrows(); i++) {
            for(int j = 0; j < covMat.GetNcols(); j++) {
                if (i != j) covMat(i, j) *= (1.0 - _shrinkage);
            }
        }

        _propCov.SetCovariance(covMat, activePars);
        SaveStepInfo();
        SetTag("");
        SetupBasicUsage();
        SetProposalFunction(_propCov);

        ChangeNIter(); 

        // 4. Tune & Run
        if (_tuneCovStep) {
            if (_tuneMode == McmcTuneMode::kMappedRhat) {
                std::cout << "--> Executing Diagnostic Tuning (R-hat + ESS)..." << std::endl;
                ExecuteTuning_MappedRhat(activePars, 10);
            } else {
                std::cout << "--> Executing Standard Acceptance Tuning..." << std::endl;
                ExecuteTuning_Acceptance(activePars, 10);
            }
        }
        
        // Official Run (Automatically merges tuning trees!)
        if (ExecutePhase4_Official() && fTreeMCMC != nullptr) {
            std::cout << "\n*** Extracting Final Posterior Covariance Matrix ***" << std::endl;
            
            // Extract the refined empirical matrix from the accepted steps.
            TMatrixDSym finalCovMat = MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kFALSE, 0.0, 0.0);
            
            // --- NEW: Calculate Diagnostics and Save to File ---
            std::cout << "\n*** Calculating Sub-Set Diagnostics (R-hat & ESS) ***" << std::endl;
            BruMappedRhat diagHelper;
            std::vector<RooArgList> paramGroups = BuildDiagnosticGroups(activePars);
            
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

            if (maxRhat != BruMappedRhat::kConvergenceFailure && maxRhat <= _rhatTarget) {
                std::cout << "--> [SUCCESS] Official chain converged perfectly!" << std::endl;
            } else {
                std::cout << "--> [WARNING] Official chain finished with sub-optimal convergence (Worst R-hat: " << maxRhat << ")." << std::endl;
            }

            if (fOutFile) {
                fOutFile->cd();
                finalCovMat.Write("PosteriorCovariance");
                std::cout << "--> Matrix successfully written to file as 'PosteriorCovariance'" << std::endl;
                
                // Write the diagnostics tree
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
    }
  
} // namespace FIT
} // namespace HS

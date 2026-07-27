#include "BruFixedCovariance.h"
#include <TFile.h>
#include <TObjString.h>
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
   BruMcmcFixedCovariance::BruMcmcFixedCovariance(const TString& covFilePath, const TString& matrixName, 
                                                   std::vector<Int_t> Niters, Int_t Nburn, Float_t norm, 
                                                   float target, float accmin, float accmax)
        : BruMcmcCovariance(Niters, Nburn, norm, target, accmin, accmax),
          fCovFilePath(covFilePath),
          fMatrixName(matrixName),
          fListName(matrixName + "_Names") 
    {
        SetNameTitle("BruMcmcFixedCovariance", "BruMcmcFixedCovariance minimiser");
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
        // ExecutePhase1_BurnIn internally calls ChangeNIter(), grabbing Niters[0]
        if (!ExecutePhase1_BurnIn(10, 1)) {
            std::cerr << "BruMcmcFixedCovariance::Run - FATAL: Phase 1 Burn-in failed." << std::endl;
            return;
        }

        // 2. Load and Align Matrix
        std::cout << "\n*** Phase 2: Loading & Aligning External Covariance Matrix ***" << std::endl;
        if (!_covReader.Load(fCovFilePath, fMatrixName, fListName)) return; 
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

        // Advance the internal iteration counter to Niters[1]
        // before launching the tuning and official covariance chains.
        ChangeNIter(); 

        // 4. Tune & Run
        // _tuneCovStep acts as your DoTune() flag.
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
            // Explicitly pass 0.0 for shrinkage/floor so the final saved matrix is pure.
            TMatrixDSym finalCovMat = MakeMcmcCovarianceMatrix(fTreeMCMC, fNumBurnInSteps, kFALSE, 0.0, 0.0);
            
            if (fOutFile) {
                fOutFile->cd();
                finalCovMat.Write("PosteriorCovariance");
                std::cout << "--> Matrix successfully written to file as 'PosteriorCovariance'" << std::endl;
            }
        }
    }
  
} // namespace FIT
} // namespace HS

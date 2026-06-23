/**
 * @file MCEventCache.cpp
 */

#include "MCEventCache.h"
#include <TVectorD.h>
#include <TEntryList.h>
#include <TDirectory.h>
#include <TMath.h>
#include <RooFitLegacy/RooCatTypeLegacy.h>
#include <iostream>

namespace bru {

    bool MCEventCache::LoadTree(TTree* tree, TString& cut,
                                const std::vector<TString>& varNames,
                                const std::vector<TString>& catNames,
                                const TString& truthPrefix,
                                const HS::FIT::WeightsConfig& wgtsConf,
                                TTree* mcGenTree) {
                                    
        if (!tree || !tree->GetEntries()) return false;

        _Nvars = varNames.size();
        _Ncats = catNames.size();
        _HasMCGenTree = (mcGenTree != nullptr);

        tree->ResetBranchAddresses();
        if (_HasMCGenTree) mcGenTree->ResetBranchAddresses();

        bool branchStatus = true;

        TVectorD MCVar(_Nvars);
        TVectorD GenVar(_Nvars);
        TVectorD MCGenVar(_Nvars);
        std::vector<Int_t> MCCat(_Ncats);
        std::vector<Int_t> GenCat(_Ncats);
        std::vector<Int_t> MCGenCat(_Ncats);

        std::vector<Int_t> GotGenVar(_Nvars, 0);
        std::vector<Int_t> GotGenCat(_Ncats, 0);

        // 1. Setup Real Variable Branches
        for (size_t i = 0; i < _Nvars; i++) {
            if (tree->GetBranch(varNames[i])) {
                tree->SetBranchStatus(varNames[i], true);
                tree->SetBranchAddress(varNames[i], &MCVar[i]);
                if (tree->GetBranch(truthPrefix + varNames[i])) {
                    tree->SetBranchStatus(truthPrefix + varNames[i], true);
                    tree->SetBranchAddress(truthPrefix + varNames[i], &GenVar[i]);
                    GotGenVar[i] = 1;
                }
            } else {
                std::cout << "Warning: Branch " << varNames[i] << " not found\n";
                if (cut.Contains(varNames[i])) {
                    TString newcut = cut;
                    newcut.Replace(newcut.Index(varNames[i]) - 2, 2, "");
                    newcut.Replace(newcut.Index(varNames[i] + ">"), (newcut.Index(varNames[i] + "<") - newcut.Index(varNames[i] + ">")) * 2 - 1, "");
                    if (newcut == TString("&&")) newcut = "";
                    if (newcut(0, 2) == "&&") newcut.Remove(0, 2);
                    cut = newcut;
                }
                branchStatus = false;
            }

            if (_HasMCGenTree && mcGenTree->GetBranch(varNames[i])) {
                mcGenTree->SetBranchStatus(varNames[i], true);
                mcGenTree->SetBranchAddress(varNames[i], &MCGenVar[i]);
            }
        }

        // 2. Setup Category Branches
        for (size_t i = 0; i < _Ncats; i++) {
            if (tree->GetBranch(catNames[i])) {
                tree->SetBranchStatus(catNames[i], true);
                tree->SetBranchAddress(catNames[i], &MCCat[i]);
                if (tree->GetBranch(truthPrefix + catNames[i])) {
                    tree->SetBranchStatus(truthPrefix + catNames[i], true);
                    tree->SetBranchAddress(truthPrefix + catNames[i], &GenCat[i]);
                    GotGenCat[i] = 1;
                }
            } else {
                std::cout << "Warning: Category Branch " << catNames[i] << " not found\n";
                branchStatus = false;
            }
        }

        // 3. Allocate Vectors
        _NTreeEntries = tree->GetEntries();
        _vecReal.resize(_NTreeEntries * _Nvars);
        _vecRealGen.resize(_NTreeEntries * _Nvars);
        _vecCat.resize(_NTreeEntries * _Ncats);
        _vecCatGen.resize(_NTreeEntries * _Ncats);

        if (_HasMCGenTree) {
            _NMCGenTreeEntries = mcGenTree->GetEntries();
            _vecRealMCGen.resize(_NMCGenTreeEntries * _Nvars);
            _vecCatMCGen.resize(_NMCGenTreeEntries * _Ncats);
        }

        // 4. Load External Weights (Isolated locally)
        Double_t idVal = 0;
        Int_t spId = -1;
        std::unique_ptr<HS::FIT::Weights> inWeights = nullptr;

        if (wgtsConf.IsValid()) {
            _EvWeights.clear();
            inWeights = std::make_unique<HS::FIT::Weights>();
            inWeights->LoadSavedDisc(wgtsConf.File(), wgtsConf.ObjName());
            if (inWeights->GetSpeciesID(wgtsConf.Species()) == -1) {
                std::cout << "ERROR: Species " << wgtsConf.Species() << " not found.\n";
            } else {
                if (tree->GetBranch(inWeights->GetIDName())) {
                    _UseEvWeights = true;
                    tree->SetBranchStatus(inWeights->GetIDName(), true);
                    tree->SetBranchAddress(inWeights->GetIDName(), &idVal);
                    _EvWeights.resize(_NTreeEntries);
                    spId = inWeights->GetSpeciesID(wgtsConf.Species());
                }
            }
        }

        // 5. Apply Cuts & Build EntryList
        tree->Draw(">>elist", cut, "entrylist");
        auto *elist = dynamic_cast<TEntryList*>(gDirectory->Get("elist"));
        tree->SetEntryList(elist);
        _NTreeEntries = elist->GetN();

        // 6. Populate Flattened Arrays
        Long64_t corrEvent = 0;
        for (Long64_t iEvent = 0; iEvent < _NTreeEntries; iEvent++) {
            Long64_t entryNumber = tree->GetEntryNumber(iEvent);
            if (entryNumber < 0) break;
            Long64_t localEntry = tree->LoadTree(entryNumber);
            if (localEntry < 0) break;
            tree->GetEntry(localEntry);

            bool removeNaNEvent = false;
            for (size_t ip = 0; ip < _Nvars; ip++) {
                if (TMath::IsNaN(MCVar[ip])) {
                    removeNaNEvent = true;
                } else {
                    _vecReal[corrEvent * _Nvars + ip] = MCVar[ip];
                    _vecRealGen[corrEvent * _Nvars + ip] = GotGenVar[ip] ? GenVar[ip] : MCVar[ip];
                }
            }
            if (removeNaNEvent) continue;

            _TreeEntryNumber.push_back(localEntry);

            if (_UseEvWeights) {
                inWeights->GetEntryBinarySearch(static_cast<Long64_t>(idVal));
                _EvWeights[corrEvent] = inWeights->GetWeight(spId);
            }

            for (size_t ip = 0; ip < _Ncats; ip++) {
                _vecCat[corrEvent * _Ncats + ip] = MCCat[ip];
                _vecCatGen[corrEvent * _Ncats + ip] = GotGenCat[ip] ? GenCat[ip] : MCCat[ip];
            }
            corrEvent++;
        }

        _NTreeEntries = corrEvent;

        tree->SetEntryList(nullptr);
        delete elist;

        // Populate MCGenTree if exists
        if (_HasMCGenTree) {
            mcGenTree->Draw(">>elistMCGen", "", "entrylistMCGen"); 
            auto* elistMCGen = dynamic_cast<TEntryList*>(gDirectory->Get("elistMCGen"));
            mcGenTree->SetEntryList(elistMCGen);
            _NMCGenTreeEntries = elistMCGen->GetN();
            
            for (Long64_t iEvent = 0; iEvent < _NMCGenTreeEntries; iEvent++) {
                Long64_t entryNumber = mcGenTree->GetEntryNumber(iEvent);
                if (entryNumber < 0) break;
                Long64_t localEntry = mcGenTree->LoadTree(entryNumber);
                if (localEntry < 0) break;
                mcGenTree->GetEntry(localEntry);
                for (size_t ip = 0; ip < _Nvars; ip++) _vecRealMCGen[iEvent * _Nvars + ip] = MCGenVar[ip];
                for (size_t ip = 0; ip < _Ncats; ip++) _vecCatMCGen[iEvent * _Ncats + ip] = MCGenCat[ip];
            }
            mcGenTree->SetEntryList(nullptr);
            delete elistMCGen;
        }

        tree->ResetBranchAddresses();  
        tree->Reset(); 
        if (_HasMCGenTree) {
            mcGenTree->ResetBranchAddresses();
            mcGenTree->Reset();
        }

        return branchStatus;
    }

} // namespace bru

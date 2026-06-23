/**
 * @file Data.h
 * @brief Classes for managing datasets, trees, and weights for fitting.
 * @details Handles the loading of ROOT trees, application of cuts, and 
 * conversion into RooDataSets for use in RooFit. Designed to be completely 
 * thread-safe and copyable for multi-process execution.
 */

#pragma once

#include "Setup.h"
#include "BootStrapper.h"
#include "FiledTree.h"
#include "Weights.h"
#include <RooAbsData.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <TString.h>
#include <TObject.h>
#include <TObjString.h>

#include <utility>
#include <memory>
#include <vector>

namespace HS {
namespace FIT {

    // Modern C++ Type Aliases
    using strings_t   = std::vector<TString>;
    using dset_uptr   = std::unique_ptr<RooDataSet>;
    using roodsets_t  = std::vector<RooDataSet*>;
    
    /**
     * @class FitData
     * @brief Base class for dataset management in BruFit.
     */
    class FitData : public TObject {
      
    public:
        FitData() = default;
        FitData(const FitData&) = default;
        FitData(FitData&&) = default;
        ~FitData() override = default;
        FitData& operator=(const FitData& other) = default;
        FitData& operator=(FitData&& other) = default;

    private:
        ClassDefOverride(HS::FIT::FitData, 1);
    };
    
    /**
     * @class DataEvents
     * @brief Concrete implementation for handling unbinned event data.
     * @details Manages the translation of raw TTrees into RooDataSets. 
     * Uses transient file pointers to ensure the class can be cleanly 
     * copied across isolated worker processes without file handle collisions.
     */
    class DataEvents : public FitData {
      
    public:
        DataEvents() = default;
        
        /** @brief Constructs a data manager pointing to specific ROOT files. */
        DataEvents(Setup& setup, TString tname, const strings_t& files);

        // Standard Rule of 5 semantics
        DataEvents(const DataEvents&) = default;
        DataEvents& operator=(const DataEvents& other) = default;
        DataEvents(DataEvents&&) = default;
        ~DataEvents() override = default;
        DataEvents& operator=(DataEvents&& other) = default;

        /** @brief Returns the total number of data files/bins currently loaded. */
        UInt_t GetN() const { return fFileNames.size(); }
        
        /**
         * @brief Opens the requested file, applies cuts, and returns a RooDataSet.
         * @param iset The index of the file/bin to load.
         * @return A unique_ptr owning the new RooDataSet.
         */
        dset_uptr Get(const UInt_t iset);
        
        /** @name Accessors */
        ///@{
        strings_t FileNames() const { return fFileNames; }
        TString FileName(UInt_t ii) const { return fFileNames[ii]; }
        TString ParentName() const { return fParentName; }
        TString ParentTreeName() const { return fTreeName; }
        void SetParentName(TString name) { fParentName = std::move(name); }
        void SetParentTreeName(TString name) { fParentTreeName = std::move(name); }
        const RooRealVar* WeightVar() const { return fWeightVar.get(); }
        ///@}

        /** @name Initialization and Configuration */
        ///@{
        void Load(Setup& setup, const TString& tname, const strings_t& files);
        void LoadSetup(Setup* setup) { fSetup = setup; }
        
        void BootStrap(Int_t N) {
            fNBoots = N;
            fBootStrap = std::make_unique<BootStrapper>(N);
        }
        void LoadBootStrap(const TString& tname, strings_t files);
        void Toys(Int_t N) { fNToys = N; }
        
        /** @brief Determines the logical bin/group index for toys and bootstraps. */
        Int_t GetDataBin(Int_t ii) {
            if (fNBoots > 0 && !fBootStrap.get()) BootStrap(fNBoots);
            if (fBootStrap.get()) return fBootStrap->GetGroup(ii);
            else if (fNToys > 0) return (int)std::round(ii / fNToys);
            return ii;
        }
        
        TString GetItemName(Int_t ii);
        
        /** @brief Configures external weights to be merged into the dataset. */
        void LoadWeights(TString wname, TString fname, TString wobj = "HSsWeights");
        ///@}
      
    private:
        HS::FIT::Setup* fSetup = nullptr;          ///<! Pointer to the active fit configuration
        strings_t fFileNames;                      ///< List of paths to the data ROOT files
        TString fTreeName;                         ///< Name of the TTree inside the ROOT files
        TString fParentName;
        TString fParentTreeName;
        
        TString fInWeightName;
        TString fInWeightFile;
        TString fInWeightObjName;
        
        std::shared_ptr<BootStrapper> fBootStrap;  ///< Shared manager for bootstrap permutations
        Int_t fNBoots = -1;
        Int_t fNToys = -1;
        
        WeightsConfig fWgtsConf;                   ///< Metadata for locating external weights

        std::shared_ptr<RooRealVar> fWeightVar;    ///< Safely shared variable representing the applied weight
      
        ClassDefOverride(HS::FIT::DataEvents, 1);
    };
    
} // namespace FIT
} // namespace HS


// ////////////////////////////////////////////////////////////////
// ///
// ///Class:               Data
// ///Description:
// ///           

// #pragma once


// #include "Setup.h"
// #include "BootStrapper.h"
// #include "FiledTree.h"
// #include "Weights.h"
// #include <RooAbsData.h>
// #include <RooDataSet.h>
// #include <RooDataHist.h>
// #include <TString.h>
// #include <TObject.h>
// #include <TObjString.h>

// #include <utility>

// #include <memory>


// namespace HS{
//   namespace FIT{

//     using strings_t = std::vector<TString>;
//     using weights_ptr = std::shared_ptr<HS::FIT::Weights>;
    
//     class FitData : public TObject {
      
//     public:
//       FitData()=default;
//       FitData(const FitData&)=default;
//       FitData(FitData&&)=default;
//       ~FitData() override =default;
//       FitData& operator=(const FitData& other)=default;
//       FitData& operator=(FitData&& other) = default;

//       //virtual RooAbsData& Get() = 0;
//     protected:
      
//     private:
//       //RooAbsData fData; //dataset to be fitted

//       ClassDefOverride(HS::FIT::FitData,1);
      
//     };//class FitData
    
//     //////////////////////////////////////////////////
//     using dset_uptr = std::unique_ptr<RooDataSet>;
//     using roodsets_t = std::vector<RooDataSet*>;
//     using filedtrees_t = std::vector<std::unique_ptr<HS::FIT::FiledTree>>;
    
//     class DataEvents  : public FitData {
      
//     public:
//       DataEvents(Setup &setup,TString tname,const strings_t& files);

//       DataEvents()=default;
//       DataEvents(const DataEvents&)=default;
//       DataEvents(DataEvents&&)=default;
//       ~DataEvents() override =default;
//       DataEvents& operator=(const DataEvents& other)=default;
//       DataEvents& operator=(DataEvents&& other) = default;

//       UInt_t GetN() const {return fFiledTrees.size();}
//       // RooAbsData& Get() final {return *(Get(0));}
//       dset_uptr Get(const UInt_t iset);
//       // RooDataSet* Get(UInt_t iset,Setup& setup) const;
//       TTree* GetTree(UInt_t ii){return fFiledTrees[ii]->Tree().get();}
//       strings_t FileNames() const {return fFileNames;}
//       TString FileName(UInt_t ii)const {return fFileNames[ii];}
//       TString ParentName() const {return fParentName;}
//       TString ParentTreeName() const {return fTreeName;}
//       void SetParentName(TString name) {fParentName=std::move(name);}
//       void SetParentTreeName(TString name) {fParentTreeName=std::move(name);}

//       void Load(Setup &setup,const TString& tname,const strings_t& files);
//       void LoadSetup(Setup *setup){fSetup=setup;}
//       void LoadBootStrap(const TString& tname,strings_t files);
//       void Reset(UInt_t ii) {fFiledTrees[ii].reset();}
      
//       void BootStrap(Int_t N){
// 	fNBoots=N;
// 	fBootStrap = std::unique_ptr<BootStrapper>(new BootStrapper{N});
//       }
//       void Toys(Int_t N){fNToys=N;}
      
//       Int_t GetDataBin(Int_t ii){
// 	if(fNBoots>0&&!fBootStrap.get())
// 	  BootStrap(fNBoots);//.recreate bootstrapper
// 	if(fBootStrap.get())
// 	  return fBootStrap->GetGroup(ii);
// 	else if(fNToys>0) //only toys if no bootstrap
// 	  return (int)std::round(ii/fNToys);
// 	//no boots or toys
// 	return ii;
//       }
//       TString GetItemName(Int_t ii);
//       void LoadWeights(TString wname,TString fname,TString wobj="HSsWeights");
//       const RooRealVar* WeightVar()const {return fWeightVar.get();}
      
//     protected:
//       void LoadWeights();

//     private:
 
//       HS::FIT::Setup *fSetup=nullptr;//!
//       strings_t fFileNames;
//       TString fTreeName;
//       TString fParentName;
//       TString fParentTreeName;
//       TString fInWeightName;
//       TString fInWeightFile;
//       TString fInWeightObjName;
      

//       filedtrees_t fFiledTrees;
//       std::shared_ptr<BootStrapper> fBootStrap;//!;

//       Int_t fNBoots=-1;
//       Int_t fNToys=-1;
      
//       weights_ptr fInWeights;//!
//       WeightsConfig fWgtsConf;

//       std::shared_ptr<RooRealVar> fWeightVar;//!
      
//       ClassDefOverride(HS::FIT::DataEvents,1);
//      };
    
//   }//namespace FIT
// }//namespace HS

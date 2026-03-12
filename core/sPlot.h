/**
 * @file sPlot.h
 * @brief Class for calculating and exporting sWeights for signal extraction.
 * @details Inherits from FitManager to perform the initial unbinned maximum 
 * likelihood fit before applying the RooStats::SPlot technique.
 */

#pragma once

#include "FitManager.h"
#include "Setup.h"
#include "Weights.h"
#include <RooStats/SPlot.h>
#include <memory>
#include <vector>

namespace HS {
namespace FIT {

    // Modern C++ Type Aliases
    using splot_uptr  = std::unique_ptr<RooStats::SPlot>;
    using weights_ptr = std::shared_ptr<Weights>;
    using weights_uptr= std::unique_ptr<Weights>;
    
    /**
     * @class sPlot
     * @brief High-level controller for calculating sWeights.
     * * The sPlot class extends the standard FitManager. It first performs a fit
     * to the data to determine the species yields. It then fixes the shape 
     * parameters of the PDFs, utilizes RooStats::SPlot to calculate per-event 
     * weights (sWeights) for each species, and exports these weights to a tree.
     */
    class sPlot : public FitManager {
      
    public:
      /** @brief Default constructor. */
        sPlot() = default;
        
        /** @brief Explicit copy constructor (bypasses unique_ptr). */
        sPlot(const sPlot& other); 
        
        sPlot(sPlot&&) = delete;
        ~sPlot() override = default;
        
        /** @brief Explicit copy assignment operator. */
        sPlot& operator=(const sPlot& other);
        
        sPlot& operator=(sPlot&& other) = delete;

      
        /**
         * @brief Executes the fit and calculates sWeights.
         * @return kTRUE if the fit and sPlot calculation succeeded, kFALSE otherwise.
         * @note This requires an unbinned dataset. Binned fits (RooDataHist) will fail.
         */
        Bool_t Run() override;
 
        /**
         * @brief Cleans up transient memory after a fit/sPlot execution.
         * Resets the underlying FitManager state as well as the local sPlot pointers.
         */
        void Reset() override {
            FitManager::Reset();
            fSPlot.reset();
            fWeights.reset();
            fWeightedFiledTree.reset();
        }

        /**
         * @brief Initializes the weights object and triggers the extraction of sWeights.
         */
        void CreateWeights();
        
        /**
         * @brief Loops over the dataset, extracts sWeights from RooStats::SPlot, 
         * and assigns them to a standard BruFit Weights object for exporting.
         */
        void ExportWeights();
        
        /**
         * @brief Merges weight files from multiple individual fits/bins.
         * @return A unique pointer to the newly merged Weights object.
         */
        weights_uptr MergeWeights();
        
        /**
         * @brief Quickly draws a variable from the internally tracked weighted tree.
         * @param var Name of the variable branch to draw.
         * @param wname Name of the weight species branch to apply.
         * @param cut Optional string cut to apply (default is "1").
         * @param opt Optional ROOT drawing options (e.g., "COLZ").
         */
        void DrawWeighted(const TString& var, const TString& wname, TString cut = "1", const TString& opt = "");

        /**
         * @brief Retrieves the active weighted TTree.
         * @return Raw pointer to the TTree, or nullptr if none is loaded.
         */
        TTree* GetWeightedTree() { 
            return fWeightedFiledTree ? fWeightedFiledTree->Tree().get() : nullptr; 
        }
        
        /** @brief Destroys the active weighted tree to free memory. */
        void DeleteWeightedTree() { fWeightedFiledTree.reset(); }
      
    protected:
        /**
         * @brief Clones the original data tree and appends the calculated sWeights.
         */
        void WeightedTree();
        
        /**
         * @brief Checks the post-fit yields. If any species collapsed to near-zero, 
         * it is safely removed from the model to prevent sPlot calculation errors.
         * @return kTRUE if only one valid species remains, kFALSE otherwise.
         */
        Bool_t ZeroYieldCheck();
      
    private:
        splot_uptr fSPlot;               ///<! Transient RooStats::SPlot calculation engine
        weights_ptr fWeights;            ///<! Transient weights container for exporting
        filed_shptr fWeightedFiledTree;  ///<! Transient open file managing the weighted TTree

        TString fSingleYield;            ///<! Caches the name of the yield if only one species exists
        std::vector<TString> fZeroYields;///<! Tracks species that collapsed to zero yield during the fit
      
        ClassDefOverride(HS::FIT::sPlot, 1);
    };
    
} // namespace FIT
} // namespace HS



// ////////////////////////////////////////////////////////////////
// ///
// ///Class:               sPlot
// ///Description:
// ///           

// #pragma once


// #include "FitManager.h"
// #include "Setup.h"
// #include "Weights.h"
// #include <RooStats/SPlot.h>
 
// namespace HS{
//   namespace FIT{

//     using splot_uptr = std::unique_ptr<RooStats::SPlot>;
//     using splot_shptr = std::shared_ptr<RooStats::SPlot>;
//     using weights_ptr =std::shared_ptr<Weights>;
//     using weights_uptr =std::unique_ptr<Weights>;
//     using tree_uptr =std::unique_ptr<TTree>;
//     using tree_shptr =std::shared_ptr<TTree>;
    
//     class sPlot  : public FitManager{
      
//     public:
//       sPlot()=default;
//       sPlot(const sPlot&)=default;
//       sPlot(sPlot&&)=delete;
//       ~sPlot() override =default;
//       sPlot& operator=(const sPlot& other) = default;
//       sPlot& operator=(sPlot&& other) = delete;
      
//       Bool_t Run() override;
 
//       void Reset() override{
// 	FitManager::Reset();
//         fSPlot.reset();
// 	fWeights.reset();
//        }

//       void CreateWeights();
//       void ExportWeights();
//       weights_uptr MergeWeights();
//       void DrawWeighted(const TString& var,const TString& wname,TString cut="1",const TString& opt="");

//       TTree* GetWeightedTree(){return fWeightedFiledTree->Tree().get();};
//       void DeleteWeightedTree(){fWeightedFiledTree.reset();};
      
//     protected:
      
//       void WeightedTree();
//       Bool_t ZeroYieldCheck();
      
//     private:
//       splot_shptr fSPlot; //!sPlot object
//       weights_ptr fWeights;//!
//       tree_shptr fWeightedTree;//!
//       filed_shptr fWeightedFiledTree;//!

//       TString fSingleYield;//!
//       std::vector<TString> fZeroYields;//!
      
//       ClassDefOverride(HS::FIT::sPlot,1);
//     };
    
//   }//namespace FIT
// }//namespace HS


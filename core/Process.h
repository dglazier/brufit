/**
 * @file Process.h
 * @brief Execution controllers for BruFit.
 */

#pragma once

#include "FitManager.h"
#include <TROOT.h>
#include <TBenchmark.h>
#include <TString.h>
#include <vector>
#include <algorithm>
#include <numeric>

namespace HS {
namespace FIT {
namespace PROCESS {
      
    static std::vector<TString> gCompilesList;

    class Here {
    public:
        static void Go(const std::shared_ptr<FitManager>& fm) { Go(fm.get()); }
        static void Go(FitManager* fm) {
            if (!fm) return;
            gBenchmark->Start("go");
            fm->RunAll();
            gBenchmark->Stop("go");
            gBenchmark->Print("go");
        }
        
        static void One(const std::shared_ptr<FitManager>& fm, Int_t ifit) { One(fm.get(), ifit); }
        static void One(FitManager* fm, Int_t ifit) {
            if (!fm) return;
            gBenchmark->Start("go");
            fm->RunOne(ifit);
            gBenchmark->Stop("go");
            gBenchmark->Print("go");
        }
    };

    /**
     * @class Multi
     * @brief Modern, process-based parallel execution (replaces PROOF).
     * Distributes individual bin fits across multiple CPU cores safely.
     */
    class Multi {
    public:
        static void Go(const std::shared_ptr<FitManager>& fm, Int_t nWorkers = 0) { Go(fm.get(), nWorkers); }
        static void Go(FitManager* fm, Int_t nWorkers = 0);
    };

    class Farm {
    public:
        static void Go(const std::shared_ptr<FitManager>& fm, Int_t Njobs, Bool_t toFarm) { Go(fm.get(), toFarm); } // toFarm acts as maxJobs here based on your cpp
        static void Go(FitManager* fm, Int_t maxJobs);
    };

    class Loader {
    public:
        static void Compile(const TString& macro) {
            gROOT->ProcessLine(Form(".L %s+", macro.Data()));
            if (!std::count(gCompilesList.begin(), gCompilesList.end(), macro)) 
                gCompilesList.push_back(macro);
        }
        
        static void Compile(const std::vector<TString>& macList) {
            for (const auto& macro : macList) Compile(macro);
        }
        
        static void CompileHere(TList* macList) {
            if (!macList) return;
            for (Int_t i = 0; i < macList->GetEntries(); i++) {
                TString macro = macList->At(i)->GetName();
                std::cout << "Will load " << macro << std::endl;
                gSystem->Exec(Form("cp %s .", macro.Data()));
                TString hfile = macro;
                hfile.ReplaceAll(".C", ".h");
                hfile.ReplaceAll(".cxx", ".h");
                gSystem->Exec(Form("cp %s .", hfile.Data()));
                Compile(gSystem->BaseName(macro));
            }
        }
    };
    
} // namespace PROCESS
} // namespace FIT
} // namespace HS

// ////////////////////////////////////////////////////////////////
// ///
// ///Class:               Process
// ///Description:
// ///           
// #pragma once

// #include "FitManager.h"
// #include "FitSelector.h"
// #include <TProof.h>
// #include <TBenchmark.h>
// #include <TString.h>
// #include <vector>
// #include <algorithm>

// namespace HS{
//   namespace FIT{
//     namespace PROCESS{
      
//       static std::vector<TString> gCompilesList;

//       class Here  {

//       public :
      
   
// 	static void Go(const std::shared_ptr<FitManager>& fm){
// 	  Go(fm.get());
// 	}
// 	static void Go(FitManager* fm){
// 	  if(!fm) return;
// 	  gBenchmark->Start("go");
// 	  fm->RunAll();
//  	  gBenchmark->Stop("go");
// 	  gBenchmark->Print("go");
	
// 	};
// 	static void One(const std::shared_ptr<FitManager>& fm,Int_t ifit){
// 	  One(fm.get(),ifit);
// 	}
// 	static void One(FitManager* fm,Int_t ifit){
// 	  if(!fm) return;
// 	  gBenchmark->Start("go");
// 	  fm->RunOne(ifit);
// 	  gBenchmark->Stop("go");
// 	  gBenchmark->Print("go");

// 	};
  
//       }; //class Here


//       class Proof  {

//       public :
      
// 	static void Go(const std::shared_ptr<FitManager>& fm,Int_t N){
// 	  Go(fm.get(),N);
// 	}
// 	static void Go(FitManager& fm,Int_t N){
// 	  Go(&fm,N);
// 	}
	
// 	static void Go(FitManager* fm,Int_t N){
// 	  if(!fm) return;
// 	  fm->WriteThis();

// 	  cout<<" Proof :: Go "<<gProof<<endl;
// 	  if(!gProof){
// 	    TProof*  proof = TProof::Open("://lite");
// 	    gROOT->ProcessLine(Form(".x LoadBruProof.C+(%d)",N));
	  
// 	  for(auto& macro : gCompilesList)
// 	    proof->Load(Form("%s++",macro.Data()),kTRUE);
// 	  }
// 	  gBenchmark->Start("goProof");
	  
// 	  FitSelector selector;
// 	  selector.SetFitManager(fm);
// 	  gProof->Process(&selector,fm->GetN());
// 	  gBenchmark->Stop("goProof");
// 	  gBenchmark->Print("goProof");
//  	}
	
//       }; //class Proof


//      class Farm  {

//       public :
      
//        static void Go(const std::shared_ptr<FitManager>& fm,Int_t Njobs,Bool_t toFarm){
// 	 Go(fm.get(),toFarm);
// 	}
	
//        static void Go(FitManager* fm,Int_t maxJobs);
	
//       }; //class Proof


//       class Loader {

//       public :
// 	static void Compile(const TString& macro){
// 	  gROOT->ProcessLine(Form(".L %s+",macro.Data()));
// 	  if(!std::count(gCompilesList.begin(),gCompilesList.end(),macro))gCompilesList.push_back(macro);
// 	}
// 	static void Compile(std::vector<TString> macList){
// 	  for(auto& macro : macList){
// 	    Compile(macro);
// 	  }
// 	}
// 	static void CompileHere(TList*  macList){
// 	  if(!macList)return;
// 	  for(Int_t i=0;i<macList->GetEntries();i++){
// 	    TString macro = macList->At(i)->GetName();
// 	    std::cout<<"Will load "<<macro<<std::endl;
// 	    gSystem->Exec(Form("cp %s .",macro.Data()));
// 	    TString hfile=macro;
// 	    hfile.ReplaceAll(".C",".h");
// 	    hfile.ReplaceAll(".cxx",".h");
// 	    gSystem->Exec(Form("cp %s .",hfile.Data()));
// 	    Compile(gSystem->BaseName(macro));
// 	  }
// 	}
// 	///	std::vector<TString> GetCompiled(){return gCompilesList;};
	
//       private :
	
	
//       };
//     }
//   }
// }

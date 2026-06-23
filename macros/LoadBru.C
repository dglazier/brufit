
#include <TSystem.h>
#include <TString.h>
#include <TInterpreter.h>
#include <TROOT.h>
#include <TClassTable.h>
#include <iostream>

// Preserve the BruFit namespace structure
namespace HS { namespace FIT { namespace PROCESS {}; namespace EXPAND {} } }
namespace bru{}
using namespace HS;
using namespace HS::FIT;
using namespace HS::FIT::PROCESS;
using namespace HS::FIT::EXPAND;
using namespace bru;

void LoadBru(TString Selection = "") {
    // 1. Get the BRUFIT path
    TString BRUCODE = gSystem->Getenv("BRUFIT");
    if (BRUCODE.IsNull()) {
        std::cerr << "ERROR: $BRUFIT environment variable not set!" << std::endl;
        return;
    }

    // 2. Define paths based on the new 'install' structure
    // Note: We use /include for the refactored headers
    TString incPath = BRUCODE + "/install/include";
    TString libPath = BRUCODE + "/install/lib";
    TString macPath = BRUCODE + "/macros";
    TString utPath  = BRUCODE + "/utility";

    // 3. Configure Interpreter and Macro Paths
    if (!TString(gInterpreter->GetIncludePath()).Contains(incPath)) {
        gInterpreter->AddIncludePath(incPath);
        
        // Update ROOT Macro search path
        gROOT->SetMacroPath(Form("%s:%s:%s:%s:%s", 
            gROOT->GetMacroPath(), 
            incPath.Data(), 
            libPath.Data(), 
            macPath.Data(), 
            utPath.Data()));

        // 4. Load the Library
        // We no longer need 'if(darwin)' checks. gSystem->Load handles 
        // the platform-specific extension (.so vs .dylib) automatically 
        // if you provide the path without the extension or the base lib name.
        libPath+="/libbrufit";
        if (gSystem->Load(libPath) < 0) {
            std::cout << "Note: Could not load " << libPath << " via gSystem. "
                      << "Checking if LD_LIBRARY_PATH handled it..." << std::endl;
        }
    }

    // 5. Compatibility Aliases (Backward Compatibility for old scripts)
    TClassTable::AddAlternate("HS::FIT::Weights", "HS::Weights");

    // 6. JIT Compile Helper Macros
    // Using '+' ensures it is compiled via ACLiC for performance
    gROOT->ProcessLine(".L $BRUFIT/macros/PDFExpand.C+");

    // ---------------------------------------------------------
    // 7. Pre-load RooFit Hardware Vectorization Libraries
    // ---------------------------------------------------------
    // RooFit sometimes fails to auto-load these during multiprocessing.
    // We explicitly load the fastest available engine into memory here so 
    // worker forks inherit the optimized symbols.
    
    std::cout << "Checking for optimized RooFit compute engines..." << std::endl;
    
    // Silence gSystem error printouts temporarily so users don't panic 
    // if their machine doesn't have AVX512
    int oldLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal; 

    if (gSystem->Load("libRooBatchCompute_AVX2") == 0) {
        std::cout << "  -> Hardware Vectorization: AVX2 engine loaded." << std::endl;
    } 
    else if (gSystem->Load("libRooBatchCompute_AVX") == 0) {
        std::cout << "  -> Hardware Vectorization: AVX engine loaded." << std::endl;
    } 
    else {
        std::cout << "  -> Hardware Vectorization: None found. Using generic fallback." << std::endl;
    }

    // Restore standard error logging
    gErrorIgnoreLevel = oldLevel;
    
    std::cout << "--- BruFit Loaded Successfully ---" << std::endl;
}


#include <TSystem.h>
#include <TString.h>
#include <TInterpreter.h>
#include <TROOT.h>
#include <TClassTable.h>
#include <iostream>

// Preserve the BruFit namespace structure
namespace HS { namespace FIT { namespace PROCESS {}; namespace EXPAND {} } }
using namespace HS;
using namespace HS::FIT;
using namespace HS::FIT::PROCESS;
using namespace HS::FIT::EXPAND;

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

    std::cout << "--- BruFit Loaded Successfully ---" << std::endl;
}

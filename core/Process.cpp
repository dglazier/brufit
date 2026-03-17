/**
 * @file Process.cpp
 * @brief Implementation of the execution controllers.
 */

#include "Process.h"
#include <TString.h>
#include <TSystem.h>
#include <ROOT/TProcessExecutor.hxx> 
#include <iostream>

// --- NEW HEADERS FOR PROGRESS BAR ---
#include <sys/mman.h> 
#include <atomic>      

namespace HS {
namespace FIT {
namespace PROCESS {
      
    // ========================================================================
    // Modern Multicore Processor (Replaces PROOF)
    // ========================================================================
    void Multi::Go(FitManager* fm, Int_t nWorkers) {
        if (!fm) return;
        
        Int_t nFits = fm->GetN();
        if (nFits == 0) return;
	
	ROOT::EnableThreadSafety();
 
        std::vector<Int_t> fitIndices(nFits);
        std::iota(fitIndices.begin(), fitIndices.end(), 0);

        // 1. ALLOCATE SHARED MEMORY
        // Because TProcessExecutor forks completely isolated processes, standard variables 
        // won't update across them. We use POSIX mmap to create a block of shared RAM 
        // to hold our thread-safe atomic counter.
        std::atomic<int>* completedFits = static_cast<std::atomic<int>*>(
            mmap(NULL, sizeof(std::atomic<int>), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0)
        );
        *completedFits = 0; // Initialize at 0

        ROOT::TProcessExecutor pool(nWorkers);

        std::cout << "\nProcess::Multi::Go starting with " << pool.GetPoolSize() 
                  << " workers for " << nFits << " isolated fits.\n" << std::endl;
                  
        gBenchmark->Start("goMulti");

        // Map distributes the lambda function across the worker processes.
        pool.Map([&](Int_t ifit) {
            
            // Reconstruct worker state safely
            FitManager workerFm(*fm);
	    // workerFm.LoadData(fm->GetDataTreeName(), fm->GetDataFileNames());
            workerFm.Data().LoadSetup(&workerFm.SetUp());
            
            // Redirect output so Minuit doesn't spam the console
            workerFm.SetRedirectOutput();
            
            // Execute the fit
            workerFm.RunOne(ifit);
            
            // 2. INCREMENT SHARED COUNTER
            int current = ++(*completedFits);
            
            // 3. TERMINAL BYPASS PROGRESS BAR
            // Since stdout is redirected to log files, we explicitly open the terminal 
            // device to draw our progress bar. (If running in a batch job without a screen, 
            // this simply fails gracefully and does nothing).
            FILE* tty = fopen("/dev/tty", "w");
            if (tty) {
                float progress = (float)current / nFits;
                int barWidth = 50;
                
                fprintf(tty, "\r["); // \r returns the carriage to the start of the line
                int pos = barWidth * progress;
                for (int i = 0; i < barWidth; ++i) {
                    if (i < pos) fprintf(tty, "=");
                    else if (i == pos) fprintf(tty, ">");
                    else fprintf(tty, " ");
                }
                fprintf(tty, "] %d/%d Fits Completed (%d%%)", current, nFits, int(progress * 100.0));
                fflush(tty);
                fclose(tty);
            }

            return 0; // Dummy return
        }, fitIndices);

        std::cout << "\n\n"; // Drop down a line so we don't overwrite the final progress bar
        
        gBenchmark->Stop("goMulti");
        gBenchmark->Print("goMulti");

        // 4. CLEAN UP SHARED MEMORY
        munmap(completedFits, sizeof(std::atomic<int>));
    }

    // ========================================================================
    // Batch Farm Processor
    // ========================================================================
    void Farm::Go(FitManager* fm, Int_t maxJobs) {
        if (!fm) return;
        fm->SetCompiledMacros(gCompilesList);
        fm->WriteThis();
        auto Njobs = fm->GetN();

        TString farmmac;
        if (gSystem->Getenv("HS_RUNMAC")) {
            std::cout << "Going to run macro " << gSystem->Getenv("HS_RUNMAC") << std::endl;
        } else {
            gSystem->Setenv("HS_RUNMAC", Form("%s/hsfit/HSFarmMac.C", gSystem->Getenv("HSCODE")));
            std::cout << "Going to run macro " << gSystem->Getenv("HS_RUNMAC") << std::endl;
        }

        TString farmrun = gSystem->Getenv("HS_FARMRUN") ? gSystem->Getenv("HS_FARMRUN") : "./pbs_run";
        TString farmsub = gSystem->Getenv("HS_FARMSUBMIT") ? gSystem->Getenv("HS_FARMSUBMIT") : "qsub";

        gSystem->Setenv("HS_LAUNCH", TString(gSystem->Getenv("PWD")));
        gSystem->Setenv("HS_OUTDIR", fm->SetUp().GetOutDir());
        
        for (Int_t i = 0; i < Njobs; i++) {
            TString JobNumber = Form("%d", i);
            std::cout << "sending JobNumber " << JobNumber << std::endl;
            gSystem->Setenv("HS_JOBNUMBER", JobNumber);
            gSystem->Setenv("HS_JOBNAME", fm->GetCurrName());
            
            if (maxJobs > 0) {
                TString njobsstring = gSystem->GetFromPipe("qstat | grep ${USER} | wc -l");
                Int_t njobs = njobsstring.Atoi();
                while (njobs > maxJobs) {
                    std::cout << "More than " << maxJobs << " running. Wait 10s." << std::endl;
                    gSystem->Sleep(10000);
                    njobsstring = gSystem->GetFromPipe("qstat | grep ${USER} | wc -l");
                    njobs = njobsstring.Atoi();
                }
                gSystem->Exec(farmsub + " " + farmrun);
            } else {
                gSystem->Exec(farmrun);
            }
        }
    }
      
} // namespace PROCESS
} // namespace FIT
} // namespace HS


// #include "Process.h"
// #include <TString.h>
// #include <TSystem.h>

// namespace HS{
//   namespace FIT{
//     namespace PROCESS{
      

//       ////////////////////////////////////////////////////////
//       ///Send jobs to Farm needs env variables
//       /// e.g. setenv HS_FARMRUN $PWD/pbs_run
//       /// e.g. setenv HS_FARMSUB qsub
//       /// e.g. setenv HS_RUNMAC FitMac.C  
//       void Farm::Go(FitManager* fm,Int_t maxJobs){

// 	if(!fm) return;
// 	fm->SetCompiledMacros(gCompilesList);
// 	fm->WriteThis();
// 	auto Njobs=fm->GetN();

// 	//Look for env variable RUNMAC for job macro 
// 	TString farmmac;
// 	if(gSystem->Getenv("HS_RUNMAC"))
// 	  std::cout<<"Going to run macro "<<gSystem->Getenv("HS_RUNMAC")<<endl;
// 	else{
// 	  gSystem->Setenv("HS_RUNMAC",Form("%s/hsfit/HSFarmMac.C",gSystem->Getenv("HSCODE")));
// 	  std::cout<<"Going to run macro "<<gSystem->Getenv("HS_RUNMAC")<<endl;
// 	}

// 	//Look for env variable FARMRUN for job submission script 
// 	TString farmrun;
// 	if(gSystem->Getenv("HS_FARMRUN"))
// 	  farmrun=gSystem->Getenv("HS_FARMRUN");
// 	else
// 	  farmrun="./pbs_run";

// 	//Look for variable FARMSUBMIT
// 	TString farmsub;
// 	if(gSystem->Getenv("HS_FARMSUBMIT"))
// 	  farmsub=gSystem->Getenv("HS_FARMSUBMIT");
// 	else
// 	  farmsub="qsub";

// 	//create a farm job for each toy requested    
// 	gSystem->Setenv("HS_LAUNCH",TString(gSystem->Getenv("PWD")));
// 	gSystem->Setenv("HS_OUTDIR",fm->SetUp().GetOutDir());
	
// 	for(Int_t i=0;i<Njobs;i++){
// 	  TString JobNumber=Form("%d",i);
// 	  cout<<"sending JobNumber "<<JobNumber<< endl;
// 	  gSystem->Setenv("HS_JOBNUMBER",JobNumber);
// 	  gSystem->Setenv("HS_JOBNAME",fm->GetCurrName());
// 	  if(maxJobs>0){
// 		TString njobsstring = gSystem->GetFromPipe("qstat | grep ${USER} | wc -l");
// 		Int_t njobs = njobsstring.Atoi();
// 		while(njobs>maxJobs){
// 			cout << "More than " << maxJobs << " running. Wait 10s." << endl;
// 			gSystem->Sleep(10000);
// 			njobsstring = gSystem->GetFromPipe("qstat | grep ${USER} | wc -l");
// 			njobs = njobsstring.Atoi();
// 		}
// 	    gSystem->Exec(farmsub+" "+farmrun);
// 	  }
// 	  else
// 	    gSystem->Exec(farmrun);

	  
//  	}
//       }//Go
      
//     }//PROCESS
//   }//FIT
// }//HS

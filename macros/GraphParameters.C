#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <RooFitResult.h>
#include <RooArgList.h>
#include <RooRealVar.h>
#include <TGraphErrors.h>
#include <iostream>
#include <memory>
#include "Bins.h"

// Usage in BruFit:
// .L GraphParameters.C
// GraphParameters("out/", "t", "ResultsBruMcmcCovariance.root", "BruMcmcCovariance")
void GraphParameters(TString DirName, TString Var, TString ResultFileName = "ResultsHSMinuit2.root", TString ResultObjName = "MinuitResult") {
    
    if (!DirName.EndsWith("/")) DirName += "/";

    std::unique_ptr<TFile> file(TFile::Open(DirName + "DataBinsConfig.root", "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open DataBinsConfig.root" << std::endl;
        return;
    }

    auto DataBins = dynamic_cast<HS::FIT::Bins*>(file->Get("HSBins"));
    if (!DataBins) {
        std::cerr << "Error: Could not find HSBins in file" << std::endl;
        return;
    }

    Int_t va = DataBins->GetAxisi(Var);
    TString AxisName = Var;

    TList* Graphs = new TList();
    Graphs->SetName("AllGraphs");

    for (Int_t ib = 0; ib < DataBins->GetN(); ib++) {
        TString redName = DataBins->GetBinName(ib);
        
        // Strip the variable we are plotting from the name so we can group the graphs
        // E.g., strips "t0.939960_" out of "Egamma8.400000_t0.939960_M2Pi0.712500_"
        Int_t idx = redName.Index(AxisName);
        if (idx != kNPOS) {
            TString remainder = redName(idx + AxisName.Sizeof() - 1, redName.Sizeof());
            Int_t firstUnderscore = remainder.First("_");
            TString axisBin = redName(idx, firstUnderscore + AxisName.Sizeof());
            redName.ReplaceAll(axisBin, ""); 
        }

        Int_t iP = DataBins->GetParti(va, DataBins->GetBinName(ib)); 
        
        // Open the specific minimizer result file
        TString resultPath = DirName + DataBins->GetBinName(ib) + "/" + ResultFileName;
        std::unique_ptr<TFile> fileR(TFile::Open(resultPath, "READ"));
        if (!fileR || fileR->IsZombie()) continue;

        // Fetch the specific minimizer result object
        RooFitResult* result = dynamic_cast<RooFitResult*>(fileR->Get(ResultObjName));
        if (!result) continue;

        const RooArgList& Pars = result->floatParsFinal();

        // Loop over parameters safely
        for (Int_t ipar = 0; ipar < Pars.getSize(); ipar++) {
            auto* parVar = dynamic_cast<RooRealVar*>(Pars.at(ipar));
            if (!parVar) continue; // Skip safely if not a RooRealVar

            // The graph name will now be exactly the prefix + parameter name
            // e.g., "Egamma8.400000_M2Pi0.712500_H_0_2_0"
            TString graphName = redName + parVar->GetName();
            TGraphErrors* graph = dynamic_cast<TGraphErrors*>(Graphs->FindObject(graphName));
            
            if (!graph) {
                graph = new TGraphErrors();
                graph->SetNameTitle(graphName, graphName);
                Graphs->Add(graph);
            }

            Int_t Npoint = graph->GetN();
            Double_t x_val = DataBins->GetAxis(va).GetBinCenter(iP + 1);
            Double_t x_err = DataBins->GetAxis(va).GetBinWidth(iP + 1) / 2.0;
            
            graph->SetPoint(Npoint, x_val, parVar->getVal());
            graph->SetPointError(Npoint, x_err, parVar->getError());
        }
    }

    TFile* fileG = TFile::Open(DirName + "ParGraphs" + Var + ".root", "RECREATE");
    Graphs->Write();
    fileG->Close();
    delete fileG; 

    std::cout << "Successfully extracted parameters to ParGraphs" << Var << ".root" << std::endl;
}

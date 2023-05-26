#define PlotAmbigs_cxx
// The class definition in PlotAmbigs.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("PlotAmbigs.C")
// root> T->Process("PlotAmbigs.C","some options")
// root> T->Process("PlotAmbigs.C+")
//


#include "PlotAmbigs.h"
#include <cmath>
#include <TH2.h>
#include <TStyle.h>

void PlotAmbigs::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void PlotAmbigs::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   std::vector<TString> names = {"S","D-","D0","D+"};
   std::vector<int> cols= {2,3,4,6};
   int ig=0;
   for(auto& gr:_graphMags){
     _multiMags->Add(&gr);
     gr.SetTitle(names[ig]);
     gr.SetMarkerStyle(20);
     gr.SetMarkerSize(1);
     gr.SetDrawOption("APE");
     gr.SetMarkerColor(cols[ig]);
     //gr1->SetLineWidth(4);
     gr.SetFillStyle(0);
     ++ig;
   }
   ig=0;
   for(auto& gr:_graphPhases){
     _multiPhases->Add(&gr);
     gr.SetTitle(names[ig]);
     gr.SetMarkerStyle(20);
     gr.SetMarkerSize(1);
     gr.SetDrawOption("APE");
     gr.SetMarkerColor(cols[ig]);
     //gr1->SetLineWidth(4);
     gr.SetFillStyle(0);
     ++ig;
 
   }

   cout<< std::setw(21)<<std::left<<"-2LogL"<<"\t"<<"|S|"<<"\t"<<"|S|err"<<"\t"<<"phi_S"<<"\t";
   cout<<"\t"<<"|D_-1|"<<"\t"<<"|D_-1|err"<<"\t"<<"phi_D-1"<<"\t"<<"phi_D-1err";
   cout<<"\t"<<"|D_0|"<<"\t"<<"|D_0|err"<<"\t"<<"phi_D0"<<"\t"<<"phi_D0err";
   cout<<"\t"<<"|D_1|"<<"\t"<<"|D_1|err"<<"\t"<<"phi_D1"<<"\t"<<"phi_D1err";
   cout<<"\t\t"<<"status";
   //cout<<"\n\n"<<"status";
   cout<<endl;
  
}

Bool_t PlotAmbigs::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);
 //   status = 1 : Covariance was made pos defined
// status = 2 : Hesse is invalid
// status = 3 : Edm is above max
// status = 4 : Reached call limit
// status = 5 : Covariance is not positive defined
   if(*NLL<0&&*fitstatus<2){
     auto a_2_1=sqrt(1-*a_0_0**a_0_0-*a_2_m1**a_2_m1-*a_2_0**a_2_0);

     
     auto d2xdy2 = [](double x){
		     return (x*x/(1-x*x));
		   };
     auto err_2_1 = sqrt(*a_0_0_err**a_0_0_err*d2xdy2(*a_0_0) + *a_2_0_err**a_2_0_err*d2xdy2(*a_2_0)  + *a_2_m1_err**a_2_m1_err*d2xdy2(*a_2_m1) );

     _graphMags[0].AddPoint(TMath::Abs(*a_0_0),*NLL);
     _graphMags[0].SetPointError(_ipoint,TMath::Abs(*a_0_0_err),0);
     
     _graphMags[1].AddPoint(TMath::Abs(*a_2_m1),*NLL);
     _graphMags[1].SetPointError(_ipoint,TMath::Abs(*a_2_m1_err),0);

    _graphMags[2].AddPoint(TMath::Abs(*a_2_0),*NLL);
    _graphMags[2].SetPointError(_ipoint,TMath::Abs(*a_2_0_err),0);

    _graphMags[3].AddPoint(sqrt(1-*a_0_0**a_0_0-*a_2_m1**a_2_m1-*a_2_0**a_2_0),*NLL);
    _graphMags[3].SetPointError(_ipoint,err_2_1,0);

    int phaseSign=1;
     if(a_2_1>0.1&&std::remainder(*aphi_2_1,2*TMath::Pi()) > 0) phaseSign=-1;
     else if(TMath::Abs(*a_2_m1)>0.1&&std::remainder(*aphi_2_1,2*TMath::Pi()) < 0)
     
     _graphPhases[0].AddPoint(phaseSign*std::remainder(*aphi_0_0,2*TMath::Pi()),*NLL);
   
     _graphPhases[1].AddPoint(phaseSign*std::remainder(*aphi_2_m1,2*TMath::Pi()),*NLL);
     if(*aphi_2_m1_err<1) _graphPhases[1].SetPointError(_ipoint,TMath::Abs(*aphi_2_m1_err),0);
 
     _graphPhases[2].AddPoint(phaseSign*std::remainder(*aphi_2_0,2*TMath::Pi()),*NLL);
     if(*aphi_2_0_err<1)_graphPhases[2].SetPointError(_ipoint,TMath::Abs(*aphi_2_0_err),0);
 
     _graphPhases[3].AddPoint(phaseSign*std::remainder(*aphi_2_1,2*TMath::Pi()),*NLL);
     if(*aphi_2_1_err<1)_graphPhases[3].SetPointError(_ipoint,TMath::Abs(*aphi_2_1_err),0);
 
 
     //text output
     cout<< std::setw(21)<<std::left<<*NLL<<"\t"<<TMath::Abs(*a_0_0)<<"\t"<<*a_0_0_err<<"\t"<<0.0;
     cout<<"\t"<<TMath::Abs(*a_2_m1)<<"\t"<<*a_2_m1_err<<"\t"<<phaseSign*std::remainder(*aphi_2_m1,2*TMath::Pi())<<"\t"<<*aphi_2_m1_err;
     cout<<"\t"<<TMath::Abs(*a_2_0)<<"\t"<<*a_2_0_err<<"\t"<<phaseSign*std::remainder(*aphi_2_0,2*TMath::Pi())<<"\t"<<*aphi_2_0_err;
     cout<<"\t"<<TMath::Abs(a_2_1)<<"\t"<<err_2_1<<"\t"<<phaseSign*std::remainder(*aphi_2_1,2*TMath::Pi())<<"\t"<<*aphi_2_1_err;
     cout<<"\t  "<<*fitstatus;
     cout<<endl;

 
     _ipoint++;
   }
   return kTRUE;
}

void PlotAmbigs::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void PlotAmbigs::Terminate()
{
  //soln 1
  std::complex<double> S(586.16,0);
  std::complex<double> D0(-662.74,64.22);
  std::complex<double> D1(106.75,-725.71);
  std::complex<double> Dm(227.74,62.56);
  //soln 2
  // std::complex<double> S(1098.19,0);
  // std::complex<double> D0(-172.68,228.61);
  // std::complex<double> D1(19.63,269.12);
  // std::complex<double> Dm(93.9 ,65.01);
  //soln 3
  // std::complex<double> S(1107.03,0);
  // std::complex<double> D0(-168.35,236.76);
  // std::complex<double> D1(46.67,-95.06);
  // std::complex<double> Dm(113.08,196.71);//????
  //soln 4
 // std::complex<double> S(584.95,0);
 //  std::complex<double> D0(-663.44 ,83.35);
 //  std::complex<double> D1(112.54,-68.24);
 //  std::complex<double> Dm(235.14,717);

  
  auto scaleNorm =TMath::Sqrt( std::abs(S)*std::abs(S)+std::abs(D0)*std::abs(D0)+std::abs(D1)*std::abs(D1)+	std::abs(Dm)*std::abs(Dm) );
  S/=scaleNorm;
  S/=1.00000000000000; //for rounding errors making a_2_2 -ve
  D0/=scaleNorm;
  D1/=scaleNorm;
  Dm/=scaleNorm;

  std::vector<int> cols= {2,3,4,6};

  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  auto can1 = new TCanvas();
  
  can1->Divide(2,1);
  //  fChain->Draw("NLL");

  can1->cd(1);
  _multiMags->SetTitle("Fit results for 10k events");
  _multiMags->GetXaxis()->SetLabelSize(0.05);
  _multiMags->GetXaxis()->SetTitleSize(0.05);
  _multiMags->GetXaxis()->SetTitle("Partial Wave Magnitude");
  _multiMags->GetXaxis()->CenterTitle();
  
  _multiMags->Draw("APE  plc");
  auto legend = gPad->BuildLegend();
  legend->SetTextSizePixels(25);

  auto nllmin = _multiMags->GetYaxis()->GetXmin();
  auto nllmax = _multiMags->GetYaxis()->GetXmax();

  TLine *Sline = new TLine(std::abs(S),nllmin,std::abs(S),nllmax);
  Sline->SetLineColor(cols[0]);
  Sline->SetLineStyle(2);  
  Sline->Draw("same");
  TLine *Dpline = new TLine(std::abs(D1),nllmin,std::abs(D1),nllmax);
  Dpline->SetLineColor(cols[3]);
  Dpline->SetLineStyle(2);  
  Dpline->Draw("same");
  TLine *D0line = new TLine(std::abs(D0),nllmin,std::abs(D0),nllmax);
  D0line->SetLineColor(cols[2]);
  D0line->SetLineStyle(2);  
  D0line->Draw("same");
  TLine *Dmline = new TLine(std::abs(Dm),nllmin,std::abs(Dm),nllmax);
  Dmline->SetLineColor(cols[1]);
  Dmline->SetLineStyle(2);  
  Dmline->Draw("same");
  
  can1->cd(2);
  _multiPhases->GetXaxis()->SetLabelSize(0.05);
  _multiPhases->GetXaxis()->SetTitleSize(0.05);
  _multiPhases->GetXaxis()->SetTitle("Partial Wave Phase");
  _multiPhases->GetXaxis()->CenterTitle();
  _multiPhases->Draw("APE  plc");
 // gPad->BuildLegend();
  TLine *SlineP = new TLine(std::arg(S),nllmin,std::arg(S),nllmax);
  SlineP->SetLineColor(cols[0]);
  SlineP->SetLineStyle(2);  
  SlineP->Draw("same");
  TLine *DplineP = new TLine(std::arg(D1),nllmin,std::arg(D1),nllmax);
  DplineP->SetLineColor(cols[3]);
  DplineP->SetLineStyle(2);  
  DplineP->Draw("same");
  TLine *D0lineP = new TLine(std::arg(D0),nllmin,std::arg(D0),nllmax);
  D0lineP->SetLineColor(cols[2]);
  D0lineP->SetLineStyle(2);  
  D0lineP->Draw("same");
  TLine *DmlineP = new TLine(std::arg(Dm),nllmin,std::arg(Dm),nllmax);
  DmlineP->SetLineColor(cols[1]);
  DmlineP->SetLineStyle(2);  
  DmlineP->Draw("same");
  
}

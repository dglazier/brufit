#!/usr/bin/env python3


import ROOT


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ${BRUFIT}/macros/LoadBru.C")

  # create the sPlot fit manager and set the output directory for fit results, plots, and weights
  RF = ROOT.sPlot()
  RF.SetUp().SetOutDir("out/")

  # define `Mmiss` as fit variable, this is also the name of the branch in the tree; set the fit range to [0, 9.5]
  RF.SetUp().LoadVariable("Mmiss[0, 9.5]")
  RF.SetUp().LoadAuxVar("Eg[0, 10]")
  # RF.SetUp().AddCut("Eg<3.2")  # cut on beam energy

  # define `fgID` as name of the event-ID variable
  # input tree should have a double branch with this name containing a unique event ID number
  RF.SetUp().SetIDBranchName("fgID")

  # define signal PDF `Signal` to be a histogram created from the `Mmiss` tree branch
  # the histogram shape is
  #   * convoluted with a Gaussian with width `smear_Sig`; initial width 0, allowed range [0, 20]
  #   * shifted by an offset `off_Sig`; initial offset 0, allowed range [-2, +2]
  #   * scaled along the x axis by a factor `scale_Sig`; initial scale 1, allowed range [0.8, 1.2]
  RF.SetUp().FactoryPDF("RooHSEventsHistPDF::Signal(Mmiss, smear_Sig[0, 0, 20], off_Sig[0, -2, 2], scale_Sig[1, 0.8, 1.2])")
  # load data from which histogram PDFs are constructed
  RF.LoadSimulated("MyModel", "SigData.root", "Signal")  # tree name, file name, PDF name
  # alternatively, define Gaussian signal PDF `Signal`
  # RF.SetUp().FactoryPDF("Gaussian::Signal(Mmiss, mean_Sig[6, 4, 7], width_Sig[0.2, 0.0001, 3])")
  RF.SetUp().LoadSpeciesPDF("Signal", 1)

  # define background PDF `BG` in the same way as above
  RF.SetUp().FactoryPDF("RooHSEventsHistPDF::BG(Mmiss, smear_Bkg[0, 0, 5], off_Bkg[0, 0, 0], scale_Bkg[1.0, 0.8, 1.2])")
  # load data from which histogram PDFs are constructed
  RF.LoadSimulated("MyModel", "BGData.root", "BG")
  # alternatively, define polynomial background PDF `BG`
  # RF.SetUp().FactoryPDF("Chebychev::BG(Mmiss, {a0_Bkg[-0.1, -1, 1], a1_Bkg[0.1, -1, 1]})")
  # alternatively, define Gaussian + polynomial background PDF `BG`
  # RF.SetUp().FactoryPDF("SUM::BG(r_Bkg[0, 0, 1] * Gaussian::BG1(Mmiss, mean_Bkg[5, 3, 7], width_Bkg[50, 0.0001, 100]), Chebychev::BG2(Mmiss, {a0_Bkg[-0.1, -1, 1], a1_Bkg[0.1, -1, 1]}))")
  RF.SetUp().LoadSpeciesPDF("BG", 1)

  # load data to be fitted
  RF.LoadData("MyModel", "Data.root")  # tree name, file name

  # perform fit and plot fit result
  # or try MCMC algorithm
  # mcmc = ROOT.BruMcmcCovariance(200, 100, 0.1, 0.23, 0.16, 0.3)
  # mcmc.TurnOffCovariance()  # BruMcmcCovariance only, do not proceed with covariance-based sampling, just perform basic stepping
  # RF.SetMinimiser(mcmc)
  ROOT.Here.Go(RF)

  # compare weighted data with true distributions
  filedTree = ROOT.FiledTree.Read("MyModel", "Data.root")
  trueTree = filedTree.Tree()
  MmissCut = "((0 <= Mmiss) && (Mmiss <= 9.5))"  # restrict distributions to the same Mmiss range that was used in the fit
  canv = ROOT.TCanvas("weightedPlots")
  canv.Divide(2, 2)
  # draw signal distributions
  canv.cd(1)
  RF.DrawWeighted("M1 >> hM1_Sig(100, 0, 10)", "Signal", MmissCut)  # data weighted with signal weights
  ROOT.gDirectory.Get("hM1_Sig").SetLineColor(ROOT.kRed + 1)
  trueTree.Draw("M1", f"(Sig == 1) && {MmissCut}", "SAME HIST")
  canv.cd(2)
  RF.DrawWeighted("M2 >> hM2_Sig(100, 0, 10)", "Signal", MmissCut)
  ROOT.gDirectory.Get("hM2_Sig").SetLineColor(ROOT.kRed + 1)
  trueTree.Draw("M2", f"(Sig == 1) && {MmissCut}", "SAME HIST")
  # draw background distributions
  canv.cd(3)
  RF.DrawWeighted("M1 >> hM1_Bkg(100, 0, 10)", "BG", MmissCut)  # data weighted with background weights
  ROOT.gDirectory.Get("hM1_Bkg").SetLineColor(ROOT.kRed + 1)
  trueTree.Draw("M1", f"(Sig == -1) && {MmissCut}", "SAME HIST")
  canv.cd(4)
  RF.DrawWeighted("M2 >> hM2_Bkg(100, 0, 10)", "BG", MmissCut)
  ROOT.gDirectory.Get("hM2_Bkg").SetLineColor(ROOT.kRed + 1)
  trueTree.Draw("M2", f"(Sig == -1) && {MmissCut}", "SAME HIST")
  canv.SaveAs(".pdf")

  # make sure weighted tree is written properly
  RF.DeleteWeightedTree()
  print("Finished successfully")

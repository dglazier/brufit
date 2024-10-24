set(BRUFIT brufit)
set(DIRNAME core)




ROOT_GENERATE_DICTIONARY(G__${BRUFIT} Weights.h FiledTree.h RooHSComplex.h RooHSComplexSumSqdTerm.h RooHSEventsPDF.h RooComponentsPDF.h RooHSEventsHistPDF.h BruEventsHistPeakPDF.h RooHSEventsHistPDF.h RooHSSphHarmonic.h RooHSDWigner.h RooHSDWignerProduct.h RelBreitWigner.h PdfParser.h PredefinedParsers.h ComponentsPdfParser.h Setup.h Binner.h Bins.h BootStrapper.h Data.h PlotResults.h MCMCPlotResults.h AutocorrPlot.h CornerPlot.h CornerFullPlot.h Minimiser.h FitManager.h sPlot.h ToyManager.h CrossSection.h RooMcmc.h AmpMinuit2.h AmpMcmc.h BruSequentialProposal.h  BruCovarianceProposal.h BruMetropolisHastings.h  BruMcmc.h HSMetropolisHastings.h Process.h FitSelector.h PhotoTwoSpin0Amps.h PhotoSpin1Spin0Amps.h  Utilities.h LINKDEF BruFitLinkDef.h)



add_library(${BRUFIT} SHARED  Weights.cpp FiledTree.cpp RooHSComplex.cpp RooHSComplexSumSqdTerm.cpp RooHSEventsPDF.cpp RooComponentsPDF.cpp  RooHSEventsHistPDF.cpp BruEventsHistPeakPDF.cpp RooHSSphHarmonic.cpp RooHSDWigner.cpp RooHSDWignerProduct.cpp RooHSEventsHistPDF.cpp RelBreitWigner.cpp PdfParser.cpp ComponentsPdfParser.cpp Setup.cpp Binner.cpp Bins.cpp  BootStrapper.cpp Data.cpp PlotResults.cpp MCMCPlotResults.cpp AutocorrPlot.cpp CornerPlot.cpp CornerFullPlot.cpp Minimiser.cpp FitManager.cpp  sPlot.cpp ToyManager.cpp CrossSection.cpp RooMcmc.cpp AmpMinuit2.cpp AmpMcmc.cpp BruSequentialProposal.cpp BruCovarianceProposal.cpp BruMetropolisHastings.cpp BruMcmc.cpp HSMetropolisHastings.cpp Process.cpp FitSelector.cpp PhotoTwoSpin0Amps.cpp PhotoSpin1Spin0Amps.cpp G__${BRUFIT}.cxx)




if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  target_link_libraries(${BRUFIT} ${ROOT_LIBRARIES} "-undefined dynamic_lookup")
else()
  target_link_libraries(${BRUFIT} ${ROOT_LIBRARIES} )
endif()


install(TARGETS ${BRUFIT}
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")

ROOT_GET_LIBRARY_OUTPUT_DIR(ROOT_LIBRARY_OUTPUT_DIR)
set(ROOT_BRUFIT_PCM ${ROOT_LIBRARY_OUTPUT_DIR}/${libprefix}${BRUFIT}_rdict.pcm)

install (FILES  ${ROOT_BRUFIT_PCM}
  DESTINATION "${CMAKE_INSTALL_LIBDIR}")

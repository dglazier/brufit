# Define the library target
set(BRUFIT_LIB brufit)

# --- 1. Dictionary Generation ---
# This macro creates G__brufit.cxx. 
# We list headers that need CINT/Cling dictionaries.
set(DICTIONARY_HEADERS
    Weights.h FiledTree.h RooHSComplex.h RooHSComplexSumSqdTerm.h 
    RooHSEventsPDF.h RooComponentsPDF.h RooHSEventsHistPDF.h 
    BruEventsHistPeakPDF.h RooHSSphHarmonic.h RooHSDWigner.h 
    RooHSDWignerProduct.h RelBreitWigner.h PdfParser.h 
    PredefinedParsers.h ComponentsPdfParser.h Setup.h 
    Binner.h Bins.h BootStrapper.h Data.h PlotResults.h 
    MCMCPlotResults.h AutocorrPlot.h CornerPlot.h 
    CornerFullPlot.h Minimiser.h FitManager.h sPlot.h 
    ToyManager.h CrossSection.h RooMcmc.h AmpMinuit2.h 
    AmpMcmc.h BruSequentialProposal.h BruCovarianceProposal.h 
    BruMetropolisHastings.h BruMcmc.h HSMetropolisHastings.h 
    Process.h FitSelector.h PhotoTwoSpin0Amps.h 
    PhotoSpin1Spin0Amps.h Utilities.h
    MCEventCache.h BruEventsPDF.h BruComponentsPDF.h BruEventsHistPDF.h
)

ROOT_GENERATE_DICTIONARY(G__${BRUFIT_LIB} ${DICTIONARY_HEADERS} LINKDEF BruFitLinkDef.h)

# --- 2. Target Definition ---
add_library(${BRUFIT_LIB} SHARED
    Weights.cpp FiledTree.cpp RooHSComplex.cpp RooHSComplexSumSqdTerm.cpp 
    RooHSEventsPDF.cpp RooComponentsPDF.cpp RooHSEventsHistPDF.cpp 
    BruEventsHistPeakPDF.cpp RooHSSphHarmonic.cpp RooHSDWigner.cpp 
    RooHSDWignerProduct.cpp RelBreitWigner.cpp PdfParser.cpp 
    ComponentsPdfParser.cpp Setup.cpp Binner.cpp Bins.cpp 
    BootStrapper.cpp Data.cpp PlotResults.cpp MCMCPlotResults.cpp 
    AutocorrPlot.cpp CornerPlot.cpp CornerFullPlot.cpp Minimiser.cpp 
    FitManager.cpp sPlot.cpp ToyManager.cpp CrossSection.cpp 
    RooMcmc.cpp AmpMinuit2.cpp AmpMcmc.cpp BruSequentialProposal.cpp 
    BruCovarianceProposal.cpp BruMetropolisHastings.cpp BruMcmc.cpp 
    HSMetropolisHastings.cpp Process.cpp FitSelector.cpp 
    PhotoTwoSpin0Amps.cpp PhotoSpin1Spin0Amps.cpp
    MCEventCache.cpp BruEventsPDF.cpp BruComponentsPDF.cpp BruEventsHistPDF.cpp
    G__${BRUFIT_LIB}.cxx
)

# --- 3. Target Properties (The "Modern" Way) ---
# This tells anything linking to brufit where to find the headers
target_include_directories(${BRUFIT_LIB} 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Link against ROOT targets
target_link_libraries(${BRUFIT_LIB} PUBLIC ROOT::RooFit ROOT::RooFitCore ROOT::RooStats ROOT::MathMore ROOT::Proof ROOT::Minuit2 ROOT::Tree ROOT::TreePlayer ROOT::Matrix ROOT::Physics)

# --- 4. Installation ---
install(TARGETS ${BRUFIT_LIB}
    EXPORT BruFitTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Grab ALL .h files in the core directory
file(GLOB ALL_HEADERS "*.h")

# Install all headers so others can develop against BruFit
install(FILES ${ALL_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install the ROOT Dictionary files (.pcm and .rootmap)
# These are generated in the library output directory
install(FILES 
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${BRUFIT_LIB}_rdict.pcm
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${BRUFIT_LIB}.rootmap
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

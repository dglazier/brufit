cmake_minimum_required(VERSION 3.16)
project(brufit VERSION 6.34.0 LANGUAGES CXX)

# --- 1. Output Directories ---
# Keep everything organized in the build folder
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# --- 2. C++ Standard ---
# BruFit R6.34 should likely use C++17 to match modern ROOT
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- 3. Dependencies ---
# Find ROOT and the specific components BruFit uses
find_package(ROOT REQUIRED COMPONENTS 
    Core 
    RIO 
    Hist 
    Graf 
    Gpad 
    Tree 
    RooFit 
    RooFitCore 
    RooStats 
    MathMore 
    Proof 
    Minuit2
)

# Use the ROOT_USE_FILE (legacy but helpful for dictionary macros)
include(${ROOT_USE_FILE})

# --- 4. Installation Setup ---
include(GNUInstallDirs)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # Optional: Set a default if BRUFIT_HOME is defined
    if(DEFINED ENV{BRUFIT})
        set(CMAKE_INSTALL_PREFIX "$ENV{BRUFIT}" CACHE PATH "..." FORCE)
    endif()
endif()

# --- 5. Subdirectories ---
add_subdirectory(core)

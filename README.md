# BruFit
## A RooFit based event based maximum likelihood fitting package 
## Installation

### 1. Prerequisites
Before building, ensure your environment is ready:
* **ROOT:** Sourced via `source /path/to/root/bin/thisroot.sh` (Bash) or `source /path/to/root/bin/thisroot.csh` (tcsh).
* **Compiler:** A C++17 compatible compiler (GCC 7+, Clang 5+, or Apple Clang).
* **CMake:** Version 3.16 or higher.

### 2. Standard Build Procedure
We use an out-of-source build. This keeps your source tree clean of temporary object files and CMake metadata, while installing the final libraries into a dedicated `install` folder.

```bash
# 1. Clone the repository and enter the directory
git clone [https://github.com/dglazier/brufit.git](https://github.com/dglazier/brufit.git)
cd brufit

# 2. Create and enter the build directory
mkdir build && cd build

# 3. Configure the project
# -DCMAKE_INSTALL_PREFIX defines where the final files go
cmake .. -DCMAKE_INSTALL_PREFIX=../install

# 4. Compile and Install
# -j$(nproc) uses all CPU cores for a faster build
cmake --build . -- -j$(nproc)
cmake --install .
```

> **Note on the `.pcm` file:** Depending on how the repository's `CMakeLists.txt` is configured, ROOT dictionary (`.pcm`) generation can sometimes lag behind the initial install step. Check if `../install/lib/libbrufit_rdict.pcm` exists. If it does not, simply run the last two commands (`cmake --build .` and `cmake --install .`) a second time.

### 3. Environment Setup & Alias
To ensure ROOT can find your new library and its dictionaries, and to enable your quick-launch alias, you need to set up your environment variables. Add the relevant block below to your `.bashrc` or `.tcshrc`.

**For Bash / Zsh:**
```bash
# Point to the top-level source directory where you cloned the repo
export BRUFIT=/path/to/brufit 

# Tell the system where the compiled libraries are
export LD_LIBRARY_PATH=$BRUFIT/install/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$BRUFIT/install/lib:$DYLD_LIBRARY_PATH # Required for macOS
export ROOT_INCLUDE_PATH=$BRUFIT/install/include:$ROOT_INCLUDE_PATH

# Create the alias to run BruFit via the macro
alias brufit='root $BRUFIT/macros/LoadBru.C'
```

**For tcsh / csh:**
```tcsh
# Point to the top-level source directory where you cloned the repo
setenv BRUFIT /path/to/brufit 

# Tell the system where the compiled libraries are
setenv LD_LIBRARY_PATH ${BRUFIT}/install/lib:$LD_LIBRARY_PATH
setenv DYLD_LIBRARY_PATH ${BRUFIT}/install/lib:$DYLD_LIBRARY_PATH # Required for macOS
setenv ROOT_INCLUDE_PATH ${BRUFIT}/install/include:$ROOT_INCLUDE_PATH

# Create the alias to run BruFit via the macro
alias brufit 'root $BRUFIT/macros/LoadBru.C'
```

### 4. Running BruFit
Once your environment variables are set and your terminal is refreshed (e.g., by running `source ~/.bashrc` or opening a new terminal), you can run the program exactly as you originally did:

```bash
brufit
```

The purpose of this package is to add to the RooFit package to allow
analysis of hadronic physics scattering reactions.

The main feature is a PDF class (RooHSEventsPDF) that allows calculation of normalisation
integrals from Monte Carlo detector simulations which allow the acceptance
of detector systems to be corrected for when extracting obervables.

In addition the (RooComponentsPDF) class provides caching of these integrals
for fast evaluation when the PDF is a sum of products

    	 A1(x_1)*B1(p_1) + A2(x_2)*B2(p_2) + ....

where x_i are measured variables and p_i are the parameters of interest.
For example

	C*cos(2*phi) +  D*sin(2*phi)
	
Has

	A1 = cos(2*phi); x_1 = phi; p_1 = C; A2=sin(2*phi); x_2=phi; p_2=D 


The idea is that more complex fits should not require more complex code
and there are components for splitting data (e.g into energy bins); running
similar fits in parallel via PROOF or a batch farm, which require minimal
extra code.

Weights can be used in the fits and can be created using the RooStats sPlot
class.

A Markov Chain Monte Carlo implementation based on Metropolis Hastings is
implemented and can provide robust (although not optimal) minimisation on
fits theat minuit may struggle to find a global minimum.



##Data

Data should be in the form of a ROOT TTree with branches that are usually double but can be int for categories, e.g. a polarisation state. If you are using weights and need an event ID branch this should also be made double so it can be read into RooFit dataset.

## Basic usage

   	 > brufit
	 root [1] FitManager fm
	 root [2] fm.SetUp().SetOutDir("out/"); //Put results files in out/
	 root [3] fm.SetUp().LoadVariable("phi[-3.1416,3.1416]"); //phi is a variable in the data tree
	 root [4] fm.SetUp().FactoryPDF("EXPR::amplitude('1+A[0,-1,1]*cos(2*phi)',phi,A)"); //Fit a cos2phi distribution
	 root [5] fm.SetUp().LoadSpeciesPDF("amplitude");//add to the total fit PDF 
	 root [6] fm.LoadData("treeName","fileName.root"); //set data (ROOT tree)
	 root [7] Here::Go(&fm); //run the fit

## Tutorials

### sPlotSimple
get the files

      cp -r $BRUFIT/tutorials/sPlotSimple .
      cd sPlotSimple


#### If running with Jupyter (Recomended)

start a notebook. Note the tutorials are written in python3 kernels.

      root --notebook  or jupyter-notebook

And open sPlotSimple.ipynb

You can also try the sPlotSimpleBins for an example of splitting the data into energy bins before making several fits.

#### If Running in ROOT interactive

First make some data

      root 'Model1.C( "Data.root" )'

Run

	brufit FitHSSimple.C

and

	brufit FitHSSimpleBins.C

### Performing fits to sinusiodal distributions with weights

      cp -r $BRUFIT/tutorials/WeightedObservable .
      cd WeightedObservable


#### If running with Jupyter (Recomended)

start a notebook. Note the tutorials are written in python3 kernels.

      root --notebook  or jupyter-notebook

And open sPlot.ipynb

Once you have found weights you can perform the weighted fit from FitWithEventsPDF.ipynb

You can also try using simulated data to give your sPlot Signal shape in sPlotWithSimulatedPDF.ipynb. And then try these weights in FitWithEventsPDF

Finally you can try splitting the Fit into Eg bins, running seperately on PROOF then plotting the result parameters as a function of Eg with FitWithComponentsPDFAndSplitBins.ipynb.

A faster more optimised method using RooComponentsPDF is given in FitWithComponentsPDF.ipynb


#### If Running in ROOT interactive

See the README in tutorials/WeightedObservable
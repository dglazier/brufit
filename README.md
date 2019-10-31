# BruFit
## A RooFit based event based maximum likelihood fitting package 

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


## Installation

### Prerequisites

ROOT with RooFit, Proof, Mathmore (if using Legendre polynomials)

### get and compile code

git clone https://github.com/dglazier/brufit.git

cd brufit

setenv BRUFIT /path/to/here (or setenv BRUFIT $PWD)

mkdir build
cd build
cmake ../
make install

alias brufit root $BRUFIT/macros/LoadBru.C

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
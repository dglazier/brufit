# Analysis of 2 Spin0 Photoproduction
-------------------------------------

The formalism of ["Moments of angular distribution and beam asymmetries in ηπ^0 photoproduction at GlueX"]( https://journals.aps.org/prd/abstract/10.1103/PhysRevD.100.054017) is implemented to analyse meson decays in linearly polarised photoproduction.

The intensity is expanded in polarised moments of spherical harmonic distributions. These moments are expressed in terms of partial wave amplitudes. Fits can be performed (or events generated) using either moments directly or (by use of RooFit Formula variables) the partial waves.

Data can be generated with given partial waves amplitudes and fit with the same or different set of partial waves. Fits can be done with Minuit and many fits with different starting values can be performed. Alternatively a MCMC algorithm can be run which gives similar, and in general more robust, results.

## Amigiuties in 2 Spin 0 Photoproduction
-----------------------------------------

This analysis was used as part of "Absence of Partial-Wave Ambiguities in Meson Photoproduction Experiments" (in preparation).

A given waveset was used to generate data, then fit many times using Minuit with different starting values. We find (as expected) only 1 true (highest likelihood) pair of solutions which are complex conjugates of each other (i.e. a trivial amiguity).

To perform this analysis you need to :

1. Generate a phase space dataset

   	    brufit GeneratePhaseSpace.C

2. Generate a data set with amplitudes

   	    brufit PhotoTwoSpin0Gen.C

3. Fit the data partial waves

       	   brufit PhotoTwoSpin0FitAmps.C

Further details can be found within these scripts for changing waves-sets and number of fits etc.

To run multiple fits you need to edit PhotoTwoSpin0FitAmps.C to use the AmpMinuit2(&config,20) minimiser, where 20 specifies the number of fits to perform.

You may then visualise the results in ROOT by running

    	root fitBruAmps/
	ResultsTreeBru->Process("PlotAmbigs.C")



## Moments Fits
-----------------------------------------
You may also prefer to just fit the moments. These should be
publishable observables themselves or may even be inverted to
give values for the amplitudes in this formalism.

3b.) Fit data moments

       	   brufit PhotoTwoSpin0FitMoments.C

To verify uncertainties bootstrapping algorithm can be applied.
In brufit this requires to iniate boostrapping with a number of
samples.

      Fitter.Data().BootStrap(10);

This must be done before the data is loaded!
IT generally makes sense to initialse the fit with the
parameters of your true fit to make sure you get the
same minima and speeden the fit.

4b.) Bootstrap moments fit

       	   brufit PhotoTwoSpin0FitMoments.C

This will produce Results file for each bootstrap fit.
Currently it is up to the user to decide what to do
with those.

e.g. hadd boot.root fitBruMomentsBoot/ResultsBoot*.root

Then take the mean and rms of each moments branch.
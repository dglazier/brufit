# Moments to Partial Waves
--------------------------

Note this tutorial may be used in conjunction with tutorials/PhotoAmps/TwoSpin0.

The concept for this technique is to take the equations for moments and partial waves given in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.100.054017
 Eqns (A9) and numerically invert them to get the partial wave amplitudes from the given moments. The inversion is done via a simple chi2 minimisation on the sum of the simultaneous equations using minuit2. A loop is performed over many randomised starting values of the partial waves to help search the full parameter space.

## Searching for ambigutites
----------------------------

One use case is given in the example RunGivenAmps.C . Where a known set of partial waves are specified from which the moments are directly calculated. The algorithm then inverts these to give the input partial waves that satisfy the simultansous equations. In case of ambigutities there may be more than 1 possible solution. The numerical values for the partial waves and the configuration of the intensity (Lmax, Mmax, reflectivities etc) are specified in the include file  ConfigureAmpsSPJune.C . Note (I) the amplitude naming convention a => +ve reflectivity, b => -ve. Then we have a_l_m specifing the reflectivity l and m. Note (II) the highest l,m +ve reflectivity wave is not a free parameter but is calculated from 1 - sum of squares of all other waves.

To run

      brufit Load.C RunGivenAmps.C


A tree will be produced which contains the numerical partial waves and moments for each minimsation.

      root resultsGivenAmpsNoH3.root
      PartialWaves->Draw("a_1_0:log_val","")

log_val is the log(chi2) for each minimisation. For the correct solution chi2->0.

Note the line to create the solver

      m2pw::EquationSolver solver{setup,0.0,{"H_3"}}; //ignore H_3


There a re a few options here, first is the possibility to give noise to the moments to simulate a more realistice experiment, rather than just a true inversion. For that you would change the value of the 0.0 to 0.05 etc. Second you may specify which polarised moments to IGNORE (by default all will be done), so in this example we will miss H3 which corresponds to circular polarised moments. Similarly we could check results for the unpolarised case by using  {"H_3","H_2","H_1"}.

## Converting Measured moments

In experiments it is usually more straightforward to measure moments than partial waves. Once this is done we may use this tool to convert the measured moments to measured partial waves.

There are 2 example here. The first allows you to type in by hand the values for each moment, having configured your intensity model (Lmax,Mmax, reflectivity etc). e.g.

       MomentHelper moments;
       moments.Set("H_0_0_0", 2.00000); //use JPAC normalisation !
       moments.Set("H_1_0_0", 0.321818);
       moments.Set("H_0_1_0", 0.0802487);
       ...

To Run

       brufit Load.C RunGivenMoments.C


Note(I) as these moments are taken from ConfigureAmpsSPJune.C the numerical partial waves found should be consistent with these. Note(II), in this case aphi_0_0 and bphi_0_0 are fixed to 0 at the top of the file.

Second if you have used brufit to extract moments from maxmimim likelhood fits you can just give the results file and it will load the given moments from there. Note you are respnsible for making sure the configuration used in this convertor is the same as when performing the fits (i.e.Lmax,Mmax, reflectivity etc) .

This example uses the moments fitting tutorial tutorials/PhotoAmps/TwoSpin0. The Results file from running PhotoTwoSpin0FitMoments.C is used. This has the following configuration,

        //Now set model options
  	//Lmax
  	config.SetLmax(2);
  	//Mmax = Lmax if not set
  	config.SetMmax(1);
  	//number of reflectivities = 1 or 2
  	config.SetNrefl(1);
 	//Only use even , S,D,... waves
  	config.SetOnlyEvenWaves();
 

Note you may want to test the convertor with 2 reflectivities or more waves, to see if it gives zero magnitudes!....

The true waveset was set as,

      //phases
      Generator.SetUp().SetParVal("aphi_0_0",0,kTRUE); //fix S real
      Generator.SetUp().SetParVal("aphi_2_-1",15.4*TMath::DegToRad(),kFALSE); //D-1
      Generator.SetUp().SetParVal("aphi_2_0",174*TMath::DegToRad(),kFALSE); //D0
      Generator.SetUp().SetParVal("aphi_2_1",-81.6*TMath::DegToRad(),kFALSE); //D+1
      //magnitudes
      Generator.SetUp().SetParVal("a_0_0",0.499,kFALSE); //S
      Generator.SetUp().SetParVal("a_2_-1",0.201,kFALSE); //D-1
      Generator.SetUp().SetParVal("a_2_0",0.567,kFALSE); //D0

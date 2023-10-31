# Example for fitting Time-Like Compton Scattering Distributions

The dilepton model is defined in Model.C and is used for both NoDetector
and Detector fits below.

    $$ I(\theta,\phi,-P)=\\&2(1+\cos^2(\theta))( B\frac{1}{\sin^2(\theta)} + T + A\frac{1}{\sin(\theta)}(ReM\cos(\phi)) $$
    
Where 	B = relative BH contribution
		T = relative TCS contribution
		A = relative INT contribution : B + T + A = 1 , so take A = 1 – B – T 
		ReM = real part of helicity ampltude
		ImM = imaginary part of helicity amplitude : ReM2  =(1- ImM2)

## No Detector

   Here we just generate a large "flat" dataset for all observables.
   We then use the "flat" dataset to generate sub data sets with given
   dilepton distribution parameters.
   Finally we fit these datasets to reproduce the parameters we put in.


Generate flat

      brufit GeneratePhaseSpace.C

Generate with model and fit

      brufit Model.C ToyMaker.C

You can edit ToyMaker to generate different sizes of data-sets (Model(toyman,10000); //create toy data with 100 events).
Or different values for our 3 parameters

      toyman.SetUp().SetParVal("BH",0.6); //Bethe Heitler contribution
      toyman.SetUp().SetParVal("TCS",0.1); //TCS contribution
      toyman.SetUp().SetParVal("ImM",0.7); //Imaginery part of M (ReM^2 =1-ImM^2)

Output in toys_test/


## Detector


For this you need to produce simulated data from some detector setup and reconstruction.
The events should be signal only and ideally truth matched.

It may help to adapt the PrepareTrees.C macro to refine your dataset so that is has :

      0<theta<pi
      -pi<phi<pi or you need to edit Model.C if you have 0<phi<2pi
      -1<pol<1

where theta and phi are the dilepton decay angles and Pol is the event circular polarisation
(=electron beam polarisation * polarisation transfer). Which should have a sign = helicity.
In addition you should have the truth values of all kinematics with the same name appended with some tag
e.g. trutheta, truphi, trupol. The tag can be changed in SimToyMaker.C.

Generate and fit the data

      brufit Model.C SimToyMaker.C

Again you can change the generated parameters in SimToyMaker.C

Or just fit data (this can be used as the basis for a real data fit, just change the file name).

      brufit Model.C FitThetaPhi.C

Note this is set to split the data into 4 t bins and fit them all.


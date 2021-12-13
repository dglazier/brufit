Here we are going to perform a fit to data using a ComponentsPDF with a given MC set of data
to calculate its normalisation integral. Then we are going to create Toy datasets based
on the resulting parameters and fit those to check the pull distributions for biases.
 model function = 1 + A*Pol*PolState*cos(2Phi) + B*Pol*PolState*sin(2Phi)


root  'Model1.C( "DataSignal.root",2 )' 

root 'Model1.C( "MC.root",0 )'


Now you can look at the FitToData.C script to see what is going on

brufit FitToData.C

Now you can look at the GenerateToys script to see what is going on

brufit GenerateToys.C

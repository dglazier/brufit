# sPlot with simple Gaussian and polynomial PDFs

You can try using the Jupyter Notebook

The fully run notebooks are saved as html files in this directory

      sPlotSimple.html

or to run with interactive ROOT session follow 

Generate a test data set using standard ROOT

       root 'Model1.C( "Data.root" )'

Perform a fit and create weights using the sPlot FitManager
Draw the signal weighted M1 distribution

       brufit FitHSSimple.C


Do the same but split the data into 4 seperate Egamma bins
and perform 4 fits. Use Proof to parallise the fits

       brufit FitHSSimpleBins.C



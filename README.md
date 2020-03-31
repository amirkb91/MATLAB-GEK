# MATLAB-GEK
Create Reduced Order Models of CFD using Gradient Enhance Kriging.
Test case is 2D NASA wall-mounted hump without a plenum: https://turbmodels.larc.nasa.gov/nasahump_val.html

Main code is in GEK_Modelling. Samples are stored inside folder and read by the main.m code.

-----------------------------------------------------------------------------------------------------------
Options for running the code:

% General options
options.platform    = platform to run on (iridis/local)
options.nfiles      = Number of files to read from samples folder
options.theta       = Theta file. If left blank found using GA
options.objective   = New sample "batch" or "verify" existing GEK prediction
options.npred       = Number of prediction points to be generated for MSE

% Options for next sample batch
options.batchnpool  = Number of pool points
options.nbatch      = Number of next sample batch points
options.batchmaxrad = Maximum exclusion radius 
options.batchtanh   = Tanh factor p. larger = more space b/w samples
options.batchxbound = New XY bounds to reduce window size
options.batchybound = Same as above
options.writebatch  = Write next sample batch to file
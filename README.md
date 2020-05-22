# MATLAB-GEK
Create Reduced Order Models of CFD using Gradient Enhance Kriging.
Test case is 2D NASA wall-mounted hump without a plenum: https://turbmodels.larc.nasa.gov/nasahump_val.html

Main code is in main.m.\
Raw contains raw data files.\
Samples contains SU2 input and output.\
Optimum_Theta contains stored theta files found using GA.\
src contains source files.

-----------------------------------------------------------------------------------------------------------
Options for running the code:

% General options\
options.platform    = platform to run on (iridis/local)\
options.nfiles      = Number of files to read from samples folder (0 = generate baseline)\
options.theta       = Theta file. If left blank found using GA\
options.objective   = New sample "batch" or "verify" existing GEK prediction\
options.npred       = Number of prediction points to be generated for MSE

% Global XY boundaries\
options.globalx     = Global X boundaries of the problem. Max if left blank\
options.globaly     = Global Y boundaries of the problem. Max if left blank

% Next sample batch\
options.nbatch      = Number of next sample batch points\
options.batchmaxrad = Maximum exclusion radius\ 
options.batchtanh   = Tanh factor p. larger = more space b/w samples\
options.batchxbound = X boundary of batch window. Equal to globalx if left blank\
options.batchybound = Y boundary of batch window. Equal to globaly if left blank\
options.writebatch  = Write next sample batch to file

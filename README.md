# MATLAB-GEK
Create Reduced Order Models of CFD using Gradient Enhance Kriging.
Test case is 2D NASA wall-mounted hump without a plenum: https://turbmodels.larc.nasa.gov/nasahump_val.html

Main code is in GEK_Modelling. Samples are stored inside folder and read by the main.m code.

-----------------------------------------------------------------------------------------------------------
Options for running the code:

platform  : Either 'local' or 'iridis'. Determines the number of nodes for parallel run. Code only runs in parellel.
nfiles    : Number of sample files to read from. All sample files located inside samples folder
theta     : Specify a theta file, inside the optimum theta folder, or leave blank for GA to calculate a new theta
objective : Either 'verify' or 'batch'. If batch code will do GEK prediction and find next batch of sample points. if verify code will verify GEK prediction by comparing velocity field to SU2 RANS solution with nominal SA values
npred     : Number of GEK prediction points
nbatch    : Number of points in next batch of samples
writebatch: If true, next batch is written to file and stored inside samples folder

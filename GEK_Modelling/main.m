%% Gradient Enhanced Kriging Modelling
% Model the product of flow velocity magnitude and angle at XY for given
% coefficients of the SA turbulence model
% Generate model from SU2 simulations

% Amir Bagheri
% Feb 2020

clear; close all; rng('shuffle');

% Change current folder to the folder of this m-file and add all to path.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end
addpath(genpath('../../SA_GEK_Hump/'));

%% Set Options for running the code

% General options
options.platform  = 'local'; % platform to run on (iridis/local)
options.nfiles    = 1; % Number of files to read from samples folder
options.theta     = 'theta01'; % theta file. If left blank found using GA
options.objective = 'batch'; % New sample "batch" or "verify" existing GEK prediction
options.npred     = 500; % number of prediction points to be generated for MSE

% Options for next sample batch
options.batchnpool  = 500; % number of pool points
options.nbatch      = 20; % number of next sample batch points
options.batchmaxrad = 0.1; % maximum exclusion radius 
options.batchtanh   = 2; % tanh factor p. larger = more space b/w samples
options.batchxbound = [0.0 1.5]; % new xy bounds to reduce window size
options.batchybound = [0.0 0.2];
options.writebatch  = false; % Write next sample batch to file

%% Run program

% Initialise the parallel run
init_parallel(options.platform);

% Read SU2 input and output from samples folder
[sample, param] = read_io(options);

% Calculate the hyperparameters theta of the Gaussian Correlation Function
tic;
[GEK.theta, GEK.ln_likelihood] = hyperparameters(sample, options.theta);
time.hyper = toc/60;

% Find Correlation matrix and Kriging mean using the hyperparameters
[GEK.R] = corrmat(sample, GEK.theta);
[GEK.mu, GEK.sighat] = kriging_mean(sample, GEK.R);

% Generate prediction points at which to predict GEK output
[pred] = predpoints(sample, param, options.objective, options.npred);

% Make GEK predictions
fprintf('\n----- Making Predictions -----\n');
tic;
[pred] = makeprediction(sample, pred, GEK);
time.prediction = toc/60;
fprintf('-Complete\n');

% Find next batch of sample points
if strcmp(options.objective, 'batch')
    fprintf('\n+++++ Generating Batch +++++\n');
    tic;
    [batch, pool, options] = nextbatch(sample, param, GEK, options);
    time.batch = toc/60;
    fprintf('-Complete\n');
else
    batch = []; pool = [];
end

% Generate plots
plotgek(sample, param, pred, batch, pool, options)

% If on IRIDIS save the workspace variables
if strcmp(options.platform, 'iridis')
    save('GEK_vars');
end

% Confirm success
fprintf('\n***** All Complete *****\n');

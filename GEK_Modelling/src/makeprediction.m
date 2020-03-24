function [pred] = makeprediction(sample, pred, GEK)
% Make GEK predictions on points to find output and MSE

fprintf('\n----- Making Predictions -----\n');

% Initialise arrays to store
predoutput = zeros(pred.npoint, 1);
predmse = zeros(pred.npoint, 1);

% de-struct for parfor loop
predpoint = pred.mapped;
theta = GEK.theta;
R = GEK.R;
mu = GEK.mu;
sighat = GEK.sighat;

% one zero matrices
o = ones(sample.npoint,1);
z = zeros(sample.npoint * sample.ndim,1);
one = [o;z];

parfor ii = 1:pred.npoint
    
    % Find correlation between prediction point and samples
    r = corrmat_pred(sample, theta, predpoint(ii,:));
    
    % Make prediction and find error
    y_p = mu + r'/R*(sample.output_aug-one*mu);
    mse = sighat*(1-r'/R*r+((1-one'/R*r)^2/(one'/R*one)));
    
    % De-normalise and store
    predoutput(ii) = y_p*sample.output_sd + sample.output_mean;
    predmse(ii) = mse/sighat;
    
end

% put results back into main pred struct
pred.output = predoutput;
pred.mse = predmse;

% Sort prediction mse from max to min
[pred.mse_sortval, pred.mse_sortindex] = sort(pred.mse,'descend');
fprintf('\n----- Predictions Complete -----\n');

end


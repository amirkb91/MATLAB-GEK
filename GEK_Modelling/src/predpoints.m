function [pred] = predpoints(sample, param, objective, npred)
% Obtain prediction points for GEK
% input options for prediction and sample
% output prediction struct

% Check if option objective is entered correctly
if ~strcmp(objective,'batch') && ~strcmp(objective,'verify')
    error('Wrong objective options specified');
end

% Create prediction points
% Need a unifrom distribution of points across the boundaries of each
% design parameter. Done using Halton sequences.

% Number of prediction points as specified by user
pred.npoint = npred;

% Halton sequence. Skip and Leap values defined here
skip = floor(rand*1e7);
% leap = nextprime(sample.ndim)*6; % integer multiple of next prime
leap = nthprime(sample.ndim+1)-1;
halton = haltonset(sample.ndim,'Skip',skip,'Leap',leap);
halton = scramble(halton,'RR2');

% Create prediction points
pred.raw = net(halton, pred.npoint);

% Map prediction points [0 1] to bounds of each parameter.
[pred.mapped] = map_samples(param, pred.raw);

if strcmp(objective,'batch')   
    % for batch objective, add original sample points to the prediction points
    % matrix. the GEK MSE at these points should be ~0    
    pred.mapped = vertcat(pred.mapped, sample.input);
    pred.npoint = pred.npoint + sample.npoint;

elseif strcmp(objective,'verify')
    % for verify objective, replace all SA values in prediction matrix with
    % nominal SA. Verify if GEK is producing same output for velocity objective
    % function than SU2 RANS at nomincal SA.
   
    % nominal SA values
    SAnom = [0.136,0.66667,0.622,0.41,0.3,2,7.1];
    SAnom = repmat(SAnom,pred.npoint,1);
    % populate prediction points
    pred.mapped = horzcat(SAnom,pred.mapped(:,param.x),pred.mapped(:,param.y));
end

end


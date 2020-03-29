function [pred] = predpoints(sample, param, objective, npred)
% Obtain prediction points for GEK
% input options for prediction and sample
% output prediction struct

% Check if option objective is entered correctly
if ~strcmp(objective,'batch') && ~strcmp(objective,'verify')
    error('Wrong objective options specified');
end

if strcmp(objective,'batch')
    % Need prediction points across all design parameters in order to
    % obtain MSE and judge location for next batch of samples.
    
    % Need a unifrom distribution of points across the boundaries of each
    % design parameter. Done using Halton sequences.
    
    % Number of prediction points (original sample points will also be added)
    pred.npoint = npred;
    
    % Halton sequence. Skip and Leap values chosen by user
    skip = floor(rand*1e7);
%     leap = nextprime(sample.ndim)*6; % integer multiple of next prime
    leap = nthprime(sample.ndim+1)-1;
    halton = haltonset(sample.ndim,'Skip',skip,'Leap',leap);
    halton = scramble(halton,'RR2');
    
    % Create prediction points
    pred.raw = net(halton, pred.npoint);
        
    % Map prediction points [0 1] to bounds of each parameter.
    [pred.mapped] = map_samples(param, pred.raw);
    
    % Add original sample points to the prediction points matrix. The GEK MSE
    % at these points should be ~0
    pred.mapped = vertcat(pred.mapped, sample.input);
    pred.npoint = pred.npoint + sample.npoint;
    
elseif strcmp(objective,'verify')
    
    % Verify if GEK is producing same output for velocity objective
    % function than SU2 RANS at nomincal SA. Need to create prediction
    % points in XY only, spaced equally in X and Y.
    
    % number of points in X-Y. Assume 5 times as many points in X than in Y.
    % 5*y^2 = npred
    len_y = ceil(sqrt(npred/5));
    len_x = len_y * 5;
    
    pred.npoint = len_x*len_y;
    % Create meshgrid
    boundary = get_boundary(param);
    xx=linspace(boundary(param.x,1),boundary(param.x,2),len_x)';
    yy=linspace(boundary(param.y,1),boundary(param.y,2),len_y)';
    [XX,YY]=meshgrid(xx,yy);
    % Change into column matrix
    XX = reshape(XX,pred.npoint,1);
    YY = reshape(YY,pred.npoint,1);
    
    % Remove points which fall inside hump
    hump_surface = load('hump_surface.mat');
    hump_surface = hump_surface.hump_surface;
    
    len = pred.npoint;
    i = 1;
    while i <= len
        if XX(i) <= 1 && XX(i) >= 0
            yhump = hump_surface(XX(i));
            if YY(i) < yhump
                % remove this point
                YY(i) = [];
                XX(i) = [];
                len = length(XX);
            end
        end
        i = i+1;
    end
    % update number of points to reflect those removed
    pred.npoint = len;
    % nominal SA values
    SAnom = [0.136,0.66667,0.622,0.41,0.3,2,7.1];
    SAnom = repmat(SAnom,pred.npoint,1);
    % populate prediction points
    pred.mapped = horzcat(SAnom,XX,YY);
end

end


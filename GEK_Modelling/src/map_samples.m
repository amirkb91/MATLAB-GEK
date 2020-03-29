function [mapped_samples] = map_samples(par,raw_samples, varargin)
%MAPTOBOUNDS Function to map samples between [0,1] to boundaries of the
%design parameters. For NASA Hump case.

% Input: par = struct with parameter names and assinged integer vals
%        raw_samples = samples between [0,1] obtained using quasi-random sets
% varargin: if required to change the xy mapping boundaries, add them here

% Output: mapped_samples = samples mapped linearly to the defined bounds

% number of parameters
npar = length(fieldnames(par));
% number of samples
nsam = length(raw_samples);
% hump surface
hump_surface = load('hump_surface.mat');
hump_surface = hump_surface.hump_surface;

% Get the boundaries of the parameters
boundary = get_boundary(par);

% over-ride x-y boundary if input is specified (used to map the batch samples)
if nargin == 4
    boundary(8,:) = varargin{1};
    boundary(9,:) = varargin{2};
end

%% Map the Samples to defined bounds

% Linear mapping of samples to limit bounds
% Map between a & b. b>a

mapped_samples = nan(size(raw_samples));
for j=1:npar
    if j ~= npar % not y, map based on fixed bounds
        a = boundary(j,1);
        b = boundary(j,2);
        for i=1:nsam
            mapped_samples(i,j) = raw_samples(i,j)*(b-a)+a;
        end
    else % y: lower bound changes based on hump coordinates
        x = mapped_samples(:,j-1);
        a = hump_surface(x);
        b = boundary(j,2);
        for i=1:nsam       
            mapped_samples(i,j) = raw_samples(i,j)*(b-a(i))+a(i);
        end
    end
end

end


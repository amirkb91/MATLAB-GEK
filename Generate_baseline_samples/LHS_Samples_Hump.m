%% Create LHS samples for SA & XY design space

% Use Latin Hypercube sampling to create a design space for the 7 SA
% coefficients and the 2 XY coefficients to run CFD on
% NASA HUMP CASE

clearvars; close all; writefile = false;


%% Create LHS Samples

nsam = 1000;     % number of samples
npar = 9;      % number of parameters
rng('shuffle');

samples.raw = lhsdesign(nsam,npar);

% halton.skip = 1e3;
% halton.leap = nextprime(npar)*6; % integer multiple of next prime
% halton.seq = haltonset(npar,'Skip',halton.skip,'Leap',halton.leap);
% halton.seq = scramble(halton.seq,'RR2');
% samples.raw = net(halton.seq, nsam);

%% Map LHS Samples onto boundaries

% Create parameter struct and set integers to each parameter
% this makes is easier to refer to each parameter later
param.cb1=1; param.sig=2; param.cb2=3; param.kar=4;
param.cw2=5; param.cw3=6; param.cv1=7; param.x=8; param.y=9;

% Find mapped samples and the boundaries
[samples.mapped] = map_samples_hump(param, samples.raw);
boundary = get_boundary(param);

%% Plot

f=figure(1); hold on;
% XY samples
plot(samples.mapped(:,end-1),samples.mapped(:,end),'*r');
axis equal;
% Geometry
load('hump_surface.mat'); % spline fitted to hump a priori. surface coords
xhump = linspace(-3,3,1000)';
yhump = hump_surface(xhump);
plot(xhump,yhump,'-','LineWidth',1.5,'Color',[0 0.4470 0.7410],'HandleVisibility','off')
plot(xhump,0.9*ones(length(xhump),1),'-','Color',[0 0.4470 0.7410],'LineWidth',1.5)
grid
% Sampling boundary
plot([boundary(param.x,1) boundary(param.x,1)],[0 boundary(param.y,2)],'--k','LineWidth',1.5,'HandleVisibility','off')
plot([boundary(param.x,2) boundary(param.x,2)],[0 boundary(param.y,2)],'--k','LineWidth',1.5,'HandleVisibility','off')
plot([boundary(param.x,1) boundary(param.x,2)],[boundary(param.y,2) boundary(param.y,2)],'--k','LineWidth',1.5)

title('Latin HyperCube XY sample space');
xlabel('x/c'); ylabel('y/c');
legend('Sample', 'Geometry boundary', 'Sampling area boundary');

%% Save samples to csv file
if writefile
    file = fopen('../samples01.dat','w');
    fprintf(file, '%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s \n', ...
        'cb1','sig','cb2','kar','cw2','cw3','cv1','X','Y');
    for i = 1:nsam
        fprintf(file, '%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f \n', ...
            samples.mapped(i,:));
    end
end




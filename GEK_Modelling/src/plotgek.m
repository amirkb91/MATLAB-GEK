function [] = plotgek(sample, param, pred, batch, objective, platform)
% Generate the plots

close all;

% Depending on objective, choose what to plot
if strcmp(objective, 'batch')
    plot_mse(sample, param, pred, batch, platform);
elseif strcmp(objective, 'verify')
    plot_vel(param, pred, platform);
end
end

%% MSE Plot
function [] = plot_mse(sample, param, pred, batch, platform)
% Plot MSE values of the prediction

% Get the boundaries of the design parameters for plotting
boundary = get_boundary(param);

fig = figure;
sgtitle(sprintf('Prediction MSE -- decluster P = %.2f',batch.p));

% MSE of prediction in X-Y space
subplot(3,1,1);
% interpolate the mse
interp = scatteredInterpolant(pred.mapped(:,param.x), ...
    pred.mapped(:,param.y),pred.mse);
% plot the mse
x = linspace(boundary(param.x,1),boundary(param.x,2),1000);
y = linspace(boundary(param.y,1),boundary(param.y,2),1000);
[X,Y] = meshgrid(x,y);
Z = interp(X,Y);
contourf(X,Y,Z,40,'LineColor','none')
axis equal; colorbar; hold on
xlabel('x/c'); ylabel('y/c')

% plot hump
hump_surface = load('hump_surface.mat');
hump_surface = hump_surface.hump_surface;
x = linspace(0,1,1000)';
y = hump_surface(x);
area(x,y,0,'FaceColor','w','HandleVisibility','off')

% plot next batch of sample points and the radius around them
% plot initial sample points
% plot GEK prediction points
plot(sample.input(:,param.x),sample.input(:,param.y),'sm');
plot(batch.point(:,param.x),batch.point(:,param.y),'*r')
% viscircles(batch.pointxy,batch.radius,'color','k','linewidth',1);
% plot(pred.mapped(:,param.x),pred.mapped(:,param.y),'.y');

%##########################################################################

% MSE of prediction in kar-cb1 space
subplot(3,1,2);
% interpolate the mse
interp = scatteredInterpolant(pred.mapped(:,param.kar), ...
    pred.mapped(:,param.cb1),pred.mse);
% plot the mse
x = linspace(boundary(param.kar,1),boundary(param.kar,2),1000);
y = linspace(boundary(param.cb1,1),boundary(param.cb1,2),1000);
[X,Y] = meshgrid(x,y);
Z = interp(X,Y);
contourf(X,Y,Z,40,'LineColor','none')
axis equal; colorbar; hold on
xlabel('kar'); ylabel('cb1')

% plot next batch of sample points
% plot initial sample points
% plot GEK prediction points
plot(sample.input(:,param.kar),sample.input(:,param.cb1),'sm');
plot(batch.point(:,param.kar),batch.point(:,param.cb1),'*r')
% plot(pred.mapped(:,param.kar),pred.mapped(:,param.cb1),'.y');

%##########################################################################

% MSE of prediction in sig-cw2 space
subplot(3,1,3);
% interpolate the mse
interp = scatteredInterpolant(pred.mapped(:,param.sig), ...
    pred.mapped(:,param.cw2),pred.mse);
% plot the mse
x = linspace(boundary(param.sig,1),boundary(param.sig,2),1000);
y = linspace(boundary(param.cw2,1),boundary(param.cw2,2),1000);
[X,Y] = meshgrid(x,y);
Z=interp(X,Y);
contourf(X,Y,Z,40,'LineColor','none')
axis equal; colorbar; hold on
xlabel('sig'); ylabel('cw2')

% plot next batch of sample points
% plot initial sample points
% plot GEK prediction points
plot(sample.input(:,param.sig),sample.input(:,param.cw2),'sm');
plot(batch.point(:,param.sig),batch.point(:,param.cw2),'*r')
% plot(pred.mapped(:,param.sig),pred.mapped(:,param.cw2),'.y');

%##########################################################################

% Save if on IRIDIS
if strcmp(platform, 'iridis')
    savefig(fig,'MSE.fig');
end

end

%% Velocity Plot
function [] = plot_vel(param, pred, platform)
% Plot value of prediction which is the velocity objective function

% Get the boundaries of the design parameters for plotting
boundary = get_boundary(param);

fig = figure;
sgtitle('Velocity Objective Function');

% output of prediction in X-Y space
subplot(3,1,1)
% interpolate the output
interp = scatteredInterpolant(pred.mapped(:,param.x), ...
    pred.mapped(:,param.y),pred.output);
% plot the output
x = linspace(boundary(param.x,1),boundary(param.x,2),1000);
y = linspace(boundary(param.y,1),boundary(param.y,2),1000);
[X,Y] = meshgrid(x,y);
Z = interp(X,Y);
contourf(X,Y,Z,30,'LineColor','none')
axis equal; colorbar; hold on; caxis([-1 1.1]);

title('GEK Prediction');
xlabel('x/c'); ylabel('y/c')

%##########################################################################

% RANS results for Nominal SA
subplot(3,1,2)

rans = load('rans.mat');
rans = rans.rans;
rans_velxinterp = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,3));
rans_velyinterp = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,4));

% Acquire velocities from interpolated function
rans_velx = rans_velxinterp(X,Y);
rans_vely = rans_velyinterp(X,Y);

% calculate objective function
rans_velmag = sqrt(rans_velx.^2 + rans_vely.^2);
rans_velang = atan2(rans_vely, rans_velx);
rans_obj = rans_velmag .* rans_velang;

contourf(X, Y, rans_obj,30,'LineColor','none')
axis equal; colorbar; hold on; caxis([-1 1.1]);

title('RANS SU2');
xlabel('x/c'); ylabel('y/c')

%##########################################################################

% Difference between RANS and GEK
subplot(3,1,3)

diff = abs(Z - rans_obj);

contourf(X,Y,diff,30,'LineColor','none')
axis equal; colorbar; hold on; caxis([0 1]);

title('|RANS - GEK|');
xlabel('x/c'); ylabel('y/c')

%##########################################################################

% plot hump
hump_surface = load('hump_surface.mat');
hump_surface = hump_surface.hump_surface;
x = linspace(0,1,1000)';
y = hump_surface(x);
for i = 1:3
    subplot(3,1,i)
    area(x,y,0,'FaceColor','w','HandleVisibility','off')
end

%##########################################################################

% Save if on IRIDIS
if strcmp(platform, 'iridis')
    savefig(fig,'Velocity.fig');
end
end
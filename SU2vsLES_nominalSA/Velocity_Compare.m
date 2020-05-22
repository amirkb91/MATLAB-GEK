%% Code to compare the velocity field between nominal SA SU2 and LES and Experiment
% Read results from the nominal steady SU2 run and compare
% Since meshes are different we need to interpolate

clear; close all
warning('off','all')

%% Read results and do interpolation
les = load('les.mat'); les = les.les;
rans = load('rans.mat'); rans = rans.rans;

% Do interpolation since meshes are different and we can't compare node to node
rans_velxint = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,3), 'linear', 'nearest');
rans_velyint = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,4), 'linear', 'nearest');
les_velxint = scatteredInterpolant(les(:,1),les(:,2),les(:,3), 'linear', 'nearest');
les_velyint = scatteredInterpolant(les(:,1),les(:,2),les(:,4), 'linear', 'nearest');

%% Create meshgrid for plotting contours
% Get hump surface
hump_surface = load('hump_surface.mat');
hump_surface = hump_surface.hump_surface;

% Number of points to plot contour
xpoint = 1000;
ypoint = 1000;

% Global limits
xbound = [0.7, 1.5];
ybound = [0.0, 0.1];
% xbound = [-1, max(les(:,1))];
% ybound = [0, max(les(:,2))];

x = linspace(xbound(1),xbound(2),xpoint)';
y = linspace(ybound(1),ybound(2),ypoint)';
[X,Y] = meshgrid(x,y);

% Restack meshgrid to remove y points inside hump
for i=1:xpoint
   if X(1,i) > 0 && X(1,i) < 1
       Y(:,i) = linspace(hump_surface(X(1,i)),ybound(2),ypoint)';
   end
end

%% Acquire velocities from interpolated function
rans_velx = rans_velxint(X,Y);
rans_vely = rans_velyint(X,Y);
les_velx  = les_velxint(X,Y);
les_vely  = les_velyint(X,Y);

% Find magnitude and angle of velocity vector
rans_velmag = sqrt(rans_velx.^2 + rans_vely.^2);
rans_velang = atan2(rans_vely, rans_velx);
les_velmag  = sqrt(les_velx.^2 + les_vely.^2);
les_velang  = atan2(les_vely, les_velx);

% Find objective function value
rans_objfunc = rans_velmag.*rans_velang;
les_objfunc  = les_velmag.*les_velang;

%% Find difference
diff_velmag  = abs(rans_velmag-les_velmag);
diff_velang  = abs(rans_velang-les_velang);
diff_objfunc = abs(rans_objfunc-les_objfunc);

%% Plot of difference
fig = figure;

% difference threshold
thresh = 0;

p{1}=subplot(3,1,1);
levels = linspace(thresh,max(diff_velmag,[],'all'),40);
contourf(X,Y,diff_velmag,levels,'LineColor','none')
axis equal; hold on;
title('Difference in Velocity Magnitude - RANS & LES')
colorbar
colormap('jet')

p{2}=subplot(3,1,2);
levels = linspace(thresh,max(diff_velang,[],'all'),40);
contourf(X,Y,diff_velang,levels,'LineColor','none')
axis equal; hold on;
title('Difference in Velocity Angle - RANS & LES')
colorbar
colormap('jet')

p{3}=subplot(3,1,3);
levels = linspace(thresh,max(diff_objfunc,[],'all'),40);
contourf(X,Y,diff_objfunc,levels,'LineColor','none')
axis equal; hold on;
title('Difference in Velocity ObjFunc - RANS & LES')
colorbar
colormap('jet')

%% Plot hump on all figures and set axis labels

xhump = linspace(xbound(1),xbound(2),1000)';
yhump = hump_surface(xhump);

for i=1:length(p)
    area(p{i},xhump,yhump,0,'FaceColor','none')
    xlabel('x/c');
    ylabel('y/c');
end

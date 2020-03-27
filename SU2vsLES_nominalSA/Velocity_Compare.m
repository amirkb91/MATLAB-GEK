%% Code to compare the velocity field between nominal SA SU2 and LES and Experiment
% Read results from the nominal steady SU2 run and compare
% Since meshes are different we need to interpolate

clear; close all
warning('off','all')

%% Read from results folder and do interpolation
% header is X/Y/velx/vely. Velocities are nondimensional wrt freestream
% rans = dlmread('/home/akb1r19/Documents/Results_CFD/NASA Hump/nominal_steady/SU2/Velocity_Field.dat','',10,0);
% les  = dlmread('/home/akb1r19/Documents/Results_CFD/NASA Hump/nominal_steady/LES/Velocity_Field.dat','',10,0);

les = load('les.mat'); les = les.les;
rans = load('rans.mat'); rans = rans.rans;

% Do interpolation since meshes are different and we can't compare node to node
rans_velx = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,3));
rans_vely = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,4));
les_velx = scatteredInterpolant(les(:,1),les(:,2),les(:,3));
les_vely = scatteredInterpolant(les(:,1),les(:,2),les(:,4));

%% Create meshgrid for plotting contours
% Number of points to plot contour
xpoint = 1000;
ypoint = 1000;
% Contour limits based on LES data (smallest mesh)
x = linspace(-1,max(les(:,1)),xpoint)';
y = linspace(0,max(les(:,2)),ypoint)';
[X,Y] = meshgrid(x,y);

%% Acquire velocities from interpolated function
rans_velx_int = rans_velx(X,Y); rans_velx_int(isnan(rans_velx_int))=0;
rans_vely_int = rans_vely(X,Y); rans_vely_int(isnan(rans_vely_int))=0;
les_velx_int  = les_velx(X,Y); les_velx_int(isnan(les_velx_int))=0;
les_vely_int  = les_vely(X,Y); les_vely_int(isnan(les_vely_int))=0;

% Find magnitude and angle of velocity vector
rans_velmag = sqrt(rans_velx_int.^2 + rans_vely_int.^2);
rans_velang = atan2(rans_vely_int, rans_velx_int);
les_velmag  = sqrt(les_velx_int.^2 + les_vely_int.^2);
les_velang  = atan2(les_vely_int, les_velx_int);

%% Find difference
diff_velmag = zeros(ypoint,xpoint);
diff_velang = zeros(ypoint,xpoint);

for i = 1:ypoint
    for j = 1:xpoint
        diff_velmag(i,j) = abs(rans_velmag(i,j)-les_velmag(i,j));
        diff_velang(i,j) = abs(rans_velang(i,j)-les_velang(i,j));
    end
end
%% Plots of SU2 and LES.
fig = cell(1);

fig{1} = figure;
contourf(X,Y,rans_velmag,40,'LineColor','none')
axis equal
title('SU2 RANS Velocity Magnitude')
colorbar

fig{2} = figure;
contourf(X,Y,les_velmag,40,'LineColor','none')
axis equal
title('LES Velocity Magnitude')
colorbar

fig{3} = figure;
contourf(X,Y,rans_velang,40,'LineColor','none')
axis equal
title('SU2 RANS Velocity Angle')
colorbar

fig{4} = figure;
contourf(X,Y,les_velang,40,'LineColor','none')
axis equal
title('LES Velocity Angle')
colorbar

% quiver to see velocity if needed
% quiver(X,Y,rans_velx_int,rans_vely_int)
%% Plot of difference
fig{5} = figure;
contourf(X,Y,diff_velmag,linspace(0,0.2,40),'LineColor','none')
axis equal
title('Difference in Velocity Magnitude - RANS & LES')
colorbar
colormap('jet')

fig{6} = figure;
contourf(X,Y,diff_velang,linspace(0,0.2,40),'LineColor','none')
axis equal
title('Difference in Velocity Angle - RANS & LES')
colorbar
colormap('jet')
%% Plot hump on all figures and set axis labels
hump_surface = load('hump_surface.mat');
hump_surface = hump_surface.hump_surface;
xhump = linspace(0,1,1000)';
yhump = hump_surface(xhump);

for i=1:length(fig)
    figure(fig{i})
    hold on
    area(xhump,yhump,0,'FaceColor','w')
    xlabel('x/c');
    ylabel('y/c');
end

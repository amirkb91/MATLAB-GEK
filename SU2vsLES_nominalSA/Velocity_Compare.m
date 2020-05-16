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
rans_velxint = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,3), 'linear', 'linear');
rans_velyint = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,4), 'linear', 'linear');
les_velxint = scatteredInterpolant(les(:,1),les(:,2),les(:,3), 'linear', 'linear');
les_velyint = scatteredInterpolant(les(:,1),les(:,2),les(:,4), 'linear', 'linear');

%% Create meshgrid for plotting contours

% Get hump surface
hump_surface = load('hump_surface.mat');
hump_surface = hump_surface.hump_surface;

% Number of points to plot contour
xpoint = 1000;
ypoint = 1000;
% Contour limits based on LES data (smallest mesh)
x = linspace(-1,max(les(:,1)),xpoint)';
y = linspace(0,max(les(:,2)),ypoint)';
[X,Y] = meshgrid(x,y);
% Restack meshgrid to remove y points inside hump
for i=1:xpoint
   if X(1,i) > 0 && X(1,i) < 1
       Y(:,i) = linspace(hump_surface(X(1,i)),max(les(:,2)),ypoint)';
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

%% Find difference
diff_velmag = zeros(ypoint,xpoint);
diff_velang = zeros(ypoint,xpoint);
diff_obj    = zeros(ypoint,xpoint);

for i = 1:ypoint
    for j = 1:xpoint
        diff_velmag(i,j) = abs(rans_velmag(i,j)-les_velmag(i,j));
        diff_velang(i,j) = abs(rans_velang(i,j)-les_velang(i,j));
    end
end

for i = 1:ypoint
    for j = 1:xpoint
        diff_obj(i,j) = abs(rans_velmag(i,j)*rans_velang(i,j)-les_velmag(i,j)*les_velang(i,j));
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

fig{7} = figure;
contourf(X,Y,diff_obj,30,'LineColor','none')
axis equal
title('Difference in Velocity ObjFunc - RANS & LES')
colorbar
colormap('jet')
xlim([0.75 1.15])
ylim([0 0.1])
%% Plot hump on all figures and set axis labels

xhump = linspace(0,1,1000)';
yhump = hump_surface(xhump);

for i=1:length(fig)
    figure(fig{i})
    hold on
    area(xhump,yhump,0,'FaceColor','none')
    xlabel('x/c');
    ylabel('y/c');
end

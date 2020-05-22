% Get the sensitivity of the objective function wrt x and y from flow.dat
% Done using central finite difference
% Replace the values inside the results files

clear;
close all;

%% Read results file
% This is generated using the python code Extract_Adjoint_results.py
resfile = 'results01.dat'; 
pyres = dlmread(resfile,'',1,0);
nfiles = size(pyres,1);

%% Get list of all folders 
dinfo = dir();
dinfo(ismember( {dinfo.name}, {'.', '..'})) = [];
dind = [dinfo.isdir]; % index
folders = dinfo(dind);

%% Create meshgrid to make interpolation of flow
% Number of points
xpoint = 1500; ypoint = 1500;

% X-Y limits
x = linspace(-1.2,2.2,xpoint)';
y = linspace(0,0.6,ypoint)';
[X,Y] = meshgrid(x,y);

%% initiate array to store new sensitivties in and set central diff step size
sens_centdiff = NaN(nfiles,2);
h = 1e-5;

%% Loop through all flow.dat files and find the gradients
fprintf('||||||||||\n');
for i = 1:nfiles

    % flow.dat file
    flowfile=strcat(folders(i).folder, '/', folders(i).name, '/flow.dat');

    % Point at which to find gradients
    point = pyres(i,1:2);

    % Read flow.dat
    rans = dlmread(flowfile,'',[3,0,44583,16]);

    % Calculate objective function at each mesh node and create interpolant
    % find u and v from flow.dat
    u = rans(:,4)./rans(:,3);
    v = rans(:,5)./rans(:,3);

    % calculate objective function and create interpolant
    objfun = sqrt(u.^2 + v.^2).*atan2(v,u);
    objfun_int = scatteredInterpolant(rans(:,1),rans(:,2),objfun, 'natural', 'none');

    % Interpolate at meshgrid points
    Z = objfun_int(X,Y);

    % Find gradients at given location using central difference
    xplus  = objfun_int(point(1)+h, point(2));
    xminus = objfun_int(point(1)-h, point(2));
    yplus  = objfun_int(point(1), point(2)+h);
    yminus = objfun_int(point(1), point(2)-h);
    if yminus < 0
        yminus = 0;
    end

    sens_centdiff(i,1) = (xplus-xminus)/(2*h);
    sens_centdiff(i,2) = (yplus-yminus)/(2*h);

    clearvars flowfile rans u v objfun objfun_int

    % display progress
    if mod(i,nfiles/10) == 0
        fprintf('|');
    end

end
fprintf('\n');

%% Write new results file with sensitivity values
% Get name of new file
name = split(resfile,'.');
newname = strcat(name{1},'_mat.',name{2});
file = fopen(newname,'w');

% Create new results matrix with new sens replacing old sens
newmat = horzcat(pyres(:,1:end-2),sens_centdiff);

% Write to file
fprintf(file, 'VARIABLES = "X", "Y", "obj_func", "sens_cb1", sens_sig", "sens_cb2","sens_kar", "sens_cw2", sens_cw3", "sens_cv1", "sens_x", "sens_y"');
fprintf(file, '\n');
for i = 1:nfiles
    fprintf(file, '%-.6E',newmat(i,:));
    fprintf(file, '\n');
end

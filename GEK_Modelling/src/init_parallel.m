function [] = init_parallel(platform)
% Set number of nodes to initialise parallel run

if strcmp(platform,'local')
    numcpu = 4;
elseif strcmp(platform,'iridis')
    numcpu = 40;
else
    error('Invalid platform name');
end

p = gcp('nocreate');
if isempty(p)
    parpool(numcpu);
end

end

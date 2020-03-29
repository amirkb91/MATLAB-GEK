function [batch] = nextbatch(sample, param, GEK, options)
% Find next batch of samples in adaptive sampling and decluster in XY
% input: sample param structs, number of batches and write to file
% output: next batch of samples

% Create points to do prediction on and determine the points with
% highest MSE. These become the samples for the next SU2 iteration.
% The highest mse value will automatically become the first sample point in
% the next SU2 iteration. After that, a radius will be assinged in x and y
% to that point, which will depend on the value of the adjoint gradient at
% that point (interpolated). If the sumabs of the adjoint gradient
% of the 7 SA parameters is low, the radius will be set to a larger value.
% The search for the next high mse point will be done outside of that
% radius. There is no need to cluster that radius with other points with
% similar XY values and different SA values, since the output is not that
% sensitive to SA anyway.

%% ****** THESE PARAMETERS CONTROL THE CLUSTERING *******
% Set the maximum possible radius allowed
batch.maxrad = options.batchmaxrad;
% Set the tanh multiplication factor. larger p = more space b/w points
batch.p = options.batchtanh;
% Set the new x-y boundaries if you want to localise the batch points
xbound_new = options.batchxbound;
ybound_new = options.batchybound;

%% Generate pool of points to extract batch from, and make predictions

pool.npoint = 1000;
% Halton sequence. Skip and Leap values chosen by user
skip = floor(rand*1e7);
leap = nthprime(sample.ndim+1)-1;
halton = haltonset(sample.ndim,'Skip',skip,'Leap',leap);
halton = scramble(halton,'RR2');
pool.raw = net(halton, pool.npoint);
    
% Map onto correct boundaries
pool.mapped = map_samples(param, pool.raw, xbound_new, ybound_new);
% Make predictions
[pool] = makeprediction(sample, pool, GEK);
% Add pool to batch struct in case needed to use later
batch.pool = pool;

%% Do interpolation of sumabs of gradients for 7 SA for all sample points
% The value of this sumabs will determine the radius around each batch
% point in x-y. This is the sum of the abs of the normalised gradients wrt
% the SA parameters as given by SU2
sumgrad_all = zeros(sample.npoint,1);
for i=1:sample.npoint
    sumgrad_all(i) = sumabs(sample.outputgrad_norm(i, param.cb1:param.cv1));
end
sumgrad_interp = ...
    scatteredInterpolant(sample.input(:,param.x), sample.input(:,param.y), sumgrad_all);

%% Populate batch points
% Taken from pool points with highest MSE

batch.npoint = options.nbatch;

% Initialise arrays
batch.pointxy = zeros(batch.npoint,2); % xy points of the batch
batch.point   = zeros(batch.npoint, sample.ndim); % full ndim points of the batch
batch.radius  = zeros(batch.npoint,1); % xy radius around each point
batch.mse     = zeros(batch.npoint,1); % mse of the prediction at these points
i = 1; % counter for qualified points
j = 0; % counter for disqualified points
disqualified = false;

while i-j <= batch.npoint
    
    if i == 1
        % The first pred point automatially gets included in batch since
        % there is no other radius yet
        batch.pointxy(i-j,:) = pool.mapped(pool.mse_sortindex(i),(param.x:param.y));
    else
        % next candidate point
        candidate = pool.mapped(pool.mse_sortindex(i),(param.x:param.y));
        % check to see if inside radius of other points, only then qualify
        for k = 1:i-j-1
            dist = pdist([candidate;batch.pointxy(k,:)]);
            if dist <= batch.radius(k)
                disqualified = true;
                break;
            else
                continue; % go to next point
            end
        end
        if ~disqualified
            % if not disqualified then assign candidate as valid point
            batch.pointxy(i-j,:) = candidate;
        else
            j = j+1; % disqualified
        end
    end
    
    if ~disqualified
        % Find closest original sample point to batch.point in X-Y space
        % and find the sumgrad at that point (COMMENTED OUT)
%         closest = dsearchn(sample.input(:,(param.x:param.y)),batch.pointxy(i-j,:));
%         sumgrad = sumabs(sample.outputgrad_norm(closest, param.cb1:param.cv1));

        % acquire the value of sumgrad at xy location of qualified batch point 
        sumgrad = sumgrad_interp(batch.pointxy(i-j,1), batch.pointxy(i-j,2));
        % Radius factor is log of the inverse of the sums, plus one so >1 always
        radius_factor = log10(1+1/sumgrad);
        % Find radius associated with new sample point
        batch.radius(i-j) = batch.maxrad * tanh(batch.p * radius_factor);
        
        % also assign full point with all 9 samples to point variable
        batch.point(i-j,:) = pool.mapped(pool.mse_sortindex(i),:);
        % include the mse of this point as well
        batch.mse(i-j,:) = pool.mse(pool.mse_sortindex(i));
    end
    
    % go to next point
    i = i+1; disqualified = false;
    
end

%% Write new batch of samples to file
if options.writebatch
    nextbatch_no = sample.nfiles + 1; % one after current sample length
    filename = sprintf('samples%02i.dat', nextbatch_no);
    file = fopen(fullfile('Samples/SU2_Input',filename),'w');
    fprintf(file, '%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s \n', ...
        'cb1','sig','cb2','kar','cw2','cw3','cv1','X','Y');
    for i = 1:batch.npoint
        fprintf(file, '%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f \n', ...
            batch.point(i,:));
    end
end
end


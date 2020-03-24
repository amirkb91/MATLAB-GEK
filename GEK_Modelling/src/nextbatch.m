function [batch] = nextbatch(sample, pred, param, nbatch, writebatch)
% Find next batch of samples in adaptive sampling and decluster in XY
% input: sample pred param structs, number of batches and write to file
% output: next batch of samples


% Number of GEK prediction points with highest MSE become the samples for
% the next SU2 iteration.
% The highest mse value will automatically become the first sample point in
% the next SU2 iteration. After that, a radius will be assinged in x and y
% to that point, which will depend on the value of the adjoint gradient at
% the closest available sample point to that point. If the adjoint gradient
% of the 7 SA parameters is low, the radius will be set to a larger value.
% The search for the next high mse point will be done outside of that
% radius. There is no need to cluster that radius with other points with
% similar XY values and different SA values, since the output is not that
% sensitive to SA anyway.
% batch points are taken from the suitable pred points on which GEK is done

%****** THESE 2 PARAMETERS CONTROL THE AMOUNT OF CLUSTERING *******
% Set the maximum possible radius allowed
batch.maxrad = 0.1;
% Set the tanh multiplication factor. larger p = more space b/w points
batch.p = 2;

% Number of points for next batch
batch.npoint = nbatch;

% Initialise arrays
batch.pointxy = zeros(batch.npoint,2); % xy points of the batch
batch.point   = zeros(batch.npoint, sample.ndim); % full ndim points of the batch
batch.radius  = zeros(batch.npoint,1); % xy radius around each point
i = 1; % counter for qualified points
j = 0; % counter for disqualified points
disqualified = false;

while i-j <= batch.npoint
    
    if i == 1
        % The first pred point automatially gets included in batch since
        % there is no other radius yet
        batch.pointxy(i-j,:) = pred.mapped(pred.mse_sortindex(i),(param.x:param.y));
    else
        % next candidate point
        candidate = pred.mapped(pred.mse_sortindex(i),(param.x:param.y));
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
        closest = dsearchn(sample.input(:,(param.x:param.y)),batch.pointxy(i-j,:));
        % Find the sum of the abs of the normalised gradients for the SA coeffs
        % at closest point
        gradsum = sumabs(sample.outputgrad_norm(closest, param.cb1:param.cv1));
        % Radius factor is log of the inverse of the sums, plus one so >1 always
        radius_factor = log10(1+1/gradsum);
        % Find radius associated with new sample point
        batch.radius(i-j) = batch.maxrad * tanh(batch.p * radius_factor);
        
        % also assign full point with all 9 samples to point variable
        batch.point(i-j,:) = pred.mapped(pred.mse_sortindex(i),:);
    end
    
    % go to next point
    i = i+1; disqualified = false;
    
end

% Write new batch of samples to file
if writebatch
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


function [optimum_theta, ln_likelihood] = hyperparameters(sample, theta_file)
% Obtain Gaussian Hyperparameter theta
% We assume Gaussian, therefore only theta is unknown and p=2.

% Check option whether to load theta or Genetic Algorithm
if ~isempty(theta_file)
    
    fprintf('\nLoading theta file %s \n', theta_file);
    optimum_theta = load(theta_file);
    optimum_theta = optimum_theta.optimum_theta;
    % find ln linkelihood for this value of theta
    R = corrmat(sample, optimum_theta); % Get correlation matrix
    [ln_likelihood] = get_lnlikelihood(R, sample);
    
else
    % Do Genetic Algorithm search for best theta to maximise ln likelihood
    % Define lower and upper search boundary of the log10(theta)
    theta_searchbound = [-4,3];
    
    % Find optimum log10(theta) which maximises the ln likelihood function
    % theta is an array with same length as number of dimension; a theta is
    % found for each dimension.
    [optimum_logtheta, ln_likelihood] = max_lnlikelihood(sample,theta_searchbound);

    % Convert log10 to real
    optimum_theta = 10.^(optimum_logtheta);   
    
end

% Print theta to screen
% list of parameters in order
paras = ["cb1";"sig";"cb2";"kar";"cw2";"cw3";"cv1";"x  ";"y  "];
[opthe, oppar]=sort(optimum_theta,'descend');
fprintf('\n');
disp('----- Optimum Theta -----')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Parameter    Theta(i)')
for i=1:sample.ndim
    line = sprintf('%s          %.4f\n',paras(oppar(i)), opthe(i));
    fprintf(line)
end
disp('~~~~~~~~~~~~~~~~~~~~~~~~~')

end


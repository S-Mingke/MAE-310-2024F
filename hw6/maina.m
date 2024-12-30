k = 2; % Regularity of the solution
l_values = [1, 2]; % Different values for l
h = linspace(0.1, 1, 10); % Mesh size values

% Calculate error estimates for each l
for l = l_values

    beta_L2 = min(k+1-0, 2*(k+1-l)); % L2-norm convergence rate
    beta_H1 = min(k+1-1, 2*(k+1-l)); % H1-norm convergence rate
    
    error_L2 = h.^beta_L2; % Error in L2-norm
    error_H1 = h.^beta_H1; % Error in H1-norm
    
    plot(h, error_H2, 'g'); 
    
    % Display results
    fprintf('For l = %d, the convergence rate in L2-norm is O(h^%f)\n', l, beta_L2);
    fprintf('For l = %d, the convergence rate in H1-norm is O(h^%f)\n', l, beta_H1);
end

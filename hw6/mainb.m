
m = 2;  % Smoothness of the solution for Bernoulli-Euler beam theory
k = 3;  % Degree of shape functions

h = linspace(0.1, 1, 10);  

% Calculate error estimates for L2-norm, H1-norm, and H2-norm
for l = 1:3  
    beta_L2 = min(k + 1 - 0, 2*(k + 1 - m)); % L2-norm convergence rate
    beta_H1 = min(k + 1 - 1, 2*(k + 1 - m)); % H1-norm convergence rate
    beta_H2 = min(k + 1 - 2, 2*(k + 1 - m)); % H2-norm convergence rate
    
    error_L2 = h.^beta_L2; % Error in L2-norm
    error_H1 = h.^beta_H1; % Error in H1-norm
    error_H2 = h.^beta_H2; % Error in H2-norm

    fprintf('For m = %d and k = %d, the convergence rates are:\n', m, k);
    fprintf('Convergence rate in L2-norm: O(h^%f)\n', beta_L2);
    fprintf('Convergence rate in H1-norm: O(h^%f)\n', beta_H1);
    fprintf('Convergence rate in H2-norm: O(h^%f)\n\n', beta_H2);
    
    % Optionally, plot error convergence for visualization
    figure;
    plot(h, error_L2, 'r', 'LineWidth', 2); hold on;
    plot(h, error_H1, 'b', 'LineWidth', 2);
    plot(h, error_H2, 'g', 'LineWidth', 2);
    xlabel('Mesh Size (h)');
    ylabel('Error');
    title(sprintf('Error Convergence for m = %d, k = %d', m, k));
    legend('L2-norm', 'H1-norm', 'H2-norm');
    grid on;
end

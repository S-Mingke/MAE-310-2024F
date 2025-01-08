% Parameters
T_x = 100;           % Applied far-field tension
R = 1;               % Radius of the hole
r_max = 5;           % Maximum radius of the plate
num_refinements = 4; % Number of mesh refinements
tol = 1e-6;          % Convergence tolerance

% Exact stress functions
sigma_rr_exact = @(r, theta) (T_x/2) * (1 - R^2./r.^2) + ...
                            (T_x/2) * (1 - 4*R^2./r.^2 + 3*R^4./r.^4) .* cos(2*theta);
sigma_theta_theta_exact = @(r, theta) (T_x/2) * (1 + R^2./r.^2) - ...
                                      (T_x/2) * (1 + 3*R^4./r.^4) .* cos(2*theta);
sigma_rtheta_exact = @(r, theta) -(T_x/2) * (1 + 2*R^2./r.^2 - 3*R^4./r.^4) .* sin(2*theta);


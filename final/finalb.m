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

% Finite Element Setup 
for refinement = 1:num_refinements
    % Generate mesh 
    [nodes, elements] = generate_mesh(refinement, R, r_max);  % Mesh generation function
    
    % Apply boundary conditions (symmetry and traction)
    boundary_conditions = apply_boundary_conditions(nodes, elements, R, r_max);  % Apply symmetry BCs and Neumann BCs
    
    % Solve the system using finite element method
    [displacements, stresses] = finite_element_solver(nodes, elements, boundary_conditions);  % FE solver
    
    % Calculate stress errors at each node (or element center)
    errors_rr = calculate_stress_error(sigma_rr_exact, nodes, stresses, 'rr');
    errors_theta_theta = calculate_stress_error(sigma_theta_theta_exact, nodes, stresses, 'theta_theta');
    errors_rtheta = calculate_stress_error(sigma_rtheta_exact, nodes, stresses, 'rtheta');
    
    % Calculate overall stress error and convergence rate
    total_error = norm(errors_rr + errors_theta_theta + errors_rtheta);
    disp(['Mesh refinement ', num2str(refinement), ' - Total Error: ', num2str(total_error)]);
    
    % Monitor convergence rate (optional)
    if refinement > 1
        convergence_rate = log(total_error_prev / total_error) / log(2);
        disp(['Convergence rate: ', num2str(convergence_rate)]);
    end
    
    % Save the previous error for convergence rate calculation
    total_error_prev = total_error;
    
    % Visualize the mesh
    figure;
    triplot(elements, nodes(:,1), nodes(:,2));  % Plot mesh
    title(['Mesh for refinement ', num2str(refinement)]);
    xlabel('r');
    ylabel('θ');
    
    % Visualize the exact stress 
    [R_grid, Theta_grid] = meshgrid(linspace(R, r_max, 50), linspace(0, pi/2, 50));
    sigma_rr_grid = sigma_rr_exact(R_grid, Theta_grid);
    
    figure;
    contourf(R_grid, Theta_grid, sigma_rr_grid, 20);  % Contour plot of σ_rr
    colorbar;
    title('Exact Stress Distribution \sigma_{rr}');
    xlabel('r');
    ylabel('\theta');
    
    % Visualize the FE stress 
    sigma_rr_FE = stresses(:, 1);  % Assuming stresses(:, 1) is σ_rr
    
    figure;
    scatter(nodes(:,1), nodes(:,2), 30, sigma_rr_FE, 'filled');
    colorbar;
    title(['FE Stress \sigma_{rr} for refinement ', num2str(refinement)]);
    xlabel('r');
    ylabel('\theta');
    
    % Visualize the error in stress
    error_plot = sqrt(errors_rr.^2 + errors_theta_theta.^2 + errors_rtheta.^2);
    
    figure;
    scatter(nodes(:,1), nodes(:,2), 30, error_plot, 'filled');
    colorbar;
    title(['Stress Error for refinement ', num2str(refinement)]);
    xlabel('r');
    ylabel('\theta');
end

% Helper Functions %%
function [nodes, elements] = generate_mesh(refinement, R, r_max)
    
    % Placeholder: Generate a simple uniform grid 
    theta = linspace(0, pi/2, 10*refinement);
    r = linspace(R, r_max, 10*refinement);
    [R_grid, Theta_grid] = meshgrid(r, theta);
    nodes = [R_grid(:), Theta_grid(:)];
    
    % Elements would be triangles, but for simplicity, let's assume quadrilaterals for this demo
    elements = delaunay(nodes(:, 1), nodes(:, 2));
end

function boundary_conditions = apply_boundary_conditions(nodes, elements, R, r_max)
    
    
    boundary_conditions = struct();
    boundary_conditions.symmetry_nodes = nodes(nodes(:,1) == R | nodes(:,2) == 0 | nodes(:,2) == pi/2, :);
    boundary_conditions.traction_nodes = nodes(nodes(:,1) == r_max, :);  % Neumann BCs at r_max
end

function [displacements, stresses] = finite_element_solver(nodes, elements, boundary_conditions)
    % Placeholder for FE solver
    % In a real implementation, you'd use an FE solver like MATLAB's PDE toolbox or write your own solver
    displacements = zeros(size(nodes, 1), 2);  % 2D displacements (r, theta)
    stresses = zeros(size(nodes, 1), 3);  % 3 stress components (rr, theta_theta, rtheta)
    
    % Placeholder for FE solution logic
    % Displacement solutions should be computed based on the FE assembly and solution process.
end

function errors = calculate_stress_error(exact_stress, nodes, stresses, stress_type)
    % Compute the error between exact stress and FE stresses
    errors = zeros(size(nodes, 1), 1);
    for i = 1:length(nodes)
        r = nodes(i, 1);
        theta = nodes(i, 2);
        
        % Calculate exact stress
        if strcmp(stress_type, 'rr')
            exact = exact_stress(r, theta);
        elseif strcmp(stress_type, 'theta_theta')
            exact = exact_stress(r, theta);
        elseif strcmp(stress_type, 'rtheta')
            exact = exact_stress(r, theta);
        end
        
        % Calculate FE stress
        fe_stress = stresses(i, find(strcmp({'rr', 'theta_theta', 'rtheta'}, stress_type)));
        
        % Compute the error (absolute difference)
        errors(i) = abs(exact - fe_stress);
    end
end

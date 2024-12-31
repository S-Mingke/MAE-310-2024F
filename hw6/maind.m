

u_analytic = @(x, y) sin(pi * x) .* sin(pi * y); % Manufactured solution

h_values = [0.5, 0.25, 0.125, 0.0625]; % Mesh sizes (h)

error_L2_quad = zeros(length(h_values), 1);
error_L2_tri = zeros(length(h_values), 1);

% Initialize arrays for H1 errors (derivatives)
error_H1_quad = zeros(length(h_values), 1);
error_H1_tri = zeros(length(h_values), 1);

% Loop over different mesh sizes to compute errors
for idx = 1:length(h_values)
    h = h_values(idx);
    
    % Generate the quadrilateral mesh and solve
    [nodes_quad, element_nodes_quad] = generate_quad_mesh(h);
    uh_quad = solve_fem(nodes_quad, element_nodes_quad, u_analytic);
    
    % Compute L2 and H1 errors for quadrilateral mesh
    error_L2_quad(idx) = compute_L2_error(nodes_quad, element_nodes_quad, uh_quad, u_analytic);
    error_H1_quad(idx) = compute_H1_error(nodes_quad, element_nodes_quad, uh_quad, u_analytic);
    
    % Generate the triangular mesh and solve
    [nodes_tri, element_nodes_tri] = generate_tri_mesh(h);
    uh_tri = solve_fem(nodes_tri, element_nodes_tri, u_analytic);
    
    % Compute L2 and H1 errors for triangular mesh
    error_L2_tri(idx) = compute_L2_error(nodes_tri, element_nodes_tri, uh_tri, u_analytic);
    error_H1_tri(idx) = compute_H1_error(nodes_tri, element_nodes_tri, uh_tri, u_analytic);
end

% Plot the log-log convergence rates for both quadrilateral and triangular elements
figure;
loglog(h_values, error_L2_quad, 'o-', 'LineWidth', 2, 'DisplayName', 'Quad L2 Error');
hold on;
loglog(h_values, error_L2_tri, 's-', 'LineWidth', 2, 'DisplayName', 'Tri L2 Error');
xlabel('Mesh Size (h)');
ylabel('L2 Error');
title('Convergence Rate for L2 Error');
legend;
grid on;

figure;
loglog(h_values, error_H1_quad, 'o-', 'LineWidth', 2, 'DisplayName', 'Quad H1 Error');
hold on;
loglog(h_values, error_H1_tri, 's-', 'LineWidth', 2, 'DisplayName', 'Tri H1 Error');
xlabel('Mesh Size (h)');
ylabel('H1 Error');
title('Convergence Rate for H1 Error');
legend;
grid on;

% Function to generate the quadrilateral mesh
function [nodes, element_nodes] = generate_quad_mesh(h)
    % Define the size of the domain
    x_max = 1;
    y_max = 1;
    
    % Number of elements per dimension
    n = round(1/h);
    
    % Generate mesh nodes
    [X, Y] = meshgrid(0:h:x_max, 0:h:y_max);
    nodes = [X(:), Y(:)];
    
    % Generate the connectivity (element_nodes)
    element_nodes = [];
    for i = 1:n
        for j = 1:n
            % Node indices for a quadrilateral
            idx1 = (i-1)*(n+1) + j;
            idx2 = idx1 + 1;
            idx3 = idx1 + (n+1);
            idx4 = idx3 + 1;
            element_nodes = [element_nodes; idx1, idx2, idx3, idx4];
        end
    end
end

% Function to generate the triangular mesh (subdivide quadrilaterals into triangles)
function [nodes, element_nodes] = generate_tri_mesh(h)
    % Define the size of the domain
    x_max = 1;
    y_max = 1;
    
    % Number of elements per dimension
    n = round(1/h);
    
    % Generate mesh nodes
    [X, Y] = meshgrid(0:h:x_max, 0:h:y_max);
    nodes = [X(:), Y(:)];
    
    % Generate the connectivity (element_nodes) for triangles
    element_nodes = [];
    for i = 1:n
        for j = 1:n
            % Node indices for a quadrilateral
            idx1 = (i-1)*(n+1) + j;
            idx2 = idx1 + 1;
            idx3 = idx1 + (n+1);
            idx4 = idx3 + 1;
            
            % Divide quadrilateral into two triangles
            element_nodes = [element_nodes; idx1, idx2, idx3];
            element_nodes = [element_nodes; idx2, idx4, idx3];
        end
    end
end

% Function to solve the FEM problem (stub for solving system, for demonstration)
function uh = solve_fem(nodes, element_nodes, u_analytic)
    % A simplified solver for FEM problem (not implemented here)
    % In a real case, this would involve assembling the stiffness matrix and solving for uh
    uh = zeros(size(nodes, 1), 1);  % Placeholder for FEM solution
end

% Compute L2 error
function error_L2 = compute_L2_error(nodes, element_nodes, uh, u_analytic)
    error_L2 = 0;
    for i = 1:size(element_nodes, 1)
        element_coords = nodes(element_nodes(i, :), :);
        u_exact = u_analytic(element_coords(:, 1), element_coords(:, 2));
        uh_values = uh(element_nodes(i, :));
        error_L2 = error_L2 + norm(uh_values - u_exact, 2);
    end
end

% Compute H1 error
function error_H1 = compute_H1_error(nodes, element_nodes, uh, u_analytic)
    error_H1 = 0;
    for i = 1:size(element_nodes, 1)
        element_coords = nodes(element_nodes(i, :), :);
        u_exact = u_analytic(element_coords(:, 1), element_coords(:, 2));
        uh_values = uh(element_nodes(i, :));
        % Compute the gradient of the error (this is a simplified version)
        error_H1 = error_H1 + norm(uh_values - u_exact, 2);  % Placeholder for H1 error
    end
end

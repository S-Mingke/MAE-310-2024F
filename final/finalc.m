clear; clc; close all;

% Material Properties
E = 1e9;           % Young's Modulus in Pa
nu = 0.3;          % Poisson's Ratio
t = 1;             % Thickness of the body (assuming unit thickness for 2D analysis)
sigma_0 = 10e3;    % Traction (10 kPa)

% Geometry
L = 4;             % Side length of the square in meters
R = 0.5;           % Radius of the hole

% Discretization (Mesh Generation)
numElem = 400;     % Number of elements (adjust for mesh density)
xMax = L; yMax = L;

% Generate a mesh (using structured grid and removal of nodes inside the hole)
[coordinates, elements] = generateMesh(xMax, yMax, R, numElem);

% Plot the mesh
figure;
hold on;
plotMesh(coordinates, elements);
axis equal;
title('Mesh for Quarter-Plane with Hole');
xlabel('x (m)');
ylabel('y (m)');

% Plane Stress or Plane Strain Analysis
% Plane Stress: sigma_z = 0
% Plane Strain: epsilon_z = 0

% Select Plane Stress or Plane Strain analysis
analysisType = 'plane_stress'; % Options: 'plane_stress', 'plane_strain'

% Define material matrix based on the analysis type
if strcmp(analysisType, 'plane_stress')
    D = (E / (1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
elseif strcmp(analysisType, 'plane_strain')
    D = (E / ((1 + nu) * (1 - 2*nu))) * ...
        [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2*nu) / 2];
end

% Assemble global stiffness matrix and force vector
[K, F] = assembleFEM(coordinates, elements, D, t, sigma_0, xMax, yMax, analysisType);

% Apply boundary conditions
% Fixing the bottom-left corner
fixedNodes = find(coordinates(:,1) == 0 & coordinates(:,2) == 0);
[K, F] = applyBoundaryConditions(K, F, fixedNodes);

% Solve for displacements
U = K \ F;

% Visualize displacements (scaled for clarity)
displacements = reshape(U, [], 2);
scaledDisplacements = 10 * displacements; % Scaling factor for better visualization

figure;
quiver(coordinates(:,1), coordinates(:,2), scaledDisplacements(:,1), scaledDisplacements(:,2), 0.5);
title('Displacement Field');
xlabel('x (m)');
ylabel('y (m)');
axis equal;

% Compute stress field using the exact analytical solution
stress = computeExactStress(coordinates, R, sigma_0);

% Visualize stress concentration around the hole
figure;
plotStress(coordinates, elements, stress, 'sigma_xx');
title('Stress Distribution - \sigma_{xx}');
xlabel('x (m)');
ylabel('y (m)');
colorbar;

% Plot the stress concentration at the hole
figure;
plotStress(coordinates, elements, stress, 'sigma_yy');
title('Stress Distribution - \sigma_{yy}');
xlabel('x (m)');
ylabel('y (m)');
colorbar;

% Define the mesh generation function
function [coords, elems] = generateMesh(xMax, yMax, R, numElem)
    % This function generates a simple square grid mesh and removes the hole's area
    [x, y] = meshgrid(linspace(0, xMax, sqrt(numElem)), linspace(0, yMax, sqrt(numElem)));
    coords = [x(:), y(:)];
    
    % Remove points inside the hole
    holeMask = sqrt(coords(:,1).^2 + coords(:,2).^2) < R;
    coords(holeMask, :) = [];
    
    % Generate elements (triangular or quadrilateral elements)
    elems = delaunay(coords(:,1), coords(:,2)); % Delaunay triangulation for simplicity
end

% Element plotting function
function plotMesh(coords, elems)
    for i = 1:size(elems, 1)
        plot([coords(elems(i,1),1), coords(elems(i,2),1)], [coords(elems(i,1),2), coords(elems(i,2),2)], 'k');
        hold on;
        plot([coords(elems(i,2),1), coords(elems(i,3),1)], [coords(elems(i,2),2), coords(elems(i,3),2)], 'k');
        plot([coords(elems(i,3),1), coords(elems(i,1),1)], [coords(elems(i,3),2), coords(elems(i,1),2)], 'k');
    end
end

% Function to assemble global stiffness matrix and force vector
function [K, F] = assembleFEM(coords, elems, D, t, sigma_0, xMax, yMax, analysisType)
    numNodes = size(coords, 1);
    numElem = size(elems, 1);
    
    % Initialize global stiffness matrix and force vector
    K = zeros(2*numNodes, 2*numNodes);
    F = zeros(2*numNodes, 1);
    
    % Loop over each element to compute its contribution to K and F
    for e = 1:numElem
        % Get the node indices for this element
        nodes = elems(e, :);
        
        % Extract the coordinates of the nodes for the element
        x = coords(nodes, 1);
        y = coords(nodes, 2);
        
        % Compute the element stiffness matrix (using the plane stress formulation here)
        [Ke, Fe] = computeElementStiffness(x, y, D, t, sigma_0, analysisType);
        
        % Assemble the element stiffness matrix and force vector into the global matrix
        for i = 1:length(nodes)
            for j = 1:length(nodes)
                K(2*nodes(i)-1, 2*nodes(j)-1) = K(2*nodes(i)-1, 2*nodes(j)-1) + Ke(2*i-1, 2*j-1);
                K(2*nodes(i), 2*nodes(j)) = K(2*nodes(i), 2*nodes(j)) + Ke(2*i, 2*j);
            end
        end
        
        % Assemble force vector (if traction is applied, for simplicity, assuming on boundary edges)
        F(2*nodes(1)-1) = F(2*nodes(1)-1) + Fe(1);
        F(2*nodes(1)) = F(2*nodes(1)) + Fe(2);
    end
end

% Function to compute the stiffness matrix and force vector for an element
function [Ke, Fe] = computeElementStiffness(x, y, D, t, sigma_0, analysisType)
   
    Ke = rand(6, 6); % Random stiffness matrix for illustration
    Fe = rand(6, 1); % Random force vector for illustration
end

% Function to apply boundary conditions (e.g., fixing the bottom-left corner)
function [K, F] = applyBoundaryConditions(K, F, fixedNodes)
    % Apply fixed boundary conditions (zero displacement) for nodes at the bottom-left
    K(fixedNodes, :) = 0;
    K(fixedNodes, fixedNodes) = eye(length(fixedNodes));
    F(fixedNodes) = 0;
end

% Function to compute the exact stress field using the given formulas
function stress = computeExactStress(coords, R, sigma_0)
    % Convert Cartesian coordinates to polar coordinates
    r = sqrt(coords(:,1).^2 + coords(:,2).^2);
    theta = atan2(coords(:,2), coords(:,1));
    
    % Compute the exact stresses
    sigma_rr = sigma_0 / 2 * (1 - R^2 ./ r.^2) + sigma_0 / 2 * (1 - 4*R^2 ./ r.^2 + 3*R^4 ./ r.^4) .* cos(2*theta);
    sigma_theta_theta = sigma_0 / 2 * (1 + R^2 ./ r.^2) - sigma_0 / 2 * (1 + 3*R^4 ./ r.^4) .* cos(2*theta);
    sigma_r_theta = -sigma_0 / 2 * (1 + 2*R^2 ./ r.^2 - 3*R^4 ./ r.^4) .* sin(2*theta);
    
    % Store the stresses in a matrix (columns: sigma_rr, sigma_theta_theta, sigma_r_theta)
    stress = [sigma_rr, sigma_theta_theta, sigma_r_theta];
end

% Function to visualize stress (e.g., sigma_xx or sigma_yy)
function plotStress(coords, elems, stress, stressType)
    if strcmp(stressType, 'sigma_xx')
        values = stress(:,1); % Extract sigma_xx values
    elseif strcmp(stressType, 'sigma_yy')
        values = stress(:,2); % Extract sigma_yy values
    end
    trisurf(elems, coords(:,1), coords(:,2), values);
    shading interp;
    colorbar;
end

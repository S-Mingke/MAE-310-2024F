

% Number of grid points along x and y
nX = 4;  % number of divisions in x-direction
nY = 4;  % number of divisions in y-direction

% Define the grid size
Lx = 10;  % Length in x-direction
Ly = 10;  % Length in y-direction

x = linspace(0, Lx, nX+1);  % x coordinates
y = linspace(0, Ly, nY+1);  % y coordinates
[X, Y] = meshgrid(x, y);    % Create mesh grid

% Flatten the node arrays to form a list of nodes
nodes = [X(:), Y(:)];

% Plot the grid of nodes
figure;
plot(nodes(:,1), nodes(:,2), 'ro');
hold on;

numNodes = size(nodes, 1);

elements_quad = [];
for i = 1:nX
    for j = 1:nY
        n1 = (i-1)*(nY+1) + j;      % Node index 1
        n2 = n1 + 1;                 % Node index 2
        n3 = n1 + (nY+1);           % Node index 3
        n4 = n3 + 1;                 % Node index 4
        
        % Create two triangles for each quadrilateral
        elements_quad = [elements_quad; n1, n2, n3];  % Triangle 1
        elements_quad = [elements_quad; n2, n4, n3];  % Triangle 2
    end
end

% Number of elements after subdivision
numElements = size(elements_quad, 1);

% Plot the mesh of triangles
for i = 1:numElements
    triNodes = nodes(elements_quad(i,:), :);
    fill(triNodes(:,1), triNodes(:,2), 'c', 'FaceAlpha', 0.3); % Draw each triangle
end
title('Mesh with Triangular Elements');
xlabel('X');
ylabel('Y');
axis equal;

% Display element connectivity (each row corresponds to one triangle)
disp('Element Connectivity (Each row is a triangle):');
disp(elements_quad);

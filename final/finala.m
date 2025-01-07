% 网格生成
Lx = 1; Ly = 1; nX = 10; nY = 10;
[x, y] = meshgrid(linspace(0, Lx, nX), linspace(0, Ly, nY));
coordinates = [x(:), y(:)];
elements = [];
for i = 1:nY-1
    for j = 1:nX-1
        node1 = (i-1)*nX + j;
        node2 = node1 + 1;
        node3 = node1 + nX;
        node4 = node3 + 1;
        elements = [elements; node1, node2, node3; node2, node4, node3];
    end
end

% 材料参数
E = 200e9; nu = 0.3; T_x = 100;
plane_stress = true;
if plane_stress
    C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
else
    C = E / ((1 + nu) * (1 - 2*nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2*nu) / 2];
end

% 边界条件
boundary_nodes = find(coordinates(:,1) == 0 | coordinates(:,2) == 0);
displacement = zeros(size(coordinates, 1), 1);
displacement(boundary_nodes) = 0;

% 应力边界条件
R = 0.2;
traction = zeros(length(boundary_nodes), 3);
for i = 1:length(boundary_nodes)
    r = sqrt(coordinates(boundary_nodes(i), 1)^2 + coordinates(boundary_nodes(i), 2)^2);
    theta = atan2(coordinates(boundary_nodes(i), 2), coordinates(boundary_nodes(i), 1));
    if r ~= 0
        sigma_rr = T_x/2 * (1 - R^2 / r^2) + T_x/2 * (1 - 4*R^2 / r^2 + 3*R^4 / r^4) * cos(2*theta);
        sigma_tt = T_x/2 * (1 + R^2 / r^2) - T_x/2 * (1 + 3*R^4 / r^4) * cos(2*theta);
        sigma_rt = -T_x/2 * (1 + 2*R^2 / r^2 - 3*R^4 / r^4) * sin(2*theta);
        traction(i, :) = [sigma_rr, sigma_tt, sigma_rt];
    else
        traction(i, :) = [0, 0, 0];
    end
end

% 假设求解有限元问题
stress = rand(size(elements, 1), 1); % 假设应力解

% 可视化
figure;
trisurf(elements, coordinates(:,1), coordinates(:,2), displacement);
title('Displacement Field'); xlabel('X'); ylabel('Y');

figure;
trisurf(elements, coordinates(:,1), coordinates(:,2), stress);
title('Stress Field'); xlabel('X'); ylabel('Y');

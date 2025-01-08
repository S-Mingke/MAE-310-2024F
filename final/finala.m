% 网格生成
Lx = 1;  % 网格的x方向长度
Ly = 1;  % 网格的y方向长度
nX = 10; % x方向的网格节点数
nY = 10; % y方向的网格节点数

% 生成均匀的节点坐标
[x, y] = meshgrid(linspace(0, Lx, nX), linspace(0, Ly, nY));
coordinates = [x(:), y(:)];

% 生成网格元素（四节点矩形元素）
elements = [];
for i = 1:nY-1
    for j = 1:nX-1
        node1 = (i-1)*nX + j;
        node2 = node1 + 1;
        node3 = node1 + nX;
        node4 = node3 + 1;
        % 将每个矩形元素分解为两个三角形
        elements = [elements; node1, node2, node3]; % 三角形1
        elements = [elements; node2, node4, node3]; % 三角形2
    end
end

% 材料参数
E = 200e9;   % 弹性模量（Pa）
nu = 0.3;    % 泊松比
T_x = 100;   % 外部加载（远场应力 Tx）

% 平面应力
plane_stress = true;
if plane_stress
    C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
else
    C = E / ((1 + nu) * (1 - 2*nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2*nu) / 2];
end

% 边界条件：应用位移（Dirichlet）或应力（Neumann）条件
boundary_nodes = findBoundaryNodes(coordinates);  % 找到边界节点
% 对于这些边界节点应用零位移（Dirichlet 边界条件）
displacement = zeros(length(coordinates), 1);
displacement(boundary_nodes) = 0;  % 施加位移边界条件

% 应力边界条件：使用解析解在边界上应用应力
R = 0.2;  % 圆孔半径
traction = zeros(length(boundary_nodes), 3);  % 3个方向上的应力
for i = 1:length(boundary_nodes)
    r = sqrt(coordinates(boundary_nodes(i), 1)^2 + coordinates(boundary_nodes(i), 2)^2);
    theta = atan2(coordinates(boundary_nodes(i), 2), coordinates(boundary_nodes(i), 1));
    
    % 使用解析公式计算应力
    if r ~= 0 
        sigma_rr = T_x/2 * (1 - R^2 / r^2) + T_x/2 * (1 - 4*R^2 / r^2 + 3*R^4 / r^4) * cos(2*theta);
        sigma_tt = T_x/2 * (1 + R^2 / r^2) - T_x/2 * (1 + 3*R^4 / r^4) * cos(2*theta);
        sigma_rt = -T_x/2 * (1 + 2*R^2 / r^2 - 3*R^4 / r^4) * cos(2*theta);
        
        % 应力边界条件
        traction(i, :) = [sigma_rr, sigma_tt, sigma_rt];
    else
        traction(i, :) = [0, 0, 0];  % 如果 r = 0，则设置为零应力
    end
end

% 求解有限元问题
[displacement, stress] = solveFEA(coordinates, elements, C, displacement, traction);

% 可视化结果：位移场和应力场
% 位移场可视化
figure;
trisurf(elements, coordinates(:,1), coordinates(:,2), displacement);
title('Displacement Field');
xlabel('X');
ylabel('Y');

% 应力场可视化
% 将元素的应力值映射到节点
node_stress = mapStressToNodes(elements, stress, length(coordinates));
figure;
trisurf(elements, coordinates(:,1), coordinates(:,2), node_stress);
title('Stress Field');
xlabel('X');
ylabel('Y');
zlabel('Stress');

% 误差计算：L2-和H1-norm
exact_stress = computeExactStress(coordinates, T_x, R);  % 使用解析解
error_L2 = computeL2Error(node_stress, exact_stress);  % L2误差
error_H1 = computeH1Error(coordinates, node_stress, exact_stress);  % H1误差
fprintf('L2 Error: %.4e\n', error_L2);
fprintf('H1 Error: %.4e\n', error_H1);

% 网格精细化与误差监测
num_refinements = 5;  % 网格精细化次数
for i = 1:num_refinements
    refined_mesh = refineMesh(coordinates, elements);  % 精细化网格
    [displacement, stress] = solveFEA(refined_mesh.Coordinates, refined_mesh.Elements, C, displacement, traction);  % 求解
    node_stress = mapStressToNodes(refined_mesh.Elements, stress, length(refined_mesh.Coordinates));  % 应力值映射
    error_L2 = computeL2Error(node_stress, exact_stress);  % L2误差
    fprintf('Refinement %d: L2 Error = %.4e\n', i, error_L2);
end



% 计算边界节点函数
function boundary_nodes = findBoundaryNodes(coordinates)
    % 假设边界节点是位于x=0或y=0的节点
    boundary_nodes = find(coordinates(:,1) == 0 | coordinates(:,2) == 0);
end

% 求解有限元问题函数
function [displacement, stress] = solveFEA(coordinates, elements, C, displacement, traction)

    displacement = rand(size(coordinates, 1), 1);  % 假设一个位移解
    stress = rand(size(elements, 1), 1);  % 假设一个应力解
end

% 计算解析应力函数
function exact_stress = computeExactStress(coordinates, T_x, R)
    % 计算解析应力
    r = sqrt(coordinates(:, 1).^2 + coordinates(:, 2).^2);
    theta = atan2(coordinates(:, 2), coordinates(:, 1));
    
    % 处理r = 0的情况，避免除以零
    sigma_rr = T_x/2 * (1 - R^2 ./ r.^2);
    sigma_tt = T_x/2 * (1 + R^2 ./ r.^2);
    sigma_rt = -T_x/2 * (1 + 2*R^2 ./ r.^2);
    
    % 在r=0时设置应力为零
    sigma_rr(r == 0) = 0;
    sigma_tt(r == 0) = 0;
    sigma_rt(r == 0) = 0;

    % 计算解析应力
    exact_stress = [sigma_rr, sigma_tt, sigma_rt];  % 包含所有应力分量
end

% L2误差计算函数
function error_L2 = computeL2Error(computed_stress, exact_stress)
    % 确保两个数组大小一致
    if length(computed_stress) ~= length(exact_stress)
        error('Stress arrays must have the same length.');
    end
    % 计算L2误差
    error_L2 = sqrt(sum((computed_stress - exact_stress).^2) / length(computed_stress));
end

% H1误差计算函数
function error_H1 = computeH1Error(coordinates, computed_stress, exact_stress)
    % 计算梯度误差，可以通过数值方法来近似计算梯度
    grad_exact = computeGradientExactStress(coordinates, exact_stress);  % 计算解析应力梯度
    error_H1 = sqrt(sum((computed_stress - grad_exact).^2) / length(computed_stress));
end

% 计算解析应力梯度（占位符）
function grad_exact = computeGradientExactStress(coordinates, exact_stress)
    grad_exact = exact_stress;  % 这里只是占位符，实际中应计算梯度
end

% 网格精细化函数
function refined_mesh = refineMesh(coordinates, elements)
    refined_mesh.Coordinates = coordinates * 1.1;  % 假设增加10%的尺度
    refined_mesh.Elements = elements;  % 假设元素数量不变
end

% 将元素应力映射到节点
function node_stress = mapStressToNodes(elements, stress, n_nodes)
    node_stress = zeros(n_nodes, 1);
    for i = 1:length(elements)
        nodes = elements(i, :);  % 获取元素的节点编号
        % 假设应力值均匀分配给每个节点
        node_stress(nodes) = node_stress(nodes) + stress(i) / 3;
    end
end

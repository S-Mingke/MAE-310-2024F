% 函数定义
u_exact = @(x) sin(pi * x);  % 精确解
u_exact_x = @(x) pi * cos(pi * x);  % 精确解的导数

% 网格划分
elements = [2, 4, 6, 8, 10, 12, 14, 16];  % 不同的元素数量
nq_values = [1, 2, 3, 4, 5, 6];  % 高斯积分点数
errors_L2 = zeros(length(elements), length(nq_values));  % 存储L2误差
errors_H1 = zeros(length(elements), length(nq_values));  % 存储H1误差
mesh_sizes = zeros(length(elements), 1);  % 存储网格大小

min_error = 1e-10;  % 设置最小误差值

for i = 1:length(elements)
    N = elements(i);  % 当前元素数量
    h = 1 / N;  % 网格大小
    mesh_sizes(i) = h;
    
    % 网格节点
    x = linspace(0, 1, N+1);  % 网格节点
    
    % 假设数值解u_h是线性插值的结果，这里用sin(pi * x)作为初始近似解
    u_h = sin(pi * x);  % 这是假设的初始解，可以根据需要调整
    
    % 对每个高斯积分点数量进行实验
    for j = 1:length(nq_values)
        nq = nq_values(j);  % 当前的高斯积分点数量
        
        % 计算L2误差
        L2_numerator = 0;  % 分子
        L2_denominator = 0;  % 分母
        for e = 1:N
            % 获取每个元素的积分区间
            xe = x(e:e+1);
            % 高斯-勒让德积分
            [quad_pts, weights] = gauss_legendre_quadrature(nq);  % 获取当前积分点和权重
            for k = 1:nq
                % 计算积分
                xi = quad_pts(k);
                w = weights(k);
                % 映射到元素坐标
                x_q = 0.5 * (1 - xi) * xe(1) + 0.5 * (1 + xi) * xe(2);
                u_exact_val = u_exact(x_q);
                u_h_val = u_h(e);  % 假设数值解是每个元素的节点值
                
                % L2误差
                L2_numerator = L2_numerator + w * (u_exact_val - u_h_val)^2;
                L2_denominator = L2_denominator + w * u_exact_val^2;
            end
        end
        errors_L2(i, j) = sqrt(L2_numerator) / sqrt(L2_denominator);
        
        % 计算H1误差
        u_h_x = pi * cos(pi * x);  % 数值解的导数
        u_exact_x_values = u_exact_x(x);  % 计算精确解的导数
        H1_numerator = 0;  % 分子
        H1_denominator = 0;  % 分母
        for e = 1:N
            % 获取每个元素的积分区间
            xe = x(e:e+1);
          
            [quad_pts, weights] = gauss_legendre_quadrature(nq); 
            for k = 1:nq
                % 计算积分
                xi = quad_pts(k);
                w = weights(k);
                % 映射到元素坐标
                x_q = 0.5 * (1 - xi) * xe(1) + 0.5 * (1 + xi) * xe(2);
                u_exact_x_val = u_exact_x(x_q);
                u_h_x_val = u_h_x(e);  % 假设数值解的导数
                
                % H1误差
                H1_numerator = H1_numerator + w * (u_h_x_val - u_exact_x_val)^2;
                H1_denominator = H1_denominator + w * u_exact_x_val^2;
            end
        end
        errors_H1(i, j) = sqrt(H1_numerator) / sqrt(H1_denominator);
    end
end

% 避免误差过小导致log-log图不可见
errors_L2(errors_L2 < min_error) = min_error;
errors_H1(errors_H1 < min_error) = min_error;

% 绘制log-log图
figure;
hold on;
for j = 1:length(nq_values)
    loglog(mesh_sizes, errors_L2(:, j), '-o', 'DisplayName', ['L2 Error (nq = ' num2str(nq_values(j)) ')']);
    loglog(mesh_sizes, errors_H1(:, j), '-s', 'DisplayName', ['H1 Error (nq = ' num2str(nq_values(j)) ')']);
end
xlabel('log(h)', 'FontSize', 12);
ylabel('log(Error)', 'FontSize', 12);
title('Relative Errors vs Mesh Size with Different Quadrature Points');
legend show;
grid on;

% 高斯-勒让德积分点和权重计算
function [quad_pts, weights] = gauss_legendre_quadrature(nq)
    % 根据积分点数目返回对应的积分点和权重
    if nq == 1
        quad_pts = 0;
        weights = 2;
    elseif nq == 2
        quad_pts = [-1/sqrt(3), 1/sqrt(3)];
        weights = [1, 1];
    elseif nq == 3
        quad_pts = [-sqrt(3/5), 0, sqrt(3/5)];
        weights = [5/9, 8/9, 5/9];
    elseif nq == 4
        quad_pts = [-sqrt((3 + 2*sqrt(6/5))/7), -sqrt((3 - 2*sqrt(6/5))/7), sqrt((3 - 2*sqrt(6/5))/7), sqrt((3 + 2*sqrt(6/5))/7)];
        weights = [(18 - sqrt(30))/36, (18 + sqrt(30))/36, (18 + sqrt(30))/36, (18 - sqrt(30))/36];
    elseif nq == 5
        quad_pts = [-sqrt(5/9), -sqrt(1/5), 0, sqrt(1/5), sqrt(5/9)];
        weights = [128/225, (322 + 13*sqrt(70))/225, 72/225, (322 - 13*sqrt(70))/225, 128/225];
    elseif nq == 6
        quad_pts = [-sqrt(1/3), -sqrt(0.2), -sqrt(0.2), 0, sqrt(0.2), sqrt(1/3)];
        weights = [5/9, 8/9, 8/9, 8/9, 8/9, 5/9];
    else
        error('Only supports 1 to 6 quadrature points.');
    end
end


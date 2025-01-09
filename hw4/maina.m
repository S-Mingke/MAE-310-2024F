% 函数定义
u_exact = @(x) sin(pi * x); 
u_exact_x = @(x) pi * cos(pi * x);  

% 网格划分
elements = [2, 4, 6, 8, 10, 12, 14, 16];  
errors_L2 = zeros(length(elements), 1);  
errors_H1 = zeros(length(elements), 1); 
mesh_sizes = zeros(length(elements), 1); 

min_error = 1e-10; 
for i = 1:length(elements)
    N = elements(i);  % 当前元素数量
    h = 1 / N;  % 网格大小
    mesh_sizes(i) = h;
    
    % 网格节点
    x = linspace(0, 1, N+1);  % 网格节点
    
    u_h = sin(pi * x) + 0.1 * randn(size(x));  % 数值解加上一些噪声
    
    % 计算L2误差
    L2_numerator = trapz(x, (u_exact(x) - u_h).^2); 
    L2_denominator = trapz(x, u_exact(x).^2); 
    errors_L2(i) = sqrt(L2_numerator) / sqrt(L2_denominator);
    
    % 计算H1误差
    u_h_x = pi * cos(pi * x) + 0.1 * randn(size(x));  
    u_exact_x_values = u_exact_x(x);  
    H1_numerator = trapz(x, (u_h_x - u_exact_x_values).^2);  
    H1_denominator = trapz(x, u_exact_x_values.^2);  
    errors_H1(i) = sqrt(H1_numerator) / sqrt(H1_denominator);
    
    % 打印调试信息
    disp(['Number of elements: ', num2str(N)]);
    disp(['L2 error: ', num2str(errors_L2(i))]);
    disp(['H1 error: ', num2str(errors_H1(i))]);
end

errors_L2(errors_L2 < min_error) = min_error;
errors_H1(errors_H1 < min_error) = min_error;

% 绘制log-log图
figure;
loglog(mesh_sizes, errors_L2, '-o', 'DisplayName', 'L2 Error');
hold on;
loglog(mesh_sizes, errors_H1, '-s', 'DisplayName', 'H1 Error');
xlabel('log(h)', 'FontSize', 12);
ylabel('log(Error)', 'FontSize', 12);
title('Relative Errors vs Mesh Size');
legend show;
grid on;

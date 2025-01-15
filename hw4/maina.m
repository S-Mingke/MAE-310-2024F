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
    N = elements(i); 
    h = 1 / N;  
    mesh_sizes(i) = h;
    
    % 网格节点
    x = linspace(0, 1, N+1); 
    
    % 系统矩阵初始化
    A = zeros(N+1, N+1);  
    F = zeros(N+1, 1);    
    
    % 装配刚度矩阵和载荷向量
    for j = 1:N
        % 单元的局部刚度矩阵
        A_local = [1, -1; -1, 1] / h; 
        % 单元的局部载荷向量
        F_local = h / 2 * [u_exact(x(j)); u_exact(x(j+1))];  
        
      
        A(j:j+1, j:j+1) = A(j:j+1, j:j+1) + A_local;
        F(j:j+1) = F(j:j+1) + F_local;
    end
    
    % 施加边界条件
    A(1,:) = 0; A(1,1) = 1; F(1) = u_exact(0);
    A(N+1,:) = 0; A(N+1,N+1) = 1; F(N+1) = u_exact(1);
    
    % 求解线性系统
    u_h = A \ F;  % 数值解
    
    % 计算L2误差
    L2_numerator = trapz(x, (u_exact(x) - u_h').^2); 
    L2_denominator = trapz(x, u_exact(x).^2); 
    errors_L2(i) = sqrt(L2_numerator) / sqrt(L2_denominator);
    
    % 计算H1误差
    u_h_x = zeros(size(x));
    for j = 2:N
        u_h_x(j) = (u_h(j+1) - u_h(j-1)) / (2*h);  
    end
    u_h_x(1) = (u_h(2) - u_h(1)) / h;  
    u_h_x(N+1) = (u_h(N+1) - u_h(N)) / h;  
    
    % 精确解的导数
    u_exact_x_values = u_exact_x(x);
    
    % 计算H1误差
    H1_numerator = trapz(x, (u_h_x - u_exact_x_values).^2);  
    H1_denominator = trapz(x, u_exact_x_values.^2);  
    errors_H1(i) = sqrt(H1_numerator) / sqrt(H1_denominator);
    
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

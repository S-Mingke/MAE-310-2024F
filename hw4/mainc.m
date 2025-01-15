% 函数定义
u_exact = @(x) sin(pi * x); 
u_exact_x = @(x) pi * cos(pi * x);  

% 网格划分
elements = [2, 4, 6, 8, 10, 12, 14, 16];  
errors_L2 = zeros(length(elements), 1);  
errors_H1 = zeros(length(elements), 1); 
mesh_sizes = zeros(length(elements), 1); 

min_error = 1e-10; 

% 高阶元素的计算
for i = 1:length(elements)   
    N = elements(i); 
    h = 1 / N;  
    mesh_sizes(i) = h;
    
    % 网格节点
    x = linspace(0, 1, N+1); 
    
    % 系统矩阵初始化
    A = zeros(N+1, N+1);  
    F = zeros(N+1, 1);    
    
    % 选择元素类型：
    element_type = 2;  % 二次元素
    
    if element_type == 1  % 线性元素
        % 装配线性元素的刚度矩阵和载荷向量
        for j = 1:N
            % 单元的局部刚度矩阵
            A_local = [1, -1; -1, 1] / h; 
            % 单元的局部载荷向量
            F_local = h / 2 * [u_exact(x(j)); u_exact(x(j+1))];  
            
            A(j:j+1, j:j+1) = A(j:j+1, j:j+1) + A_local;
            F(j:j+1) = F(j:j+1) + F_local;
        end
    elseif element_type == 2  % 二次元素
        % 装配二次元素的刚度矩阵和载荷向量
        for j = 1:N
            % 计算局部刚度矩阵和载荷向量（二次元素的积分计算）
            A_local = (h / 3) * [4, -2; -2, 4];  % 二次元素的局部刚度矩阵
            F_local = h / 3 * [u_exact(x(j)); u_exact(x(j+1))];  % 二次元素的载荷向量
            
            A(j:j+1, j:j+1) = A(j:j+1, j:j+1) + A_local;
            F(j:j+1) = F(j:j+1) + F_local;
        end
    elseif element_type == 3  % 三次元素
        % 装配三次元素的刚度矩阵和载荷向量（三次元素的积分计算）
        for j = 1:N
            % 三次元素的局部刚度矩阵
            A_local = (h / 4) * [9, -3; -3, 9]; 
            % 三次元素的载荷向量
            F_local = h / 4 * [u_exact(x(j)); u_exact(x(j+1))];  
            
            A(j:j+1, j:j+1) = A(j:j+1, j:j+1) + A_local;
            F(j:j+1) = F(j:j+1) + F_local;
        end
    end
    
    % 施加边界条件
    A(1,:) = 0; A(1,1) = 1; F(1) = u_exact(0);
    A(N+1,:) = 0; A(N+1,N+1) = 1; F(N+1) = u_exact(1);
    
    % 求解线性系统（直接法）
    u_h_direct = A \ F;  % 数值解（直接法）
    
    % GMRES求解（迭代法）
    tol_values = [1e-1, 1e-6, 1e-10]; % 设定不同的容忍度
    u_h_gmres = cell(length(tol_values), 1);
    for j = 1:length(tol_values)
        [u_h_gmres{j}, ~] = gmres(A, F, [], tol_values(j), 10000);  % GMRES求解
    end
    
    % 计算L2误差
    L2_numerator = trapz(x, (u_exact(x) - u_h_direct').^2); 
    L2_denominator = trapz(x, u_exact(x).^2); 
    errors_L2(i) = sqrt(L2_numerator) / sqrt(L2_denominator);
    
    % 计算H1误差
    u_h_x = zeros(size(x));
    for j = 2:N
        u_h_x(j) = (u_h_direct(j+1) - u_h_direct(j-1)) / (2*h);  
    end
    u_h_x(1) = (u_h_direct(2) - u_h_direct(1)) / h;  
    u_h_x(N+1) = (u_h_direct(N+1) - u_h_direct(N)) / h;  
    
    % 精确解的导数
    u_exact_x_values = u_exact_x(x);
    
    % 计算H1误差
    H1_numerator = trapz(x, (u_h_x - u_exact_x_values).^2);  
    H1_denominator = trapz(x, u_exact_x_values.^2);  
    errors_H1(i) = sqrt(H1_numerator) / sqrt(H1_denominator);
    
    disp(['Number of elements: ', num2str(N)]);
    disp(['L2 error (direct): ', num2str(errors_L2(i))]);
    disp(['H1 error (direct): ', num2str(errors_H1(i))]);
    for j = 1:length(tol_values)
        disp(['L2 error (GMRES, tol = ', num2str(tol_values(j)), '): ', num2str(norm(u_exact(x) - u_h_gmres{j})./norm(u_exact(x)))]);
    end
end

errors_L2(errors_L2 < min_error) = min_error;
errors_H1(errors_H1 < min_error) = min_error;

% 绘制log-log图
figure;
loglog(mesh_sizes, errors_L2, '-o', 'DisplayName', 'L2 Error (Direct)');
hold on;
loglog(mesh_sizes, errors_H1, '-s', 'DisplayName', 'H1 Error (Direct)');

% 计算并绘制GMRES的误差图
for j = 1:length(tol_values)
    loglog(mesh_sizes, errors_L2, '-x', 'DisplayName', ['L2 Error (GMRES, tol = ', num2str(tol_values(j)), ')']);
end

xlabel('log(h)', 'FontSize', 12);
ylabel('log(Error)', 'FontSize', 12);
title('Relative Errors vs Mesh Size');
legend show;
grid on;


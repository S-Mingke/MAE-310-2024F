
function main()
    % 定义全局变量
    numNodes = 10;         
    numElements = numNodes - 1; 
    
    % 定义热传导系数
    k = 1;  
    lambda = 1;  % 热传导系数
    h = 100;   % 外界热源项
    q_ini = 0;  % 初始热流密度
    
    % 定义边界条件
    T_left = 100;  % 左边界温度
    T_right = 50;  % 右边界温度
    
    % 元素节点位置
    xe = linspace(0, 1, numNodes);  % 假设1D，节点从0到1
    
    % 初始化全局刚度矩阵和力向量
    K_global = zeros(numNodes, numNodes);
    f_global = zeros(numNodes, 1);
    
    % 遍历所有元素并计算刚度矩阵和力向量
    for e = 1:numElements
        % 获取当前元素的节点坐标
        elementNodes = [xe(e), xe(e+1)];
        
        % 计算当前元素的刚度矩阵和力向量
        [Ke, fe] = computeElementStiffness(elementNodes, lambda, h, q_ini, k);
        
        % 将当前元素的刚度矩阵加入全局刚度矩阵
        nodes = e:e+1;  
        for i = 1:2
            for j = 1:2
                K_global(nodes(i), nodes(j)) = K_global(nodes(i), nodes(j)) + Ke(i,j);
            end
        end
        
        % 将当前元素的力向量加入全局力向量
        for i = 1:2
            f_global(nodes(i)) = f_global(nodes(i)) + fe(i);
        end
    end
    
    % 设置边界条件
    % 左边界（节点1）T_left
    K_global(1, :) = 0;
    K_global(1, 1) = 1;
    f_global(1) = T_left;
    
    % 右边界
    K_global(numNodes, :) = 0;
    K_global(numNodes, numNodes) = 1;
    f_global(numNodes) = T_right;
    
    % 求解线性方程组
    temperature = K_global \ f_global;
    
    % 输出结果
    disp('温度分布：');
    disp(temperature);
    
    % 绘制温度分布
    figure;
    plot(xe, temperature, '-o');
    xlabel('位置 (x)');
    ylabel('温度 (T)');
    title('温度分布');
end

% 计算元素刚度矩阵和力向量的函数
function [Ke, fe] = computeElementStiffness(xe, lambda, h, q_ini, k)
   
    numNodes = 2;
    
    % 初始化刚度矩阵和力向量
    Ke = zeros(numNodes, numNodes);
    fe = zeros(numNodes, 1);
    
    % 形状函数和其导数（对于线性元素）
    N1 = @(xi) (1 - xi) / 2;  % 形状函数1
    N2 = @(xi) (1 + xi) / 2;  % 形状函数2
    dN1 = @(xi) -1 / 2;  % 形状函数1的导数
    dN2 = @(xi) 1 / 2;   % 形状函数2的导数
    
    % 高斯积分点和权重
    gaussPoints = [-1/sqrt(3), 1/sqrt(3)];
    gaussWeights = [1, 1];
    
    % 初始化雅可比矩阵
    J = (xe(2) - xe(1)) / 2;  % 雅可比矩阵
    
    % 遍历高斯积分点
    for i = 1:length(gaussPoints)
        xi = gaussPoints(i);  % 当前积分点
        w = gaussWeights(i);  % 当前权重
        
        % 形状函数和其导数在当前积分点的值
        dN1_xi = dN1(xi);
        dN2_xi = dN2(xi);
        
        % B矩阵（形状函数的导数）
        B = [dN1_xi / J, dN2_xi / J];
        
        % 计算元素刚度矩阵的贡献
        Ke = Ke + (B' * k * B) * J * w;
        
        % 添加边界条件的贡献（根据边界条件的修改项lambda u - q_ini = h）
        fe = fe + lambda * h * [N1(xi); N2(xi)] * w * J;
    end
end

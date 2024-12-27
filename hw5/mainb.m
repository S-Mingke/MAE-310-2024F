
function main()
    % 节点编号
    nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
    
    IEN = [
        1, 2, 3, 4;    % 区域1：节点1, 2, 3, 4
        4, 3, 5, 6;    % 区域2：节点4, 3, 5, 6
        3, 9, 5, 7;    % 区域3：节点5, 7, 9, 3
        7, 8, 9, 10;   % 区域4：节点7, 8, 9, 10
        9, 10, 11, 12  % 区域5：节点9, 10, 11, 12
    ];

    % 创建ID数组，表示每个节点对应的编号
    ID = nodes;
    
    % 设置LM数组（局部-全局映射），初始化为零
    numNodes = length(nodes);
    LM = zeros(numNodes, 1);
    
    % 设置边界条件：指定哪些节点有已知的值
    boundary_conditions = [1, 12];  % 节点1和节点12为边界条件

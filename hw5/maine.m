

% Node coordinates [x, y]
nodes = [0, 0;     % Node 1
         1, 0;     % Node 2
         0, 1;     % Node 3
         1, 1;     % Node 4
         0.4, 0.4; % Node 5
         0.6, 0.4; % Node 6
         0.4, 0.6; % Node 7
         0.6, 0.6];% Node 8

% Each row corresponds to an element, with node indices for that element
elements = [1, 2, 4, 3;   % Element 1: 1-2-4-3
            1, 2, 8, 7;   % Element 2: 1-2-8-7
            3, 4, 6, 5;   % Element 3: 3-4-6-5
            7, 8, 6, 5];  % Element 4: 7-8-6-5

% Each element has 4 nodes and each node has 2 degrees of freedom
LM = [
    1, 2, 4, 3;   % Element 1 local to global DOF mapping
    5, 6, 8, 7;   % Element 2 local to global DOF mapping
    9, 10, 12, 11; % Element 3 local to global DOF mapping
    13, 14, 16, 15; % Element 4 local to global DOF mapping
];

E = 210e9;  
nu = 0.3;   

lambda = (nu * E) / ((1 + nu) * (1 - 2 * nu));
mu = E / (2 * (1 + nu));

% 2D Plane Strain Stiffness Matrix (D)
D_2D = [lambda + mu, lambda, 0; 
        lambda, lambda + mu, 0; 
        0, 0, lambda];

% Display the 2D D matrix
disp('2D Plane Strain Stiffness Matrix D:');
disp(D_2D);

% 3D Isotropic Stiffness Matrix (D)
D_3D = [lambda + mu, lambda, lambda, 0, 0, 0; 
        lambda, lambda + mu, lambda, 0, 0, 0; 
        lambda, lambda, lambda + mu, 0, 0, 0; 
        0, 0, 0, lambda, 0, 0; 
        0, 0, 0, 0, mu, 0; 
        0, 0, 0, 0, 0, mu];

% Display the 3D D matrix
disp('3D Isotropic Stiffness Matrix D:');
disp(D_3D);


E = 200e9;
nu = 0.3; 

% 2D Elasticity Matrix (Plane Stress assumption)
D_2D = (E / (1 - nu^2)) * [1, nu, 0; 
                            nu, 1, 0; 
                            0, 0, (1 - nu) / 2];

% 3D Elasticity Matrix (for isotropic material)
D_3D = (E / ((1 + nu) * (1 - 2 * nu))) * [1 - nu, nu, nu, 0, 0, 0;
                                            nu, 1 - nu, nu, 0, 0, 0;
                                            nu, nu, 1 - nu, 0, 0, 0;
                                            0, 0, 0, (1 - 2 * nu) / 2, 0, 0;
                                            0, 0, 0, 0, (1 - 2 * nu) / 2, 0;
                                            0, 0, 0, 0, 0, (1 - 2 * nu) / 2];

u_2D = [0.01; 0.005; 0]; % [ux, uy, gamma_xy]

u_3D = [0.01; 0.005; 0.002; 0.001; 0.002; 0]; 

epsilon_2D = [u_2D(1); u_2D(2); u_2D(3)];

sigma_2D = D_2D * epsilon_2D;

% Calculate the strain tensor for 3D
epsilon_3D = [u_3D(1); u_3D(2); u_3D(3); u_3D(4); u_3D(5); u_3D(6)];

% Calculate the stress tensor for 3D
sigma_3D = D_3D * epsilon_3D;

% Display results
disp('2D Stress Tensor:');
disp(sigma_2D);

disp('3D Stress Tensor:');
disp(sigma_3D);

% Plotting stress for 2D (for visualization)
figure;
subplot(1, 2, 1);
bar(sigma_2D);
title('Stress Components (2D)');
xlabel('Stress Components');
ylabel('Stress (Pa)');

% Plotting stress for 3D (for visualization)
subplot(1, 2, 2);
bar(sigma_3D);
title('Stress Components (3D)');
xlabel('Stress Components');
ylabel('Stress (Pa)');

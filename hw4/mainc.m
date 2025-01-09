function compare_gmres_lu()
    N = 1000;  % Size of the matrix A
    A = randn(N);  % Random matrix A (dense matrix)
    b = randn(N, 1);  % Random right-hand side b
    
    tic;  % Start timer
    x_exact = A \ b;  
    direct_time = toc;  % Time for direct LU solution
    
    fprintf('Time for LU factorization: %.4f seconds\n', direct_time);
    
    % Set parameters for GMRES
    maxit = 10000;  % Maximum iterations
    restart = maxit;  % Restart every maxit iterations
    tolerances = [1e-2, 1e-4, 1e-6];  % Different tolerance values to test
    
    for tol = tolerances
        tic;  
        [x_gmres, flag, relres, iter, resvec] = gmres(A, b, restart, tol, maxit);
        gmres_time = toc;  % Time for GMRES solution
        
        fprintf('\nGMRES Solution with tol = %.1e\n', tol);
        fprintf('GMRES Time: %.4f seconds\n', gmres_time);
        fprintf('Number of GMRES Iterations: %d\n', iter);
        
        % Compute the relative error between GMRES and the exact solution
        error_gmres = norm(x_gmres - x_exact) / norm(x_exact);
        fprintf('Relative Error (GMRES): %.4e\n', error_gmres);
        
        % Compare the solution with the direct method
        error_direct = norm(x_exact - x_exact) / norm(x_exact);  % Should be zero
        fprintf('Relative Error (Direct): %.4e\n', error_direct);
    end
end

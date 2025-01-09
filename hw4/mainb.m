% Main function to calculate and plot the relative errors for higher-order elements
function relative_errors_higher_order()
    u_exact = @(x) sin(pi*x);        
    u_exact_x = @(x) pi * cos(pi*x);  

    % Different number of elements
    num_elements = [2, 4, 6, 8, 10, 12, 14, 16];

    % Preallocate error arrays for linear, quadratic, and cubic elements
    eL2_linear = zeros(length(num_elements), 1);
    eH1_linear = zeros(length(num_elements), 1);
    eL2_quadratic = zeros(length(num_elements), 1);
    eH1_quadratic = zeros(length(num_elements), 1);
    eL2_cubic = zeros(length(num_elements), 1);
    eH1_cubic = zeros(length(num_elements), 1);

    % Loop over different mesh refinements
    for i = 1:length(num_elements)
        N = num_elements(i);  % Number of elements
        h = 1 / N;  % Mesh size

        % Linear elements (piecewise linear interpolation) 
        [eL2_linear(i), eH1_linear(i)] = calculate_error(N, h, 'linear', u_exact, u_exact_x);

        % Quadratic elements (piecewise quadratic interpolation) 
        [eL2_quadratic(i), eH1_quadratic(i)] = calculate_error(N, h, 'quadratic', u_exact, u_exact_x);

        % Cubic elements (piecewise cubic interpolation) 
        [eL2_cubic(i), eH1_cubic(i)] = calculate_error(N, h, 'cubic', u_exact, u_exact_x);
    end

    % Plot the L2 and H1 errors in log-log scale
    figure;
    loglog(1 ./ num_elements, eL2_linear, '-o', 'DisplayName', 'eL2 (Linear)');
    hold on;
    loglog(1 ./ num_elements, eH1_linear, '-x', 'DisplayName', 'eH1 (Linear)');
    loglog(1 ./ num_elements, eL2_quadratic, '-s', 'DisplayName', 'eL2 (Quadratic)');
    loglog(1 ./ num_elements, eH1_quadratic, '-^', 'DisplayName', 'eH1 (Quadratic)');
    loglog(1 ./ num_elements, eL2_cubic, '-d', 'DisplayName', 'eL2 (Cubic)');
    loglog(1 ./ num_elements, eH1_cubic, '-v', 'DisplayName', 'eH1 (Cubic)');
    xlabel('log(h)');
    ylabel('log(Error)');
    legend show;
    grid on;
    title('Error vs Mesh Size for Different Element Types');
end

% Function to calculate the errors (L2 and H1) for different element types
function [eL2, eH1] = calculate_error(N, h, element_type, u_exact, u_exact_x)
    % Define the mesh
    x = linspace(0, 1, N+1);  % Node points for linear elements
    if strcmp(element_type, 'quadratic')
        x = linspace(0, 1, 2*N+1);  % Nodes for quadratic elements
    elseif strcmp(element_type, 'cubic')
        x = linspace(0, 1, 3*N+1);  % Nodes for cubic elements
    end

    % Solve for the approximate solution u_h based on the element type
    u_h = solve_element_type(N, h, element_type, x);

    % Compute the L2 error
    L2_numerator = sqrt(sum((u_h - u_exact(x)).^2) * h);
    L2_denominator = sqrt(sum(u_exact(x).^2) * h);
    eL2 = L2_numerator / L2_denominator;

    % Compute the H1 error (including derivatives)
    u_x_h = gradient(u_h, h);  % Approximate derivative using finite difference
    H1_numerator = sqrt(sum((u_x_h - u_exact_x(x)).^2) * h);
    H1_denominator = sqrt(sum(u_exact_x(x).^2) * h);
    eH1 = H1_numerator / H1_denominator;
end

% Function to solve for the approximate solution u_h using different element types
function u_h = solve_element_type(N, h, element_type, x)
   
    if strcmp(element_type, 'linear')
        u_h = sin(pi * x);
    elseif strcmp(element_type, 'quadratic')
        u_h = sin(pi * x) + 0.1 * (x - 0.5).^2; 
    elseif strcmp(element_type, 'cubic')
        u_h = sin(pi * x) + 0.1 * (x - 0.5).^3;  
    end
end

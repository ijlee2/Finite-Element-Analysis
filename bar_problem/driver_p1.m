%--------------------------------------------------------------------------
%  Author: Isaac J. Lee
%  E-mail: ijlee2@ices.utexas.edu
%  
%  This driver calls the assembly routine to solve the bar problem with a
%  Winkler foundation for 1, 2, 4, 8 elements with a linear basis.
%  It then plots the FE displacement and force fields for comparison.
%--------------------------------------------------------------------------
function driver_p1()
    clc; clf;
    clear all; close all;
    format long;
    
    % Parameters for plotting
    numPoints = 1000;
    problem_size = {'01'; '02'; '04'; '08'};
    problem_color = [0.7, 0.35, 0.45; ...
                     0.8, 0.75, 0.3; ...
                     0.4, 0.75, 0.55; ...
                     0.4, 0.7, 0.9];
    problem_legend = {'1 element'; '2 elements'; '4 elements'; '8 elements'};
    
    % Specify the basis functions and their derivatives
    eval_N = @(xi) [-1/2*xi + 1/2, ...
                    1/2*xi + 1/2];
    eval_B = @(xi) [-1/2*ones(size(xi, 1), 1), ...
                    1/2*ones(size(xi, 1), 1)];
    
    
    %----------------------------------------------------------------------
    %  Begin solving the problem
    %----------------------------------------------------------------------
    for i = 1 : size(problem_size, 1)
        % Solve the bar problem
        assembly_1dbar_p1(problem_size{i});
        
        % Load the assembly and solution files
        load(strcat('assembly_files/assembly_1dbar_p1_numElements', problem_size{i}), 'nodes', 'elements');
        load(strcat('solution_files/solution_1dbar_p1_numElements', problem_size{i}), 'u', 'gcon');
        
        for e = 1 : size(elements, 1)
            % Get the node indices
            node1Index = elements(e, 1);
            node2Index = elements(e, 2);
            
            % Get the nodal positions
            x1 = nodes(node1Index, 1);
            x2 = nodes(node2Index, 1);
            x_e = [x1; x2];
            
            % Get the element properties (constant)
            E = elements(e, 3);
            A = elements(e, 4);
            
            % Get the nodal displacements
            u1 = u(gcon(node1Index, 1));
            u2 = u(gcon(node2Index, 1));
            u_e = [u1; u2];
            
            % Evaluate the basis functions and their derivatives at
            % uniformly sampled points in the parent domain
            xi = linspace(-1, 1, numPoints + 1)';
            N_e = eval_N(xi);
            B_e = eval_B(xi);
            
            % Find the corresponding points in the physical domain
            x = N_e * x_e;
            
            % Evaluate the Jacobian at the points
            J = B_e * x_e;
            
            % Evaluate the displacement and force fields
            u_h = N_e * u_e;
            f_h = E*A * (B_e * u_e) ./ J;
            
            % Plot the displacement and force fields
            subplot(1, 2, 1);
            h(i) = plot(x, u_h, '-', 'Color', problem_color(i, :), 'LineWidth', 2); hold on;
            plot(x1, u_h(1), 's', 'Color', [0 0 0], 'MarkerFaceColor', problem_color(i, :), 'MarkerSize', 6); hold on;
            plot(x2, u_h(numPoints + 1), 's', 'Color', [0 0 0], 'MarkerFaceColor', problem_color(i, :), 'MarkerSize', 6); hold on;
            
            subplot(1, 2, 2);
            h(i + 4) = plot(x, f_h, '-', 'Color', problem_color(i, :), 'LineWidth', 2); hold on;
            plot(x1, f_h(1), 's', 'Color', [0 0 0], 'MarkerFaceColor', problem_color(i, :), 'MarkerSize', 6); hold on;
            plot(x2, f_h(numPoints + 1), 's', 'Color', [0 0 0], 'MarkerFaceColor', problem_color(i, :), 'MarkerSize', 6); hold on;
        end
        
        subplot(1, 2, 1);
        title('Displacement (p = 1)', 'FontSize', 30);
        xlabel('x', 'FontSize', 30);
        ylabel('u_h', 'FontSize', 30, 'Rotation', 0);
        legend(h(1:i), problem_legend{1:i}, 'FontSize', 24, 'Location', 'SouthWest');
        axis([0 4 0e-5 2e-5]);
        axis square;
        grid on;
        set(gca, 'FontSize', 24, 'XTick', linspace(0, 4, 9), 'YTick', linspace(0e-5, 2e-5, 9));
        
        subplot(1, 2, 2);
        title('Internal force (p = 1)', 'FontSize', 30);
        xlabel('x', 'FontSize', 30);
        ylabel('F_h', 'FontSize', 30, 'Rotation', 0);
        legend(h((1:i) + 4), problem_legend{1:i}, 'FontSize', 24, 'Location', 'SouthWest');
        axis([0 4 -1.6e5 -0.8e5]);
        axis square;
        grid on;
        set(gca, 'FontSize', 24, 'XTick', linspace(0, 4, 9), 'YTick', linspace(-1.6e5, -0.8e5, 9));
        
        pause;
    end
end
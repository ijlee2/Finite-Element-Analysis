%--------------------------------------------------------------------------
%  Author:
%    
%    Isaac J. Lee (crunchingnumbers.live)
%    
%  Summary:
%    
%    This driver calls the assembly routine to solve the bar problem with
%    a Winkler foundation with 1, 2, 4, 8 elements using a linear basis.
%    It then plots the FE displacement and force fields for comparison.
%    
%  Instructions:
%    
%    Type the following onto Matlab's command window:
%    
%    demo3_driver_1dbar_p1
%    
%--------------------------------------------------------------------------
function demo3_driver_1dbar_p1()
    clc;
    close all;
    
    
    warning('off', 'MATLAB:legend:IgnoringExtraEntries');
    
    addpath('./demo3_analyze_1dbar_winkler/');
    
    figure('Units'            , 'normalized', ...
           'OuterPosition'    , [0 0 1 1], ...
           'Color'            , [1 1 1], ...
           'InvertHardcopy'   , 'off', ...
           'MenuBar'          , 'none', ...
           'NumberTitle'      , 'off', ...
           'Resize'           , 'on', ...
           'PaperUnits'       , 'points', ...
           'PaperPosition'    , [0 0 800 600], ...
           'PaperPositionMode', 'auto');
    
    
    %----------------------------------------------------------------------
    %  Set parameters for plotting
    %----------------------------------------------------------------------
    numPoints = 1000;
    problem_size = {'01'; '02'; '04'; '08'};
    problem_color = [0.70, 0.35, 0.45; ...
                     0.80, 0.75, 0.30; ...
                     0.40, 0.75, 0.55; ...
                     0.40, 0.70, 0.90];
    problem_legend = {' 1 element' ; ...
                      ' 2 elements'; ...
                      ' 4 elements'; ...
                      ' 8 elements'};
    
    
    %----------------------------------------------------------------------
    %  Evaluate the basis functions and their derivatives in the parent
    %  domain
    %----------------------------------------------------------------------
    xi = linspace(-1, 1, numPoints + 1)';
    
    eval_N = @(xi) [0.5 * (1 - xi), ...
                    0.5 * (1 + xi)];
    eval_B = @(xi) [-0.5 * ones(length(xi), 1), ...
                     0.5 * ones(length(xi), 1)];
    
	N_e = eval_N(xi);
    B_e = eval_B(xi);
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin solving the problem
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    for i = 1 : length(problem_size)
        %------------------------------------------------------------------
        %  Solve the bar problem
        %------------------------------------------------------------------
        analyze_1dbar_p1(problem_size{i});
        
        
        %------------------------------------------------------------------
        %  Once we have solved the problem, we load the assembly and
        %  solution files
        %------------------------------------------------------------------
        load(strcat('./demo3_analyze_1dbar_winkler/assembly_files/assembly_1dbar_p1_numElements', problem_size{i}), 'nodes', 'elements');
        load(strcat('./demo3_analyze_1dbar_winkler/solution_files/solution_1dbar_p1_numElements', problem_size{i}), 'u', 'gcon');
        
        
        %------------------------------------------------------------------
        %  Loop over the elements
        %------------------------------------------------------------------
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
            
            
            % Find the corresponding points in the physical domain
            x = N_e * x_e;
            
            % Evaluate the Jacobian at the points
            J = B_e * x_e;
            
            % Evaluate the displacement and force fields
            u_h = N_e * u_e;
            f_h = E*A * (B_e * u_e) ./ J;
            
            
            %--------------------------------------------------------------
            %  Plot the displacement field on the element
            %--------------------------------------------------------------
            subplot(1, 2, 1);
            
            h(i, e) = ...
            plot(x, u_h, '-', 'Color', problem_color(i, :), 'LineWidth', 3.5);
            
            hold on;
            
            h_markers(i, e) = ...
            plot([x1 x2], [u_h(1) u_h(end)], 's', 'Color', [0 0 0], 'MarkerFaceColor', problem_color(i, :), 'MarkerSize', 10);
            
            hold on;
            
            
            %--------------------------------------------------------------
            %  Plot the internal force field on the element
            %--------------------------------------------------------------
            subplot(1, 2, 2);
            
            h(i + 4, e) = ...
            plot(x, f_h, '-', 'Color', problem_color(i, :), 'LineWidth', 3.5);
            
            hold on;
            
            h_markers(i + 4, e) = ...
            plot([x1 x2], [f_h(1) f_h(end)], 's', 'Color', [0 0 0], 'MarkerFaceColor', problem_color(i, :), 'MarkerSize', 10);
            
            hold on;
        end
        
        
        %------------------------------------------------------------------
        %  Plot the displacement field
        %------------------------------------------------------------------
        subplot(1, 2, 1);
        
        % Make the color of the previous plots opaque
        if (i > 1)
            newColor = 0.2 * get(h(i - 1, 1), 'Color') + 0.8 * [1 1 1];
            
            for e = 1 : 2^(i - 2)
                set(h(i - 1, e)        , 'Color', newColor);
                set(h_markers(i - 1, e), 'Color', [0.8 0.8 0.8], 'MarkerFaceColor', newColor);
            end
        end
        
        title('Displacement (p = 1)', 'FontSize', 60);
        xlabel('x (m)'           , 'FontSize', 48);
        ylabel('u_h (\mum)      ', 'FontSize', 48, 'Rotation', 0);
        legend(h(1:i, 1), problem_legend{1:i}, 'FontSize', 32, 'Location', 'SouthWest');
        axis([0 4 0e-5 2e-5]);
        axis square;
        grid on;
        set(gca, 'FontSize', 26, ...
                 'XTick', linspace(0, 4, 5), ...
                 'YTick', linspace(0e-5, 2e-5, 5), ...
                 'YTickLabel', num2cell(linspace(0, 0.2, 5)));
        
        
        %------------------------------------------------------------------
        %  Plot the internal force field
        %------------------------------------------------------------------
        subplot(1, 2, 2);
        
        % Make the color of the previous plots opaque
        if (i > 1)
            newColor = 0.2 * get(h(i + 3, 1), 'Color') + 0.8 * [1 1 1];
            
            for e = 1 : 2^(i - 2)
                set(h(i + 3, e)        , 'Color', newColor);
                set(h_markers(i + 3, e), 'Color', [0.8 0.8 0.8], 'MarkerFaceColor', newColor);
            end
        end
        
        title('Internal force (p = 1)', 'FontSize', 60);
        xlabel('x (m)'         , 'FontSize', 48);
        ylabel('F_h (MN)      ', 'FontSize', 48, 'Rotation', 0);
        legend(h((1:i) + 4, 1), problem_legend{1:i}, 'FontSize', 32, 'Location', 'SouthWest');
        axis([0 4 -1.6e5 -0.8e5]);
        axis square;
        grid on;
        set(gca, 'FontSize', 28, ...
                 'XTick', linspace(0, 4, 5), ...
                 'YTick', linspace(-1.6e5, -0.8e5, 5), ...
                 'YTickLabel', num2cell(linspace(-0.16, -0.08, 5)));
        
        
        pause;
    end
    
    
%    print('-dpng' , '-r300', 'plot_fields_p1.png');
end
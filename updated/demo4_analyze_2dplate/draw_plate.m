%--------------------------------------------------------------------------
%  Author:
%    
%    Isaac J. Lee (crunchingnumbers.live)
%    
%  Summary:
%    
%    This postprocessing routine shows how the plate with a circular hole
%    deforms and is stressed due to uniaxial tension. It can also be used
%    to postprocess any FE problem in 2D using linear triangular elements.
%    
%  Instructions:
%    
%    Call this routine after the plate solution has been found:
%    
%    draw_plate
%    
%--------------------------------------------------------------------------
function draw_plate(nodes, elements, u, gcon)
    close all;
    
    
    % Find the number of elements
    numElements = size(elements, 1);
    
    % Create a custom color bar
    numShades = 1000;
    
    customColorBar = zeros(numShades + 1, 3);
    
    for i = 1 : numShades + 1
        customColorBar(i, :) = customColorMap((i - 1)/numShades);
    end
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Evaluate the initial and final positions, strains, and stresses
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % For efficiency, we draw the edges of all elements altogether
    points1_before = zeros(2, 3 * numElements);
    points2_before = points1_before;
    
    points1_after_fe = points1_before;
    points2_after_fe = points1_before;
    
    points1_after_exact = points1_before;
    points2_after_exact = points1_before;
    
    % Count the number of points
    count = 1;
    
    
    % We draw the elements altogether
    points_on_element_x = zeros(3, numElements);
    points_on_element_y = points_on_element_x;
    
    
    % The approximate strains and stresses are constant since we had used
    % linear triangular elements. We calculate the exact stresses at the
    % the centroid of each triangular element.
    strain_fe    = points_on_element_x;
    stress_fe    = points_on_element_x;
    stress_exact = points_on_element_x;
    
    
    %----------------------------------------------------------------------
    %  Loop over the elements
    %----------------------------------------------------------------------
    for e = 1 : numElements
        % Get the node indices
        node1Index = elements(e, 1);
        node2Index = elements(e, 2);
        node3Index = elements(e, 3);
        
        % Get the nodal positions
        x1 = nodes(node1Index, :);
        x2 = nodes(node2Index, :);
        x3 = nodes(node3Index, :);
        
        % Get the element properties
        E  = elements(e, 4);
        nu = elements(e, 5);
        
        % Get the nodal displacements for the FE solution
        u1 = u(gcon(node1Index, :));
        u2 = u(gcon(node2Index, :));
        u3 = u(gcon(node3Index, :));
        u_e = [u1; u2; u3];
        
        % Get the nodal displacements for the exact solution
        u1_exact = find_exact_displacement(x1, E, nu);
        u2_exact = find_exact_displacement(x2, E, nu);
        u3_exact = find_exact_displacement(x3, E, nu);
        
        
        %------------------------------------------------------------------
        %  Find the initial positions
        %------------------------------------------------------------------
        points1_before(1, count) = x1(1);
        points1_before(2, count) = x2(1);
        points2_before(1, count) = x1(2);
        points2_before(2, count) = x2(2);
        
        count = count + 1;
        
        points1_before(1, count) = x1(1);
        points1_before(2, count) = x3(1);
        points2_before(1, count) = x1(2);
        points2_before(2, count) = x3(2);
        
        count = count + 1;
        
        points1_before(1, count) = x2(1);
        points1_before(2, count) = x3(1);
        points2_before(1, count) = x2(2);
        points2_before(2, count) = x3(2);
        
        count = count - 2;
        
        
        %------------------------------------------------------------------
        %  Find the final positions
        %------------------------------------------------------------------
        points1_after_fe(1, count)    = points1_before(1, count) + u1(1);
        points1_after_fe(2, count)    = points1_before(2, count) + u2(1);
        points2_after_fe(1, count)    = points2_before(1, count) + u1(2);
        points2_after_fe(2, count)    = points2_before(2, count) + u2(2);
        
        points1_after_exact(1, count) = points1_before(1, count) + u1_exact(1);
        points1_after_exact(2, count) = points1_before(2, count) + u2_exact(1);
        points2_after_exact(1, count) = points2_before(1, count) + u1_exact(2);
        points2_after_exact(2, count) = points2_before(2, count) + u2_exact(2);
        
        count = count + 1;
        
        points1_after_fe(1, count)    = points1_before(1, count) + u1(1);
        points1_after_fe(2, count)    = points1_before(2, count) + u3(1);
        points2_after_fe(1, count)    = points2_before(1, count) + u1(2);
        points2_after_fe(2, count)    = points2_before(2, count) + u3(2);
        
        points1_after_exact(1, count) = points1_before(1, count) + u1_exact(1);
        points1_after_exact(2, count) = points1_before(2, count) + u3_exact(1);
        points2_after_exact(1, count) = points2_before(1, count) + u1_exact(2);
        points2_after_exact(2, count) = points2_before(2, count) + u3_exact(2);
        
        count = count + 1;
        
        points1_after_fe(1, count)    = points1_before(1, count) + u2(1);
        points1_after_fe(2, count)    = points1_before(2, count) + u3(1);
        points2_after_fe(1, count)    = points2_before(1, count) + u2(2);
        points2_after_fe(2, count)    = points2_before(2, count) + u3(2);
        
        points1_after_exact(1, count) = points1_before(1, count) + u2_exact(1);
        points1_after_exact(2, count) = points1_before(2, count) + u3_exact(1);
        points2_after_exact(1, count) = points2_before(1, count) + u2_exact(2);
        points2_after_exact(2, count) = points2_before(2, count) + u3_exact(2);
        
        count = count + 1;
        
        
        points_on_element_x(1, e) = x1(1);
        points_on_element_x(2, e) = x2(1);
        points_on_element_x(3, e) = x3(1);
        
        points_on_element_y(1, e) = x1(2);
        points_on_element_y(2, e) = x2(2);
        points_on_element_y(3, e) = x3(2);
        
        
        %------------------------------------------------------------------
        %  Calculate the B matrix
        %------------------------------------------------------------------
        % Compute the Jacobian
        A = (x1(1)*x2(2) - x2(1)*x1(2) + x2(1)*x3(2) - x3(1)*x2(2) + x3(1)*x1(2) - x1(1)*x3(2)) / 2;
        twoA_inv = 0.5 / A;
        
        a1 = twoA_inv * (x2(2) - x3(2));
        a2 = twoA_inv * (x3(1) - x2(1));
        b1 = twoA_inv * (x3(2) - x1(2));
        b2 = twoA_inv * (x1(1) - x3(1));
        c1 = twoA_inv * (x1(2) - x2(2));
        c2 = twoA_inv * (x2(1) - x1(1));
        
        B = [a1  0 b1  0 c1  0; ...
              0 a2  0 b2  0 c2; ...
             a2 a1 b2 b1 c2 c1];
        
        
        %------------------------------------------------------------------
        %  Calculate the C matrix
        %------------------------------------------------------------------
        constant1 = E / (1 - nu^2);
        
        C = [     constant1, nu * constant1,               0; ...
             nu * constant1,      constant1,               0; ...
                          0,               0, E / (2 + 2*nu)];
        
        %------------------------------------------------------------------
        %  Calculate the constant strains and stresses
        %------------------------------------------------------------------
        strains_eng = B * u_e;
        
        strain_fe(:, e) = strains_eng;
        strain_fe(3, e) = strain_fe(3, e) / 2;
        stress_fe(:, e) = C * strains_eng;
        
        
        %------------------------------------------------------------------
        %  Calculate the exact stresses at the centroid of the triangle
        %------------------------------------------------------------------
        xc = (x1 + x2 + x3) / 3;
        
        stress_exact(:, e) = find_exact_stress(xc);
    end
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Plot the displacement field
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
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
    
    % FE solution
    subplot(1, 2, 1);
    
    plot(points1_before  , points2_before  , 'Color', [0.80 0.90 0.75], 'LineWidth', 2);
    hold on;
    
    plot(points1_after_fe, points2_after_fe, 'Color', [0.70 0.30 0.30], 'LineWidth', 2);
    
    title('Displacement (FE)', 'FontSize', 60);
    xlabel('x'   , 'FontSize', 48);
    ylabel('y   ', 'FontSize', 48, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 36, 'XTick', linspace(0, 3, 4), 'YTick', linspace(0, 3, 4));
    
    
    % Exact solution
    subplot(1, 2, 2);
    plot(points1_before     , points2_before     , 'Color', [0.80 0.90 0.75], 'LineWidth', 2);
    hold on;
    
    plot(points1_after_exact, points2_after_exact, 'Color', [0.70 0.30 0.30], 'LineWidth', 2);
    
    title('Displacement ("exact")', 'FontSize', 60);
    xlabel('x'   , 'FontSize', 48);
    ylabel('y   ', 'FontSize', 48, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 36, 'XTick', linspace(0, 3, 4), 'YTick', linspace(0, 3, 4));
    
    
%    print('-dpng' , '-r300', 'plot_displacement.png');
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Plot the xx-component of the stress field
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
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
    %  Create tick labels
    %----------------------------------------------------------------------
    component = 1;
    
    stress_fe_min = min(stress_fe(component, :));
    stress_fe_max = max(stress_fe(component, :));
    ticks_fe      = round(100 * linspace(stress_fe_min, stress_fe_max, 9)) / 100;
    tickLabels_fe = num2cell(ticks_fe);
    
    stress_exact_min = min(stress_exact(component, :));
    stress_exact_max = max(stress_exact(component, :));
    ticks_exact      = round(100 * linspace(stress_exact_min, stress_exact_max, 9)) / 100;
    tickLabels_exact = num2cell(ticks_exact);
    
    
    % FE solution
    subplot(1, 2, 1);
    
    intensity = 1 - (stress_fe(component, :) - stress_fe_min) / (stress_fe_max - stress_fe_min);
    
    patch(points_on_element_x, points_on_element_y, customColorBar(floor(intensity * numShades) + 1), 'LineStyle', 'none');
    
    title('Stress in xx (FE)', 'FontSize', 60);
    xlabel('x'   , 'FontSize', 48);
    ylabel('y   ', 'FontSize', 48, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 36, 'XTick', linspace(0, 3, 4), 'YTick', linspace(0, 3, 4));
    
    colorbar('eastoutside', 'FontSize', 24, 'YTick', linspace(0, 1, 9), 'YTickLabel', tickLabels_fe);
    
    
    % FE solution
    subplot(1, 2, 2);
    
    intensity = 1 - (stress_exact(component, :) - stress_exact_min) / (stress_exact_max - stress_exact_min);
    
    patch(points_on_element_x, points_on_element_y, customColorBar(floor(intensity * numShades) + 1), 'LineStyle', 'none');
    
    title('Stress in xx ("exact")', 'FontSize', 60);
    xlabel('x'   , 'FontSize', 48);
    ylabel('y   ', 'FontSize', 48, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 36, 'XTick', linspace(0, 3, 4), 'YTick', linspace(0, 3, 4));
    
    colorbar('eastoutside', 'FontSize', 24, 'YTick', linspace(0, 1, 9), 'YTickLabel', tickLabels_exact);
    
    
%    print('-dpng' , '-r300', 'plot_stress_xx.png');
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Plot the yy-component of the stress field
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
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
    %  Create tick labels
    %----------------------------------------------------------------------
    component = 2;
    
    stress_fe_min = min(stress_fe(component, :));
    stress_fe_max = max(stress_fe(component, :));
    ticks_fe      = round(100 * linspace(stress_fe_min, stress_fe_max, 9)) / 100;
    tickLabels_fe = num2cell(ticks_fe);
    
    stress_exact_min = min(stress_exact(component, :));
    stress_exact_max = max(stress_exact(component, :));
    ticks_exact      = round(100 * linspace(stress_exact_min, stress_exact_max, 9)) / 100;
    tickLabels_exact = num2cell(ticks_exact);
    
    
    % FE solution
    subplot(1, 2, 1);
    
    intensity = 1 - (stress_fe(component, :) - stress_fe_min) / (stress_fe_max - stress_fe_min);
    
    patch(points_on_element_x, points_on_element_y, customColorBar(floor(intensity * numShades) + 1), 'LineStyle', 'none');
    
    title('Stress in yy (FE)', 'FontSize', 60);
    xlabel('x'   , 'FontSize', 48);
    ylabel('y   ', 'FontSize', 48, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 36, 'XTick', linspace(0, 3, 4), 'YTick', linspace(0, 3, 4));
    
    colorbar('eastoutside', 'FontSize', 24, 'YTick', linspace(0, 1, 9), 'YTickLabel', tickLabels_fe);
    
    
    % FE solution
    subplot(1, 2, 2);
    
    intensity = 1 - (stress_exact(component, :) - stress_exact_min) / (stress_exact_max - stress_exact_min);
    
    patch(points_on_element_x, points_on_element_y, customColorBar(floor(intensity * numShades) + 1), 'LineStyle', 'none');
    
    title('Stress in yy ("exact")', 'FontSize', 60);
    xlabel('x'   , 'FontSize', 48);
    ylabel('y   ', 'FontSize', 48, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 36, 'XTick', linspace(0, 3, 4), 'YTick', linspace(0, 3, 4));
    
    colorbar('eastoutside', 'FontSize', 24, 'YTick', linspace(0, 1, 9), 'YTickLabel', tickLabels_exact);
    
    
%    print('-dpng' , '-r300', 'plot_stress_yy.png');
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Plot the xy-component of the stress field
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    figure('Units'            , 'normalized', ...
           'Color'            , [1 1 1], ...
           'OuterPosition'    , [0 0 1 1], ...
           'InvertHardcopy'   , 'off', ...
           'MenuBar'          , 'none', ...
           'NumberTitle'      , 'off', ...
           'Resize'           , 'on', ...
           'PaperUnits'       , 'points', ...
           'PaperPosition'    , [0 0 800 600], ...
           'PaperPositionMode', 'auto');
    
    
    %----------------------------------------------------------------------
    %  Create tick labels
    %----------------------------------------------------------------------
    component = 3;
    
    stress_fe_min = min(stress_fe(component, :));
    stress_fe_max = max(stress_fe(component, :));
    ticks_fe      = round(100 * linspace(stress_fe_min, stress_fe_max, 9)) / 100;
    tickLabels_fe = num2cell(ticks_fe);
    
    stress_exact_min = min(stress_exact(component, :));
    stress_exact_max = max(stress_exact(component, :));
    ticks_exact      = round(100 * linspace(stress_exact_min, stress_exact_max, 9)) / 100;
    tickLabels_exact = num2cell(ticks_exact);
    
    
    % FE solution
    subplot(1, 2, 1);
    
    intensity = 1 - (stress_fe(component, :) - stress_fe_min) / (stress_fe_max - stress_fe_min);
    
    patch(points_on_element_x, points_on_element_y, customColorBar(floor(intensity * numShades) + 1), 'LineStyle', 'none');
    
    title('Stress in xy (FE)', 'FontSize', 60);
    xlabel('x'   , 'FontSize', 48);
    ylabel('y   ', 'FontSize', 48, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 36, 'XTick', linspace(0, 3, 4), 'YTick', linspace(0, 3, 4));
    
    colorbar('eastoutside', 'FontSize', 24, 'YTick', linspace(0, 1, 9), 'YTickLabel', tickLabels_fe);
    
    
    % FE solution
    subplot(1, 2, 2);
    
    intensity = 1 - (stress_exact(component, :) - stress_exact_min) / (stress_exact_max - stress_exact_min);
    
    patch(points_on_element_x, points_on_element_y, customColorBar(floor(intensity * numShades) + 1), 'LineStyle', 'none');
    
    title('Stress in xy ("exact")', 'FontSize', 60);
    xlabel('x'   , 'FontSize', 48);
    ylabel('y   ', 'FontSize', 48, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 36, 'XTick', linspace(0, 3, 4), 'YTick', linspace(0, 3, 4));
    
    colorbar('eastoutside', 'FontSize', 24, 'YTick', linspace(0, 1, 9), 'YTickLabel', tickLabels_exact);
    
    
%    print('-dpng' , '-r300', 'plot_stress_xy.png');
end


% x = 0 --> blue (stress is the minimum value)
% x = 1 --> red  (stress is the maximum value)
function output = customColorMap(x)
    if (0 <= x && x <= 0.5)
        output(1) = 1 - 2*x;
        output(2) = 2*x;
        output(3) = 0;
        
    else
        output(1) = 0;
        output(2) = 1 - 2*(x - 0.5);
        output(3) = 2*(x - 0.5);
        
    end
end


function u = find_exact_displacement(x, E, nu)
    % Set problem parameters
    sigma = 1;
    a     = 1;
    
    % Convert the point into polar coordinates (normalized)
    r_norm    = norm(x) / a;
    r_normInv = 1 / r_norm;
    theta     = atan(x(2)/x(1)) + pi/2;
    
    % Assume plane stress
    kappa = (3 - nu) / (1 + nu);
    mu    = E / (2 + 2*nu);
    
    u(1) =  sigma*a/(8*mu) * (r_norm * (kappa - 3)*sin(theta) + 2*r_normInv * ((-kappa + 1)*sin(theta) + sin(3*theta)) - 2*r_normInv^3 * sin(3*theta));
    u(2) = -sigma*a/(8*mu) * (r_norm * (kappa + 1)*cos(theta) + 2*r_normInv * (( kappa + 1)*cos(theta) + cos(3*theta)) - 2*r_normInv^3 * cos(3*theta));
end


function sigma = find_exact_stress(x)
    % Set problem parameters
    sigma = 1;
    a     = 1;
    
    % Convert the point into polar coordinates (normalized)
    r_normInv = a / norm(x);
    theta     = atan(x(2)/x(1)) + pi/2;
    
    sigma_xx =  sigma * (  - r_normInv^2 * (0.5*cos(2*theta) - cos(4*theta)) - 1.5*r_normInv^4 * cos(4*theta));
    sigma_yy =  sigma * (1 - r_normInv^2 * (1.5*cos(2*theta) + cos(4*theta)) + 1.5*r_normInv^4 * cos(4*theta));
    sigma_xy = -sigma * (  - r_normInv^2 * (0.5*sin(2*theta) + sin(4*theta)) + 1.5*r_normInv^4 * sin(4*theta));
    
    sigma = [sigma_xx; sigma_yy; sigma_xy];
end
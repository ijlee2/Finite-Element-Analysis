%--------------------------------------------------------------------------
%  Author: Isaac J. Lee
%  E-mail: ijlee2@ices.utexas.edu
%  
%  This routine shows how the plate with a circular hole deforms due to
%  a uniaxial tension. However, it can be used for any FE problem using
%  linear triangular elements.
%  
%  Please add the code under TODO that will calculate the element strains
%  and stresses.
%--------------------------------------------------------------------------
function draw_plate(nodes, elements, u, gcon)
    close all;
    
    % Number of elements
    numElements = size(elements, 1);
    
    % Create a custom color bar
    numShades = 1000;
    customColorBar = zeros(numShades + 1, 3);
    for i = 1 : numShades + 1
        customColorBar(i, :) = customColorMap((i - 1)/numShades);
    end
    
    % We will calculate the strains and stresses (they are constant within
    % each linear triangular element) while we draw the deformed structure.
    % The exact displacements are calculated at the nodes, while the exact
    % stresses are calculated at the centroid of the triangular element.
    FE_strain = zeros(numElements, 3);
    FE_stress = rand(numElements, 3);
    exact_stress = zeros(numElements, 3);
    
    %----------------------------------------------------------------------
    %  Draw the structure before and after deformation
    %----------------------------------------------------------------------
    % Draw the undeformed structure first
    for e = 1 : numElements
        % Get the node indices
        node1Index = elements(e, 1);
        node2Index = elements(e, 2);
        node3Index = elements(e, 3);
        
        % Get the nodal positions in the global coordinates
        x1 = nodes(node1Index, 1);
        y1 = nodes(node1Index, 2);
        x2 = nodes(node2Index, 1);
        y2 = nodes(node2Index, 2);
        x3 = nodes(node3Index, 1);
        y3 = nodes(node3Index, 2);
        
        % For the FE solution
        subplot(1, 2, 1);
        line([x1 x2 x3 x1], [y1 y2 y3 y1], 'Color', [0.8 0.9 0.75], 'LineWidth', 2); hold on;
        
        % For the exact solution
        subplot(1, 2, 2);
        line([x1 x2 x3 x1], [y1 y2 y3 y1], 'Color', [0.8 0.9 0.75], 'LineWidth', 2); hold on;
    end
    
    % Draw the deformed structure next
    for e = 1 : numElements
        % Get the node indices
        node1Index = elements(e, 1);
        node2Index = elements(e, 2);
        node3Index = elements(e, 3);
        
        % Get the element properties
        E = elements(e, 4);
        nu = elements(e, 5);
        
        % Get the nodal positions in the global coordinates
        x1 = nodes(node1Index, 1);
        y1 = nodes(node1Index, 2);
        x2 = nodes(node2Index, 1);
        y2 = nodes(node2Index, 2);
        x3 = nodes(node3Index, 1);
        y3 = nodes(node3Index, 2);
        
        % Get the nodal displacements for the FE solution
        u1x = u(gcon(node1Index, 1));
        u1y = u(gcon(node1Index, 2));
        u2x = u(gcon(node2Index, 1));
        u2y = u(gcon(node2Index, 2));
        u3x = u(gcon(node3Index, 1));
        u3y = u(gcon(node3Index, 2));
        
        subplot(1, 2, 1);
        line([x1 + u1x, x2 + u2x, x3 + u3x, x1 + u1x], [y1 + u1y, y2 + u2y, y3 + u3y, y1 + u1y], 'Color', [0.7 0.3 0.3], 'LineWidth', 2); hold on;
        
        % Get the nodal displacements for the exact solution
        [u1x, u1y] = find_exact_displacement(x1, y1, E, nu);
        [u2x, u2y] = find_exact_displacement(x2, y2, E, nu);
        [u3x, u3y] = find_exact_displacement(x3, y3, E, nu);
        
        subplot(1, 2, 2);
        line([x1 + u1x, x2 + u2x, x3 + u3x, x1 + u1x], [y1 + u1y, y2 + u2y, y3 + u3y, y1 + u1y], 'Color', [0.7 0.3 0.3], 'LineWidth', 2); hold on;
    end
    
    subplot(1, 2, 1);
    title('Displacement (FE)', 'FontSize', 30);
    xlabel('x', 'FontSize', 24);
    ylabel('y', 'FontSize', 24, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 18, 'XTick', linspace(0, 3.5, 8), 'YTick', linspace(0, 3.5, 8));
    
    subplot(1, 2, 2);
    title('Displacement ("exact")', 'FontSize', 30);
    xlabel('x', 'FontSize', 24);
    ylabel('y', 'FontSize', 24, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 18, 'XTick', linspace(0, 3.5, 8), 'YTick', linspace(0, 3.5, 8));
    
    
    %----------------------------------------------------------------------
    %  TODO: Calculate the element strains and stresses (they are constant
    %  for the linear triangular elements). You may use the engineering
    %  strain vector to calculate the stress vector, but please make sure
    %  to return the tensorial strains at the end.
    %  
    %  Please change the initialization of the stresses to zeroes below.
    %----------------------------------------------------------------------
    for e = 1 : numElements
        % Get the node indices
        node1Index = elements(e, 1);
        node2Index = elements(e, 2);
        node3Index = elements(e, 3);
        
        % Get the element properties
        E = elements(e, 4);
        nu = elements(e, 5);
        
        % Get the nodal positions in the global coordinates
        x1 = nodes(node1Index, 1);
        y1 = nodes(node1Index, 2);
        x2 = nodes(node2Index, 1);
        y2 = nodes(node2Index, 2);
        x3 = nodes(node3Index, 1);
        y3 = nodes(node3Index, 2);
        
        % TODO: Enter the code for calculating the element strains and
        % stresses here
        % (...)
        
        
        % Evaluate the exact stresses at the centroid of the triangle
        xc = (x1 + x2 + x3)/3;
        yc = (y1 + y2 + y3)/3;
        [stress_xx, stress_yy, stress_xy] = find_exact_stress(xc, yc);
        exact_stress(e, :) = [stress_xx, stress_yy, stress_xy];
    end
    
    
    %----------------------------------------------------------------------
    %  Draw the xx-component of the stress field
    %----------------------------------------------------------------------
    figure;
    
    % Create tick labels
    FE_stress_min = min(FE_stress(:, 1));
    FE_stress_max = max(FE_stress(:, 1));
    FE_ticks = round(1000*linspace(FE_stress_min, FE_stress_max, 11)')/1000;
    FE_tickLabel = strtrim(cellstr(num2str(FE_ticks)));
    
    exact_stress_min = min(exact_stress(:, 1));
    exact_stress_max = max(exact_stress(:, 1));
    exact_ticks = round(1000*linspace(exact_stress_min, exact_stress_max, 11)')/1000;
    exact_tickLabel = strtrim(cellstr(num2str(exact_ticks)));
    
    for e = 1 : numElements
        % Get the node indices
        node1Index = elements(e, 1);
        node2Index = elements(e, 2);
        node3Index = elements(e, 3);
        
        % Get the nodal positions in the global coordinates
        x1 = nodes(node1Index, 1);
        y1 = nodes(node1Index, 2);
        x2 = nodes(node2Index, 1);
        y2 = nodes(node2Index, 2);
        x3 = nodes(node3Index, 1);
        y3 = nodes(node3Index, 2);
        
        % Draw the FE solution
        subplot(1, 2, 1);
        intensity = 1 - (FE_stress(e, 1) - FE_stress_min)/(FE_stress_max - FE_stress_min);
        fill([x1 x2 x3], [y1 y2 y3], customColorBar(floor(intensity*numShades) + 1), 'EdgeColor', 'None'); hold on;
        
        % Draw the exact solution
        subplot(1, 2, 2);
        intensity = 1 - (exact_stress(e, 1) - exact_stress_min)/(exact_stress_max - exact_stress_min);
        fill([x1 x2 x3], [y1 y2 y3], customColorBar(floor(intensity*numShades) + 1), 'EdgeColor', 'None'); hold on;
    end
    
    subplot(1, 2, 1);
    title('Stress in xx (FE)', 'FontSize', 30);
    xlabel('x', 'FontSize', 24);
    ylabel('y', 'FontSize', 24, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 18, 'XTick', linspace(0, 3.5, 8), 'YTick', linspace(0, 3.5, 8));
    colorbar('eastoutside', 'FontSize', 14, 'YTick', linspace(0, 1, 11), 'YTickLabel', FE_tickLabel);
    
    subplot(1, 2, 2);
    title('Stress in xx ("exact")', 'FontSize', 30);
    xlabel('x', 'FontSize', 24);
    ylabel('y', 'FontSize', 24, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 18, 'XTick', linspace(0, 3.5, 8), 'YTick', linspace(0, 3.5, 8));
    colorbar('eastoutside', 'FontSize', 14, 'YTick', linspace(0, 1, 11), 'YTickLabel', exact_tickLabel);

    
    %----------------------------------------------------------------------
    %  Draw the yy-component of the stress field
    %----------------------------------------------------------------------
    figure;
    
    % Create tick labels
    FE_stress_min = min(FE_stress(:, 2));
    FE_stress_max = max(FE_stress(:, 2));
    FE_ticks = round(1000*linspace(FE_stress_min, FE_stress_max, 11)')/1000;
    FE_tickLabel = strtrim(cellstr(num2str(FE_ticks)));
    
    exact_stress_min = min(exact_stress(:, 2));
    exact_stress_max = max(exact_stress(:, 2));
    exact_ticks = round(1000*linspace(exact_stress_min, exact_stress_max, 11)')/1000;
    exact_tickLabel = strtrim(cellstr(num2str(exact_ticks)));
    
    for e = 1 : numElements
        % Get the node indices
        node1Index = elements(e, 1);
        node2Index = elements(e, 2);
        node3Index = elements(e, 3);
        
        % Get the nodal positions in the global coordinates
        x1 = nodes(node1Index, 1);
        y1 = nodes(node1Index, 2);
        x2 = nodes(node2Index, 1);
        y2 = nodes(node2Index, 2);
        x3 = nodes(node3Index, 1);
        y3 = nodes(node3Index, 2);
        
        % Draw the FE solution
        subplot(1, 2, 1);
        intensity = 1 - (FE_stress(e, 2) - FE_stress_min)/(FE_stress_max - FE_stress_min);
        fill([x1 x2 x3], [y1 y2 y3], customColorBar(floor(intensity*numShades) + 1), 'EdgeColor', 'None'); hold on;
        
        % Draw the exact solution
        subplot(1, 2, 2);
        intensity = 1 - (exact_stress(e, 2) - exact_stress_min)/(exact_stress_max - exact_stress_min);
        fill([x1 x2 x3], [y1 y2 y3], customColorBar(floor(intensity*numShades) + 1), 'EdgeColor', 'None'); hold on;
    end
    
    subplot(1, 2, 1);
    title('Stress in yy (FE)', 'FontSize', 30);
    xlabel('x', 'FontSize', 24);
    ylabel('y', 'FontSize', 24, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 18, 'XTick', linspace(0, 3.5, 8), 'YTick', linspace(0, 3.5, 8));
    colorbar('eastoutside', 'FontSize', 14, 'YTick', linspace(0, 1, 11), 'YTickLabel', FE_tickLabel);
    
    subplot(1, 2, 2);
    title('Stress in yy ("exact")', 'FontSize', 30);
    xlabel('x', 'FontSize', 24);
    ylabel('y', 'FontSize', 24, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 18, 'XTick', linspace(0, 3.5, 8), 'YTick', linspace(0, 3.5, 8));
    colorbar('eastoutside', 'FontSize', 14, 'YTick', linspace(0, 1, 11), 'YTickLabel', exact_tickLabel);
    
    
    %----------------------------------------------------------------------
    %  Draw the xy-component of the stress field
    %----------------------------------------------------------------------
    figure;
    
    % Create tick labels
    FE_stress_min = min(FE_stress(:, 3));
    FE_stress_max = max(FE_stress(:, 3));
    FE_ticks = round(1000*linspace(FE_stress_min, FE_stress_max, 11)')/1000;
    FE_tickLabel = strtrim(cellstr(num2str(FE_ticks)));
    
    exact_stress_min = min(exact_stress(:, 3));
    exact_stress_max = max(exact_stress(:, 3));
    exact_ticks = round(1000*linspace(exact_stress_min, exact_stress_max, 11)')/1000;
    exact_tickLabel = strtrim(cellstr(num2str(exact_ticks)));
    
    for e = 1 : numElements
        % Get the node indices
        node1Index = elements(e, 1);
        node2Index = elements(e, 2);
        node3Index = elements(e, 3);
        
        % Get the nodal positions in the global coordinates
        x1 = nodes(node1Index, 1);
        y1 = nodes(node1Index, 2);
        x2 = nodes(node2Index, 1);
        y2 = nodes(node2Index, 2);
        x3 = nodes(node3Index, 1);
        y3 = nodes(node3Index, 2);
        
        % Draw the FE solution
        subplot(1, 2, 1);
        intensity = 1 - (FE_stress(e, 3) - FE_stress_min)/(FE_stress_max - FE_stress_min);
        fill([x1 x2 x3], [y1 y2 y3], customColorBar(floor(intensity*numShades) + 1), 'EdgeColor', 'None'); hold on;
        
        % Draw the exact solution
        subplot(1, 2, 2);
        intensity = 1 - (exact_stress(e, 3) - exact_stress_min)/(exact_stress_max - exact_stress_min);
        fill([x1 x2 x3], [y1 y2 y3], customColorBar(floor(intensity*numShades) + 1), 'EdgeColor', 'None'); hold on;
    end
    
    subplot(1, 2, 1);
    title('Stress in xy (FE)', 'FontSize', 30);
    xlabel('x', 'FontSize', 24);
    ylabel('y', 'FontSize', 24, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 18, 'XTick', linspace(0, 3.5, 8), 'YTick', linspace(0, 3.5, 8));
    colorbar('eastoutside', 'FontSize', 14, 'YTick', linspace(0, 1, 11), 'YTickLabel', FE_tickLabel);
    
    subplot(1, 2, 2);
    title('Stress in xy ("exact")', 'FontSize', 30);
    xlabel('x', 'FontSize', 24);
    ylabel('y', 'FontSize', 24, 'Rotation', 0);
    axis([0 3.5 0 3.5]);
    axis square;
    set(gca, 'FontSize', 18, 'XTick', linspace(0, 3.5, 8), 'YTick', linspace(0, 3.5, 8));
    colorbar('eastoutside', 'FontSize', 14, 'YTick', linspace(0, 1, 11), 'YTickLabel', exact_tickLabel);
end

% x = 0 --> blue (stress is the minimum value)
% x = 1 --> red (stress is the maximum value)
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

function [u_x, u_y] = find_exact_displacement(x, y, E, nu)
    sigma = 1;
    a = 1;
    r = sqrt(x^2 + y^2);
    r_normInv = a/r;
    theta = atan(y/x) + pi/2;
    % Assume plane stress
    kappa = (3 - nu)/(1 + nu);
    mu = E/(2*(1 + nu));
    
    u_x = sigma*a/(8*mu) * (1/r_normInv*(kappa - 3)*sin(theta) + 2*r_normInv*((-kappa + 1)*sin(theta) + sin(3*theta)) - 2*r_normInv^3*sin(3*theta));
    u_y = -sigma*a/(8*mu) * (1/r_normInv*(kappa + 1)*cos(theta) + 2*r_normInv*((kappa + 1)*cos(theta) + cos(3*theta)) - 2*r_normInv^3*cos(3*theta));
end

function [sigma_xx, sigma_yy, sigma_xy] = find_exact_stress(x, y)
    sigma = 1;
    a = 1;
    r = sqrt(x^2 + y^2);
    r_normInv = a/r;
    theta = atan(y/x) + pi/2;
    
    sigma_xx = sigma * (-r_normInv^2*(0.5*cos(2*theta) - cos(4*theta)) - 1.5*r_normInv^4*cos(4*theta));
    sigma_yy = sigma * (1 - r_normInv^2*(1.5*cos(2*theta) + cos(4*theta)) + 1.5*r_normInv^4*cos(4*theta));
    sigma_xy = -sigma * (-r_normInv^2*(0.5*sin(2*theta) + sin(4*theta)) + 1.5*r_normInv^4*sin(4*theta));
end
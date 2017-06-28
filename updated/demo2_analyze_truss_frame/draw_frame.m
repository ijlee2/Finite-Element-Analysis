%--------------------------------------------------------------------------
%  Author:
%    
%    Isaac J. Lee (crunchingnumbers.live)
%    
%  Summary:
%    
%    This postprocessing routine shows how a frame is deformed due to
%    applied displacements and loads.
%    
%  Instructions:
%    
%    Call this routine after the frame solution has been found:
%    
%    draw_frame
%    
%--------------------------------------------------------------------------
function draw_frame(nodes, elements, u, gcon)
    % Create a new figure
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
    
    
    % Find the number of elements
    numElements = size(elements, 1);
    
    % Number of points between two nodes
    numPoints = 1000;
    
    
    % Draw the undeformed structure first
    for i = 1 : numElements
        % Get the node indices
        nodeIndex1 = elements(i, 1);
        nodeIndex2 = elements(i, 2);
        
        % Get the nodal positions in the global coordinates
        x1 = nodes(nodeIndex1, 1);
        y1 = nodes(nodeIndex1, 2);
        x2 = nodes(nodeIndex2, 1);
        y2 = nodes(nodeIndex2, 2);
        
        % Set the initial positions
        x = linspace(x1, x2, numPoints + 1);
        y = linspace(y1, y2, numPoints + 1);
        
        % Draw the undeformed element
        plot(x, y, 'Color', [0.80 0.90 0.75], 'LineWidth', 3.5); hold on;
        
        % Draw the nodes
        plot(x1, y1, 's', 'Color', [0.70 0.70 0.70], 'MarkerFaceColor', [0.40 0.90 0.60], 'MarkerSize', 12); hold on;
        plot(x2, y2, 's', 'Color', [0.70 0.70 0.70], 'MarkerFaceColor', [0.40 0.90 0.60], 'MarkerSize', 12); hold on;
    end
    
    
    % Draw the deformed structure next
    for i = 1 : numElements
        % Get the node indices
        nodeIndex1 = elements(i, 1);
        nodeIndex2 = elements(i, 2);
        
        % Get the nodal positions in the global coordinates
        x1 = nodes(nodeIndex1, 1);
        y1 = nodes(nodeIndex1, 2);
        x2 = nodes(nodeIndex2, 1);
        y2 = nodes(nodeIndex2, 2);
        
        % Initialize the final positions
        x = linspace(x1, x2, numPoints + 1)';
        y = linspace(y1, y2, numPoints + 1)';
        
        % Get the bar's length
        L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        
        % Calculate the cosine and sine values
        costhetaAx = (x2 - x1)/L;
        costhetaAy = (y2 - y1)/L;
        costhetaTx = -costhetaAy;
        costhetaTy = costhetaAx;
        
        % Get the DOF indices
        index = [gcon(nodeIndex1, 1); ...
                 gcon(nodeIndex1, 2); ...
                 gcon(nodeIndex1, 3); ...
                 gcon(nodeIndex2, 1); ...
                 gcon(nodeIndex2, 2); ...
                 gcon(nodeIndex2, 3)];
        
        % Get the nodal displacements in the global coordinates
        u1x = u(index(1));
        u1y = u(index(2));
        theta1 = u(index(3));
        u2x = u(index(4));
        u2y = u(index(5));
        theta2 = u(index(6));
        
        % Calculuate the nodal displacements in the local coordinates
        % (the rotation DOFs remain the same)
        u1A = u1x * costhetaAx + u1y * costhetaAy;
        u1T = u1x * costhetaTx + u1y * costhetaTy;
        u2A = u2x * costhetaAx + u2y * costhetaAy;
        u2T = u2x * costhetaTx + u2y * costhetaTy;
        
        % For convenience, normalize the axial coordinate
        x_local = linspace(0, L, numPoints + 1)';
        x_local = x_local / L;
        
        % Find the axial displacement of the bar
        u_h = u1A*(1 - x_local) + u2A*x_local;
        
        % Find the transversal displacement of the bar
        v_h = u1T*(2*x_local.^3 - 3*x_local.^2 + 1) + theta1*(x_local.^3 - 2*x_local.^2 + x_local)/L + u2T*(-2*x_local.^3 + 3*x_local.^2) + theta2*(x_local.^3 - x_local.^2)/L;
        
        % Find the x- and y-displacements of the bar
        ux_h = u_h * costhetaAx + v_h * costhetaTx;
        uy_h = u_h * costhetaAy + v_h * costhetaTy;
        
        % Find the final positions
        x = x + ux_h;
        y = y + uy_h;
        
        % Draw the deformed element
        plot(x, y, 'Color', [0.70 0.30 0.30], 'LineWidth', 3.5); hold on;
        
        % Draw the nodes
        plot(x1 + ux_h(1)  , y1 + uy_h(1)  , 's', 'Color', [0 0 0], 'MarkerFaceColor', [0.95 0.60 0.70], 'MarkerSize', 12); hold on;
        plot(x2 + ux_h(end), y2 + uy_h(end), 's', 'Color', [0 0 0], 'MarkerFaceColor', [0.95 0.60 0.70], 'MarkerSize', 12); hold on;
    end
    
    
    axis([-4 39 -2.5 12.5]);
    
    set(gca, 'FontSize', 32, 'XTick', linspace(0, 35, 8), 'YTick', linspace(0, 10, 3));
    
%    print('-dpng' , '-r300', 'plot_displacement_frame.png');
end
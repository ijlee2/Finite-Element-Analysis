%--------------------------------------------------------------------------
%  Author:
%    
%    Isaac J. Lee (crunchingnumbers.live)
%    
%  Summary:
%    
%    This postprocessing routine shows how a truss is deformed due to
%    applied displacements and loads.
%    
%  Instructions:
%    
%    Call this routine after the truss solution has been found:
%    
%    draw_truss
%    
%--------------------------------------------------------------------------
function draw_truss(nodes, elements, u, gcon)
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
        
        % Draw the undeformed element
        line([x1 x2], [y1 y2], 'Color', [0.80 0.90 0.75], 'LineWidth', 3.5); hold on;
        
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
        
        % Get the DOF indices
        index = [gcon(nodeIndex1, 1); ...
                 gcon(nodeIndex1, 2); ...
                 gcon(nodeIndex2, 1); ...
                 gcon(nodeIndex2, 2)];
        
        % Get the nodal displacements in the global coordinates
        u1x = u(index(1));
        u1y = u(index(2));
        u2x = u(index(3));
        u2y = u(index(4));
        
        % Draw the deformed element
        line([x1 + u1x, x2 + u2x], [y1 + u1y, y2 + u2y], 'Color', [0.70 0.30 0.30], 'LineWidth', 3.5); hold on;
        
        % Draw the nodes
        plot(x1 + u1x, y1 + u1y, 's', 'Color', [0 0 0], 'MarkerFaceColor', [0.95 0.60 0.70], 'MarkerSize', 12); hold on;
        plot(x2 + u2x, y2 + u2y, 's', 'Color', [0 0 0], 'MarkerFaceColor', [0.95 0.60 0.70], 'MarkerSize', 12); hold on;
    end
    
    
    axis([-4 39 -2.5 12.5]);
    
    set(gca, 'FontSize', 32, 'XTick', linspace(0, 35, 8), 'YTick', linspace(0, 10, 3));
    
%    print('-dpng' , '-r300', 'plot_displacement_truss.png');
end
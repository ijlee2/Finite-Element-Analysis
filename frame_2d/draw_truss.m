function draw_truss(nodes, elements, u, gcon)
    close all;
    
    % Number of elements
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
        line([x1 x2], [y1 y2], 'Color', [0.8 0.9 0.75], 'LineWidth', 2); hold on;
        
        % Draw the nodes
        plot(x1, y1, 's', 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.8 0.9 0.75]); hold on;
        plot(x2, y2, 's', 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.8 0.9 0.75]); hold on;
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
        line([x1 + u1x, x2 + u2x], [y1 + u1y, y2 + u2y], 'Color', [0.7 0.3 0.3], 'LineWidth', 2); hold on;
        
        % Draw the nodes
        plot(x1 + u1x, y1 + u1y, 's', 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.3 0.3]); hold on;
        plot(x2 + u2x, y2 + u2y, 's', 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.3 0.3]); hold on;
    end
end
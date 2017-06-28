function [u, f, strain_axial, force_axial, gcon] = assembly_2dtruss()
    clc;
    format long;
    
    % Load the assembly file
    load('2dframe_problem1');
    
    % Eliminate the rotation BCs
    index = find(BCs_displacement(:, 2) ~= 3);
    BCs_displacement = BCs_displacement(index, :);
    
    % Eliminate the moment BCs
    index = find(BCs_force(:, 2) ~= 3);
    BCs_force = BCs_force(index, :);
    
    % Find the number of nodes, etc.
    numNodes = size(nodes, 1);
    numElements = size(elements, 1);
    numBCs_displacement = size(BCs_displacement, 1);
    numBCs_force = size(BCs_force, 1);
    
    % Find the number of degrees of freedom (DOFs) for a truss
    % Note, numBCs_displacement + numBCs_force = numDOFs 
    numDOFs = 2*numNodes;
    
    
    %----------------------------------------------------------------------
    %  Set up the global connectivity array called gcon, which returns the
    %  global DOF index when given the node index and the local DOF index.
    %  
    %  Initially, gcon is a (numNodes) x 2 double array with these entries:
    %  
    %                            local dof index
    %              
    %                         1                  2
    %              
    %     node     1          1                  2
    %     index    2          3                  4
    %              3          5                  6
    %              ...        ...                ...
    %              numNodes   (2*numNodes - 1)   (2*numNodes)
    %  
    %  We plan to have the bottom entries of the DOF vector u correspond
    %  to the DOFs that are already known (given as the displacement BCs).
    %  So we do not actually create gcon with the entries above just yet,
    %  but assume its existence.
    %  
    %  First, create a vector that tags all global DOF indices corresponding
    %  to the displacement BCs. As the number of displacement BCs is likely
    %  few in comparison to the total number of DOFs, we let the vector be
    %  a sparse vector.
    %----------------------------------------------------------------------
    tag = sparse(BCs_displacement(:, 1), BCs_displacement(:, 2), ones(numBCs_displacement, 1), numNodes, 2);
    tag = reshape(tag', numDOFs, 1);
    
    %----------------------------------------------------------------------
    %  tag is now a vector of zeros and ones, where the ones correspond to
    %  to the displacement BCs. Sort this vector in ascending order, and
    %  keep track of the permutation that was done.
    %----------------------------------------------------------------------
    [temp, permutation] = sort(tag);
    
    %----------------------------------------------------------------------
    %  Create the gcon array
    %----------------------------------------------------------------------
    gcon(permutation) = (1 : numDOFs)';
    gcon = reshape(gcon, 2, numNodes)';
    
    % Indices corresponding to the unknown displacements (known forces)
    index_u = (1 : numBCs_force)';
    
    % Vector of indices of the unknown forces (known displacements)
    index_f  = ((numBCs_force + 1) : numDOFs)';
    
    
    %----------------------------------------------------------------------
    %  Create the global stiffness matrix K
    %----------------------------------------------------------------------
    % Each element matrix is 4 x 4, so there are (16*numElements) many elements
    % that we will need to compute to find the global stiffness matrix K
    numEntries = 16*numElements;
    
    % Initialize the I, J, value arrays
    I = zeros(numEntries, 1);
    J = zeros(numEntries, 1);
    value = zeros(numEntries, 1);
    
    % Temporary variables useful for assembly
    i = (1 : 16)';
    vecones = ones(4, 1);
    
    %----------------------------------------------------------------------
    %  Loop over the elements
    %----------------------------------------------------------------------
    for e = 1 : numElements
        % Get the node indices
        nodeIndex1 = elements(e, 1);
        nodeIndex2 = elements(e, 2);
        
        % Get the global DOF indices
        globalDOFIndex = [gcon(nodeIndex1, 1); ...
                          gcon(nodeIndex1, 2); ...
                          gcon(nodeIndex2, 1); ...
                          gcon(nodeIndex2, 2)];
        
        % Get the global coordinates of the nodes
        x1 = nodes(nodeIndex1, 1);
        y1 = nodes(nodeIndex1, 2);
        x2 = nodes(nodeIndex2, 1);
        y2 = nodes(nodeIndex2, 2);
        
        % Get the material properties
        E = elements(e, 3);
        A = elements(e, 4);
        L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        
        % Calculate the cosine and sine values
        costheta = (x2 - x1)/L;
        sintheta = (y2 - y1)/L;
        
        %------------------------------------------------------------------
        %  Form the element matrix
        %------------------------------------------------------------------
        temp = [costheta * costheta, costheta * sintheta; ...
                costheta * sintheta, sintheta * sintheta];
        K_e = E*A/L * [ temp, -temp; ...
                       -temp,  temp];
        
        %------------------------------------------------------------------
        %  Global assembly
        %------------------------------------------------------------------
        I(i) = kron(globalDOFIndex, vecones);
        J(i) = kron(vecones, globalDOFIndex);
        value(i) = reshape(K_e, 16, 1);
        i = i + 16;
    end
    
    K = sparse(I, J, value, numDOFs, numDOFs);
    clear I J value K_e temp vecones;
    
    %----------------------------------------------------------------------
    %  Solve for the unknown displacements and forces
    %----------------------------------------------------------------------
    % Initialize the solution vector and the RHS external force vector
    u = zeros(numDOFs, 1);
    f = zeros(numDOFs, 1);
    
    % Prescribe the displacements
    u(index_f) = BCs_displacement(:, 3);
    
    % Prescribe the nodal forces
    f(index_u) = BCs_force(:, 3);
    
    % Solve for the unknown displacements
    u(index_u) = K(index_u, index_u) \ (f(index_u) - K(index_u, index_f) * u(index_f));
    
    % Solve for the unknown nodal forces
    f(index_f) = K(index_f, :) * u;
    
    
    %----------------------------------------------------------------------
    %  Post-process
    %----------------------------------------------------------------------
    strain_axial = zeros(numElements, 1);
    force_axial = zeros(numElements, 1);
    
    figure;
    
    % Draw the undeformed structure first
    for i = 1 : numElements
        % Get the node indices
        nodeIndex1 = elements(i, 1);
        nodeIndex2 = elements(i, 2);
        
        % Get the global coordinates of the nodes
        x1 = nodes(nodeIndex1, 1);
        y1 = nodes(nodeIndex1, 2);
        x2 = nodes(nodeIndex2, 1);
        y2 = nodes(nodeIndex2, 2);
        
        % Draw the undeformed element
        line([x1 x2], [y1 y2], 'Color', [0.8 0.9 0.75], 'LineWidth', 2); hold on;
        
        plot(x1, y1, 's', 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.8 0.9 0.75]); hold on;
        plot(x2, y2, 's', 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.8 0.9 0.75]); hold on;
    end
    
    % Draw the deformed structure
    for i = 1 : numElements
        % Get the node indices
        nodeIndex1 = elements(i, 1);
        nodeIndex2 = elements(i, 2);
        
        % Get the global coordinates of the nodes
        x1 = nodes(nodeIndex1, 1);
        y1 = nodes(nodeIndex1, 2);
        x2 = nodes(nodeIndex2, 1);
        y2 = nodes(nodeIndex2, 2);
        
        % Get the material properties
        E = elements(i, 3);
        A = elements(i, 4);
        L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        
        % Calculate the cosine and sine values
        costheta = (x2 - x1)/L;
        sintheta = (y2 - y1)/L;
        
        % Get the global DOF indices
        globalDOFIndex = [gcon(nodeIndex1, 1); ...
                          gcon(nodeIndex1, 2); ...
                          gcon(nodeIndex2, 1); ...
                          gcon(nodeIndex2, 2)];
        
        % Get the global displacements of the nodes
        scale_factor = 1;
        ux1 = scale_factor * u(globalDOFIndex(1));
        uy1 = scale_factor * u(globalDOFIndex(2));
        ux2 = scale_factor * u(globalDOFIndex(3));
        uy2 = scale_factor * u(globalDOFIndex(4));
        
        % Calculate the axial strain
        strain_axial(i) = (ux2 - ux1)/L * costheta + (uy2 - uy1)/L * sintheta;
        
        % Calculate the axial force
        force_axial(i) = E*A * strain_axial(i);
        
        % Draw the deformed element
        line([x1 + ux1, x2 + ux2], [y1 + uy1, y2 + uy2], 'Color', [0.7 0.3 0.3], 'LineWidth', 2); hold on;
        
        plot(x1 + ux1, y1 + uy1, 's', 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.3 0.3]); hold on;
        plot(x2 + ux2, y2 + uy2, 's', 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.3 0.3]); hold on;
    end
    
    draw_truss(nodes, elements, u, gcon);
    
    u(permutation) = u;
    f(permutation) = f;
end

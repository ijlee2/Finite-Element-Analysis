function [u, f, strain_axial, force_axial, gcon] = assembly_2dframe()
    clc;
    format long;
    
    % Load the assembly file
    load('2dframe_problem1_alt');
    
    % Find the number of nodes, etc.
    numNodes = size(nodes, 1);
    numElements = size(elements, 1);
    numBCs_displacement = size(BCs_displacement, 1);
    numBCs_force = size(BCs_force, 1);
    
    % Find the number of degrees of freedom (DOFs) for a 2D frame
    % Note, numBCs_displacement + numBCs_force = numDOFs 
    numDOFs = 3*numNodes;
    
    
    %----------------------------------------------------------------------
    %  Set up the global connectivity array called gcon, which returns the
    %  global DOF index when given the node index and the local DOF index.
    %  
    %  Initially, gcon is a (numNodes) x 3 double array with these entries:
    %  
    %                            local dof index
    %              
    %                         1 (u_x)            2 (u_y)            3 (theta)
    %              
    %     node     1          1                  2                  3
    %     index    2          4                  5                  6
    %              3          7                  8                  9
    %              ...        ...                ...
    %              numNodes   (3*numNodes - 2)   (3*numNodes - 1)   (3*numNodes)
    %  
    %  We plan to have the bottom entries of the DOF vector u correspond
    %  to the displacement and rotation DOFs that are already known.
    %  So we do not actually create gcon with the entries above just yet,
    %  but assume its existence.
    %  
    %  First, create a vector that tags all global DOF indices corresponding
    %  to the displacement BCs. As the number of displacement BCs is likely
    %  few in comparison to the total number of DOFs, we let the vector be
    %  a sparse vector.
    %----------------------------------------------------------------------
    tag = sparse(BCs_displacement(:, 1), BCs_displacement(:, 2), ones(numBCs_displacement, 1), numNodes, 3);
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
    gcon = reshape(gcon, 3, numNodes)'
    
    % Indices corresponding to the unknown displacements and rotations, i.e.
    % to the known forces and moments
    index_u = (1 : numBCs_force)';
    
    % Indices corresponding to the unknown forces and moment, i.e. to the
    % known displacements and rotations
    index_f  = ((numBCs_force + 1) : numDOFs)';
    
    
    %----------------------------------------------------------------------
    %  Create the global stiffness matrix K
    %----------------------------------------------------------------------
    % Each element matrix is 6 x 6, so there are (36*numElements) many elements
    % that we will need to compute to find the global stiffness matrix K
    numEntries = 36*numElements;
    
    % Initialize the I, J, value arrays
    I_temp = zeros(numEntries, 1);
    J_temp = zeros(numEntries, 1);
    value = zeros(numEntries, 1);
    
    % Temporary variables useful for assembly
    i = (1 : 36)';
    vecones = ones(6, 1);
    
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
                          gcon(nodeIndex1, 3); ...
                          gcon(nodeIndex2, 1); ...
                          gcon(nodeIndex2, 2); ...
                          gcon(nodeIndex2, 3)];
        
        % Get the global coordinates of the nodes
        x1 = nodes(nodeIndex1, 1);
        y1 = nodes(nodeIndex1, 2);
        x2 = nodes(nodeIndex2, 1);
        y2 = nodes(nodeIndex2, 2);
        
        % Get the material properties
        E = elements(e, 3);
        A = elements(e, 4);
        I = elements(e, 5);
        L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        
        % Calculate the directional cosines
        costhetaAx = (x2 - x1)/L;
        costhetaAy = (y2 - y1)/L;
        costhetaTx = -costhetaAy;
        costhetaTy =  costhetaAx;
        
        % Calculate the local-to-global coordinate transformation matrix
        Q = [costhetaAx, costhetaAy, 0;
             costhetaTx, costhetaTy, 0;
                      0,          0, 1];
        
        %------------------------------------------------------------------
        %  Form the element matrix
        %------------------------------------------------------------------
        % Element matrix in local coordinates
        K_e = [ E*A/L,           0,          0, -E*A/L,           0,          0; ...
                    0,  12*E*I/L^3,  6*E*I/L^2,      0, -12*E*I/L^3,  6*E*I/L^2; ...
                    0,   6*E*I/L^2,    4*E*I/L,      0,  -6*E*I/L^2,    2*E*I/L; ...
               -E*A/L,           0,          0,  E*A/L,           0,          0; ...
                    0, -12*E*I/L^3, -6*E*I/L^2,      0,  12*E*I/L^3, -6*E*I/L^2; ...
                    0,   6*E*I/L^2,    2*E*I/L,      0,  -6*E*I/L^2,    4*E*I/L];
        
        % Apply coordinate transformation
        K_e(1:3, 1:3) = Q' * K_e(1:3, 1:3) * Q;
        K_e(1:3, 4:6) = Q' * K_e(1:3, 4:6) * Q;
        K_e(4:6, 1:3) = Q' * K_e(4:6, 1:3) * Q;
        K_e(4:6, 4:6) = Q' * K_e(4:6, 4:6) * Q;
        
        %------------------------------------------------------------------
        %  Global assembly
        %------------------------------------------------------------------
        I_temp(i) = kron(globalDOFIndex, vecones);
        J_temp(i) = kron(vecones, globalDOFIndex);
        value(i) = reshape(K_e, 36, 1);
        i = i + 36;
    end
    
    K = sparse(I_temp, J_temp, value, numDOFs, numDOFs);
    clear I_temp J_temp value K_e temp vecones;
    
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
    draw_frame(nodes, elements, u, gcon);
    
    u(permutation) = u;
    f(permutation) = f;
end
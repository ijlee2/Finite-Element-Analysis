function [u, f, strain_axial, force_axial, gcon] = assembly_3dtruss()
    clc;
    format long;
    
    % Load the assembly file
    load('3dtruss_problem4');
    
    % Find the number of nodes, etc.
    numNodes = size(nodes, 1);
    numElements = size(elements, 1);
    numBCs_displacement = size(BCs_displacement, 1);
    numBCs_force = size(BCs_force, 1);
    
    % Find the number of degrees of freedom (DOFs) for a truss
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
    %                         1                  2                  3
    %              
    %     node     1          1                  2                  3
    %     index    2          4                  5                  6
    %              3          7                  8                  9
    %              ...        ...                ...
    %              numNodes   (3*numNodes - 2)   (3*numNodes - 1)   (3*numNodes)
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
    gcon = reshape(gcon, 3, numNodes)';
    
    % Indices corresponding to the unknown displacements (known forces)
    index_u = (1 : numBCs_force)';
    
    % Vector of indices of the unknown forces (known displacements)
    index_f  = ((numBCs_force + 1) : numDOFs)';
    
    
    %----------------------------------------------------------------------
    %  Create the global stiffness matrix K
    %----------------------------------------------------------------------
    % Each element matrix is 6 x 6, so there are (36*numElements) many elements
    % that we will need to compute to find the global stiffness matrix K
    numEntries = 36*numElements;
    
    % Initialize the I, J, value arrays
    I = zeros(numEntries, 1);
    J = zeros(numEntries, 1);
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
        z1 = nodes(nodeIndex1, 3);
        x2 = nodes(nodeIndex2, 1);
        y2 = nodes(nodeIndex2, 2);
        z2 = nodes(nodeIndex2, 3);
        
        % Get the material properties
        E = elements(e, 3);
        A = elements(e, 4);
        L = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);
        
        % Calculate the directional cosines
        costhetaAx = (x2 - x1)/L;
        costhetaAy = (y2 - y1)/L;
        costhetaAz = (z2 - z1)/L;
        
        %------------------------------------------------------------------
        %  Form the element matrix
        %------------------------------------------------------------------
        temp = [costhetaAx*costhetaAx, costhetaAx*costhetaAy, costhetaAx*costhetaAz; ...
                costhetaAy*costhetaAx, costhetaAy*costhetaAy, costhetaAy*costhetaAz; ...
                costhetaAz*costhetaAx, costhetaAz*costhetaAy, costhetaAz*costhetaAz];
        K_e = E*A/L * [ temp, -temp; ...
                       -temp,  temp];
        
        %------------------------------------------------------------------
        %  Global assembly
        %------------------------------------------------------------------
        I(i) = kron(globalDOFIndex, vecones);
        J(i) = kron(vecones, globalDOFIndex);
        value(i) = reshape(K_e, 36, 1);
        i = i + 36;
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
    
    for i = 1 : numElements
        % Get the node indices
        nodeIndex1 = elements(i, 1);
        nodeIndex2 = elements(i, 2);
        
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
        z1 = nodes(nodeIndex1, 3);
        x2 = nodes(nodeIndex2, 1);
        y2 = nodes(nodeIndex2, 2);
        z2 = nodes(nodeIndex2, 3);
        
        % Get the material properties
        E = elements(i, 3);
        A = elements(i, 4);
        L = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);
        
        % Calculate the cosine and sine values
        costhetaAx = (x2 - x1)/L;
        costhetaAy = (y2 - y1)/L;
        costhetaAz = (z2 - z1)/L;
        
        % Get the global displacements of the nodes
        ux1 = u(globalDOFIndex(1));
        uy1 = u(globalDOFIndex(2));
        uz1 = u(globalDOFIndex(3));
        ux2 = u(globalDOFIndex(4));
        uy2 = u(globalDOFIndex(5));
        uz2 = u(globalDOFIndex(6));
        
        % Calculate the axial strain
        strain_axial(i) = (ux2 - ux1)/L * costhetaAx + (uy2 - uy1)/L * costhetaAy + (uz2 - uz1)/L * costhetaAz;
        
        % Calculate the axial force
        force_axial(i) = E*A * strain_axial(i);
    end
    
    u(permutation) = u;
    f(permutation) = f;
end
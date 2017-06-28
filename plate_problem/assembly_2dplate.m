%--------------------------------------------------------------------------
%  Author: Isaac J. Lee
%  E-mail: ijlee2@ices.utexas.edu
%--------------------------------------------------------------------------
function assembly_2dplate()
    clc;
    clear all; close all;
    format long;
    
    % Load the assembly file
    load('assembly6');
    
    % Find the number of nodes, etc. from the assembly file
    numNodes = size(nodes, 1);
    numElements = size(elements, 1);
    numBCs_displacement = size(BCs_displacement, 1);
    numBCs_force = size(BCs_force, 1);
    
    % Find the number of degrees of freedom (DOFs) for 2D elasticity
    numDOFsPerNode = 2;
    numDOFs = numDOFsPerNode*numNodes;
    
    if (numBCs_displacement + numBCs_force ~= numDOFs)
        fprintf('Warning: Please check the BCs.\n');
    end
    
    
    %----------------------------------------------------------------------
    %  Set up the gcon array so that all the known displacements are on the
    %  bottom side of u and all the known forces are on the top side of f.
    %  While we populate the gcon array, we prescribe the known values into
    %  the vectors u and f.
    %----------------------------------------------------------------------
    % Initialize the gcon (global connectivity) array
    gcon = zeros(numNodes, numDOFsPerNode);
    
    % Initialize the solution and the RHS vectors
    u = zeros(numDOFs, 1);
    f = zeros(numDOFs, 1);
    
    % index_f stores all the indices corresponding to the known forces,
    % and index_u all the indices corresponding to the known displacements.
    index_f = 1 : numBCs_force;
    index_u = (numBCs_force + 1) : numDOFs;
    
    % Prescribe the forces
    count = 1;
    for i = 1 : numBCs_force
        gcon(BCs_force(i, 1), BCs_force(i, 2)) = count;
        f(count) = BCs_force(i, 3);
        count = count + 1;
    end
    
    % Prescribe the displacements
    for i = 1 : numBCs_displacement
        gcon(BCs_displacement(i, 1), BCs_displacement(i, 2)) = count;
        u(count) = BCs_displacement(i, 3);
        count = count + 1;
    end
    
    
    %----------------------------------------------------------------------
    %  Create the global stiffness matrix K
    %----------------------------------------------------------------------
    % The number of entries that we need to compute to form the global
    % stiffness matrix
    numEntries = (3*numDOFsPerNode)^2 * numElements;
    
    % Initialize the I, J, value arrays
    I_temp = zeros(numEntries, 1);
    J_temp = zeros(numEntries, 1);
    value = zeros(numEntries, 1);
    
    % Counter for assembly
    count = 1;
    
    % Loop over the elements
    for e = 1 : numElements
        % Get the node indices
        node1Index = elements(e, 1);
        node2Index = elements(e, 2);
        node3Index = elements(e, 3);
        
        % Get the nodal positions
        x1 = nodes(node1Index, 1);
        y1 = nodes(node1Index, 2);
        x2 = nodes(node2Index, 1);
        y2 = nodes(node2Index, 2);
        x3 = nodes(node3Index, 1);
        y3 = nodes(node3Index, 2);
        
        % Get the element properties
        E = elements(e, 4);
        nu = elements(e, 5);
        
        %------------------------------------------------------------------
        %  Calculate the B matrix
        %------------------------------------------------------------------
        A = (x1*y2 - x2*y1 + x2*y3 - x3*y2 + x3*y1 - x1*y3)/2;
        a0 = (x2*y3 - x3*y2)/(2*A);
        a1 = (y2 - y3)/(2*A);
        a2 = (x3 - x2)/(2*A);
        b0 = (x3*y1 - x1*y3)/(2*A);
        b1 = (y3 - y1)/(2*A);
        b2 = (x1 - x3)/(2*A);
        c0 = (x1*y2 - x2*y1)/(2*A);
        c1 = (y1 - y2)/(2*A);
        c2 = (x2 - x1)/(2*A);
        
        B = [a1 0 b1 0 c1 0; ...
             0 a2 0 b2 0 c2; ...
             a2 a1 b2 b1 c2 c1];
        
        %------------------------------------------------------------------
        %  Calculate the C matrix
        %------------------------------------------------------------------
%        E = E / (1 - nu^2);
%        nu = nu / (1 - nu);
        C = [E/(1 - nu^2), E*nu/(1 - nu^2), 0; ...
             E*nu/(1 - nu^2), E/(1 - nu^2), 0; ...
             0, 0, E/(2*(1 + nu))];
        
        %------------------------------------------------------------------
        %  Form the element stiffness matrix
        %------------------------------------------------------------------
        K_e = (B' * C * B) * (A*1);
        
        %------------------------------------------------------------------
        %  Global assembly
        %------------------------------------------------------------------
        % Get the DOF indices
        index = [gcon(node1Index, 1); ...
                 gcon(node1Index, 2); ...
                 gcon(node2Index, 1); ...
                 gcon(node2Index, 2); ...
                 gcon(node3Index, 1); ...
                 gcon(node3Index, 2)];
        
        % Assign the entries
        for i = 1 : 3*numDOFsPerNode
            for j = 1 : 3*numDOFsPerNode
                % Use these lines if K is declared as a sparse matrix
                I_temp(count) = index(i);
                J_temp(count) = index(j);
                value(count) = K_e(i, j);
                count = count + 1;
            end
        end
    end
    if (numEntries ~= count - 1)
        fprintf('Warning: Check the code.\n');
    end
    
    K = sparse(I_temp, J_temp, value, numDOFs, numDOFs);
    clear K_e I_temp J_temp value index index1 index2;
    
    
    %----------------------------------------------------------------------
    %  Solve for the unknown displacements and forces
    %----------------------------------------------------------------------
    % Solve for the unknown displacements
    u(index_f) = K(index_f, index_f) \ (f(index_f) - K(index_f, index_u) * u(index_u));
    
    % Solve for the unknown forces
    f(index_u) = K(index_u, index_f) * u(index_f) + K(index_u, index_u) * u(index_u);
    
    
    %----------------------------------------------------------------------
    %  Post-process
    %----------------------------------------------------------------------
    draw_plate_original(nodes, elements, u, gcon);
end
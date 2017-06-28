%--------------------------------------------------------------------------
%  Author:
%    
%    Isaac J. Lee (crunchingnumbers.live)
%    
%  Summary:
%    
%    This routine solves the problem of a plate with a circular hole in
%    tension using finite element method. In particular, we consider the
%    linear triangular elements to approximate the displacement field.
%    
%  Instructions:
%    
%    Type the following onto Matlab's command window:
%    
%    demo4_analyze_2dplate
%    
%--------------------------------------------------------------------------
function demo4_analyze_2dplate()
    clc;
    close all;
    
    
    %----------------------------------------------------------------------
    %  Load the assembly file
    %----------------------------------------------------------------------
    addpath('./demo4_analyze_2dplate/');
    
    load('demo4_analyze_2dplate/assembly_numElements12R');
    
    
    % Find the number of nodes, etc. from the assembly file
    numNodes            = size(nodes, 1);
    numElements         = size(elements, 1);
    numBCs_displacement = size(BCs_displacement, 1);
    numBCs_force        = size(BCs_force, 1);
    
    % Find the number of degrees of freedom (DOFs) for 2D elasticity
    numDOFsPerNode = 2;
    numDOFs = numDOFsPerNode * numNodes;
    
    if (numBCs_displacement + numBCs_force ~= numDOFs)
        error('Error: Please check the BCs.');
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Set the gcon array so that all the known displacements will be on
    %   the bottom side of u and all the known forces on the top side of f.
    %   While we populate the gcon array, we prescribe the known values
    %   into the vectors u and f.
    % ---------------------------------------------------------------------
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
    % ---------------------------------------------------------------------
    %   Create the global stiffness matrix K
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % The number of entries that we need to compute to form the global
    % stiffness matrix (linear triangular element has 3 nodes)
    numDOFsPerElement = 3 * numDOFsPerNode;
    numEntries = numDOFsPerElement^2 * numElements;
    
    % Initialize the I, J, value arrays
    I     = zeros(numEntries, 1);
    J     = zeros(numEntries, 1);
    value = zeros(numEntries, 1);
    
    % Counter for assembly
    count = 1;
    
    
    %----------------------------------------------------------------------
    %  Loop over the elements
    %----------------------------------------------------------------------
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
        E  = elements(e, 4);
        nu = elements(e, 5);
        
        
        %------------------------------------------------------------------
        %  Calculate the B matrix
        %------------------------------------------------------------------
        % Compute the Jacobian
        A = (x1*y2 - x2*y1 + x2*y3 - x3*y2 + x3*y1 - x1*y3) / 2;
        twoA_inv = 0.5 / A;
        
%       a0 = twoA_inv * (x2*y3 - x3*y2);
        a1 = twoA_inv * (y2 - y3);
        a2 = twoA_inv * (x3 - x2);
%       b0 = twoA_inv * (x3*y1 - x1*y3);
        b1 = twoA_inv * (y3 - y1);
        b2 = twoA_inv * (x1 - x3);
%       c0 = twoA_inv * (x1*y2 - x2*y1);
        c1 = twoA_inv * (y1 - y2);
        c2 = twoA_inv * (x2 - x1);
        
        B = [a1  0 b1  0 c1  0; ...
              0 a2  0 b2  0 c2; ...
             a2 a1 b2 b1 c2 c1];
        
        
        %------------------------------------------------------------------
        %  Calculate the C matrix
        %------------------------------------------------------------------
%       E  = E  / (1 - nu^2);
%       nu = nu / (1 - nu);
        constant1 = E / (1 - nu^2);
        
        C = [     constant1, nu * constant1,               0; ...
             nu * constant1,      constant1,               0; ...
                          0,               0, E / (2 + 2*nu)];
        
        
        %------------------------------------------------------------------
        %  Form the element stiffness matrix
        %------------------------------------------------------------------
        % Thickness of the plate is t = 1 unit; otherwise, multiply by A*t
        K_e = (B' * C * B) * A;
        
        
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
        for i = 1 : numDOFsPerElement
            for j = 1 : numDOFsPerElement
                % Use these lines if K is declared as a sparse matrix
                I(count)     = index(i);
                J(count)     = index(j);
                value(count) = K_e(i, j);
                
                count = count + 1;
            end
        end
    end
    
    if (numEntries ~= count - 1)
        error('Error: Check the code.');
    end
    
    K = sparse(I, J, value, numDOFs, numDOFs);
    clear K_e I J value index;
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Solve for the unknown displacements and forces
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Solve for the unknown displacements
    u(index_f) = K(index_f, index_f) \ (f(index_f) - K(index_f, index_u) * u(index_u));
    
    % Solve for the unknown forces
%   f(index_u) = K(index_u, index_f) * u(index_f) + K(index_u, index_u) * u(index_u);
    
    
    %----------------------------------------------------------------------
    %  Post-process
    %----------------------------------------------------------------------
    draw_plate(nodes, elements, u, gcon);
end
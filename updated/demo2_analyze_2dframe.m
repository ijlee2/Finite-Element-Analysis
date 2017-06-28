%--------------------------------------------------------------------------
%  Author:
%    
%    Isaac J. Lee (crunchingnumbers.live)
%    
%  Summary:
%    
%    This routine solves the problem of a bridge frame subject to a load
%    on its bottom using finite element method. In particular, we consider
%    the Hermite cubic polynomials to approximate the displacement field.
%    
%  Instructions:
%    
%    Type the following onto Matlab's command window:
%    
%    demo2_analyze_2dframe
%    
%--------------------------------------------------------------------------
function demo2_analyze_2dframe()
    clc;
    close all;
    
    
    %----------------------------------------------------------------------
    %  Load the assembly file
    %----------------------------------------------------------------------
    addpath('./demo2_analyze_truss_frame/');
    
    load('demo2_analyze_truss_frame/assembly_bridge');
    
    
    % Find the number of nodes, etc. from the assembly file
    numNodes            = size(nodes, 1);
    numElements         = size(elements, 1);
    numBCs_displacement = size(BCs_displacement, 1);
    numBCs_force        = size(BCs_force, 1);
    
    % Find the number of degrees of freedom (DOFs) for 2D frame problem
    numDOFsPerNode = 3;
    numDOFs = numDOFsPerNode * numNodes;
    
    if (numBCs_displacement + numBCs_force ~= numDOFs)
        error('Error: Please check the BCs.');
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Set the global connectivity array called gcon, which returns the
    %   global DOF index when given the node index and the local DOF index.
    %   
    %   Had no displacement or rotation BCs been specified, gcon would
    %   contain these entries:
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
    %   We set gcon so that all the known displacements will be on the
    %   bottom side of u and all the known forces on the top side of f.
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
    % stiffness matrix (Hermite cubic polynomial element has 2 nodes)
    numDOFsPerElement = 2 * numDOFsPerNode;
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
        
        % Get the nodal positions
        x1 = nodes(node1Index, 1);
        y1 = nodes(node1Index, 2);
        x2 = nodes(node2Index, 1);
        y2 = nodes(node2Index, 2);
        
        % Get the element properties
        E  = elements(e, 3);
        A  = elements(e, 4);
        Ii = elements(e, 5);
        L  = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        
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
        K_e = [ E*A/L,            0,           0, -E*A/L,            0,           0; ...
                    0,  12*E*Ii/L^3,  6*E*Ii/L^2,      0, -12*E*Ii/L^3,  6*E*Ii/L^2; ...
                    0,   6*E*Ii/L^2,    4*E*Ii/L,      0,  -6*E*Ii/L^2,    2*E*Ii/L; ...
               -E*A/L,            0,           0,  E*A/L,            0,           0; ...
                    0, -12*E*Ii/L^3, -6*E*Ii/L^2,      0,  12*E*Ii/L^3, -6*E*Ii/L^2; ...
                    0,   6*E*Ii/L^2,    2*E*Ii/L,      0,  -6*E*Ii/L^2,    4*E*Ii/L];
        
        % Apply coordinate transformation
        K_e(1:3, 1:3) = Q' * K_e(1:3, 1:3) * Q;
        K_e(1:3, 4:6) = Q' * K_e(1:3, 4:6) * Q;
        K_e(4:6, 1:3) = Q' * K_e(4:6, 1:3) * Q;
        K_e(4:6, 4:6) = Q' * K_e(4:6, 4:6) * Q;
        
        
        %------------------------------------------------------------------
        %  Global assembly
        %------------------------------------------------------------------
        % Get the DOF indices
        index = [gcon(node1Index, 1); ...
                 gcon(node1Index, 2); ...
                 gcon(node1Index, 3); ...
                 gcon(node2Index, 1); ...
                 gcon(node2Index, 2); ...
                 gcon(node2Index, 3)];
        
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
    draw_frame(nodes, elements, u, gcon);
end
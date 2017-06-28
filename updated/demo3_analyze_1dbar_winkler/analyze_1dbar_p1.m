%--------------------------------------------------------------------------
%  Author:
%    
%    Isaac J. Lee (crunchingnumbers.live)
%    
%  Summary:
%    
%    This routine solves the bar (pile) problem with a Winkler foundation
%    using a linear polynomial basis.
%    
%  Instructions:
%    
%    Call this routine from the driver file demo3_driver_1dbar_p1:
%    
%    analyze_1dbar_p1(problem_size)
%    
%    where problem_size is '01', '02', '04', or '08'.
%    
%--------------------------------------------------------------------------
function analyze_1dbar_p1(problem_size)
    %----------------------------------------------------------------------
    %  Load the assembly file
    %----------------------------------------------------------------------
    load(strcat('./demo3_analyze_1dbar_winkler/assembly_files/assembly_1dbar_p1_numElements', problem_size));
    
    
    % Find the number of nodes, etc. from the assembly file
    numNodes            = size(nodes, 1);
    numElements         = size(elements, 1);
    numBCs_displacement = size(BCs_displacement, 1);
    numBCs_force        = size(BCs_force, 1);
    
    % Set the number of degrees of freedom (DOFs) for 1D elasticity
    numDOFsPerNode = 1;
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
    %   Set the quadrature rule (3-point Gauss quadrature) and evaluate
    %   the basis functions and their derivatives at the quadrature points
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    z = [-sqrt(3/5); 0; sqrt(3/5)];
    w = [5/9; 8/9; 5/9];
    numQuadraturePoints = length(z);
    
    eval_N = @(xi) [0.5 * (1 - xi), ...
                    0.5 * (1 + xi)];
    eval_B = @(xi) [-0.5 * ones(length(xi), 1), ...
                     0.5 * ones(length(xi), 1)];
    
    % Evaluate at the quadrature points
    N = eval_N(z);
    B = eval_B(z);
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Create the stiffness matrix K and the load vector f
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Initialize the stiffness matrix
    K = zeros(numDOFs);
    
    numNodesPerElement = 2;
    numDOFsPerElement = numNodesPerElement * numDOFsPerNode;
    
    
    %----------------------------------------------------------------------
    %  Specify the distributed load
    %----------------------------------------------------------------------
    P_0 = abs(BCs_force(1, 3));
    L   = nodes(numNodes, 1);
    eval_distributed_load = @(x) P_0 / L;
    
    
    %----------------------------------------------------------------------
    %  Loop over the elements
    %----------------------------------------------------------------------
    for e = 1 : numElements
        % Get the node indices
        node1Index = elements(e, 1);
        node2Index = elements(e, 2);
        
        % Get the nodal positions
        x1 = nodes(node1Index, 1);
        x2 = nodes(node2Index, 1);
        x_e = [x1; x2];
        
        % Get the element properties (constant)
        E = elements(e, 3);
        A = elements(e, 4);
        k = E*A / L^2;
        
        
        % Evaluate the elasticity matrix (constant)
        C = [E];
        
        % Initialize the element stiffness matrix and element load vector
        K_e = zeros(numNodesPerElement);
        f_e = zeros(numNodesPerElement, 1);
        
        
        %------------------------------------------------------------------
        %  Loop over the quadrature points
        %------------------------------------------------------------------
        for i = 1 : numQuadraturePoints
            % Get the values of the basis functions and their derivatives
            % at the quadrature point
            N_e = N(i, :);
            B_e = B(i, :);
            
            % Find the point x in the physical domain that corresponds to
            % the quadrature point
            x = N_e * x_e;
            
            % Evaluate the Jacobian at the quadrature point
            J = B_e * x_e;
            
            % Evaluate the distributed load at the quadrature point
            distributed_load = eval_distributed_load(x);
            
            
            %--------------------------------------------------------------
            %  Form the element stiffness matrix
            %--------------------------------------------------------------
            K_e = K_e + w(i) * ((A / J) * (B_e' * C * B_e) + (k * J) * (N_e' * N_e));
            
            
            %--------------------------------------------------------------
            %  Form the element load vector
            %--------------------------------------------------------------
            f_e = f_e + w(i) * ((distributed_load * J) * N_e');
        end
        
        
        %------------------------------------------------------------------
        %  Global assembly
        %------------------------------------------------------------------
        % Get the global DOF indices
        index = [gcon(node1Index, 1); ...
                 gcon(node2Index, 1)];
        
        % Assign the entries
        for i = 1 : numDOFsPerElement
            for j = 1 : numDOFsPerElement
                K(index(i), index(j)) = K(index(i), index(j)) + K_e(i, j);
            end
            
            f(index(i)) = f(index(i)) + f_e(i);
        end
    end
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Solve for the unknown displacements and forces
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Solve for the unknown displacements
    u(index_f) = K(index_f, index_f) \ (f(index_f) - K(index_f, index_u) * u(index_u));
    
    % Solve for the unknown forces
    f(index_u) = K(index_u, index_f) * u(index_f) + K(index_u, index_u) * u(index_u);
    
    
    save(strcat('./demo3_analyze_1dbar_winkler/solution_files/solution_1dbar_p1_numElements', problem_size), 'u', 'f', 'gcon', '-v7.3');
end
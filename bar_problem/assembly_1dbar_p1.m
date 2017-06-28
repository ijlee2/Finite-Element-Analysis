%--------------------------------------------------------------------------
%  Author: Isaac J. Lee
%  E-mail: ijlee2@ices.utexas.edu
%  
%  This routine solves the bar (pile) problem with a Winkler foundation
%  using a linear polynomial basis.
%--------------------------------------------------------------------------
function assembly_1dbar_p1(problem_size)
    % Load the assembly file
    load(strcat('assembly_files/assembly_1dbar_p1_numElements', problem_size));
    
    % Find the number of nodes, etc. from the assembly file
    numNodes = size(nodes, 1);
    numElements = size(elements, 1);
    numBCs_displacement = size(BCs_displacement, 1);
    numBCs_force = size(BCs_force, 1);
    
    % Set the number of degrees of freedom (DOFs) for 1D elasticity
    numDOFsPerNode = 1;
    numDOFs = numDOFsPerNode * numNodes;
    
    % Initialize the stiffness matrix, the solution and RHS vectors
    K = zeros(numDOFs, numDOFs);
    u = zeros(numDOFs, 1);
    f = zeros(numDOFs, 1);
    
    
    %----------------------------------------------------------------------
    %  Set the gcon array
    %  
    %  We create gcon so that the known displacements are stored on the
    %  bottom side of u and the known forces are stored on the top of f.
    %  While we populate the gcon array, we prescribe the known values
    %  into the vectors u and f.
    %----------------------------------------------------------------------
    gcon = zeros(numNodes, numDOFsPerNode);
    
    % index_f is a vector of indices for the known forces, and
    % index_u is a vector of indices for the known displacements
    index_f = 1 : numBCs_force;
    index_u = (numBCs_force + 1) : numDOFs;
    
    % Counter to keep track of how many DOFs we have encountered so far
    count = 1;
    
    % Prescribe the forces
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
    %  Set the quadrature rule (3-point Gauss quadrature)
    %----------------------------------------------------------------------
    z = [-sqrt(3/5); 0; sqrt(3/5)];
    w = [5/9; 8/9; 5/9];
    numQuadraturePoints = size(z, 1);
    
    % Specify the basis functions and their derivatives
    eval_N = @(xi) [-1/2*xi + 1/2, ...
                    1/2*xi + 1/2];
    eval_B = @(xi) [-1/2*ones(size(xi, 1), 1), ...
                    1/2*ones(size(xi, 1), 1)];
    
    % Evaluate them at the quadrature points
    N = eval_N(z);
    B = eval_B(z);
    
    % Specify the distributed load
    P_0 = abs(BCs_force(1, 3));
    L = nodes(numNodes, 1);
    eval_distributed_load = @(x) P_0 / L;
    
    
    %----------------------------------------------------------------------
    %  Create the stiffness matrix K and the load vector f
    %----------------------------------------------------------------------
    numNodesPerElement = 2;
    numDOFsPerElement = numNodesPerElement * numDOFsPerNode;
    
    % Loop over the elements
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
        k = E*A/L^2;
        
        % Evaluate the elasticity matrix (constant)
        C = [E];
        
        % Initialize the element stiffness matrix and element load vector
        K_e = zeros(numNodesPerElement, numNodesPerElement);
        f_e = zeros(numNodesPerElement, 1);
        
        % Loop over the quadrature points
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
            
            % Form the element stiffness matrix
            K_e = K_e + w(i) * (A * (B_e' * C * B_e) * 1/J + k * (N_e' * N_e) * J);
            
            % Form the element load vector
            f_e = f_e + w(i) * (distributed_load * N_e' * J);
        end
        
        % Get the global DOF indices
        index = [gcon(node1Index, 1); ...
                 gcon(node2Index, 1)];
        
        % Global assembly
        for i = 1 : numDOFsPerElement
            for j = 1 : numDOFsPerElement
                K(index(i), index(j)) = K(index(i), index(j)) + K_e(i, j);
            end
            
            f(index(i)) = f(index(i)) + f_e(i);
        end
    end
    
    
    %----------------------------------------------------------------------
    %  Solve for the unknown displacements and forces
    %----------------------------------------------------------------------
    % Solve for the unknown displacements
    u(index_f) = K(index_f, index_f) \ (f(index_f) - K(index_f, index_u) * u(index_u));
    
    % Solve for the unknown forces
    f(index_u) = K(index_u, index_f) * u(index_f) + K(index_u, index_u) * u(index_u);
    
    save(strcat('solution_files/solution_1dbar_p1_numElements', problem_size), 'u', 'f', 'gcon', '-v6');
end
%--------------------------------------------------------------------------
%  Author: Isaac J. Lee
%  E-mail: ijlee2@ices.utexas.edu
%  
%  This routine reads the assembly files provided by Dr. Landis, and
%  converts them into the .mat file format that we had used for truss
%  and frame problems. The .mat file contains four arrays called nodes,
%  elements, BCs_displacement, and BCs_force.
%  
%  The .bin files are assumed to be in the same directory as this routine,
%  and the .mat file will be created in this directory.
%  
%  To run this routine, type into Matlab's command window,
%      convertAssemblyFile(problem)
%  where
%      problem = 1, 2, 3, or 4
%--------------------------------------------------------------------------
function convertAssemblyFile(problem)
    clc;
    
    switch problem
        case 1
            fileName = '6';
        case 2
            fileName = '12';
        case 3
            fileName = '24';
        case 4
            fileName = 'R';
    end
    
    
    %----------------------------------------------------------------------
    %  Create the nodes array
    %  (x, y)
    %----------------------------------------------------------------------
    fileID = fopen(strcat('nodes', fileName), 'r');
    numNodes = str2num(fgetl(fileID));
    numDOFsPerNode = 2;
    numDOFs = numDOFsPerNode*numNodes;
    
    % Initialize the nodes array
    nodes = zeros(numNodes, 2);
    
    for i = 1 : numNodes
        temp = str2num(fgetl(fileID));
        
        % x-coordinate
        nodes(i, 1) = temp(2);
        
        % y-coordinate
        nodes(i, 2) = temp(3);
    end
    
    fclose(fileID);
    
    
    %----------------------------------------------------------------------
    %  Create the elements array
    %  (node 1 index, node 2 index, node 3 index, E, nu)
    %----------------------------------------------------------------------
    fileID = fopen(strcat('elements', fileName), 'r');
    temp = str2num(fgetl(fileID));
    numElements = temp(1);
%   E = temp(2);
    E = 50;
    nu = temp(3);
    
    % Initialize the elements array
    elements = zeros(numElements, 5);
    
    for i = 1 : numElements
        temp = str2num(fgetl(fileID));
        
        % node 1 index
        elements(i, 1) = temp(2);
        
        % node 2 index
        elements(i, 2) = temp(3);
        
        % node 3 index
        elements(i, 3) = temp(4);
    end
    
    % Young's modulus
    elements(:, 4) = E * ones(numElements, 1);
    
    % Poisson's ratio
    elements(:, 5) = nu * ones(numElements, 1);
    
    fclose(fileID);
    
    
    %----------------------------------------------------------------------
    %  Create the BCs_displacement array
    %  (node index, dof index, value)
    %----------------------------------------------------------------------
    fileID = fopen(strcat('displacements', fileName), 'r');
    numBCs_displacement = str2num(fgetl(fileID));
    
    % Initialize the BCs_displacement array
    BCs_displacement = zeros(numBCs_displacement, 3);
    
    for i = 1 : numBCs_displacement
        temp = str2num(fgetl(fileID));
        
        % node index
        BCs_displacement(i, 1) = temp(1);
        
        % dof index
        BCs_displacement(i, 2) = temp(2);
        
        % value
        BCs_displacement(i, 3) = temp(3);
    end
    
    fclose(fileID);
    
    
    %----------------------------------------------------------------------
    %  Create the BCs_force array
    %  (node index, dof index, value)
    %----------------------------------------------------------------------
    fileID = fopen(strcat('forces', fileName), 'r');
    numBCs_force = str2num(fgetl(fileID));
    
    % Initialize the BCs_force array
    BCs_force = zeros(numDOFs - numBCs_displacement, 3);
    
    % The forces.bin file only lists the nodal forces with nonzero values.
    % However, it is easier if we list all nodal forces, and set to zero
    % if there was no applied force at a node in a given direction.
    count = 1;
    for i = 1 : numNodes
        for j = 1 : numDOFsPerNode
            % Consider applying a zero nodal force only if we did not
            % prescribe a displacement BC already
            if (sum(ismember(BCs_displacement(:, [1 2]), [i j], 'rows')) == 0)
                BCs_force(count, 1) = i;
                BCs_force(count, 2) = j;
                BCs_force(count, 3) = 0;
                
                count = count + 1;
            end
        end
    end
    
    if ((numDOFs - numBCs_displacement) ~= (count - 1))
        fprintf('Warning: Please check this part of the code.\n');
    end
    
    % Prescribe the nodal forces with nonzero values
    for i = 1 : numBCs_force
        temp = str2num(fgetl(fileID));
        
        % Find the right row
        count = find(ismember(BCs_force(:, [1 2]), [temp(1) temp(2)], 'rows') == 1);
        
        BCs_force(count, 3) = temp(3);
    end
    
    fclose(fileID);
    
    
    %----------------------------------------------------------------------
    %  Create the .mat assembly file
    %----------------------------------------------------------------------
    save(strcat('assembly', fileName, '.mat'), 'nodes', 'elements', 'BCs_displacement', 'BCs_force', '-v6');
end
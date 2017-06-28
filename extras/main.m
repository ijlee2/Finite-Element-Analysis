%--------------------------------------------------------------------------
%  Author: Isaac J. Lee
%  E-mail: ijlee2@ices.utexas.edu
%  
%  This routine creates a GUI (graphical user interface) for analyzing
%  a 2D truss problem. The routine will be updated to also handle a 2D
%  frame problem.
%  
%  The primary outputs of this routine are the information needed for a
%  global assembly of the stiffness matrix and the external force vector.
%  The students are to provide the assembly and postprocessing routines
%  to tie into this GUI.
%  
%  To run this program, type into the Matlab command window:
%      main
%--------------------------------------------------------------------------
function main()
    % Close other Matlab programs, clear the screen and existing variables
    clc; clf;
    close all;
    clear all;
    
    % When displaying a number on Matlab's command window (for debugging),
    % show all digits after the decimal mark
    format long;
    
    %----------------------------------------------------------------------
    %  Declare global variables
    %  
    %  By declaring a variable to be global, we will be able to use it in
    %  a routine that sits "outside," without having to pass the variable
    %  as an input to that routine.
    %  
    %  workspaceName, workspaceDirectory
    %     -- Name of the workspace, and the directory under which the
    %        workspace is saved
    %  
    %  handle_gui, handle_prompt
    %     -- Handles (pointers) to the GUI and the prompt window
    %  
    %  nodes, elements, BCs
    %     -- Arrays containing information about the nodes, elements, and
    %        boundary conditions.
    %        
    %        The two columns of nodes correspond to,
    %      
    %            [x_coordinate, y_coordinate],
    %        
    %        the four columns of elements,
    %        
    %            [start node index, end node index, Young's modulus,
    %             cross-sectional area],
    %        
    %        and the four columns of boundary conditions,
    %        
    %            [node index, boundary condition type, x-component value,
    %             y-component value].
    %        
    %        Please note, we do not store the information of to which node
    %        a row of the nodes array corresponds, because this program is
    %        written such that the node indices are always consecutive
    %        integers, starting with the number 1. The same goes for the
    %        elements and BCs arrays.
    %        
    %        As a result, we can easily check the number of nodes that we
    %        have by calling Matlab's size routine:
    %        
    %            size(nodes, 1)
    %        
    %        The same goes for the number of elements and that of BCs.
    %  
    %  numNodes, numElements, numBCs
    %     -- The number of nodes, elements, and BCs that we have
    %  
    %  option_displayNodeIndex, option_displayElementIndex,
    %  option_displayBCs, option_displayGridAxes
    %     -- Strings (with a value of 'off' or 'on') indicating whether to
    %        display node indices, etc. on screen.
    %----------------------------------------------------------------------
    global workspaceName workspaceDirectory;
    global handle_gui handle_prompt;
    global nodes elements BCs;
    global numNodes numElements numBCs;
    global option_displayNodeIndex option_displayElementIndex option_displayBCs option_displayGridAxes;
    
    %----------------------------------------------------------------------
    %  Initialize the GUI
    %  
    %  We will create a GUI with two sides (they are often called panels).
    %  The left panel will display the nodes, elements, and loads, whereas
    %  the right panel will allow us to add or remove them.
    %  
    %  Please note that how the GUI actually looks depends on the computer
    %  and possibly the version of Matlab installed. As a result, things
    %  may look out of place and will need to be fixed by the students.
    %  Unfortunately, there is no easy fix to this problem.
    %----------------------------------------------------------------------
    % Set the name and directory of the workspace
    workspaceName = 'untitled';
    workspaceDirectory = '';
    
    % Set the handles
    handle_gui = ...
    figure('Name', ['2D Analysis of a Truss - ', workspaceName], ...
           'Units', 'pixels', ...
           'Position', [0 0 1000 720], ...
           'Color', [0.91 0.91 0.9], ...
           'DockControls', 'off', ...
           'MenuBar', 'none', ...
           'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'ToolBar', 'none');
    movegui('northwest');
    handle_prompt = [];
    
    % Set the nodes, elements, BCs arrays
    nodes = [];
    elements = [];
    BCs = [];
    numNodes = 0;
    numElements = 0;
    numBCs = 0;
    
    % Set the options
    option_displayNodeIndex = 'on';
    option_displayElementIndex = 'off';
    option_displayBCs = 'on';
    option_displayGridAxes = 'on';
    
    % Display the GUI
    refreshGUI();
end


%--------------------------------------------------------------------------
%  This routine refreshes the GUI screen. It is to be called whenever the
%  user makes a change by adding or removing a node, element, or BC.
%--------------------------------------------------------------------------
function refreshGUI()
    global handle_gui handle_prompt;
    global nodes elements BCs;
    global numNodes numElements numBCs;
    global option_displayNodeIndex option_displayElementIndex option_displayBCs option_displayGridAxes;
    
    % Set the current window to be the GUI, and clear the screen
    figure(handle_gui);
    clf;
    
    % Copyright information (PLEASE DO NOT REMOVE OR CHANGE THIS)
    uicontrol('Style', 'text', ...
              'String', [char(169), ' 2015 Isaac J. Lee (ijlee2@ices.utexas.edu)'], ...
              'Position', [15 0 470 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0 0 0], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    % Create a menu bar
    menu(1) = uimenu(handle_gui, 'Label', 'File');
    menu(2) = uimenu(handle_gui, 'Label', 'Tools');
    
    % Menu items
    menu_file(1) = uimenu(menu(1), 'Label', 'New workspace', 'Accelerator', 'n', 'Callback', @newWorkspace);
    menu_file(2) = uimenu(menu(1), 'Label', 'Open workspace', 'Accelerator', 'o', 'Callback', @openWorkspace);
    menu_file(3) = uimenu(menu(1), 'Label', 'Save workspace', 'Accelerator', 's', 'Callback', @saveWorkspace);
    menu_file(4) = uimenu(menu(1), 'Label', 'Export', 'Separator', 'on');
    menu_file(5) = uimenu(menu(1), 'Label', 'Close program', 'Accelerator', 'w', 'Separator', 'on', 'Callback', @closeProgram);
    menu_tools(1) = uimenu(menu(2), 'Label', 'Analyze structure', 'Accelerator', 'a', 'Callback', @featureNotAvailableYet);
    menu_tools(2) = uimenu(menu(2), 'Label', 'Postprocess', 'Accelerator', 'p', 'Callback', @featureNotAvailableYet);
    menu_tools(3) = uimenu(menu(2), 'Label', 'Display options', 'Separator', 'on');
    
    % Submenu items
    uimenu(menu_file(4), 'Label', 'Assembly file', 'Callback', @exportAssemblyFile);
    uimenu(menu_file(4), 'Label', 'Screenshot', 'Callback', @exportScreenshot);
    uimenu(menu_tools(3), 'Label', 'Node index', 'Accelerator', '1', 'Checked', option_displayNodeIndex, 'Callback', {@toggleOption, 'displayNodeIndex'});
    uimenu(menu_tools(3), 'Label', 'Element index', 'Accelerator', '2', 'Checked', option_displayElementIndex, 'Callback', {@toggleOption, 'displayElementIndex'});
    uimenu(menu_tools(3), 'Label', 'Boundary conditions', 'Accelerator', '3', 'Checked', option_displayBCs, 'Callback', {@toggleOption, 'displayBCs'});
    uimenu(menu_tools(3), 'Label', 'Grid', 'Accelerator', '4', 'Checked', option_displayGridAxes, 'Callback', {@toggleOption, 'displayGridAxes'});
    
    %----------------------------------------------------------------------
    %  Set up the left panel
    %----------------------------------------------------------------------
    subplot('Position', [0.08 0.1 0.36 0.8]);
    
    % Set the plot range
    if (numNodes == 0)
        xMin = -4; xMax = 4;
        yMin = -4; yMax = 4;
    else
        xMin = nodes(1, 1); xMax = nodes(1, 1);
        yMin = nodes(1, 2); yMax = nodes(1, 2);
        
        % Although not necessary, we create these arrays, because we will
        % need the nodal coordinates when we draw other objects later
        nodes_x = nodes(:, 1);
        nodes_y = nodes(:, 2);
        
        % Based on the current minimum and maximum x- and y- coordinates,
        % we try to extend the range so that the nodes and their indices
        % can be seen on screen
        nodes_xMin = min(nodes_x); nodes_xMax = max(nodes_x);
        nodes_yMin = min(nodes_y); nodes_yMax = max(nodes_y);
        
        padding_x = 0.05*max([xMax - nodes_xMin, nodes_xMax - xMin, nodes_xMax - nodes_xMin, 4]);
        padding_y = 0.05*max([yMax - nodes_yMin, nodes_yMax - yMin, nodes_yMax - nodes_yMin, 4]);
        
        xMin = min(xMin, floor(nodes_xMin - padding_x));
        xMax = max(xMax, ceil(nodes_xMax + padding_x));
        yMin = min(yMin, floor(nodes_yMin - padding_y));
        yMax = max(yMax, ceil(nodes_yMax + padding_y));
    end
    
    % Find how large the plot range is in the x- and y-directions
    plotSize_x = xMax - xMin;
    plotSize_y = yMax - yMin;
    
    % Draw the elements
    for i = 1 : numElements
        % Find the coordinates of the nodes that make up this element
        x1 = nodes_x(elements(i, 1)); y1 = nodes_y(elements(i, 1));
        x2 = nodes_x(elements(i, 2)); y2 = nodes_y(elements(i, 2));
        
        line([x1 x2], [y1 y2], 'Color', [0.7 0.3 0.3], 'LineWidth', 2); hold on;
        
        % Display the element index if requested
        if (strcmp(option_displayElementIndex, 'on'))
            if (abs(x2 - x1) >= abs(y2 - y1))
                padding_x = 0;
                padding_y = 0.03*plotSize_y;
            else
                padding_x = 0.02*plotSize_x;
                padding_y = 0;
            end
            
            text((x1 + x2)/2 + padding_x, (y1 + y2)/2 + padding_y, num2str(i), 'Color', [0.4 0.7 0.2], 'FontWeight', 'bold');
        end
    end
    
    % Draw the BCs if requested
    if (strcmp(option_displayBCs, 'on'))
        for i = 1 : numBCs
            % Find the coordinates of the node
            x1 = nodes_x(BCs(i, 1));
            y1 = nodes_y(BCs(i, 1));
            
            % Get the boundary condition type
            BCtype = BCs(i, 2);
                    
            % These are the components of a vector that points in the same
            % direction as the BC vector, but with a reasonable length such
            % that the arrow representing the vector can be seen on screen
            Delta_x = BCs(i, 3);
            Delta_y = BCs(i, 4);
            
            % Note, the nodes are colored dark green, while the elements
            % are colored dark red. We encourage the user to associate
            % the displacement BCs with the nodes and the force BCs with
            % the elements by coloring the corresponding arrows light green
            % or light red.
            if (BCtype == 1 || BCtype == 3)
                arrow_x_color = [0.4 0.9 0.6];
            else
                arrow_x_color = [0.95 0.6 0.7];
            end
            if (BCtype == 1 || BCtype == 4)
                arrow_y_color = [0.4 0.9 0.6];
            else
                arrow_y_color = [0.95 0.6 0.7];
            end
            
            % If both components are zero, just draw a square around the node
            if (Delta_x == 0 && Delta_y == 0)
                rectangle('Position', [x1 - plotSize_x/64, y1 - plotSize_y/64, plotSize_x/64, plotSize_y/32], 'EdgeColor', arrow_x_color, 'FaceColor', arrow_x_color); hold on;
                rectangle('Position', [x1,                 y1 - plotSize_y/64, plotSize_x/64, plotSize_y/32], 'EdgeColor', arrow_y_color, 'FaceColor', arrow_y_color); hold on;
                
            % For a pure Dirichlet or pure Neumann BC, we draw a single arrow
            elseif (BCtype == 1 || BCtype == 2)
                % Sign change for the arrow head
                if (Delta_x >= 0)
                    arrowhead_sign = 1;
                else
                    arrowhead_sign = -1;
                end
                arrow_length = norm([Delta_x Delta_y]);
                arrow_angle = arrowhead_sign * atan(Delta_y/abs(Delta_x));
                
                % Find the coordinates of the arrow head
                arrowhead_x = x1 + plotSize_x/18 * (Delta_x/arrow_length);
                arrowhead_y = y1 + plotSize_y/18 * (Delta_y/arrow_length);
                
                % If the arrow is horizontal or vertical, we shift the two ends
                % of the arrow head by a bit
                shift_x = 0;
                shift_y = 0;
                if (Delta_x == 0)
                    shift_x = -plotSize_x/400;
                elseif (Delta_y == 0)
                    shift_y = plotSize_y/400;
                end
                
                % Draw the arrow
                line([x1 arrowhead_x], [y1 arrowhead_y], 'Color', arrow_x_color, 'LineWidth', 3.5); hold on;
                line(arrowhead_x + shift_x - arrowhead_sign * [0, plotSize_x/36 * cos(arrow_angle + 0.4)], ...
                     arrowhead_y + shift_y - arrowhead_sign * [0, plotSize_y/36 * sin(arrow_angle + 0.4)], 'Color', arrow_x_color, 'LineWidth', 3); hold on;
                line(arrowhead_x + shift_x - arrowhead_sign * [0, plotSize_x/36 * cos(arrow_angle - 0.4)], ...
                     arrowhead_y + shift_y - arrowhead_sign * [0, plotSize_y/36 * sin(arrow_angle - 0.4)], 'Color', arrow_x_color, 'LineWidth', 3); hold on;
                
            % For a mixed BC, we draw two separate arrows
            else
                if (Delta_x == 0)
                    rectangle('Position', [x1 - plotSize_x/64, y1 - plotSize_y/64, plotSize_x/32, plotSize_y/32], 'EdgeColor', arrow_x_color, 'FaceColor', arrow_x_color); hold on;
                else
                    % Sign change for the arrow head
                    arrowhead_sign = sign(Delta_x);
                    
                    % Find the coordinates of the arrow head
                    arrowhead_x = x1 + arrowhead_sign * plotSize_x/18;
                    arrowhead_y = y1;
                    shift_y = plotSize_y/400;
                    
                    % Draw the arrow (horizontal)
                    line([x1 arrowhead_x], [y1 arrowhead_y], 'Color', arrow_x_color, 'LineWidth', 3.5); hold on;
                    line(arrowhead_x           - arrowhead_sign * [0, plotSize_x/36 * cos(0.4)], ...
                         arrowhead_y + shift_y - arrowhead_sign * [0, plotSize_y/36 * sin(0.4)], 'Color', arrow_x_color, 'LineWidth', 3); hold on;
                    line(arrowhead_x           - arrowhead_sign * [0, plotSize_x/36 * cos(-0.4)], ...
                         arrowhead_y + shift_y - arrowhead_sign * [0, plotSize_y/36 * sin(-0.4)], 'Color', arrow_x_color, 'LineWidth', 3); hold on;
                end
                
                if (Delta_y == 0)
                    rectangle('Position', [x1 - plotSize_x/64, y1 - plotSize_y/64, plotSize_x/32, plotSize_y/32], 'EdgeColor', arrow_y_color, 'FaceColor', arrow_y_color); hold on;
                else
                    % Sign change for the arrow head
                    arrowhead_sign = sign(Delta_y);
                    
                    % Find the coordinates of the arrow head
                    arrowhead_x = x1;
                    arrowhead_y = y1 + arrowhead_sign * plotSize_y/18;
                    shift_x = -plotSize_x/400;
                    
                    % Draw the arrow (vertical)
                    line([x1 arrowhead_x], [y1 arrowhead_y], 'Color', arrow_y_color, 'LineWidth', 3.5); hold on;
                    line(arrowhead_x + shift_x - [0, plotSize_x/36 * cos(arrowhead_sign * pi/2 + 0.4)], ...
                         arrowhead_y           - [0, plotSize_y/36 * sin(arrowhead_sign * pi/2 + 0.4)], 'Color', arrow_y_color, 'LineWidth', 3); hold on;
                    line(arrowhead_x + shift_x - [0, plotSize_x/36 * cos(arrowhead_sign * pi/2 - 0.4)], ...
                         arrowhead_y           - [0, plotSize_y/36 * sin(arrowhead_sign * pi/2 - 0.4)], 'Color', arrow_y_color, 'LineWidth', 3); hold on;
                end
            end
        end
    end
    
    % Draw the nodes
    for i = 1 : numNodes
        % Find the coordinates of the node
        x1 = nodes_x(i);
        y1 = nodes_y(i);
        
        plot(x1, y1, 's', 'Color', [0 0 0], 'MarkerFaceColor', [0.4 0.8 0.3]); hold on;
        
        % Display the node index if requested
        if (strcmp(option_displayNodeIndex, 'on') == 1)
            padding_x = 0.02*plotSize_x;
            padding_y = 0.03*plotSize_y;
            
            text(x1 + padding_x, y1 + padding_y, num2str(i), 'Color', [0.5 0.4 0.1], 'FontWeight', 'bold');
        end
    end
    
    % Label the x- and y-axes
    xlabel('x', 'FontSize', 18);
    ylabel('y', 'FontSize', 18, 'Rotation', 0);
    
    % Create ticks as visual aids
    set(gca, 'FontSize', 13, ...
             'XTick', linspace(xMin, xMax, 9), 'XTickLabel', sprintf('%.1f|', linspace(xMin, xMax, 9)), ...
             'YTick', linspace(yMin, yMax, 9), 'YTickLabel', sprintf('%.1f|', linspace(yMin, yMax, 9)));
    
    % Draw the x- and y-axes and a grid if requested
    if (strcmp(option_displayGridAxes, 'on'))
        grid on;
    end
    
    % Set the plot axes
    axis([xMin xMax yMin yMax]);
    axis square;
    
    %----------------------------------------------------------------------
    %  Set up the right panel
    %----------------------------------------------------------------------
    % Nodes (keypoints) section
    uicontrol('Style', 'text', ...
              'String', '  Nodes (keypoints)', ...
              'Position', [530 660 400 25], ...
              'BackgroundColor', [0.2 0.6 0.3], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'text', ...
              'String', ' Index       Coordinates (x, y)', ...
              'Position', [530 625 400 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'listbox', ...
              'String', listNodes(), ...
              'Position', [530 538 401 90], ...
              'BackgroundColor', [0.9 0.92 0.9]);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Add', ...
              'Position', [530 498 120 25], ...
              'FontSize', 10, ...
              'Callback', @addNode);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Edit', ...
              'Position', [670 498 120 25], ...
              'FontSize', 10, ...
              'Callback', @editNode);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Remove', ...
              'Position', [810 498 120 25], ...
              'FontSize', 10, ...
              'Callback', @removeNode);
    
    % Elements (lines) section
    uicontrol('Style', 'text', ...
              'String', '  Elements (lines)', ...
              'Position', [530 430 400 25], ...
              'BackgroundColor', [0.2 0.6 0.3], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'text', ...
              'String', ' Index       Node connectivity      Material properties (E, A)', ...
              'Position', [530 395 400 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'listbox', ...
              'String', listElements(), ...
              'Position', [530 308 401 90], ...
              'BackgroundColor', [0.9 0.92 0.9]);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Add', ...
              'Position', [670 268 120 25], ...
              'FontSize', 10, ...
              'Callback', @addElement);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Remove', ...
              'Position', [810 268 120 25], ...
              'FontSize', 10, ...
              'Callback', @removeElement);
    
    % Boundary conditions section
    uicontrol('Style', 'text', ...
              'String', '  Boundary conditions', ...
              'Position', [530 200 400 25], ...
              'BackgroundColor', [0.2 0.6 0.3], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'text', ...
              'String', ' Node index    BC type (x, y)        Component values (x, y)', ...
              'Position', [530 165 400 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'listbox', ...
              'String', listBCs(), ...
              'Position', [530 78 401 90], ...
              'BackgroundColor', [0.9 0.92 0.9]);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Add', ...
              'Position', [670 38 120 25], ...
              'FontSize', 10, ...
              'Callback', @addBC);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Remove', ...
              'Position', [810 38 120 25], ...
              'FontSize', 10, ...
              'Callback', @removeBC);
    
    % If refreshGUI was called from a prompt window, switch back to it
    if (~isempty(handle_prompt))
        figure(handle_prompt);
    end
end


%--------------------------------------------------------------------------
%  These routines specify what to do when the menu items are selected.
%--------------------------------------------------------------------------
function newWorkspace(hObject, eventdata)
    global workspaceName workspaceDirectory;
    global handle_gui handle_prompt;
    global nodes elements BCs;
    global numNodes numElements numBCs;
    global option_displayNodeIndex option_displayElementIndex option_displayBCs option_displayGridAxes;
    
    % If there is a prompt window already, close it
    if (ishandle(handle_prompt))
        close_prompt();
    end
    
    % Set to the default values
    workspaceName = 'untitled';
    workspaceDirectory = '';
    nodes = [];
    elements = [];
    BCs = [];
    numNodes = 0;
    numElements = 0;
    numBCs = 0;
    option_displayNodeIndex = 'on';
    option_displayElementIndex = 'off';
    option_displayBCs = 'on';
    option_displayGridAxes = 'on';
    
    % Change the workspace name
    set(handle_gui, 'Name', ['2D Analysis of a Truss - ', workspaceName]);
    
    % Display the GUI
    refreshGUI();
end


function openWorkspace(hObject, eventdata)
    global workspaceName workspaceDirectory;
    global handle_gui handle_prompt;
    global nodes elements BCs;
    global numNodes numElements numBCs;
    global option_displayNodeIndex option_displayElementIndex option_displayBCs option_displayGridAxes;
    
    % Prompt window for opening a workspace
    [fileName, workspaceDirectory] = uigetfile({'*.mat', 'Workspace for ASE 321K (*.mat)'}, 'Please select the workspace file to open');
    
    if (~isempty(strcat(workspaceDirectory, fileName)))
        load(strcat(workspaceDirectory, fileName), 'nodes', 'elements', 'BCs', 'options');
        
        % Change the workspace name
        workspaceName = fileName(1 : (end - 4));
        set(handle_gui, 'Name', ['2D Analysis of a Truss - ', workspaceName]);
        
        % Set the number of nodes, elements, BCs
        numNodes = size(nodes, 1);
        numElements = size(elements, 1);
        numBCs = size(BCs, 1);
        
        % Set the options
        option_displayNodeIndex = options{1};
        option_displayElementIndex = options{2};
        option_displayBCs = options{3};
        option_displayGridAxes = options{4};
        
        % If there is a prompt window already, close it
        if (ishandle(handle_prompt))
            close_prompt();
        end
        
        % Display the GUI
        refreshGUI();
    end
end


function saveWorkspace(hObject, eventdata)
    global workspaceName workspaceDirectory;
    global handle_gui;
    global nodes elements BCs;
    global option_displayNodeIndex option_displayElementIndex option_displayBCs option_displayGridAxes;
    
    options = {option_displayNodeIndex;
               option_displayElementIndex;
               option_displayBCs;
               option_displayGridAxes};
    
    % Prompt window for saving the workspace
    [fileName, workspaceDirectory] = uiputfile({'*.mat', 'Workspace for ASE 321K (*.mat)'}, 'Save the workspace as', workspaceName);
    save(strcat(workspaceDirectory, fileName), 'nodes', 'elements', 'BCs', 'options', '-v6');
    
    % Change the workspace name
    workspaceName = fileName(1 : (end - 4));
    set(handle_gui, 'Name', ['2D Analysis of a Truss - ', workspaceName]);
end


function exportAssemblyFile(hObject, eventdata)
    global workspaceName workspaceDirectory;
    global handle_gui handle_prompt;
    global nodes elements BCs;
    global numNodes numElements numBCs;
    global option_displayNodeIndex option_displayElementIndex option_displayBCs option_displayGridAxes;
    
    % Find the current position of the GUI window
    gui_position = get(handle_gui, 'Position');
    gui_outerPosition = get(handle_gui, 'OuterPosition');
    
    % Check that there are at least two nodes and an element
    if (numNodes < 2 || numElements < 1)
        handle_prompt = ...
        figure('Name', 'Error', ...
               'Position', [0 0 340 80], ...
               'Color', [1 0.8 0.25], ...
               'MenuBar', 'none', ...
               'NumberTitle', 'off', ...
               'Resize', 'off', ...
               'ToolBar', 'none');
        movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
        
        uicontrol('Style', 'text', ...
                  'String', sprintf('Umm... a structure must have at least two nodes\nand an element.'), ...
                  'Position', [0 22 340 50], ...
                  'BackgroundColor', [1 0.8 0.25], ...
                  'ForegroundColor', [0.8 0.1 0.1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'center');
        
        uicontrol('Style', 'pushbutton', ...
                  'String', 'OK', ...
                  'Position', [120 10 100 20], ...
                  'FontSize', 9, ...
                  'Callback', @close_prompt);
        
        return;
    end
        
    % Check that the nodes and elements form a connected graph
    adjacencyMatrix = eye(numNodes);
    for i = 1 : numElements
        % Find the node indices
        index1 = elements(i, 1);
        index2 = elements(i, 2);
        
        % Set the corresponding entries in the adjacency matrix to 1
        adjacencyMatrix(index1, index2) = 1;
        adjacencyMatrix(index2, index1) = 1;
    end
    [p, q, r, s] = dmperm(adjacencyMatrix);
    
    % By typing in p(r(1):r(2) - 1), p(r(2):r(3) - 1), etc., we can find
    % the nodes of each disjoint graph. This means, we can use the size
    % of the vector r to determine whether the nodes and elements form
    % a connected graph
    if (size(r, 2) ~= 2)
        handle_prompt = ...
        figure('Name', 'Error', ...
               'Position', [0 0 340 80], ...
               'Color', [1 0.8 0.25], ...
               'MenuBar', 'none', ...
               'NumberTitle', 'off', ...
               'Resize', 'off', ...
               'ToolBar', 'none');
        movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
        
        uicontrol('Style', 'text', ...
                  'String', sprintf('Umm... the structure has disjoint members.\nPlease add elements to connect these members.'), ...
                  'Position', [0 22 340 50], ...
                  'BackgroundColor', [1 0.8 0.25], ...
                  'ForegroundColor', [0.8 0.1 0.1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'center');
        
        uicontrol('Style', 'pushbutton', ...
                  'String', 'OK', ...
                  'Position', [120 10 100 20], ...
                  'FontSize', 9, ...
                  'Callback', @close_prompt);
        
        return;
    end
    
    % If there are nodes without a BC specified, assume that they are
    % traction-free (i.e. zero force)
    if (numBCs < numNodes)
        % Find the nodes without a BC
        if (numBCs > 0)
            nodesWithoutBCs = setdiff((1 : numNodes)', BCs(:, 1));
        else
            nodesWithoutBCs = (1 : numNodes)';
        end
        
        % Set zero force BC to these nodes
        BCs = [BCs; ...
               [nodesWithoutBCs, ones(numNodes - numBCs, 1) * [2 0 0]]];
        
        % Update the number of BCs
        numBCs = numNodes;
        
        % Sort the BCs in ascending node index
        if (numBCs > 1)
            [temp, permutation] = sort(BCs(:, 1));
            BCs = BCs(permutation, :);
        end
    end
    
    % Initialize
    BCs_displacement  = [];
    BCs_force = [];
    
    for i = 1 : numBCs
        switch BCs(i, 2)
            % Dirichlet (displacement in x and y)
            case 1
                BCs_displacement  = [BCs_displacement; ...
                                     BCs(i, 1), 1, BCs(i, 3); ...
                                     BCs(i, 1), 2, BCs(i, 4)];
            % Neumann (force in x and y)
            case 2
                BCs_force         = [BCs_force; ...
                                     BCs(i, 1), 1, BCs(i, 3); ...
                                     BCs(i, 1), 2, BCs(i, 4)];
            % Mixed 1 (displacement in x, force in y)
            case 3
                BCs_displacement  = [BCs_displacement; ...
                                     BCs(i, 1), 1, BCs(i, 3)];
                BCs_force         = [BCs_force; ...
                                     BCs(i, 1), 2, BCs(i, 4)];
            % Mixed 2 (force in x, displacement in y)
            case 4
                BCs_force         = [BCs_force; ...
                                     BCs(i, 1), 1, BCs(i, 3)];
                BCs_displacement  = [BCs_displacement; ...
                                     BCs(i, 1), 2, BCs(i, 4)];
        end
    end
    
    options = {option_displayNodeIndex;
               option_displayElementIndex;
               option_displayBCs;
               option_displayGridAxes};
    
    % Save the files
    save(strcat(workspaceDirectory, workspaceName, '.mat'), 'nodes', 'elements', 'BCs', 'options', '-v6');
    save(strcat(workspaceDirectory, workspaceName, '_assembly.mat'), 'nodes', 'elements', 'BCs_displacement', 'BCs_force', '-v6');
    
    %----------------------------------------------------------------------
    %  Create text files in case the students want to use these in Matlab,
    %  or want to open the files in another program. A whitespace will be
    %  used as the delimiter.
    %----------------------------------------------------------------------
    % Create a text file for nodes
    fileID = fopen(strcat(workspaceDirectory, workspaceName, '_assembly_nodes.txt'), 'w');
    fprintf(fileID, '%d\n', numNodes);
    for i = 1 : numNodes
        fprintf(fileID, '%.15f %.15f\n', nodes(i, 1), nodes(i, 2));
    end
    fclose(fileID);
    
    % Create a text file for elements
    fileID = fopen(strcat(workspaceDirectory, workspaceName, '_assembly_elements.txt'), 'w');
    fprintf(fileID, '%d\n', numElements);
    for i = 1 : numElements
        fprintf(fileID, '%d %d %.15f %.15f\n', elements(i, 1), elements(i, 2), elements(i, 3), elements(i, 4));
    end
    fclose(fileID);
    
    % Create a text file for displacement BCs
    fileID = fopen(strcat(workspaceDirectory, workspaceName, '_assembly_BCs_displacement.txt'), 'w');
    fprintf(fileID, '%d\n', size(BCs_displacement, 1));
    for i = 1 : size(BCs_displacement, 1)
        fprintf(fileID, '%d %d %.15f\n', BCs_displacement(i, 1), BCs_displacement(i, 2), BCs_displacement(i, 3));
    end
    fclose(fileID);
    
    % Create a text file for force BCs
    fileID = fopen(strcat(workspaceDirectory, workspaceName, '_assembly_BCs_force.txt'), 'w');
    fprintf(fileID, '%d\n', size(BCs_force, 1));
    for i = 1 : size(BCs_force, 1)
        fprintf(fileID, '%d %d %.15f\n', BCs_force(i, 1), BCs_force(i, 2), BCs_force(i, 3));
    end
    fclose(fileID);
    
    % Display the GUI
    refreshGUI();
    
    % Feedback message
    handle_prompt = ...
    figure('Name', 'Success', ...
           'Position', [0 0 340 80], ...
           'Color', [0.2 0.6 0.3], ...
           'MenuBar', 'none', ...
           'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'ToolBar', 'none');
    movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
    
    uicontrol('Style', 'text', ...
              'String', sprintf('The assembly file has been created.\nPlease check %s_assembly.mat.', workspaceName), ...
              'Position', [0 22 340 50], ...
              'BackgroundColor', [0.2 0.6 0.3], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'center');
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'OK', ...
              'Position', [120 10 100 20], ...
              'FontSize', 9, ...
              'Callback', @close_prompt);
end


function exportScreenshot(hObject, eventdata)
    global workspaceName workspaceDirectory;
    global handle_gui handle_prompt;
    
    set(handle_gui, 'PaperUnits', 'points', 'PaperPosition', [0, 0, 400, 400]);
    print('-dpng', '-r300', strcat(workspaceDirectory, workspaceName, '_plot.png'));
    
    % Find the current position of the GUI window
    gui_position = get(handle_gui, 'Position');
    gui_outerPosition = get(handle_gui, 'OuterPosition');
    
    % Feedback message
    handle_prompt = ...
    figure('Name', 'Success', ...
           'Position', [0 0 340 80], ...
           'Color', [0.2 0.6 0.3], ...
           'MenuBar', 'none', ...
           'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'ToolBar', 'none');
    movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
    
    uicontrol('Style', 'text', ...
              'String', sprintf('The screenshot has been created.\nPlease check %s_plot.png.', workspaceName), ...
              'Position', [0 22 340 50], ...
              'BackgroundColor', [0.2 0.6 0.3], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'center');
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'OK', ...
              'Position', [120 10 100 20], ...
              'FontSize', 9, ...
              'Callback', @close_prompt);
end


function closeProgram(hObject, eventdata)
    close all;
end


function toggleOption(hObject, eventdata, option)
    global option_displayNodeIndex option_displayElementIndex option_displayBCs option_displayGridAxes;
    
    switch option
        case 'displayNodeIndex'
            if (strcmp(option_displayNodeIndex, 'off'))
                option_displayNodeIndex = 'on';
            else
                option_displayNodeIndex = 'off';
            end
            
        case 'displayElementIndex'
            if (strcmp(option_displayElementIndex, 'off'))
                option_displayElementIndex = 'on';
            else
                option_displayElementIndex = 'off';
            end
            
        case 'displayBCs'
            if (strcmp(option_displayBCs, 'off'))
                option_displayBCs = 'on';
            else
                option_displayBCs = 'off';
            end
            
        case 'displayGridAxes'
            if (strcmp(option_displayGridAxes, 'off'))
                option_displayGridAxes = 'on';
            else
                option_displayGridAxes = 'off';
            end
            
    end
    
    refreshGUI();
end


function output = listNodes()
    global nodes numNodes;
    
    output = {};
    
    for i = 1 : numNodes
        x = nodes(i, 1);
        y = nodes(i, 2);
        
        numSpaces = max(12 - 2*floor(log10(i)), 0);
        
        output{i} = sprintf(' %d %s (%.3f, %.3f)\n', i, repmat(' ', [1 numSpaces]), x, y);
    end
end


function output = listElements()
    global elements numElements;
    
    output = {};
    
    for i = 1 : numElements
        index1 = elements(i, 1);
        index2 = elements(i, 2);
        E = elements(i, 3);
        A = elements(i, 4);
        
        numSpaces = [max(12 - 2*floor(log10(i)), 0); ...
                     max(24 - 2*floor(log10(index1)) - 2*floor(log10(index2)), 0)];
        
        output{i} = sprintf(' %d %s (%d, %d) %s %5.3f, %5.3f\n', i, repmat(' ', [1 numSpaces(1)]), index1, index2, repmat(' ', [1 numSpaces(2)]), E, A);
    end
end


function output = listBCs()
    global BCs numBCs;
    
    output = {};
    
    for i = 1 : numBCs
        index = BCs(i, 1);
        BCtype = BCs(i, 2);
        x = BCs(i, 3);
        y = BCs(i, 4);
        
        switch BCtype
            case 1
                BCtype_str = 'Disp., Disp.';
            case 2
                BCtype_str = 'Force, Force';
            case 3
                BCtype_str = 'Disp., Force';
            case 4
                BCtype_str = 'Force, Disp.';
        end
        
        numSpaces = [max(18 - 2*floor(log10(i)), 0); ...
                     max(20 - size(BCtype_str, 2), 0)];
        
        output{i} = sprintf(' %d %s %s %s %5.3f, %5.3f\n', index, repmat(' ', [1 numSpaces(1)]), BCtype_str, repmat(' ', [1 numSpaces(2)]), x, y);
    end
end


function output = listIndices(type)
    global nodes elements BCs;
    global numNodes numElements numBCs;
    
    % Initialize a string array of nodes
    output = {};
    
    switch type
        case 'nodes'
            for i = 1 : numNodes;
                output{i} = sprintf('%d', i);
            end
            
        case 'elements'
            for i = 1 : numElements;
                output{i} = sprintf('%d', i);
            end
            
        case 'nodesWithBCs'
            for i = 1 : numBCs
                output{i} = sprintf('%d', BCs(i, 1));
            end
            
        case 'nodesWithoutBCs'
            if (numBCs > 0)
                temp = setdiff((1 : numNodes)', BCs(:, 1));
            else
                temp = (1 : numNodes)';
            end
            
            for i = 1 : size(temp, 1);
                output{i} = sprintf('%d', temp(i));
            end
            
    end
end


function featureNotAvailableYet(hObject, eventdata)
    global handle_gui handle_prompt;
    
    % Find the current position of the GUI window
    gui_position = get(handle_gui, 'Position');
    gui_outerPosition = get(handle_gui, 'OuterPosition');
    
    % Feedback message
    handle_prompt = ...
    figure('Name', 'Error', ...
           'Position', [0 0 340 80], ...
           'Color', [1 0.8 0.25], ...
           'MenuBar', 'none', ...
           'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'ToolBar', 'none');
    movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
    
    uicontrol('Style', 'text', ...
              'String', sprintf('Sorry, this feature is not available yet.\nPlease try it again later.'), ...
              'Position', [0 22 340 50], ...
              'BackgroundColor', [1 0.8 0.25], ...
              'ForegroundColor', [0.8 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'center');
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'OK', ...
              'Position', [120 10 100 20], ...
              'FontSize', 9, ...
              'Callback', @close_prompt);
end


%--------------------------------------------------------------------------
%  These routines specify what the prompt window looks like. The prompt
%  is created to the right of the GUI window, so that it is easier to
%  navigate between the two windows. We will specify later what to do
%  when the user presses the ok button.
%--------------------------------------------------------------------------
function addNode(hObject, eventdata)
    global handle_gui handle_prompt handle_edit;
    
    % Find the current position of the GUI window
    gui_position = get(handle_gui, 'Position');
    gui_outerPosition = get(handle_gui, 'OuterPosition');
    
    %----------------------------------------------------------------------
    %  Initialize the prompt window
    %----------------------------------------------------------------------
    % If there is a prompt window already, close it
    if (ishandle(handle_prompt))
        close_prompt();
    end
    
    handle_prompt = ...
    figure('Name', 'Prompt window', ...
           'Position', [0 0 300 745], ...
           'Color', [0.91 0.91 0.9], ...
           'MenuBar', 'none', ...
           'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'ToolBar', 'none');
    movegui(handle_prompt, [gui_position(1) + gui_outerPosition(3), gui_position(2)]);
    
    %----------------------------------------------------------------------
    %  Set up the panel
    %----------------------------------------------------------------------
    uicontrol('Style', 'text', ...
              'String', sprintf('  To add a node,\n  please enter these information.'), ...
              'Position', [15 645 271 60], ...
              'BackgroundColor', [0.3 0.4 0.7], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'text', ...
              'String', ' x-coordinate:', ...
              'Position', [15 580 276 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(1) = ...
    uicontrol('Style', 'edit', ...
              'Position', [18 550 248 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'text', ...
              'String', ' y-coordinate:', ...
              'Position', [15 495 282 20], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(2) = ...
    uicontrol('Style', 'edit', ...
              'Position', [18 465 248 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'OK', ...
              'Position', [20 375 120 25], ...
              'FontSize', 10, ...
              'Callback', @ok_addNode);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Return', ...
              'Position', [160 375 120 25], ...
              'FontSize', 10, ...
              'Callback', @close_prompt);
    
end


function editNode(hObject, eventdata)
    global handle_gui handle_prompt handle_edit;
    global nodes;
    global numNodes;
    
    % Find the current position of the GUI window
    gui_position = get(handle_gui, 'Position');
    gui_outerPosition = get(handle_gui, 'OuterPosition');
    
    %----------------------------------------------------------------------
    %  Initialize the prompt window
    %----------------------------------------------------------------------
    % If there is a prompt window already, close it
    if (ishandle(handle_prompt))
        close_prompt();
    end
    
    if (numNodes == 0)
        handle_prompt = ...
        figure('Name', 'Error', ...
               'Position', [0 0 340 80], ...
               'Color', [1 0.8 0.25], ...
               'MenuBar', 'none', ...
               'NumberTitle', 'off', ...
               'Resize', 'off', ...
               'ToolBar', 'none');
        movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
        
        uicontrol('Style', 'text', ...
                  'String', sprintf('Umm... there are no nodes.\nTo edit a node, please add a node first.'), ...
                  'Position', [0 22 340 50], ...
                  'BackgroundColor', [1 0.8 0.25], ...
                  'ForegroundColor', [0.8 0.1 0.1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'center');
        
        uicontrol('Style', 'pushbutton', ...
                  'String', 'OK', ...
                  'Position', [120 10 100 20], ...
                  'FontSize', 9, ...
                  'Callback', @close_prompt);
        
        return;
    end
    
    handle_prompt = ...
    figure('Name', 'Prompt window', ...
           'Position', [0 0 300 745], ...
           'Color', [0.91 0.91 0.9], ...
           'MenuBar', 'none', ...
           'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'ToolBar', 'none');
    movegui(handle_prompt, [gui_position(1) + gui_outerPosition(3), gui_position(2)]);
    
    %----------------------------------------------------------------------
    %  Set up the panel
    %----------------------------------------------------------------------
    uicontrol('Style', 'text', ...
              'String', sprintf('  To edit a node,\n  please enter these information.'), ...
              'Position', [15 645 271 60], ...
              'BackgroundColor', [0.3 0.4 0.7], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'text', ...
              'String', ' Index of the node:', ...
              'Position', [15 580 276 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(1) = ...
    uicontrol('Style', 'popupmenu', ...
              'String', listIndices('nodes'), ...
              'Position', [19 550 89 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'text', ...
              'String', ' x-coordinate:', ...
              'Position', [15 495 276 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(2) = ...
    uicontrol('Style', 'edit', ...
              'Position', [18 465 248 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'text', ...
              'String', ' y-coordinate:', ...
              'Position', [15 410 282 20], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(3) = ...
    uicontrol('Style', 'edit', ...
              'Position', [18 380 248 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'OK', ...
              'Position', [20 290 120 25], ...
              'FontSize', 10, ...
              'Callback', @ok_editNode);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Return', ...
              'Position', [160 290 120 25], ...
              'FontSize', 10, ...
              'Callback', @close_prompt);
    
end


function removeNode(hObject, eventdata)
    global handle_gui handle_prompt handle_edit;
    global nodes;
    global numNodes;
    
    % Find the current position of the GUI window
    gui_position = get(handle_gui, 'Position');
    gui_outerPosition = get(handle_gui, 'OuterPosition');
    
    %----------------------------------------------------------------------
    %  Initialize the prompt window
    %----------------------------------------------------------------------
    % If there is a prompt window already, close it
    if (ishandle(handle_prompt))
        close_prompt();
    end
    
    if (numNodes == 0)
        handle_prompt = ...
        figure('Name', 'Error', ...
               'Position', [0 0 340 80], ...
               'Color', [1 0.8 0.25], ...
               'MenuBar', 'none', ...
               'NumberTitle', 'off', ...
               'Resize', 'off', ...
               'ToolBar', 'none');
        movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
        
        uicontrol('Style', 'text', ...
                  'String', sprintf('Umm... there are no nodes.\nPlease add a node first.'), ...
                  'Position', [0 22 340 50], ...
                  'BackgroundColor', [1 0.8 0.25], ...
                  'ForegroundColor', [0.8 0.1 0.1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'center');
        
        uicontrol('Style', 'pushbutton', ...
                  'String', 'OK', ...
                  'Position', [120 10 100 20], ...
                  'FontSize', 9, ...
                  'Callback', @close_prompt);
        
        return;
    end
    
    handle_prompt = ...
    figure('Name', 'Prompt window', ...
           'Position', [0 0 300 745], ...
           'Color', [0.91 0.91 0.9], ...
           'MenuBar', 'none', ...
           'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'ToolBar', 'none');
    movegui(handle_prompt, [gui_position(1) + gui_outerPosition(3), gui_position(2)]);
    
    %----------------------------------------------------------------------
    %  Set up the panel
    %----------------------------------------------------------------------
    uicontrol('Style', 'text', ...
              'String', sprintf('  To remove a node,\n  please enter these information.'), ...
              'Position', [15 645 271 60], ...
              'BackgroundColor', [0.3 0.4 0.7], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'text', ...
              'String', ' Index of the node:', ...
              'Position', [15 580 276 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(1) = ...
    uicontrol('Style', 'popupmenu', ...
              'String', listIndices('nodes'), ...
              'Position', [19 550 89 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'OK', ...
              'Position', [20 460 120 25], ...
              'FontSize', 10, ...
              'Callback', @ok_removeNode);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Return', ...
              'Position', [160 460 120 25], ...
              'FontSize', 10, ...
              'Callback', @close_prompt);
    
    % Set the focus to the first input box
    uicontrol(handle_edit(1));
end


function addElement(hObject, eventdata)
    global handle_gui handle_prompt handle_edit;
    global nodes elements;
    global numNodes numElements;
    
    % Find the current position of the GUI window
    gui_position = get(handle_gui, 'Position');
    gui_outerPosition = get(handle_gui, 'OuterPosition');
    
    %----------------------------------------------------------------------
    %  Initialize the prompt window
    %----------------------------------------------------------------------
    % If there is a prompt window already, close it
    if (ishandle(handle_prompt))
        close_prompt();
    end
    
    if (numNodes < 2)
        handle_prompt = ...
        figure('Name', 'Error', ...
               'Position', [0 0 340 80], ...
               'Color', [1 0.8 0.25], ...
               'MenuBar', 'none', ...
               'NumberTitle', 'off', ...
               'Resize', 'off', ...
               'ToolBar', 'none');
        movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
        
        if (numNodes == 0)
            uicontrol('Style', 'text', ...
                      'String', sprintf('Umm... there are no nodes.\nTo add an element, you need at least two nodes.'), ...
                      'Position', [0 22 340 50], ...
                      'BackgroundColor', [1 0.8 0.25], ...
                      'ForegroundColor', [0.8 0.1 0.1], ...
                      'FontSize', 10, ...
                      'HorizontalAlignment', 'center');
        else
            uicontrol('Style', 'text', ...
                      'String', sprintf('Umm... there is only one node.\nTo add an element, you need at least two nodes.'), ...
                      'Position', [0 22 340 50], ...
                      'BackgroundColor', [1 0.8 0.25], ...
                      'ForegroundColor', [0.8 0.1 0.1], ...
                      'FontSize', 10, ...
                      'HorizontalAlignment', 'center');
        end
        
        uicontrol('Style', 'pushbutton', ...
                  'String', 'OK', ...
                  'Position', [120 10 100 20], ...
                  'FontSize', 9, ...
                  'Callback', @close_prompt);
        
        return;
        
    elseif (numElements == nchoosek(numNodes, 2))
        handle_prompt = ...
        figure('Name', 'Error', ...
               'Position', [0 0 340 80], ...
               'Color', [1 0.8 0.25], ...
               'MenuBar', 'none', ...
               'NumberTitle', 'off', ...
               'Resize', 'off', ...
               'ToolBar', 'none');
        movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
        
        uicontrol('Style', 'text', ...
                  'String', sprintf('There are no more elements that you can create.\nPlease add a node or remove an element first.'), ...
                  'Position', [0 22 340 50], ...
                  'BackgroundColor', [1 0.8 0.25], ...
                  'ForegroundColor', [0.8 0.1 0.1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'center');
        
        uicontrol('Style', 'pushbutton', ...
                  'String', 'OK', ...
                  'Position', [120 10 100 20], ...
                  'FontSize', 9, ...
                  'Callback', @close_prompt);
        
        return;
    end
    
    handle_prompt = ...
    figure('Name', 'Prompt window', ...
           'Position', [0 0 300 745], ...
           'Color', [0.91 0.91 0.9], ...
           'MenuBar', 'none', ...
           'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'ToolBar', 'none');
    movegui(handle_prompt, [gui_position(1) + gui_outerPosition(3), gui_position(2)]);
    
    %----------------------------------------------------------------------
    %  Set up the panel
    %----------------------------------------------------------------------
    uicontrol('Style', 'text', ...
              'String', sprintf('  To add an element,\n  please enter these information.'), ...
              'Position', [15 645 271 60], ...
              'BackgroundColor', [0.3 0.4 0.7], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'text', ...
              'String', ' Indices of the two nodes:', ...
              'Position', [15 580 276 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(1) = ...
    uicontrol('Style', 'popupmenu', ...
              'String', listIndices('nodes'), ...
              'Value', 1, ...
              'Position', [19 550 89 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'text', ...
              'String', ',', ...
              'Position', [115 540 25 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(2) = ...
    uicontrol('Style', 'popupmenu', ...
              'String', listIndices('nodes'), ...
              'Value', 2, ...
              'Position', [141 550 89 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'text', ...
              'String', ' Young''s modulus (E):', ...
              'Position', [15 495 282 20], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(3) = ...
    uicontrol('Style', 'edit', ...
              'Position', [18 465 248 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'text', ...
              'String', ' Cross-sectional area (A):', ...
              'Position', [15 410 282 20], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(4) = ...
    uicontrol('Style', 'edit', ...
              'Position', [18 380 248 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'OK', ...
              'Position', [20 290 120 25], ...
              'FontSize', 10, ...
              'Callback', @ok_addElement);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Return', ...
              'Position', [160 290 120 25], ...
              'FontSize', 10, ...
              'Callback', @close_prompt);
end


function removeElement(hObject, eventdata)
    global handle_gui handle_prompt handle_edit;
    global elements;
    global numElements;
    
    % Find the current position of the GUI window
    gui_position = get(handle_gui, 'Position');
    gui_outerPosition = get(handle_gui, 'OuterPosition');
    
    %----------------------------------------------------------------------
    %  Initialize the prompt window
    %----------------------------------------------------------------------
    % If there is a prompt window already, close it
    if (ishandle(handle_prompt))
        close_prompt();
    end
    
    if (numElements == 0)
        handle_prompt = ...
        figure('Name', 'Error', ...
               'Position', [0 0 340 80], ...
               'Color', [1 0.8 0.25], ...
               'MenuBar', 'none', ...
               'NumberTitle', 'off', ...
               'Resize', 'off', ...
               'ToolBar', 'none');
        movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
        
        uicontrol('Style', 'text', ...
                  'String', sprintf('Umm... there are no elements.\nPlease add an element first.'), ...
                  'Position', [0 22 340 50], ...
                  'BackgroundColor', [1 0.8 0.25], ...
                  'ForegroundColor', [0.8 0.1 0.1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'center');
        
        uicontrol('Style', 'pushbutton', ...
                  'String', 'OK', ...
                  'Position', [120 10 100 20], ...
                  'FontSize', 9, ...
                  'Callback', @close_prompt);
        
        return;
    end
    
    handle_prompt = ...
    figure('Name', 'Prompt window', ...
           'Position', [0 0 300 745], ...
           'Color', [0.91 0.91 0.9], ...
           'MenuBar', 'none', ...
           'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'ToolBar', 'none');
    movegui(handle_prompt, [gui_position(1) + gui_outerPosition(3), gui_position(2)]);
    
    %----------------------------------------------------------------------
    %  Set up the panel
    %----------------------------------------------------------------------
    uicontrol('Style', 'text', ...
              'String', sprintf('  To remove an element,\n  please enter these information.'), ...
              'Position', [15 645 271 60], ...
              'BackgroundColor', [0.3 0.4 0.7], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'text', ...
              'String', ' Index of the element:', ...
              'Position', [15 580 276 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(1) = ...
    uicontrol('Style', 'popupmenu', ...
              'String', listIndices('elements'), ...
              'Position', [19 550 89 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'OK', ...
              'Position', [20 460 120 25], ...
              'FontSize', 10, ...
              'Callback', @ok_removeElement);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Return', ...
              'Position', [160 460 120 25], ...
              'FontSize', 10, ...
              'Callback', @close_prompt);
    
    % Set the focus to the first input box
    uicontrol(handle_edit(1));
end


function addBC(hObject, eventdata)
    global handle_gui handle_prompt handle_edit;
    global nodes BCs;
    global numNodes numBCs;
    
    % Find the current position of the GUI window
    gui_position = get(handle_gui, 'Position');
    gui_outerPosition = get(handle_gui, 'OuterPosition');
    
    %----------------------------------------------------------------------
    %  Initialize the prompt window
    %----------------------------------------------------------------------
    % If there is a prompt window already, close it
    if (ishandle(handle_prompt))
        close_prompt();
    end
    
    if (numNodes == 0)
        handle_prompt = ...
        figure('Name', 'Error', ...
               'Position', [0 0 340 80], ...
               'Color', [1 0.8 0.25], ...
               'MenuBar', 'none', ...
               'NumberTitle', 'off', ...
               'Resize', 'off', ...
               'ToolBar', 'none');
        movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
        
        uicontrol('Style', 'text', ...
                  'String', sprintf('Umm... there are no nodes.\nTo add a BC, please add a node first.'), ...
                  'Position', [0 22 340 50], ...
                  'BackgroundColor', [1 0.8 0.25], ...
                  'ForegroundColor', [0.8 0.1 0.1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'center');
        
        uicontrol('Style', 'pushbutton', ...
                  'String', 'OK', ...
                  'Position', [120 10 100 20], ...
                  'FontSize', 9, ...
                  'Callback', @close_prompt);
        
        return;
    elseif (numBCs == numNodes)
        handle_prompt = ...
        figure('Name', 'Error', ...
               'Position', [0 0 340 80], ...
               'Color', [1 0.8 0.25], ...
               'MenuBar', 'none', ...
               'NumberTitle', 'off', ...
               'Resize', 'off', ...
               'ToolBar', 'none');
        movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
        
        uicontrol('Style', 'text', ...
                  'String', sprintf('All the nodes have a BC specified.\nPlease remove a boundary condition first.'), ...
                  'Position', [0 22 340 50], ...
                  'BackgroundColor', [1 0.8 0.25], ...
                  'ForegroundColor', [0.8 0.1 0.1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'center');
        
        uicontrol('Style', 'pushbutton', ...
                  'String', 'OK', ...
                  'Position', [120 10 100 20], ...
                  'FontSize', 9, ...
                  'Callback', @close_prompt);
        
        return;
    end
    
    handle_prompt = ...
    figure('Name', 'Prompt window', ...
           'Position', [0 0 300 745], ...
           'Color', [0.91 0.91 0.9], ...
           'MenuBar', 'none', ...
           'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'ToolBar', 'none');
    movegui(handle_prompt, [gui_position(1) + gui_outerPosition(3), gui_position(2)]);
    
    %----------------------------------------------------------------------
    %  Set up the panel
    %----------------------------------------------------------------------
    uicontrol('Style', 'text', ...
              'String', sprintf('  To add a boundary condition,\n  please enter these information.'), ...
              'Position', [15 645 271 60], ...
              'BackgroundColor', [0.3 0.4 0.7], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'text', ...
              'String', ' Index of the node:', ...
              'Position', [15 580 276 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(1) = ...
    uicontrol('Style', 'popupmenu', ...
              'String', listIndices('nodesWithoutBCs'), ...
              'Value', 1, ...
              'Position', [19 550 89 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'text', ...
              'String', ' Type of boundary condition (x, y):', ...
              'Position', [15 495 282 20], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(2) = ...
    uicontrol('Style', 'popupmenu', ...
              'String', {'Displacement, Displacement'; 'Force, Force'; 'Displacement, Force'; 'Force, Displacement'}, ...
              'Position', [18 465 249 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'text', ...
              'String', ' Value of the x-component:', ...
              'Position', [15 410 282 20], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(3) = ...
    uicontrol('Style', 'edit', ...
              'Position', [18 380 248 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'text', ...
              'String', ' Value of the y-component:', ...
              'Position', [15 325 282 20], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(4) = ...
    uicontrol('Style', 'edit', ...
              'Position', [18 295 248 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'OK', ...
              'Position', [20 205 120 25], ...
              'FontSize', 10, ...
              'Callback', @ok_addBC);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Return', ...
              'Position', [160 205 120 25], ...
              'FontSize', 10, ...
              'Callback', @close_prompt);
end


function removeBC(hObject, eventdata)
    global handle_gui handle_prompt handle_edit;
    global BCs;
    global numBCs;
    
    % Find the current position of the GUI window
    gui_position = get(handle_gui, 'Position');
    gui_outerPosition = get(handle_gui, 'OuterPosition');
    
    %----------------------------------------------------------------------
    %  Initialize the prompt window
    %----------------------------------------------------------------------
    % If there is a prompt window already, close it
    if (ishandle(handle_prompt))
        close_prompt();
    end
    
    if (numBCs == 0)
        handle_prompt = ...
        figure('Name', 'Error', ...
               'Position', [0 0 340 80], ...
               'Color', [1 0.8 0.25], ...
               'MenuBar', 'none', ...
               'NumberTitle', 'off', ...
               'Resize', 'off', ...
               'ToolBar', 'none');
        movegui(handle_prompt, [gui_position(1) + (gui_position(3) - 340)/2, gui_position(2) + (gui_position(4) - 55)/2]);
        
        uicontrol('Style', 'text', ...
                  'String', sprintf('Umm... there are no boundary conditions.\nPlease add a boundary condition first.'), ...
                  'Position', [0 22 340 50], ...
                  'BackgroundColor', [1 0.8 0.25], ...
                  'ForegroundColor', [0.8 0.1 0.1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'center');
        
        uicontrol('Style', 'pushbutton', ...
                  'String', 'OK', ...
                  'Position', [120 10 100 20], ...
                  'FontSize', 9, ...
                  'Callback', @close_prompt);
        
        return;
    end
    
    handle_prompt = ...
    figure('Name', 'Prompt window', ...
           'Position', [0 0 300 745], ...
           'Color', [0.91 0.91 0.9], ...
           'MenuBar', 'none', ...
           'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'ToolBar', 'none');
    movegui(handle_prompt, [gui_position(1) + gui_outerPosition(3), gui_position(2)]);
    
    %----------------------------------------------------------------------
    %  Set up the panel
    %----------------------------------------------------------------------
    uicontrol('Style', 'text', ...
              'String', sprintf('  To remove a boundary condition,\n  please enter these information.'), ...
              'Position', [15 645 271 60], ...
              'BackgroundColor', [0.3 0.4 0.7], ...
              'ForegroundColor', [0.8 0.9 0.5], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'text', ...
              'String', ' Index of the node:', ...
              'Position', [15 580 276 25], ...
              'BackgroundColor', [0.91 0.91 0.9], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left');
    
    handle_edit(1) = ...
    uicontrol('Style', 'popupmenu', ...
              'String', listIndices('nodesWithBCs'), ...
              'Position', [19 550 89 25], ...
              'BackgroundColor', [0.9 0.92 0.95], ...
              'ForegroundColor', [0.1 0.1 0.1], ...
              'FontSize', 10, ...
              'HorizontalAlignment', 'left', ...
              'Enable', 'Inactive', ...
              'ButtonDownFcn', @changeColor);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'OK', ...
              'Position', [20 460 120 25], ...
              'FontSize', 10, ...
              'Callback', @ok_removeBC);
    
    uicontrol('Style', 'pushbutton', ...
              'String', 'Return', ...
              'Position', [160 460 120 25], ...
              'FontSize', 10, ...
              'Callback', @close_prompt);
end


%--------------------------------------------------------------------------
%  This routine allows us to change the color of the input boxes depending
%  on whether a box has been selected (is in focus).
%--------------------------------------------------------------------------
function changeColor(hObject, eventdata)
    global handle_edit;
    
    numHandles = size(handle_edit, 2);
    
    % Highlight the input box if selected
    index = find(hObject == handle_edit);
    set(handle_edit(index), 'BackgroundColor', [0.8 0.7 0.3], 'Enable', 'on');
    
    for i = 1 : numHandles
        % Un-highlight all the other input boxes
        if (i ~= index)
            set(handle_edit(i), 'BackgroundColor', [0.9 0.92 0.95], 'Enable', 'off');
        end
    end
    
    % Set the focus to the object
    uicontrol(hObject);
end


%--------------------------------------------------------------------------
%  Close the prompt window
%--------------------------------------------------------------------------
function close_prompt(hObject, eventdata)
    global handle_prompt handle_edit;
    
    close(handle_prompt);
    
    handle_prompt = [];
    handle_edit = [];
end


%--------------------------------------------------------------------------
%  These routines specify what to do with the user's input. They check for
%  the validity of the input (as much as possible), process the input, and
%  provide a feedback to the user, indicating success or failure.
%  
%  The checks that are common throughout the routines are done by calling
%  the corresponding checkinput routine. The checks that are unique to the
%  routine are specified within.
%--------------------------------------------------------------------------
function ok_addNode(hObject, eventdata)
    global handle_edit;
    global nodes elements;
    global numNodes numElements;
    
    % Retrieve the inputs
    x_str = get(handle_edit(1), 'String');
    y_str = get(handle_edit(2), 'String');
    x = str2num(x_str);
    y = str2num(y_str);
    
    %----------------------------------------------------------------------
    %  Check for validity
    %----------------------------------------------------------------------
    error_flag = 0;
    error_message = '';
    
    % Check the x-coordinate
    [flag, message] = checkinput_isReal('x-coordinate', x_str);
    error_flag = max(error_flag, flag);
    error_message = strcat(error_message, message);
    
    % Check if the x-coordinate is too big (this causes a problem with plotting)
    if (error_flag == 0)
        if (abs(x) >= 1e10)
            error_flag = 1;
            error_message = strcat(error_message, '  Error in the x-coordinate:\n   The coordinate is too big for plotting.\n   Please rescale it by changing the unit.\n\n');
        end
    end
    
    % Check the y-coordinate
    [flag, message] = checkinput_isReal('y-coordinate', y_str);
    error_flag = max(error_flag, flag);
    error_message = strcat(error_message, message);
    
    % Check if the y-coordinate is too big (this causes a problem with plotting)
    if (error_flag == 0)
        if (abs(y) >= 1e10)
            error_flag = 1;
            error_message = strcat(error_message, '  Error in the y-coordinate:\n   The coordinate is too big for plotting.\n   Please rescale it by changing the unit.\n\n');
        end
    end
    
    % Check if there is already a node with the given coordinates
    if (error_flag == 0)
        if (numNodes > 0)
            % If there is a match, the corresponding row in temp will be [1 1].
            % Otherwise, a row will be [0 0], [0 1], or [1 0].
            temp = (ones(numNodes, 1) * [x y] == nodes);
            
            if (sum(temp(:, 1) .* temp(:, 2)) > 0)
                error_flag = 1;
                error_message = strcat(error_message, '  Error in the node:\n   A node with the given coordinates\n   already exists.\n\n');
            end
        end
    end
    
    % Finally, check if the node lies on top of an element
    if (error_flag == 0)
        for i = 1 : numElements
            % Get the coordinates of the two nodes connecting the element
            x1 = nodes(elements(i, 1), 1);
            y1 = nodes(elements(i, 1), 2);
            x2 = nodes(elements(i, 2), 1);
            y2 = nodes(elements(i, 2), 2);
            
            % The coordinates of the given node
            x3 = x;
            y3 = y;
            
            % Check if the three points are collinear (warning: a false
            % positive is possible)
            if (abs((y2 - y1)*(x3 - x2) - (y3 - y2)*(x2 - x1)) < 1e-15)
                % Check if (x3, y3) is the middle point
                if ((min(x1, x2) <= x3 && x3 <= max(x1, x2)) && (min(y1, y2) <= y3 && y3 <= max(y1, y2)))
                    error_flag = 1;
                    error_message = strcat(error_message, '  Error in the node:\n   A node cannot lie on top of an element.\n\n');
                end
            end
            
            if (error_flag == 1)
                break;
            end
        end
    end
    
    %----------------------------------------------------------------------
    %  Respond accordingly
    %----------------------------------------------------------------------
    if (error_flag == 0)
        % Increment the number of nodes
        numNodes = numNodes + 1;
        
        % Add the node to the nodes array
        nodes(numNodes, 1) = x;
        nodes(numNodes, 2) = y;
        
        % Refresh the GUI
        refreshGUI();
        
        % Display feedback
        uicontrol('Style', 'text', ...
                  'String', sprintf('  Success:\n   The node has been added.'), ...
                  'Position', [15 25 271 160], ...
                  'BackgroundColor', [0.2 0.6 0.3], ...
                  'ForegroundColor', [0.8 0.9 0.5], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'left');
    else
        % Display feedback
        uicontrol('Style', 'text', ...
                  'String', sprintf(error_message), ...
                  'Position', [15 25 271 160], ...
                  'BackgroundColor', [0.85 0.2 0.25], ...
                  'ForegroundColor', [0.3 0.88 1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'left');
    end
end


function ok_editNode(hObject, eventdata)
    global handle_edit;
    global nodes elements;
    global numNodes numElements;
    
    % Retrieve the inputs
    index = get(handle_edit(1), 'Value');
    x_str = get(handle_edit(2), 'String');
    y_str = get(handle_edit(3), 'String');
    x = str2num(x_str);
    y = str2num(y_str);
    
    % Find all nodes that are not this node
    otherNodes = find(index ~= (1 : numNodes)');
    numOtherNodes = numNodes - 1;
    
    % Find all elements that are not connected to this node
    if (numElements > 0)
        % Check if an element consists of the given node.
        % If so, the corresponding row in temp will be [0 1] or [1 0].
        % Otherwise, a row will be [0 0].
        temp = (index * ones(numElements, 2) == elements(:, [1 2]));
        
        otherElements = elements(find(sum(temp, 2) == 0), [1 2]);
    else
        otherElements = [];
    end
    numOtherElements = size(otherElements, 1);
    
    %----------------------------------------------------------------------
    %  Check for validity
    %----------------------------------------------------------------------
    error_flag = 0;
    error_message = '';
    
    % Check the x-coordinate
    [flag, message] = checkinput_isReal('x-coordinate', x_str);
    error_flag = max(error_flag, flag);
    error_message = strcat(error_message, message);
    
    % Check if the x-coordinate is too big (this causes a problem with plotting)
    if (error_flag == 0)
        if (abs(x) >= 1e10)
            error_flag = 1;
            error_message = strcat(error_message, '  Error in the x-coordinate:\n   The coordinate is too big for plotting.\n   Please rescale it by changing the unit.\n\n');
        end
    end
    
    % Check the y-coordinate
    [flag, message] = checkinput_isReal('y-coordinate', y_str);
    error_flag = max(error_flag, flag);
    error_message = strcat(error_message, message);
    
    % Check if the y-coordinate is too big (this causes a problem with plotting)
    if (error_flag == 0)
        if (abs(y) >= 1e10)
            error_flag = 1;
            error_message = strcat(error_message, '  Error in the y-coordinate:\n   The coordinate is too big for plotting.\n   Please rescale it by changing the unit.\n\n');
        end
    end
    
    % Check if there is already a node with the given coordinates
    if (error_flag == 0)
        if (numOtherNodes > 0)
            temp = (ones(numOtherNodes, 1) * [x, y] == nodes(otherNodes, :));
            
            if (sum(temp(:, 1) .* temp(:, 2)) > 0)
                error_flag = 1;
                error_message = strcat(error_message, '  Error in the node:\n   A node cannot lie on top of another\n   node.\n\n');
            end
        end
    end
    
    % Finally, check if the node lies on top of an element
    if (error_flag == 0)
        for i = 1 : numOtherElements
            % Get the coordinates of the two nodes connecting the element
            x1 = nodes(otherElements(i, 1), 1);
            y1 = nodes(otherElements(i, 1), 2);
            x2 = nodes(otherElements(i, 2), 1);
            y2 = nodes(otherElements(i, 2), 2);
            
            % The coordinates of the given node
            x3 = x;
            y3 = y;
            
            % Check if the three points are collinear (warning: a false
            % positive is possible)
            if (abs((y2 - y1)*(x3 - x2) - (y3 - y2)*(x2 - x1)) < 1e-15)
                % Check if (x3, y3) is the middle point
                if ((min(x1, x2) <= x3 && x3 <= max(x1, x2)) && (min(y1, y2) <= y3 && y3 <= max(y1, y2)))
                    error_flag = 1;
                    error_message = strcat(error_message, '  Error in the node:\n   A node cannot lie on top of an element.\n\n');
                end
            end
            
            if (error_flag == 1)
                break;
            end
        end
    end
    
    %----------------------------------------------------------------------
    %  Respond accordingly
    %----------------------------------------------------------------------
    if (error_flag == 0)
        % Add the node to the nodes array
        nodes(index, 1) = x;
        nodes(index, 2) = y;
        
        % Refresh the GUI
        refreshGUI();
        
        % Display feedback
        uicontrol('Style', 'text', ...
                  'String', sprintf('  Success:\n   The node has been updated.'), ...
                  'Position', [15 25 271 160], ...
                  'BackgroundColor', [0.2 0.6 0.3], ...
                  'ForegroundColor', [0.8 0.9 0.5], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'left');
    else
        % Display feedback
        uicontrol('Style', 'text', ...
                  'String', sprintf(error_message), ...
                  'Position', [15 25 271 160], ...
                  'BackgroundColor', [0.85 0.2 0.25], ...
                  'ForegroundColor', [0.3 0.88 1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'left');
    end
end


function ok_removeNode(hObject, eventdata)
    global handle_edit;
    global nodes elements BCs;
    global numNodes numElements numBCs;
    
    % Retrieve the inputs
    index = get(handle_edit(1), 'Value');
    
    %----------------------------------------------------------------------
    %  Respond accordingly (no possibility of a user input error)
    %----------------------------------------------------------------------
    % Keep the nodes that do not match the given node index
    nodes = nodes(find(index ~= (1 : numNodes)'), :);
    
    % Check if there are elements that need to be removed along with the
    % node
    if (numElements > 0)
        % If there are, the corresponding row in temp will be [0 1] or [1 0].
        % Otherwise, a row will be [0 0].
        temp = (index == elements(:, [1 2]));
        
        % The elements to keep correspond to the rows in temp that have
        % a column sum of 0
        elements = elements(find(sum(temp, 2) == 0), :);
        
        % Don't forget to reassign the nodes in the elements. Because
        % the node indices are assumed to be consecutive integers,
        % we just need to subtract the node index by 1 if it was
        % previously larger than the input index
        elements(:, [1 2]) = elements(:, [1 2]) - (elements(:, [1 2]) > index);
    end
    
    % Check if there are any BCs that need to be removed
    if (numBCs > 0)
        BCs = BCs(find(index ~= BCs(:, 1)), :);
        
        % Don't forget to reassign the nodes in the BCs
        BCs(:, 1) = BCs(:, 1) - (BCs(:, 1) > index);
    end
    
    % Update the number of nodes, elements, and BCs
    numNodes = numNodes - 1;
    numElements = size(elements, 1);
    numBCs = size(BCs, 1);
    
    % Refresh the GUI
    refreshGUI();
    
    % If there are no more nodes that can be removed, close the prompt window
    if (numNodes == 0)
        close_prompt();
    else
        % Refresh the prompt window
        set(handle_edit(1), 'String', listIndices('nodes'), 'Value', min(index, numNodes));
        
        % Display feedback
        uicontrol('Style', 'text', ...
                  'String', sprintf('  Success:\n   The node has been removed.'), ...
                  'Position', [15 25 271 160], ...
                  'BackgroundColor', [0.2 0.6 0.3], ...
                  'ForegroundColor', [0.8 0.9 0.5], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'left');
    end
end


function ok_addElement(hObject, eventdata)
    global handle_edit;
    global nodes elements;
    global numNodes numElements;
    
    % Retrieve the inputs
    index1 = get(handle_edit(1), 'Value');
    index2 = get(handle_edit(2), 'Value');
    E_str = get(handle_edit(3), 'String');
    A_str = get(handle_edit(4), 'String');
    
    %----------------------------------------------------------------------
    %  Check for validity
    %----------------------------------------------------------------------
    error_message = '';
    error_flag = 0;
    
    if (index1 == index2)
        error_flag = 1;
        error_message = strcat(error_message, '  Error in the node indices:\n   The two nodes must be distinct.\n\n');
    elseif (numElements > 0)
        % Check if there is already an element connected by the nodes
        % index1 and index2. If so, the corresponding row in temp will
        % be [1 1]. Otherwise, a row will be [0 0], [0 1], or [1 0].
        temp1 = (ones(numElements, 1) * [index1 index2] == elements(:, [1 2]));
        temp2 = (ones(numElements, 1) * [index2 index1] == elements(:, [1 2]));
        
        if (sum(temp1(:, 1) .* temp1(:, 2)) > 0 || sum(temp2(:, 1) .* temp2(:, 2)) > 0)
            error_flag = 1;
            error_message = strcat(error_message, '  Error in the node indices:\n   An element with the given node indices\n   already exists.\n\n');
        end
    end
    
    % Finally, check if the element lies on top of a node
    if (error_flag == 0)
        % Get the coordinates of the two nodes
        x1 = nodes(index1, 1);
        y1 = nodes(index1, 2);
        x2 = nodes(index2, 1);
        y2 = nodes(index2, 2);
        
        for i = 1 : numNodes
            if (i ~= index1 && i ~= index2)
                % Get the coordinates of the third node
                x3 = nodes(i, 1);
                y3 = nodes(i, 2);
                
                % Check if the three points are collinear (warning: a false
                % positive is possible)
                if (abs((y2 - y1)*(x3 - x2) - (y3 - y2)*(x2 - x1)) < 1e-15)
                    % Check if (x3, y3) is the middle point
                    if ((min(x1, x2) <= x3 && x3 <= max(x1, x2)) && (min(y1, y2) <= y3 && y3 <= max(y1, y2)))
                        error_flag = 1;
                        error_message = strcat(error_message, '  Error in the node:\n   A node cannot lie on top of an element.\n\n');
                    end
                end
                
                if (error_flag == 1)
                    break;
                end
            end
        end
    end
    
    % Check the Young's modulus
    [flag, message] = checkinput_isRealPositive('Young''s modulus', E_str);
    error_flag = max(error_flag, flag);
    error_message = strcat(error_message, message);
    
    % Check the cross-sectional area
    [flag, message] = checkinput_isRealPositive('cross-sectional area', A_str);
    error_flag = max(error_flag, flag);
    error_message = strcat(error_message, message);
    
    %----------------------------------------------------------------------
    %  Respond accordingly
    %----------------------------------------------------------------------
    if (error_flag == 0)
        % Increment the number of elements
        numElements = numElements + 1;
        
        % Add the node to the nodes array
        elements(numElements, 1) = index1;
        elements(numElements, 2) = index2;
        elements(numElements, 3) = str2num(E_str);
        elements(numElements, 4) = str2num(A_str);
        
        % Refresh the GUI
        refreshGUI();
        
        % If there are no more elements that can be added, close the prompt
        % window
        if (numElements == nchoosek(numNodes, 2))
            close_prompt();
        else
            % Display feedback
            uicontrol('Style', 'text', ...
                      'String', sprintf('  Success:\n   The element has been added.'), ...
                      'Position', [15 25 271 160], ...
                      'BackgroundColor', [0.2 0.6 0.3], ...
                      'ForegroundColor', [0.8 0.9 0.5], ...
                      'FontSize', 10, ...
                      'HorizontalAlignment', 'left');
        end
    else
        % Display feedback
        uicontrol('Style', 'text', ...
                  'String', sprintf(error_message), ...
                  'Position', [15 25 271 160], ...
                  'BackgroundColor', [0.85 0.2 0.25], ...
                  'ForegroundColor', [0.3 0.88 1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'left');
    end
end


function ok_removeElement(hObject, eventdata)
    global handle_edit;
    global elements;
    global numElements;
    
    % Retrieve the user's inputs
    index = get(handle_edit(1), 'Value');
    
    %----------------------------------------------------------------------
    %  Respond accordingly (no possibility of a user input error)
    %----------------------------------------------------------------------
    % Keep the nodes that do not match the given element index
    elements = elements(find(index ~= (1 : numElements)'), :);
    
    % Update the number of elements
    numElements = numElements - 1;
    
    % Refresh the GUI
    refreshGUI();
    
    % If there are no more elements that can be removed, close the prompt
    % window
    if (numElements == 0)
        close_prompt();
    else
        % Refresh the prompt window
        set(handle_edit(1), 'String', listIndices('elements'), 'Value', min(index, numElements));
        
        % Display feedback
        uicontrol('Style', 'text', ...
                  'String', sprintf('  Success:\n   The element has been removed.'), ...
                  'Position', [15 25 271 160], ...
                  'BackgroundColor', [0.2 0.6 0.3], ...
                  'ForegroundColor', [0.8 0.9 0.5], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'left');
    end
end


function ok_addBC(hObject, eventdata)
    global handle_edit;
    global nodes BCs;
    global numNodes numBCs;
    
    % Find the nodes without a BC specified
    nodesWithoutBCs = get(handle_edit(1), 'String');
    
    % Retrieve the inputs
    index = str2num(nodesWithoutBCs{get(handle_edit(1), 'Value')});
    BCtype = get(handle_edit(2), 'Value');
    x_str = get(handle_edit(3), 'String');
    y_str = get(handle_edit(4), 'String');
    x = str2num(x_str);
    y = str2num(y_str);
    
    %----------------------------------------------------------------------
    %  Check for validity
    %----------------------------------------------------------------------
    error_message = '';
    error_flag = 0;
    
    % Check the x-component value
    [flag, message] = checkinput_isReal('x-component value', x_str);
    error_flag = max(error_flag, flag);
    error_message = strcat(error_message, message);
    
    % Check the y-component value
    [flag, message] = checkinput_isReal('y-component value', y_str);
    error_flag = max(error_flag, flag);
    error_message = strcat(error_message, message);
    
    %----------------------------------------------------------------------
    %  Respond accordingly
    %----------------------------------------------------------------------
    if (error_flag == 0)
        % Increment the number of BCs
        numBCs = numBCs + 1;
        
        % Add the node to the nodes array
        BCs(numBCs, 1) = index;
        BCs(numBCs, 2) = BCtype;
        BCs(numBCs, 3) = x;
        BCs(numBCs, 4) = y;
        
        % Sort the BCs in ascending node index
        if (numBCs > 1)
            [temp, permutation] = sort(BCs(:, 1));
            BCs = BCs(permutation, :);
        end
        
        % Refresh the GUI
        refreshGUI();
        
        % If there are no more boundary conditions that can be added, close
        % the prompt window
        if (numBCs == numNodes)
            close_prompt();
        else
            % Refresh the prompt window
            set(handle_edit(1), 'String', listIndices('nodesWithoutBCs'), 'Value', 1);
            
            % Display feedback
            uicontrol('Style', 'text', ...
                      'String', sprintf('  Success:\n   The boundary condition has been added\n   to the node.'), ...
                      'Position', [15 25 271 160], ...
                      'BackgroundColor', [0.2 0.6 0.3], ...
                      'ForegroundColor', [0.8 0.9 0.5], ...
                      'FontSize', 10, ...
                      'HorizontalAlignment', 'left');
        end
    else
        % Display feedback
        uicontrol('Style', 'text', ...
                  'String', sprintf(error_message), ...
                  'Position', [15 25 271 160], ...
                  'BackgroundColor', [0.85 0.2 0.25], ...
                  'ForegroundColor', [0.3 0.88 1], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'left');
    end
end


function ok_removeBC(hObject, eventdata)
    global handle_edit;
    global BCs;
    global numBCs;
    
    % Find the nodes without a BC specified
    nodesWithBCs = get(handle_edit(1), 'String');
    
    % Find the row that corresponds to the node to be removed
    index = find(str2num(nodesWithBCs{get(handle_edit(1), 'Value')}) == BCs(:, 1));
    
    %----------------------------------------------------------------------
    %  Respond accordingly (no possibility of a user input error)
    %----------------------------------------------------------------------
    % Keep the boundary conditions for all the other nodes
    BCs = BCs(find(index ~= (1 : numBCs)'), :);
    
    % Update the number of BCs
    numBCs = numBCs - 1;
    
    % Refresh the GUI
    refreshGUI();
    
    % If there are no more elements, close the prompt window
    if (numBCs == 0)
        close_prompt();
    else
        % Refresh the prompt window
        set(handle_edit(1), 'String', listIndices('nodesWithBCs'), 'Value', 1);
        
        % Display feedback
        uicontrol('Style', 'text', ...
                  'String', sprintf('  Success:\n   The boundary condition at the node\n   has been removed.'), ...
                  'Position', [15 25 271 160], ...
                  'BackgroundColor', [0.2 0.6 0.3], ...
                  'ForegroundColor', [0.8 0.9 0.5], ...
                  'FontSize', 10, ...
                  'HorizontalAlignment', 'left');
    end
end


%--------------------------------------------------------------------------
%  These routines check an input against various criteria.
%--------------------------------------------------------------------------
function [flag, message] = checkinput_isIndex(name, value_str)
    % Convert the input from string to a number
    value = str2num(value_str);
    
    % Assume that there is an error
    flag = 1;
    message = ['  Error in the ', name, ':\n   '];
    
    % Check these criteria
    if (isempty(value_str))
        message = [message, 'Please enter the index.\n\n'];
    elseif (isempty(value))
        message = [message, 'The input must be a positive integer.\n\n'];
    elseif (~isscalar(value))
        message = [message, 'The input must be a scalar, and cannot \n   be a vector or a matrix.\n\n'];
    elseif (~isreal(value))
        message = [message, 'The input must be a positive integer.\n\n'];
    elseif (isinf(value) || isnan(value))
        message = [message, 'The input must be a positive integer.\n\n'];
    elseif ((mod(value, 1) ~= 0) || (value <= 0))
        message = [message, 'The input must be a positive integer.\n\n'];
    else
        flag = 0;
        message = '';
    end
end


function [flag, message] = checkinput_isReal(name, value_str)
    % Convert the input from string to a number
    value = str2num(value_str);
    
    % Assume that there is an error
    flag = 1;
    message = ['  Error in the ', name, ':\n   '];
    
    % Check these criteria
    if (isempty(value_str))
        message = [message, 'Please enter the ', name, '.\n\n'];
    elseif (isempty(value))
        message = [message, 'The input must be a real number.\n\n'];
    elseif (~isscalar(value))
        message = [message, 'The input must be a scalar, and cannot\n   be a vector or a matrix.\n\n'];
    elseif (~isreal(value))
        message = [message, 'The input must be a real number.\n\n'];
    elseif (isinf(value) || isnan(value))
        message = [message, 'The input must be a real number.\n\n'];
    else
        flag = 0;
        message = '';
    end
end


function [flag, message] = checkinput_isRealPositive(name, value_str)
    % Convert the input from string to a number
    value = str2num(value_str);
    
    % Assume that there is an error
    flag = 1;
    message = ['  Error in the ', name, ':\n   '];
    
    % Check these criteria
    if (isempty(value_str))
        message = [message, 'Please enter the ', name, '.\n\n'];
    elseif (isempty(value))
        message = [message, 'The input must be a positive number.\n\n'];
    elseif (~isscalar(value))
        message = [message, 'The input must be a scalar, and cannot\n   be a vector or a matrix.\n\n'];
    elseif (value <= 0)
        message = [message, 'The input must be a positive number.\n\n'];
    elseif (~isreal(value))
        message = [message, 'The input must be a positive number.\n\n'];
    elseif (isinf(value) || isnan(value))
        message = [message, 'The input must be a positive number.\n\n'];
    else
        flag = 0;
        message = '';
    end
end
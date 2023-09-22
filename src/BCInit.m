classdef BCInit < handle
    properties
        inputReader
        mesh
    end
    
    properties (Access = private)
        utils
        elements
        meshFileName
        initialdofs
        loadVector
        integrationFunction
    end
    
    methods
    function obj = BCInit(inputReader, mesh)
        obj.inputReader_ = inputReader;
        obj.mesh_ = mesh;
        obj.elements_ = Elements(); % Assuming Elements is a class in your project
        obj.loadVector_ = zeros(2 * mesh.getNumAllNodes(), 1);
        
        % Initialize initialdofs_ based on the number of nodes
        obj.initialdofs_ = zeros(2 * mesh.getNumAllNodes(), 1);
        
        % Call boundaryConditions to perform any necessary setup
        obj.boundaryConditions();
    end
        
        function boundaryConditions(obj)
            % Implement the boundaryConditions method here
        end
        
        function createLoadVector(obj)
            % Implement the createLoadVector method here
        end
        
        function initialDof = getInitialDof(obj, i)
            initialDof = obj.initialdofs(i);
        end
        
        function allInitialDof = getAllInitialDof(obj)
            allInitialDof = obj.initialdofs;
        end
        
        function loadVec = getLoadVector(obj)
            loadVec = obj.loadVector;
        end
        
        function F_q = CteSurfBC(obj, natcoords, coords, value, element)

            % Define variables
            nodesperelement_etype = obj.mesh_.getNumNodesForElement(element);
            nodes_per_element = nodesperelement_etype(1); % Number of nodes
            etype = nodesperelement_etype(2); % Element type
            shapeFunctions = zeros(nodes_per_element, 1);          % Shape functions as a column vector
            shapeFunctionDerivatives = zeros(nodes_per_element, 2); % Shape function derivatives
            
            % Extract natural coordinates
            xi = 0;
            eta = 0;
            if size(natcoords, 1) >= 2 && size(natcoords, 2) >= 1
                xi = natcoords(1, 1);  % Extract the first element (a)
                eta = natcoords(2, 1); % Extract the second element (b)
                % Rest of your code...
            else
                % Handle the case when natcoords doesn't have the expected dimensions.
                disp('Wrong natural coordinates dimension.');
                % You may want to print an error message or take appropriate action.
                return;
            end
            
            % Calculate shape functions and their derivatives
            obj.mesh_.selectShapeFunctionsAndDerivatives(etype, xi, eta, -1, shapeFunctions, shapeFunctionDerivatives);
        
            % Calculate Jacobian matrix JM
            JM = shapeFunctionDerivatives' * coords';
        
            % Calculate the determinant of the Jacobian
            detJ = det(JM);
        
            % Calculate F_q (heat flow) using your equations
            F_q = detJ * (shapeFunctions * value);
        
            % Return the heat flow as a column vector
            return;

        end

        function result = gaussIntegrationBC(obj, dimension, order, elementTag, bcvalue, func)
            % Check for invalid dimension or order
            if dimension < 1 || order < 1
                fprintf('Invalid dimension or order for Gauss integration.\n');
                result = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                return;
            end

            % Get the number of nodes per element and element type from the mesh
            nodesperelement_etype = obj.mesh_.getNumNodesForElement(elementTag);
            nodes_per_element = nodesperelement_etype(1); % Number of nodes
            etype = nodesperelement_etype(2); % Element type

            % Initialize variables
            weights = [];
            gaussPoints = [];
            nodeTags_el = obj.mesh_.getElementInfo(elementTag);
            coordinates = obj.mesh_.getCoordinates(nodeTags_el);
            coordinates_tr = obj.utils_.calculate_inverse_T3(coordinates);
            coordinates_tr = coordinates_tr * coordinates;

            % Initialize coordinates_tr_XY as an empty matrix
            coordinates_tr_XY = [];

            % Define tolerance
            tolerance = 1e-6;

            % Iterate through rows of coordinates_tr and copy rows that do not meet the tolerance condition
            for i = 1:size(coordinates_tr, 1)
                removeRow = true; % Assume we want to remove the row by default
                firstValue = coordinates_tr(i, 1); % Store the first value in the row

                % Check if all values in the row are within the tolerance
                for j = 2:size(coordinates_tr, 2)
                    if abs(coordinates_tr(i, j) - firstValue) > tolerance
                        removeRow = false; % Row contains a different value outside the tolerance, so we keep it
                        break;
                    end
                end

                if ~removeRow
                    % Copy the row to coordinates_tr_XY
                    if isempty(coordinates_tr_XY)
                        coordinates_tr_XY = coordinates_tr(i, :);
                    else
                        coordinates_tr_XY = [coordinates_tr_XY; coordinates_tr(i, :)];
                    end
                end
            end

            if isempty(coordinates_tr_XY)
                % Handle the case when all rows were removed, resulting in an empty matrix.
                fprintf('Warning: All rows removed. coordinates_tr_XY is empty.\n');
                % You may want to print a warning message or take appropriate action.
            end

            % Get the Gauss weights and points for the specified order
            obj.utils_.getGaussWeightsAndPoints(order, weights, gaussPoints);

            if isempty(weights) || isempty(gaussPoints)
                fprintf('Invalid order for Gauss integration.\n');
                result = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                return;
            end

            if size(weights, 1) ~= size(gaussPoints, 1)
                fprintf('Weights and Gauss points have mismatched dimensions.\n');
                result = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                return;
            end

            % Initialize result to zero matrix of appropriate dimensions
            result = zeros(dimension, dimension);

            if dimension == 1
                % 1D integration using a single loop
                natcoords = zeros(1, 1);
                f = zeros(nodes_per_element, 1);
                for i = 1:size(weights, 1)
                    f = func(natcoords, coordinates_tr_XY, bcvalue, elementTag) * weights(i);
                    result = result + f;
                end
            elseif dimension == 2
                % 2D integration using a double loop
                natcoords = zeros(2, 1);
                f = zeros(nodes_per_element, 1);
                for i = 1:size(weights, 1)
                    for j = 1:size(weights, 1)
                        natcoords(1) = gaussPoints(i);
                        natcoords(2) = gaussPoints(j);
                        f = func(natcoords, coordinates_tr_XY, bcvalue, elementTag) * (weights(i) * weights(j));
                        result = result + f;
                    end
                end
            elseif dimension == 3
                if order == 14
                    % Special case for 3D integration with order 14
                    result = weights * weights' .* weights * weights' .* weights * weights' .* func(gaussPoints, coordinates_tr_XY, bcvalue, elementTag);
                else
                    % Generic 3D integration using a triple loop
                    natcoords = zeros(3, 1);
                    f = zeros(nodes_per_element, 1);
                    for i = 1:size(weights, 1)
                        for j = 1:size(weights, 1)
                            for k = 1:size(weights, 1)
                                natcoords(1) = gaussPoints(i);
                                natcoords(2) = gaussPoints(j);
                                natcoords(3) = gaussPoints(k);
                                f = func(natcoords, coordinates_tr_XY, bcvalue, elementTag) * (weights(i) * weights(j) * weights(k));
                                result = result + f;
                            end
                        end
                    end
                end
            else
                fprintf('Invalid dimension for Gauss integration.\n');
                result = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
            end

        end    

    end
end

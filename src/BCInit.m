classdef BCInit < handle
    properties
        inputReader_
        mesh_
        elements_
        loadVector_
        initialdofs_
    end
    
    properties (Access = private)
        integrationFunction
        numNodes
        numDofs
        dofspernode
    end
    
    methods
    function obj = BCInit(inputReader, mesh)
        obj.inputReader_ = inputReader;
        obj.mesh_ = mesh;
        obj.elements_ = Elements(); % Assuming Elements is a class in your project
        numNodes=size(mesh.data.NODE);
        numNodes=numNodes(2);
        dofspernode=2;
        numDofs=numNodes*dofspernode;
        obj.loadVector_ = zeros(numDofs, 1);
        
        % Initialize initialdofs_ based on the number of nodes
        obj.initialdofs_ = zeros(numDofs, 1);
        
        % Call boundaryConditions to perform any necessary setup
        obj.boundaryConditions();
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function T3 = calculate_T3(nodes)
            % Check if the input matrix is 3x(nodes per element)
            if size(nodes, 1) ~= 3
                error('Input matrix must be 3x(nodes per element)');
            end
        
            % Extract node coordinates
            V1 = nodes(:, 1);
            V2 = nodes(:, 2);
            V3 = nodes(:, 3);
        
            % Calculate the unit normal vector XN
            V12 = V2 - V1;
            V31 = V3 - V1;
            XN = V12 / norm(V12);
            lx = XN(1);
            mx = XN(2);
            nx = XN(3);
        
            % Calculate the unit normal vector ZN
            ZN = cross(V12, V31) / norm(cross(V12, V31));
            lz = ZN(1);
            mz = ZN(2);
            nz = ZN(3);
        
            % Calculate ly, my, and ny
            ly = mz * nx - nz * mx;
            my = nz * lx - lz * nx;
            ny = lz * mx - mz * lx;
        
            % Create and return the transformation matrix T3
            T3 = [lx, mx, nx;
                  ly, my, ny;
                  lz, mz, nz];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function result = gaussIntegrationBC(dimension, order, elementTag, mesh, bcvalue, func)
            if dimension < 1 || order < 1
                fprintf('Invalid dimension or order for Gauss integration.\n');
                result = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                return;
            end
            
            mesh.getElementInfo(elementTag, nodeTags_el);
            element_nodes = mesh.data.ELEMENTS{elementTag};
            coordinates = zeros(3,length(element_nodes));
            for element = element_nodes
                coordinates(:,element)=mesh.data.NODE{element};
            end
            iT3 = Utils.calculate_inverse_T3(coordinates);
            coordinates_tr = iT3 * coordinates;
            
            % Create a new matrix to store the result without the last row
            coordinates_tr_XY = [];
            tolerance = 1e-6; % Define your tolerance here
            
            % Iterate through the rows of coordinates_tr and copy rows that do not meet the tolerance condition
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
                        coordinates_tr_XY = vertcat(coordinates_tr_XY, coordinates_tr(i, :));
                    end
                end
            end
            
            if isempty(coordinates_tr_XY)
                % Handle the case when all rows were removed, resulting in an empty matrix.
                fprintf('Warning: All rows removed. coordinates_tr_XY is empty.\n');
                % You may want to print a warning message or take appropriate action.
            end
           
            [weights, gaussPoints] = getGaussWeightsAndPoints(order);
            
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
            
            if dimension == 1
                % 1D integration using a single loop.
                natcoords = zeros(1, 1);
                for i = 1:size(weights, 1)
                    natcoords(1) = gaussPoints(i);
                    % Explicitly use the element-wise multiplication .* for arrays
                    f = func(natcoords, coordinates_tr_XY, bcvalue, elementTag, mesh) .* weights(i);
                    result = result + f;
                end
            elseif dimension == 2
                % 2D integration using a double loop.
                natcoords = zeros(2, 1);
                for i = 1:size(weights, 1)
                    for j = 1:size(weights, 1)
                        natcoords(1) = gaussPoints(i);
                        natcoords(2) = gaussPoints(j);
                        % Explicitly use the element-wise multiplication .* for arrays
                        f = func(natcoords, coordinates_tr_XY, bcvalue, elementTag, mesh) .* (weights(i) * weights(j));
                        result = result + f;
                    end
                end
            elseif dimension == 3
                if order == 14
                    % Special case for 3D integration with order 14.
                    result = weights * weights' .* weights * weights' .* weights * weights' .* func(gaussPoints, coordinates_tr_XY, bcvalue, elementTag, mesh);
                else
                    % Generic 3D integration using a triple loop.
                    natcoords = zeros(3, 1);
                    for i = 1:size(weights, 1)
                        for j = 1:size(weights, 1)
                            for k = 1:size(weights, 1)
                                natcoords(1) = gaussPoints(i);
                                natcoords(2) = gaussPoints(j);
                                natcoords(3) = gaussPoints(k);
                                % Explicitly use the element-wise multiplication .* for arrays
                                f = func(natcoords, coordinates_tr_XY, bcvalue, elementTag, mesh) .* (weights(i) * weights(j) * weights(k));
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function F_q = CteSurfBC(natcoords, coords, value, element, mesh)
            % Define variables
            etype=mesh.ElementTypes{element};
        
            % Extract natural coordinates
            xi = 0;
            eta = 0;
            if size(natcoords, 1) >= 2 && size(natcoords, 2) >= 1
                xi = natcoords(1);  % Extracts the first element (a)
                eta = natcoords(2); % Extracts the second element (b)
            else
                % Handle the case when natcoords doesn't have the expected dimensions.
                fprintf('Wrong natural coordinates dimension.\n');
                % You may want to print an error message or take appropriate action.
                return;
            end
        
            % Calculate shape functions and their derivatives
            [shapeFunctions,shapeFunctionDerivatives]=mesh.selectShapeFunctionsAndDerivatives(etype, xi, eta, -1);
        
            % Calculate Jacobian matrix JM
            JM = shapeFunctionDerivatives' * coords';
        
            % Calculate the determinant of the Jacobian
            detJ = det(JM);
        
            % Calculate F_q (heat flow) using your equations
            F_q = detJ * (shapeFunctions * value);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function boundaryConditions(obj)
            fprintf('#BCInit::boundaryConditions\n');
            
            % Get the boundary conditions from the InputReader
            boundaryConditions = obj.inputReader_.bctype;
            
            % Iterate through each boundary condition
            for k = 1:length(boundaryConditions)
                bc = boundaryConditions{k};
                boundaryName = bc;
                surfaceName = obj.inputReader_.bcloc{k};
                value = obj.inputReader_.bcval(k);
                
                fprintf('Boundary Name: %s, Surface Name: %s, Value: %f\n', boundaryName, surfaceName, value);
                
                if strcmp(boundaryName, 'Temperature')
                    % Assuming 'getNodesForPhysicalGroup' returns a vector of unsigned long long integers
                    nodes=obj.mesh_.data.NSET{strcmp(obj.mesh_.data.NodalSelectionNames, surfaceName)};
                    
                    % Loop through the nodes and set the corresponding values in initialdofs_ (loadVector_)
                    for node = nodes
                        if (node*2) < length(obj.initialdofs_)
                            obj.initialdofs_(node*2) = value;
                            obj.setFixedof(node*2);
                        else
                            % Handle the case where the node index is out of bounds.
                            % This could be an error condition depending on your application.
                            fprintf('Node index out of bounds: %d\n', node);
                        end
                    end
                elseif strcmp(boundaryName, 'Voltage')
                    % Assuming 'getNodesForPhysicalGroup' returns a vector of unsigned long long integers
                    nodes=obj.mesh_.data.NSET{strcmp(obj.mesh_.data.NodalSelectionNames, surfaceName)};
                    
                    % Loop through the nodes and set the corresponding values in initialdofs_ (loadVector_)
                    for node = nodes
                        if (node*2+1) < length(obj.initialdofs_)
                            obj.initialdofs_(node*2+1) = value;
                            obj.mesh_.setFixedof(node*2+1);
                        else
                            % Handle the case where the node index is out of bounds.
                            % This could be an error condition depending on your application.
                            fprintf('Node index out of bounds: %d\n', node);
                        end
                    end
                elseif strcmp(boundaryName, 'heat_n')
                    % Assuming 'getElementsAndNodeTagsForPhysicalGroup' returns a cell array of vectors
                    elementindexVector=obj.mesh_.data.ELSET{strcmp(obj.mesh_.data.ElementSelectionNames, surfaceName)};
                    
                    % Debugging: Print the sizes of elements and element_nodes vectors
                    fprintf('Size of ''elements'' vector: %d\n', length(elementindexVector));
                    
                    fprintf('HEAT INTEGRATION.\n');
                    % Initialize temporary load vectors for each worker
                    numWorkers = numel(gcp); % Get the number of workers
                    tempLoadVectors = cell(1, numWorkers);
                    
                    % Create a function handle for gaussIntegrationBC and CteSurfBC
                    %gaussIntegrationBCFun = @(element) gaussIntegrationBC(2, 3, element, value, CteSurfBC);
                    gaussIntegrationBCFun = @(element) gaussIntegrationBC(2, 3, element, value, CteSurfBC,obj.mesh_);
                    
                    % Loop through elements in parallel
                    parfor workerIdx = 1:numWorkers
                        % Initialize temporary load vector for this worker
                        tempLoadVector = zeros(size(obj.loadVector_));
                        
                        % Calculate elements for this worker
                        workerElements = elementindexVector(workerIdx:numWorkers:end);
                        
                        % Loop through elements for this worker
                        for elementIdx = workerElements
                            % Calculate element_load_vector for the current element
                            element_load_vector = gaussIntegrationBCFun(elementIdx);
                            
                            % Loop through element nodes and update temporary loadVector_
                            for i = (elementIdx - 1) * nodes_per_element + 1 : elementIdx * nodes_per_element
                                if i <= length(element_nodes)
                                    node = element_nodes(i);
                                    if (node * 2) <= length(tempLoadVector) % debugging
                                        tempLoadVector(node * 2) = tempLoadVector(node * 2) + element_load_vector(i - (elementIdx - 1) * nodes_per_element);
                                    else
                                        fprintf('Error: Node index out of bounds!\n');
                                    end
                                else
                                    fprintf('Error: Element node index out of bounds!\n');
                                end
                            end
                        end
                        
                        % Store the temporary load vector in the cell array
                        tempLoadVectors{workerIdx} = tempLoadVector;
                    end
                    
                    % Accumulate the results from all workers into the final loadVector_
                    for workerIdx = 1:numWorkers
                        obj.loadVector_ = obj.loadVector_ + tempLoadVectors{workerIdx};
                    end
                    fprintf('HEAT INTEGRATION FINISHED.\n');
                end
            end
            
            obj.mesh_.SetFreedofsIdx(); % store the index of the dofs that are not fixed in order
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end
end

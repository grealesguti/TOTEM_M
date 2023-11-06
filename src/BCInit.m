classdef BCInit < handle
    properties
        inputReader_
        mesh_
        elements_
        loadVector_
        initialdofs_
        dofs_fixed_
        numDofs
        numNodes
        dofspernode
        dofs_free_
        initialdofs_mech_
        dofs_fixed_mech_
        loadVector_mech
        numDofs_mech
        dofs_free_mech
        dim
    end
    
    properties (Access = private)
        integrationFunction
    end
    
    methods
    function obj = BCInit(inputReader, mesh)
        obj.elements_ = Elements(); % Assuming Elements is a class in your project

        obj.numNodes=length(mesh.data.NODE);
        obj.dofspernode=2;
        obj.numDofs=obj.numNodes*obj.dofspernode;

        % thermoelectricity
        obj.loadVector_ = zeros(obj.numDofs, 1);
        obj.dofs_fixed_=zeros(length(mesh.data.NODE)*2,1);
        obj.initialdofs_ = zeros(obj.numDofs, 1);

        %thermomech
        meshelements=mesh.retrieveElementalSelection(inputReader.MeshEntityName);
            etype=mesh.data.ElementTypes{meshelements(1)};
            obj.dim = mesh.retrieveelementdimension(etype); 
            obj.numDofs_mech=obj.numNodes*obj.dim;
        obj.loadVector_mech = zeros(obj.numDofs_mech, 1);
        obj.dofs_fixed_mech_=zeros(length(mesh.data.NODE)*obj.dim,1);
        obj.initialdofs_mech_ = zeros(obj.numDofs_mech, 1);

        % Call boundaryConditions to perform any necessary setup
        obj.boundaryConditions(mesh,inputReader);
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function T3 = calculate_T3(~,nodes,dim)
            % Check if the input matrix is 3x(nodes per element)
            if size(nodes, 1) ~= 3
                error('Input matrix must be 3x(nodes per element)');
            end
            if dim>1
            % Extract node coordinates
            V1 = nodes(:, 1);
            V2 = nodes(:, 2);
            V3 = nodes(:, 3);
            elseif dim==1
            V1 = nodes(:, 1);
            V2 = nodes(:, 2);
                % find which axis is common for both
                V3=[];
                for i=1:3
                    axis=V1(i)-V2(i);
                    if axis==0
                        V3 =(V1+V2)/2;
                        Vnadd=zeros(3,1);
                        Vnadd(axis)=V1(axis);
                        V3=V3+Vnadd;
                    end
                end
                if V3==[]
                   fprintf('Warning: For 1D T3, no V3 was created.\n');
                end
            end
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
        function result = gaussIntegrationBC(obj,geometric_dimension, order, elementTag, bcvalue, mesh)
            if geometric_dimension < 1 || order < 1
                fprintf('Invalid dimension or order for Gauss integration.\n');
                result = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                return;
            end
            
            element_nodes = mesh.data.ELEMENTS{elementTag};
            coordinates = zeros(3,length(element_nodes));
            for i=1: length(element_nodes)
                element= element_nodes(i);
                node_coordinates=mesh.data.NODE{element};
                coordinates(:,i)=node_coordinates';
            end
            etype=mesh.data.ElementTypes{elementTag};
            dim = mesh.retrieveelementdimension(etype);
            T3=obj.calculate_T3(coordinates,dim);
            coordinates_tr = inv(T3) \ coordinates;
            result= zeros(length(element_nodes),1);
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
            
            if geometric_dimension == 1
                % 1D integration using a single loop.
                natcoords = zeros(1, 1);
                for i = 1:size(weights, 1)
                    natcoords(1) = gaussPoints(i);
                    % Explicitly use the element-wise multiplication .* for arrays
                    f = obj.CteSurfBC(natcoords, coordinates_tr_XY, bcvalue, elementTag,mesh) .* weights(i);
                    result = result + f;
                end
            elseif geometric_dimension == 2
                % 2D integration using a double loop.
                natcoords = zeros(2, 1);
                for i = 1:size(weights, 1)
                    for j = 1:size(weights, 1)
                        natcoords(1) = gaussPoints(i);
                        natcoords(2) = gaussPoints(j);
                        % Explicitly use the element-wise multiplication .* for arrays
                        f = obj.CteSurfBC(natcoords, coordinates_tr_XY, bcvalue, elementTag,mesh) .* (weights(i) * weights(j));
                        result = result + f;
                    end
                end
            elseif geometric_dimension == 3
                if order == 14
                    % Special case for 3D integration with order 14.
                    result = weights * weights' .* weights * weights' .* weights * weights' .* obj.CteSurfBC(gaussPoints, coordinates_tr_XY, bcvalue, elementTag,mesh);
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
                                f = obj.CteSurfBC(natcoords, coordinates_tr_XY, bcvalue, elementTag,mesh) .* (weights(i) * weights(j) * weights(k));
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
        function F_q = CteSurfBC(~,natcoords, coords, value, element,mesh)
            % Define variables
            etype=mesh.data.ElementTypes{element};
            dim = mesh.retrieveelementdimension(etype);
            % Extract natural coordinates
            if dim==2
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
            elseif dim==1
                xi = 0;
                if size(natcoords, 1) >= 2 && size(natcoords, 2) >= 1
                    xi = natcoords(1);  % Extracts the first element (a)
                else
                    % Handle the case when natcoords doesn't have the expected dimensions.
                    fprintf('Wrong natural coordinates dimension.\n');
                    % You may want to print an error message or take appropriate action.
                    return;
                end
            
                % Calculate shape functions and their derivatives
                [shapeFunctions,shapeFunctionDerivatives]=mesh.selectShapeFunctionsAndDerivatives(etype, xi, -1, -1);
            
                % Calculate Jacobian matrix JM
                JM = shapeFunctionDerivatives' * coords';

            end
            % Calculate the determinant of the Jacobian
            detJ = det(JM);
        
            % Calculate F_q (heat flow) using your equations
            F_q = detJ * (shapeFunctions * value);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function boundaryConditions(obj, mesh, inputReader)
            fprintf('#BCInit::boundaryConditions\n');
            
            % Get the boundary conditions from the InputReader
            boundaryConditions = inputReader.bctype;
            
            % Iterate through each boundary condition
            for k = 1:length(boundaryConditions)
                bc = boundaryConditions{k};
                boundaryName = bc;
                surfaceName = inputReader.bcloc{k};
                value = inputReader.bcval(k);
                
                fprintf('Boundary Name: %s, Surface Name: %s, Value: %f\n', boundaryName, surfaceName, value);
                
                if strcmp(boundaryName, 'Temperature')
                    if all(strcmp(mesh.data.NodalSelectionNames, surfaceName) == 0)
                        warning(append('BoundaryName ',surfaceName,' non-existent. Check your input data.'));
                        error('Terminating the code due to the condition.');
                    end
                    % Assuming 'getNodesForPhysicalGroup' returns a vector of unsigned long long integers
                    nodes=mesh.data.NSET{strcmp(mesh.data.NodalSelectionNames, surfaceName)};
                    
                    % Loop through the nodes and set the corresponding values in initialdofs_ (loadVector_)
                    for node = nodes
                        if (node*2) < length(obj.initialdofs_)
                            obj.initialdofs_(node*2-1) = value;
                            obj.dofs_fixed_(node*2-1)=1;
                        else
                            % Handle the case where the node index is out of bounds.
                            % This could be an error condition depending on your application.
                            fprintf('Node index out of bounds: %d\n', node);
                        end
                    end
                elseif strcmp(boundaryName, 'Voltage')
                    % Assuming 'getNodesForPhysicalGroup' returns a vector of unsigned long long integers
                    if all(strcmp(mesh.data.NodalSelectionNames, surfaceName) == 0)
                        warning(append('BoundaryName ',surfaceName,' non-existent. Check your input data.'));
                        error('Terminating the code due to the condition.');
                    end
                    nodes=mesh.data.NSET{strcmp(mesh.data.NodalSelectionNames, surfaceName)};
                    
                    % Loop through the nodes and set the corresponding values in initialdofs_ (loadVector_)
                    for node = nodes
                        if (node*2+1) < length(obj.initialdofs_)
                            obj.initialdofs_(node*2) = value;
                            obj.dofs_fixed_(node*2)=1;
                        else
                            % Handle the case where the node index is out of bounds.
                            % This could be an error condition depending on your application.
                            fprintf('Node index out of bounds: %d\n', node);
                        end
                    end
                elseif strcmp(boundaryName, 'FreeDofs_Voltage')
                    % Assuming 'getNodesForPhysicalGroup' returns a vector of unsigned long long integers
                    if all(strcmp(mesh.data.NodalSelectionNames, surfaceName) == 0)
                        warning(append('BoundaryName ',surfaceName,' non-existent. Check your input data.'));
                        error('Terminating the code due to the condition.');
                    end
                    nodes=mesh.data.NSET{strcmp(mesh.data.NodalSelectionNames, surfaceName)};
                    
                    % Loop through the nodes and set the corresponding values in initialdofs_ (loadVector_)
                    for node = nodes
                        if (node*2+1) < length(obj.initialdofs_)
                            obj.initialdofs_(node*2) = value;
                            obj.dofs_fixed_(node*2)=0;
                        else
                            % Handle the case where the node index is out of bounds.
                            % This could be an error condition depending on your application.
                            fprintf('Node index out of bounds: %d\n', node);
                        end
                    end
                elseif startsWith(boundaryName, 'Displacement')
                    if all(strcmp(mesh.data.NodalSelectionNames, surfaceName) == 0)
                        warning(append('BoundaryName ',surfaceName,' non-existent. Check your input data.'));
                        error('Terminating the code due to the condition.');
                    end
                    % Assuming 'getNodesForPhysicalGroup' returns a vector of unsigned long long integers
                    nodes=mesh.data.NSET{strcmp(mesh.data.NodalSelectionNames, surfaceName)};
                    
                    % Loop through the nodes and set the corresponding values in initialdofs_ (loadVector_)
                    direction=0;
                    if strcmp(boundaryName(end),'x')
                        direction=1;
                    elseif strcmp(boundaryName(end),'y')
                        direction=2;
                    elseif strcmp(boundaryName(end),'z')
                        direction=3;
                    end
                    for node = nodes
                        if (node*obj.dim) < length(obj.initialdofs_mech_)
                            obj.initialdofs_mech_((node-1)*obj.dim+direction) = value;
                            obj.dofs_fixed_mech_ ((node-1)*obj.dim+direction) = 1;
                        else
                            % Handle the case where the node index is out of bounds.
                            % This could be an error condition depending on your application.
                            fprintf('Node index out of bounds: %d\n', node);
                        end
                    end
                elseif strcmp(boundaryName, 'heat_n')
                    %% Heat integration

                    % Assuming 'getElementsAndNodeTagsForPhysicalGroup' returns a cell array of vectors
                    elementindexVector=mesh.retrieveElementalSelection(surfaceName);
                    %data.ELSET{strcmp(mesh.data.ElementSelectionNames, surfaceName)};
                    
                    % Debugging: Print the sizes of elements and element_nodes vectors
                    fprintf('Size of ''elements'' vector: %d\n', length(elementindexVector));
                    
                    fprintf('HEAT INTEGRATION.\n');
                    % Initialize temporary load vectors for each worker
                    %numWorkers = numel(gcp); % Get the number of workers
                    numWorkers=1;
                    tempLoadVectors = cell(1, numWorkers);
                    etype=mesh.data.ElementTypes{elementindexVector(1)};
                    dim = mesh.retrieveelementdimension(etype);
                    % Create a function handle for gaussIntegrationBC and CteSurfBC
                    gaussIntegrationBCFun = @(element) obj.gaussIntegrationBC(dim, 5, element, value, mesh);
                    %gaussIntegrationBCFun = @(element) obj.gaussIntegrationBC(2, 3, element, value);
                    
                    % Loop through elements in parallel
                    for workerIdx = 1:numWorkers
                        % Initialize temporary load vector for this worker
                        tempLoadVector = zeros(size(obj.loadVector_));
                        
                        % Calculate elements for this worker
                        workerElements = elementindexVector(workerIdx:numWorkers:end);
                        
                        % Loop through elements for this worker
                        for elementIdx = workerElements
                            % Calculate element_load_vector for the current element
                            element_load_vector = gaussIntegrationBCFun(elementIdx);
                            element_nodes = mesh.data.ELEMENTS{elementIdx};
                            obj.loadVector_(element_nodes*2-1) = obj.loadVector_(element_nodes*2-1) + element_load_vector;
                        end
                    end
                    
                    fprintf('HEAT INTEGRATION FINISHED.\n');
               
                elseif startsWith(boundaryName, 'Pressure')
                    %% Force integration
                    % Assuming 'getElementsAndNodeTagsForPhysicalGroup' returns a cell array of vectors
                    elementindexVector=mesh.retrieveElementalSelection(surfaceName);
                    %data.ELSET{strcmp(mesh.data.ElementSelectionNames, surfaceName)};
                    
                    % Debugging: Print the sizes of elements and element_nodes vectors
                    fprintf('Size of ''elements'' vector: %d\n', length(elementindexVector));
                    
                    fprintf('Force X INTEGRATION.\n');
                    % Initialize temporary load vectors for each worker
                    %numWorkers = numel(gcp); % Get the number of workers
                    numWorkers=1;
                    tempLoadVectors = cell(1, numWorkers);
                    etype=mesh.data.ElementTypes{elementindexVector(1)};
                    dim = mesh.retrieveelementdimension(etype);
                    % Create a function handle for gaussIntegrationBC and CteSurfBC
                    gaussIntegrationBCFun = @(element) obj.gaussIntegrationBC(dim, 3, element, value, mesh);
                    %gaussIntegrationBCFun = @(element) obj.gaussIntegrationBC(2, 3, element, value);
                    % Loop through the nodes and set the corresponding values in initialdofs_ (loadVector_)
                    direction=0;
                    if strcmp(boundaryName(end),'x')
                        direction=1;
                    elseif strcmp(boundaryName(end),'y')
                        direction=2;
                    elseif strcmp(boundaryName(end),'z')
                        direction=3;
                    end
                    % Loop through elements in parallel
                    for workerIdx = 1:numWorkers
                        % Initialize temporary load vector for this worker
                        tempLoadVector = zeros(size(obj.loadVector_));
                        
                        % Calculate elements for this worker
                        workerElements = elementindexVector(workerIdx:numWorkers:end);
                        
                        % Loop through elements for this worker
                        for elementIdx = workerElements
                            % Calculate element_load_vector for the current element
                            element_load_vector = gaussIntegrationBCFun(elementIdx);
                            obj.loadVector_mech((mesh.data.ELEMENTS{elementIdx}-1)*3+direction) = obj.loadVector_mech(mesh.data.ELEMENTS{elementIdx}*3+direction) + element_load_vector;
                        end
                    end
                    
                    fprintf('HEAT INTEGRATION FINISHED.\n');
                end
            end
            obj.dofs_free_=find(obj.dofs_fixed_ == 0);
            obj.dofs_free_mech=find(obj.dofs_fixed_mech_ == 0);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end
end

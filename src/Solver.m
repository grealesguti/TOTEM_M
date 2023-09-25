classdef Solver < handle
    properties
        inputReader_
        mesh_
        bcinit_
        TO
        tolerance
        max_iterations
        soldofs
        loadVector
    end
    
    properties (Hidden)
        elements
        utils
        meshFileName
        freedofidxs
    end
    
    methods
        function obj = Solver(inputReader, mesh, bcinit)
                obj.inputReader_ = inputReader;
                obj.mesh_ = mesh;
                obj.bcinit_ = bcinit;
                obj.elements = Elements(); % Initialize Elements here
        
                % Get the number of nodes from the mesh
                numNodes = obj.mesh.getNumAllNodes();
        
                % Initialize the loadVector_ member with a size double the number of nodes
                obj.loadVector = obj.bcinit.getloadVector();
                obj.soldofs = zeros(1, 2 * numNodes);
                obj.soldofs(1:numel(obj.bcinit.getAllInitialDof())) = obj.bcinit.getAllInitialDof();
                obj.freedofidxs = obj.mesh.GetFreedofsIdx();
        
                % Initialize the thermoelectricityintegrationFunction_ using a function handle
                obj.thermoelectricityintegrationFunction = @(natcoords, coords, dofs, elementTag) obj.thermoelectricityintegration(natcoords, coords, dofs, elementTag);
        
                obj.max_iterations=10;
                obj.tolerance=1e-6;
                fprintf('### SOLVER Initialized.\n');
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [KJs,Ra]=Assembly(obj)

            mesh_elements_name = obj.reader.MeshEntityName;
            mesh_elements = obj.mesh.getElementsFromName (mesh_elements_name);
            total_number_of_elements = length(mesh_elements);
            total_number_of_nodes = length(obj.mesh.NODE);
            total_number_of_dofs = total_number_of_nodes*2;
            dofs_per_element = 20*2;
            % FIXME:
                % define nnv and nele and improve their naming
            Ra=zeros(total_number_of_dofs,1);%KJa=sparse(ndof,ndof);
            KJvnz=zeros(dofs_per_element,total_number_of_elements);
            KJvc=zeros(dofs_per_element,total_number_of_elements);
            KJvr=zeros(dofs_per_element,total_number_of_elements);
            for  i=1:total_number_of_elements

                % recover each element tag
                element_Tag=elements[i];

                % initialization to zeros for each element
                Rs=zeros(ndof,1);%KJs=sparse(ndof,ndof);

                % FIXME
                % modify the function to run with the gauss
                % previoously defined
                Ke,Re = gaussIntegrationK(3, 5, elementTag, mesh, func);
                % assembly in global residual and jacobian matrix in sparse format
                Rs(doforder,1)=Re(:,1);      
                Ra=Ra+Rs;
                
                KJvnz(:,i)=reshape(Ke,nnv,1);
                KJvc(:,i)=reshape(repmat(doforder,20*2,1),nnv,1);
                KJvr(:,i)=reshape(repmat(doforder,20*2,1)',nnv,1);

            end
            KJvnzr=reshape(KJvnz,total_number_of_elements*dofs_per_element,1);
            KJvcr=reshape(KJvc,total_number_of_elements*dofs_per_element,1);
            KJvrr=reshape(KJvr,total_number_of_elements*dofs_per_element,1); % is this needed??
            KJs=sparse(KJvrr,KJvcr,KJvnzr,total_number_of_dofs,total_number_of_dofs);
            Ra=Ra-obj.load_vector;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [K,R] = gaussIntegrationK(dimension, order, elementTag, mesh, dofs, func)
            if dimension < 1 || order < 1
                fprintf('Invalid dimension or order for Gauss integration.\n');
                K = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                R = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                return;
            end
           
            [weights, gaussPoints] = getGaussWeightsAndPoints(order);
            
            if isempty(weights) || isempty(gaussPoints)
                fprintf('Invalid order for Gauss integration.\n');
                K = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                R = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                return;
            end
            
            if size(weights, 1) ~= size(gaussPoints, 1)
                fprintf('Weights and Gauss points have mismatched dimensions.\n');
                K = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                R = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                return;
            end

            number_of_nodes = length(mesh.data.ELEMENTS{elementTag});
            dof_per_node = 2;
            K=zeros(number_of_nodes*dof_per_node,number_of_nodes*dof_per_node);
            R=zeros(number_of_nodes*dof_per_node,1);
            
            if dimension == 1
                % 1D integration using a single loop.
                natcoords = zeros(1, 1);
                for i = 1:size(weights, 1)
                    natcoords(1) = gaussPoints(i);
                    % Explicitly use the element-wise multiplication .* for arrays
                    Ke,Re = func(natcoords, coordinates_tr_XY, bcvalue, elementTag, mesh) .* weights(i);
                    K = K + Ke;
                    R = R + Re;
                end
            elseif dimension == 2
                % 2D integration using a double loop.
                natcoords = zeros(2, 1);
                for i = 1:size(weights, 1)
                    for j = 1:size(weights, 1)
                        natcoords(1) = gaussPoints(i);
                        natcoords(2) = gaussPoints(j);
                        % Explicitly use the element-wise multiplication .* for arrays
                        Ke,Re = func(natcoords, coordinates_tr_XY, bcvalue, elementTag, mesh) .* (weights(i) * weights(j));
                        K = K + Ke;
                        R = R + Re;
                    end
                end
            elseif dimension == 3
                if order == 14
                    % Special case for 3D integration with order 14.
                    Ke,Re = weights * weights' .* weights * weights' .* weights * weights' .* func(gaussPoints, coordinates_tr_XY, bcvalue, elementTag, mesh);
                    K = K + Ke;
                    R = R + Re;
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
                                Ke,Re = func(natcoords, coordinates_tr_XY, bcvalue, elementTag, mesh) .* (weights(i) * weights(j) * weights(k));
                                K = K + Ke;
                                R = R + Re;
                            end
                        end
                    end
                end
            else
                fprintf('Invalid dimension for Gauss integration.\n');
                K = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                R = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [KJ, R] = thermoelectricityintegrationFunction(elementTag,mesh,xi,eta,zeta,etype,reader)
            %% Thermoelectricity Simulation
            % This function calculates thermoelectric properties using finite element analysis.
            % Inputs:
            %
            % Outputs:
            %   - KJ: Jacobian matrix
            %   - R: Residues
            
            % Extract coordinates for the current order
            element_nodes = mesh.data.ELEMENTS{elementTag};
            element_coordinates=zeros(3,length(element_nodes));
            Tee=zeros(1,element_nodes);
            Vee=zeros(1,element_nodes);
            for i = 1: length(element_nodes)
                element_coordinates(:,i)=mesh.data.NODE{element_nodes(i)};
                Tee(i)=dofs(element_nodes*2);
                Vee(i)=dofs(element_nodes*2+1);
            end
            element_material_index=mesh.elements_material(elementTag);
            
           
                % FIXME: Calculate shape functions
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, xi, eta, zeta);

                JM = dShape * element_coordinates;
                %Jacinv = inv(JM);
                DN = Jacinv \ dShape; % FIXME and check it is the same!
            
                % FIXME, calculate from all dofs input
                Th = N' * Tee;
            
                % FIXME: Calculate material properties
                De = reader.getmaterialproperty(element_material_index,'electricalconductivity');
                Dk = reader.getmaterialproperty(element_material_index,'thermalconductivity');
                Da = reader.getmaterialproperty(element_material_index,'Seebeck');
                %if TO.isTO
                %    De=De;
                %    Dk=Dk;
                %    Da=Da;
                %end

                Dde=0;Dda=0;Ddk=0;
            
                detJ = det(JM);
            
                % Calculate current density and heat flux
                je = -De * DN * Vee - Da * De * DN * Tee;
                qe = Da * (N' * Tee) * je - Dk * DN * Tee;
            
                % Calculate derivatives
                djdt = -Da * De * DN - Dda * De * DN * Tee * N' - Dde * (DN * Ve + Da * DN * Tee) * N';
                djdv = -De * DN;
                dqdt = Da * Th * djdt + Da * je * N' - Dk * DN + Dda * Th * je * N' - Ddk * DN * Tee * N';
                dqdv = -Da * De * Th * DN;
            
                % Update residues
                RT = detJ*(-(DN' * qe) + (N * je') * (DN * Vee));
                RV = detJ*(-DN' * je);
            
                % Update stiffness matrices
                K11 = detJ*(DN' * dqdt - N * (djdt' * DN * Vee)');
                K12 = detJ*(DN' * dqdv - N * (djdv' * DN * Vee)' - N * (je' * DN));
                K21 = detJ*(DN' * djdt);
                K22 = detJ*(DN' * djdv);
            
                % Construct the elemental Jacobian matrix and Residues
                KJ = [K11, K12; K21, K22];
                R = [RT(:, 1); RV(:, 1)];
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function solution = solveSparseSystem(obj, K, R)
            solution=K(obj.mesh_.freedofs,obj.mesh_.freedofs)\R(obj.mesh_.freedofs);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function residual = runNewtonRaphson(obj)
            % Initialize some parameters and initial guess
            
            for iter = 1:obj.max_iterations
                % Assembly the system matrix
                [Ksparse,Rsparse] = Assembly();
                
                if (residual < obj.tolerance && iter>1)
                    % Converged, return both the solution and the final residual
                    fprintf('### NR. CONVERGED\n');
                    return;
                end
                
                % Solve the system using the tangential matrix and residual: KT*dU=R->dU
                delta_degreesoffreedom = solveSparseSystem(Ksparse(obj.freedofidxs,obj.freedofidxs),Rsparse(obj.freedofidxs));
                
                % Update solution
                for i = 1:length(obj.freedofidxs)
                    obj.soldofs(obj.freedofidxs(i)) = obj.soldofs(obj.freedofidxs(i)) + delta_degreesoffreedom(i);
                end
                
                % Calculate the residual (error)
                residual = norm(Rsparse); % Calculate the norm
                fprintf('### NR. Iteration %d residual %f\n', iter, residual);
              
            end
        
            % If we reach here, the Newton-Raphson method did not converge
            error('Newton-Raphson did not converge: iteration limit');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function residual = runModifiedNewtonRaphsonSolver(obj)
            % Initialize some parameters and initial guess
                % Assembly the system matrix
            [Ksparse,Rsparse] = Assembly();
            for iter = 1:obj.max_iterations
                % Assembly the system matrix
                
                if (residual < obj.tolerance && iter>1)
                    % Converged, return both the solution and the final residual
                    fprintf('### NR. CONVERGED\n');
                    return;
                end
                
                % Solve the system using the tangential matrix and residual: KT*dU=R->dU
                delta_degreesoffreedom = solveSparseSystem(Ksparse(obj.freedofidxs,obj.freedofidxs),Rsparse(obj.freedofidxs));
                
                % Update solution
                for i = 1:length(obj.freedofidxs)
                    obj.soldofs(obj.freedofidxs(i)) = obj.soldofs(obj.freedofidxs(i)) + delta_degreesoffreedom(i);
                end
                
                % Calculate the residual (error)
                residual = norm(Rsparse); % Calculate the norm
                fprintf('### NR. Iteration %d residual %f\n', iter, residual);
              
            end
        
            % If we reach here, the Newton-Raphson method did not converge
            error('Newton-Raphson did not converge: iteration limit');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function result = runNewtonRaphsonWithUniformIncrements(obj, totalLoadVector, numUniformIncrements)
            % Implement your code here
        end
    end
    
    methods (Access = private)
        function result = thermoelectricityintegration(obj, natcoords, coords, dofs, elementTag)
            % Implement your code here
        end
    end
end

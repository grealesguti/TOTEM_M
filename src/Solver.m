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
        numNodes
        thermoelectricityintegrationFunctionFun
        KT
        Residual
    end
    
    properties (Hidden)
        utils
        meshFileName
        freedofidxs
    end
    
    methods
        function obj = Solver( mesh, bcinit)
        
                % Get the number of nodes from the mesh
                obj.numNodes = length(mesh.data.NODE);
        
                % Initialize the loadVector_ member with a size double the number of nodes
                obj.loadVector = bcinit.loadVector_;
                obj.soldofs = bcinit.initialdofs_;
                
                %[obj.KT,obj.Residual]=obj.Assembly(inputReader,mesh,bcinit);

                %obj.SolveLinearSystemInParallel(bcinit);
                
                obj.max_iterations=10;
                obj.tolerance=1e-6;
                fprintf('### SOLVER Initialized.\n');
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [KJs, Ra] = Assembly(obj, reader, mesh)
            mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);
            total_number_of_elements = length(mesh_elements);
            total_number_of_nodes = length(mesh.data.NODE);
            dofs_per_element = 2;
            total_number_of_dofs = total_number_of_nodes * dofs_per_element;
            
            [node_el, etype_element] = mesh.retrievemeshtype(reader);
            dofs_per_element = (node_el * dofs_per_element);
            
            Ra = zeros(total_number_of_dofs, 1);
            KJvnz = zeros(dofs_per_element^2, total_number_of_elements);
            KJvc = zeros(dofs_per_element^2, total_number_of_elements);
            KJvr = zeros(dofs_per_element^2, total_number_of_elements);
            
            initialdofs = obj.soldofs;
            
            % Initialize the parallel pool with the desired number of workers
            %numWorkers = 4; % Adjust the number of workers as needed
            %pool = parpool(numWorkers);
            
            for i = 1:total_number_of_elements
                % Create a separate variable for each parallel iteration
                Rs = zeros(dofs_per_element, 1);

                % Recover each element tag
                elementTag = mesh_elements(i);
                
                % Compute element stiffness matrix and residual
                [Ke, Re, element_dof_indexes] = obj.gaussIntegrationK(3, 14, elementTag, mesh, initialdofs, reader, etype_element);
                
                % Assembly in global residual
                Rs(element_dof_indexes, 1) = Re(:, 1);
                
                % Accumulate residuals and stiffness matrix contributions
                Ra = Ra + Rs;
                KJvnz(:, i) = reshape(Ke, dofs_per_element^2, 1);
                repmat_idx=repmat(element_dof_indexes,1, dofs_per_element);
                KJvr(:, i) = reshape(repmat_idx, dofs_per_element^2, 1);
                KJvc(:, i) = reshape(repmat_idx', dofs_per_element^2, 1);               
            end
            
            % Clean up parallel pool
            %delete(pool);
            
            % Reshape and construct stiffness matrix
            KJvnzr = reshape(KJvnz, total_number_of_elements * dofs_per_element^2, 1);
            KJvcr = reshape(KJvc, total_number_of_elements * dofs_per_element^2, 1);
            KJvrr = reshape(KJvr, total_number_of_elements * dofs_per_element^2, 1);
            
            KJs = sparse(KJvrr, KJvcr, KJvnzr, total_number_of_dofs, total_number_of_dofs);
            
            % Subtract external load
            Ra = Ra - obj.loadVector;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [K,R,element_dof_indexes] = gaussIntegrationK(obj,dimension, order, elementTag, mesh, initialdofs,reader,etype)
            
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% Ininitalization of elemental integration variables %%%%%
            element_nodes = mesh.data.ELEMENTS{elementTag};
            number_of_nodes = length(element_nodes);
            element_coordinates=zeros(3,number_of_nodes);
            Tee=zeros(1,number_of_nodes);
            Vee=zeros(1,number_of_nodes);
            element_dof_indexes=zeros(number_of_nodes*2,1);
            for i=1:number_of_nodes
                element_dof_indexes(i)=element_nodes(i)*2-1;
                element_dof_indexes(number_of_nodes+i)=element_nodes(i)*2;
                element_coordinates(:,i)=mesh.data.NODE{element_nodes(i)};
                Tee(i)=initialdofs(element_nodes(i)*2-1);
                Vee(i)=initialdofs(element_nodes(i)*2);
            end

            element_material_index=mesh.elements_material(elementTag);

            dof_per_node = 2;
            K=zeros(number_of_nodes*dof_per_node,number_of_nodes*dof_per_node);
            R=zeros(number_of_nodes*dof_per_node,1);

            integrationFunction = @(natcoords) obj.thermoelectricityintegrationFunction(natcoords, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype, mesh.elements_density(elementTag));      
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
             if dimension == 1
                % 1D integration using a single loop.
                natcoords = zeros(1, 1);
                for i = 1:size(weights, 1)
                    natcoords(1) = gaussPoints(i);
                    % Explicitly use the element-wise multiplication .* for arrays
                    [Ke,Re] = integrationFunction(natcoords) ;
                    Ke=Ke.* weights(i);
                    Re=Re.* weights(i);
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
                        Ke=Ke.* (weights(i) * weights(j));
                        Re=Re.* (weights(i) * weights(j));
                        K = K + Ke;
                        R = R + Re;
                    end
                end
            elseif dimension == 3
                if order == 14
                    % Special case for 3D integration with order 14.
                    %[Ke,Re] = integrationFunction(natcoords) ;
                    %[Ke,Re] = weights * weights' .* weights * weights' .* weights * weights' .* integrationFunction(gaussPoints);
                    for i=1:14
                        natcoords(1) = gaussPoints(1,i);
                        natcoords(2) = gaussPoints(2,i);
                        natcoords(3) = gaussPoints(3,i);
                        [Ke,Re] = integrationFunction(natcoords) ;
                        Ke=Ke .* (weights(i));
                        Re=Re .* (weights(i));
                        K = K + Ke;
                        R = R + Re;
                    end

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
                                [Ke,Re] = integrationFunction(natcoords) ;
                                Ke=Ke .* (weights(i) * weights(j) * weights(k));
                                Re=Re .* (weights(i) * weights(j) * weights(k));
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
        function [KJ, R] = thermoelectricityintegrationFunction(~,natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)
            %% Thermoelectricity Simulation
            % This function calculates thermoelectric properties using finite element analysis.
            % Inputs:
            %
            % Outputs:
            %   - KJ: Jacobian matrix
            %   - R: Residues
            
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));

                JM = dShape' * element_coordinates';
                %Jacinv = inv(JM);
                DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
                %DN1 = JM \ dShape'; % FIXME and check it is the same!

                % FIXME, calculate from all dofs input
                Th = N * Tee';
            
                % FIXME: Calculate material properties
                Dep = reader.getmaterialproperty(element_material_index,'ElectricalConductivity');
                Dkp = reader.getmaterialproperty(element_material_index,'ThermalConductivity');
                Dap = reader.getmaterialproperty(element_material_index,'Seebeck');
                
                [De,Dde]=CalculateMaterialProperties(Dep,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
                [Da,Dda]=CalculateMaterialProperties(Dap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
                [Dk,Ddk]=CalculateMaterialProperties(Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));

                Vee=Vee';
                Tee=Tee';
                detJ = det(JM);
            
                % Calculate current density and heat flux
                je = -De * DN * Vee - Da * De * DN * Tee;
                qe = Da * (N * Tee) * je - Dk * DN * Tee;
            
                % Calculate derivatives
                djdt = -Da * De * DN - Dda * De * DN * Tee * N - Dde * (DN * Vee + Da * DN * Tee) * N;
                djdv = -De * DN;
                dqdt = Da * Th * djdt + Da * je * N - Dk * DN + Dda * Th * je * N - Ddk * DN * Tee * N;
                dqdv = -Da * De * Th * DN;
            
                % Update residues
                RT = detJ*(-(DN' * qe) + (N' * je') * (DN * Vee));
                RV = detJ*(-DN' * je);
            
                % Update stiffness matrices
                K11 = detJ*(DN' * dqdt - N' * (djdt' * DN * Vee)');
                K12 = detJ*(DN' * dqdv - N' * (djdv' * DN * Vee)' - N' * (je' * DN));
                K21 = detJ*(DN' * djdt);
                K22 = detJ*(DN' * djdv);
                %K12=detJ*(DN0*dqdv-N*(djdv'*DN0'*Vee)'-N*(je'*DN0'));

                % Construct the elemental Jacobian matrix and Residues
                KJ = [K11, K12; K21, K22];
                R = [RT(:, 1); RV(:, 1)];
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function SolveLinearSystemInParallel(obj, bcinit)
            % Extract the necessary variables
            dofs_free = bcinit.dofs_free_;
        
            KT_Distr = distributed(obj.KT(dofs_free, dofs_free)); 
            R_Distr=distributed(obj.Residual(dofs_free));
            dU =( KT_Distr ) \ ( R_Distr );  % Calculation of step
            obj.soldofs(dofs_free)=obj.soldofs(dofs_free)+dU;

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function residual = runNewtonRaphson(obj,reader, mesh,bcinit)
            % Initialize some parameters and initial guess
            residual = obj.tolerance*1000;
            for iter = 1:obj.max_iterations
                % Assembly the system matrix
                [obj.KT,obj.Residual] = obj.Assembly(reader, mesh);
                
                if (residual < obj.tolerance && iter>1)
                    % Converged, return both the solution and the final residual
                    fprintf('### NR. CONVERGED\n');
                    return;
                end
                
                % Solve the system using the tangential matrix and residual: KT*dU=R->dU
               obj.SolveLinearSystemInParallel(bcinit)
             
                % Calculate the residual (error)
                dofs_free = bcinit.dofs_free_;
                residual = norm(obj.Residual(dofs_free)); % Calculate the norm
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

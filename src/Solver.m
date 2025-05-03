classdef Solver < handle
    properties
        inputReader_
        mesh_
        bcinit_
        TO
        tolerance
        max_iterations
        soldofs
        soldofs_mech
        loadVector
        loadVector_copy
        numNodes
        thermoelectricityintegrationFunctionFun
        KT
        Residual
        Residual_mech
        loadVector_mech
        KStiff
        KUT
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
                obj.loadVector_copy=obj.loadVector;
                obj.loadVector_mech = bcinit.loadVector_mech;
                obj.soldofs = bcinit.initialdofs_;
                obj.soldofs_mech = bcinit.initialdofs_mech_;
                %[obj.KT,obj.Residual]=obj.Assembly(inputReader,mesh,bcinit);

                %obj.SolveLinearSystemInParallel(bcinit);
                obj.max_iterations=20;
                obj.tolerance=1e-8;
                fprintf('### SOLVER Initialized.\n');
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [KJs, Ral, Ra] = Assembly_Thermoelectricity(obj, reader, mesh)
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
            etype=mesh.data.ElementTypes{mesh_elements(1)};
            dim = mesh.retrieveelementdimension(etype);            
            for i = 1:total_number_of_elements
                % Create a separate variable for each parallel iteration
                Rs = zeros(total_number_of_dofs, 1);

                % Recover each element tag
                elementTag = mesh_elements(i);

                % Compute element stiffness matrix and residual
                [Ke, Re, element_dof_indexes] = obj.gaussIntegrationK(dim, reader.GI_order, elementTag, mesh, initialdofs, reader, etype_element,'Thermoelectricity');
                
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
            Ral = Ra - obj.loadVector;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [KStiffness, KTemp_Load,KTheta] = Assembly_DecoupledThermoMech(obj, reader, mesh)
            mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);
            total_number_of_elements = length(mesh_elements);
            total_number_of_nodes = length(mesh.data.NODE);
            etype=mesh.data.ElementTypes{mesh_elements(1)};
            dim = mesh.retrieveelementdimension(etype); 
            dofs_per_node = dim;
            total_number_of_dofs = total_number_of_nodes * dofs_per_node;
            
            [node_el, etype_element] = mesh.retrievemeshtype(reader);
            dofs_per_element = (node_el * dofs_per_node);
            
            KJvnz = zeros(dofs_per_element*node_el, total_number_of_elements);
            KJvc = zeros(dofs_per_element*node_el, total_number_of_elements);
            KJvr = zeros(dofs_per_element*node_el, total_number_of_elements);
            KJvnzT = zeros(dofs_per_element*node_el, total_number_of_elements);
            KJvcT = zeros(dofs_per_element*node_el, total_number_of_elements);
            KJvrT = zeros(dofs_per_element*node_el, total_number_of_elements);            
            initialdofs_TV = obj.soldofs;

            parfor i = 1:total_number_of_elements

                % Recover each element tag
                elementTag = mesh_elements(i);
                
                % Compute element stiffness matrix and residual
                [Ke, Ke1, element_dof_indexes_M,element_dof_indexes_T] = obj.gaussIntegrationK(dim,  reader.GI_order, elementTag, mesh, initialdofs_TV, reader, etype_element,'DecoupledThermoMechanical_Load');
                
                % Accumulate residuals and stiffness matrix contributions
                %KJvnz(:, i) = reshape(Ke, dofs_per_element*node_el, 1);
                %repmat_idx_M=repmat(element_dof_indexes_M,1, node_el);
                %repmat_idx_T=repmat(element_dof_indexes_T',dofs_per_element, 1);
                %KJvr(:, i) = reshape(repmat_idx_M, dofs_per_element*node_el, 1);
                %KJvc(:, i) = reshape(repmat_idx_T', dofs_per_element*node_el, 1);   

                nnv=dim*node_el*node_el;
                KJvnz(:,i)=reshape(Ke,nnv,1);
                KJvc(:,i)=reshape(repmat(element_dof_indexes_T,node_el*dim,1),nnv,1);
                KJvr(:,i)=reshape(repmat(element_dof_indexes_M,node_el,1)',nnv,1);

                KJvnzT(:,i)=reshape(Ke1,nnv,1);
                KJvcT(:,i)=reshape(repmat(element_dof_indexes_T,node_el*dim,1),nnv,1);
                KJvrT(:,i)=reshape(repmat(element_dof_indexes_M,node_el,1)',nnv,1);
            end
            KJvnzr=reshape(KJvnz,total_number_of_elements*dofs_per_element*node_el,1);
            KJvcr=reshape(KJvc,total_number_of_elements*dofs_per_element*node_el,1);
            KJvrr=reshape(KJvr,total_number_of_elements*dofs_per_element*node_el,1);
            KTemp_Load=sparse(KJvrr,KJvcr,KJvnzr,total_number_of_dofs,total_number_of_nodes);

            KJvnzrT=reshape(KJvnzT,total_number_of_elements*dofs_per_element*node_el,1);
            KJvcrT=reshape(KJvcT,total_number_of_elements*dofs_per_element*node_el,1);
            KJvrrT=reshape(KJvrT,total_number_of_elements*dofs_per_element*node_el,1);
            KTheta=sparse(KJvrrT,KJvcrT,KJvnzrT,total_number_of_dofs,total_number_of_nodes);

           % [F,KTemp_Load]=linkMT_TOIso_Test(reader,mesh,obj)
            
            % Initialize the parallel pool with the desired number of workers
            %numWorkers = 4; % Adjust the number of workers as needed
            %pool = parpool(numWorkers);
            KJvnz = zeros(dofs_per_element^2, total_number_of_elements);
            KJvc = zeros(dofs_per_element^2, total_number_of_elements);
            KJvr = zeros(dofs_per_element^2, total_number_of_elements);
            parfor i = 1:total_number_of_elements
                % Recover each element tag
                elementTag = mesh_elements(i);
                
                % Compute element stiffness matrix and residual
                [Ke, empt, element_dof_indexes] = obj.gaussIntegrationK(dim,  reader.GI_order, elementTag, mesh, initialdofs_TV, reader, etype_element,'Mechanical');
                
                % Accumulate residuals and stiffness matrix contributions
                KJvnz(:, i) = reshape(Ke, dofs_per_element^2, 1);
                repmat_idx=repmat(element_dof_indexes,1, dofs_per_element);
                KJvr(:, i) = reshape(repmat_idx, 1,  dofs_per_element^2);
                KJvc(:, i) = reshape(repmat_idx', dofs_per_element^2, 1);               
            end
            
            % Clean up parallel pool
            %delete(pool);
            
            % Reshape and construct stiffness matrix
            KJvnzr = reshape(KJvnz, total_number_of_elements * dofs_per_element^2, 1);
            KJvcr = reshape(KJvc, total_number_of_elements * dofs_per_element^2, 1);
            KJvrr = reshape(KJvr, total_number_of_elements * dofs_per_element^2, 1);
            
            KStiffness = sparse(KJvrr, KJvcr, KJvnzr, total_number_of_dofs, total_number_of_dofs);
                        %Temperature_solution = obj.soldofs(1:2:end);
                        %obj.KUT=KThermalLoad;
                        %obj.loadVector_mech=obj.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
            %Ra=KStiffness*obj.soldofs_mech - KTheta*(obj.soldofs(1:2:end)-str2double(reader.T0));
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function [K,R,element_dof_indexes,empt] = gaussIntegrationK(obj,dimension, order, elementTag, mesh, initialdofs,reader,etype, physics)
            empt=[];
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
            if strcmp(physics,'Thermoelectricity')
                element_nodes = mesh.data.ELEMENTS{elementTag};
                number_of_nodes = length(element_nodes);
                element_coordinates=zeros(3,number_of_nodes);
                dof_per_node=2;
                Tee=zeros(1,number_of_nodes);
                Vee=zeros(1,number_of_nodes);
                element_dof_indexes=zeros(number_of_nodes*dof_per_node,1);
                for i=1:number_of_nodes
                    element_dof_indexes(i)=element_nodes(i)*dof_per_node-1;
                    element_dof_indexes(number_of_nodes+i)=element_nodes(i)*dof_per_node;
                    element_coordinates(:,i)=mesh.data.NODE{element_nodes(i)};
                    Tee(i)=initialdofs(element_nodes(i)*dof_per_node-1);
                    Vee(i)=initialdofs(element_nodes(i)*dof_per_node);
                end
    
                element_material_index=mesh.elements_material(elementTag);
                    
% Validate material index
if element_material_index < 1 || element_material_index > numel(reader.MaterialProperties)
    error('Invalid material index (%d) for element %d. Material index must be a positive integer within bounds (1 to %d).', ...
          element_material_index, elementTag, numel(reader.MaterialProperties));
end

                K=zeros(number_of_nodes*dof_per_node,number_of_nodes*dof_per_node);
                R=zeros(number_of_nodes*dof_per_node,1);
    
                integrationFunction = @(natcoords) obj.thermoelectricityintegrationFunction(natcoords, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype, mesh.elements_density(elementTag));      
            elseif strcmp(physics,'DecoupledThermoMechanical_Load')
                element_nodes = mesh.data.ELEMENTS{elementTag};
                number_of_nodes = length(element_nodes);
                element_coordinates=zeros(3,number_of_nodes);
                [dof_per_node] = mesh.retrieveelementdimension(etype);

                Tee=zeros(1,number_of_nodes);
                element_dof_indexes=zeros(number_of_nodes*dof_per_node,1);
                element_dof_indexes_T=zeros(number_of_nodes,1);
                for i=1:number_of_nodes
                    for j=1:dof_per_node
                        element_dof_indexes((i-1)*dof_per_node+j)=(element_nodes(i)-1)*dof_per_node+j;
                    end
                    element_coordinates(:,i)=mesh.data.NODE{element_nodes(i)};
                    Tee(i)=initialdofs(element_nodes(i)*2-1);
                end
    
                element_dof_indexes_T=zeros(1,number_of_nodes);
                element_dof_indexes=zeros(1,dof_per_node*number_of_nodes);
                for i=1:number_of_nodes
                    for j=1:dof_per_node
                        element_dof_indexes((i-1)*dof_per_node+j)=(element_nodes(i)-1)*dof_per_node+j;
                    end
                    element_dof_indexes_T(i)=element_nodes(i);
                end

                element_material_index=mesh.elements_material(elementTag);
    
                K=zeros(number_of_nodes*dof_per_node,number_of_nodes);
                R=zeros(number_of_nodes*dof_per_node,number_of_nodes);
    
                integrationFunction = @(natcoords) obj.thermalloadintegrationFunction(natcoords, element_coordinates, Tee, [], element_material_index, reader, mesh, etype, mesh.elements_density(elementTag));   
            elseif strcmp(physics,'Mechanical')
                element_nodes = mesh.data.ELEMENTS{elementTag};
                number_of_nodes = length(element_nodes);
                element_coordinates=zeros(3,number_of_nodes);
                [dof_per_node] = mesh.retrieveelementdimension(etype);
                Tee=zeros(1,number_of_nodes);
                Uee=zeros(1,number_of_nodes*dof_per_node);
                element_dof_indexes=zeros(number_of_nodes*dof_per_node,1);
                for i=1:number_of_nodes
                    for j=1:dof_per_node
                        element_dof_indexes((i-1)*dof_per_node+j)=(element_nodes(i)-1)*dof_per_node+j;
                    end
                    element_coordinates(:,i)=mesh.data.NODE{element_nodes(i)};
                    Tee(i)=initialdofs(element_nodes(i)*2-1);
                end
    
                element_material_index=mesh.elements_material(elementTag);
    
                K=zeros(number_of_nodes*dof_per_node,number_of_nodes*dof_per_node);
                R=[];
    
                integrationFunction = @(natcoords) obj.stiffnessintegrationFunction(natcoords, element_coordinates, Tee, [], element_material_index, reader, mesh, etype, mesh.elements_density(elementTag));   
           
            end
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
                        [Ke,Re] = integrationFunction(natcoords) ;
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
             if strcmp(physics,'DecoupledThermoMechanical_Load')
                 empt=element_dof_indexes_T;
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
            [dim] = mesh.retrieveelementdimension(etype);
            if dim==2
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), 0);
                element_coordinates=element_coordinates(1:2,:);
                N=N';
            else
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));
            end
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
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_ElectricalConductivity'),reader.getmaterialproperty(element_material_index,'Tmax_ElectricalConductivity')];
                [De,Dde]=CalculateMaterialProperties(1e-6,Dep,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_Seebeck'),reader.getmaterialproperty(element_material_index,'Tmax_Seebeck')];
                [Da,Dda]=CalculateMaterialProperties(1e-6,Dap,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_ThermalConductivity'),reader.getmaterialproperty(element_material_index,'Tmax_ThermalConductivity')];
                [Dk,Ddk]=CalculateMaterialProperties(reader.kmin,Dkp,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));

                if (reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity')>1)
                    found=reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity');
                end

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

                % Check if there are any NaN values
                %if any(isnan(KJ))
                %    % Raise a warning
                %    warning('Array KJ contains NaN values.');
                %elseif any(isnan(R))
                %    % Raise a warning
                %    warning('Array R contains NaN values.');                    
                %end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F_T, F_TR] = thermalloadintegrationFunction(~,natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)
            %% Thermoelectricity Simulation
            % This function calculates thermoelectric properties using finite element analysis.
            % Inputs:
            %
            % Outputs:
            %   - KJ: Jacobian matrix
            %   - R: Residues
            [dim] = mesh.retrieveelementdimension(etype);
            if dim==2
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), 0);
                element_coordinates=element_coordinates(1:2,:);
                N=N';
            else
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));
            end
                Number_of_Nodes = length(element_coordinates(1,:));
                JM = dShape' * element_coordinates';
                %Jacinv = inv(JM);
                DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
                %DN1 = JM \ dShape'; % FIXME and check it is the same!
                
                % FIXME, calculate from all dofs input
                Th = N * Tee';
            
                % FIXME: Calculate material properties
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_YoungModulus'),reader.getmaterialproperty(element_material_index,'Tmax_YoungModulus')];
                Dalpha_x = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_x');
                [Dax,Daxdt]=CalculateMaterialProperties(1e-6,Dalpha_x,Tmat,1,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
                Dalpha_y = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_y');
                [Day,Daydt]=CalculateMaterialProperties(1e-6,Dalpha_y,Tmat,1,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
                Dalpha_z = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_z');
                [Daz,Dazdt]=CalculateMaterialProperties(1e-6,Dalpha_z,Tmat,1,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
                DEp = reader.getmaterialproperty(element_material_index,'YoungModulus');
                nu = reader.getmaterialproperty(element_material_index,'PoissonRatio');
                [DE,DdE]=CalculateMaterialProperties(1e-6,DEp,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
                %[Dalpha,Ddalpha]=CalculateMaterialProperties(1e-6,Dalphap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalExpansionCoefficient'));
                
                
                detJ = det(JM);
            if dim==2
               alphav=zeros(3,1);
               alphav(1:2,1)=[Dax,Day];
              alphavdt=zeros(3,1);
                alphavdt(1:2,1)=[Daxdt,Daydt];
                C = DE / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2]; % plane stress
                % C = DE / (1 + nu) / (1 - 2 * nu) * [1 - nu, nu, 0; nu, 1
                % - nu, 0; 0, 0, (1 - 2 * nu) / 2]; % plane strain

                Bi = zeros(3, 2);
                B = zeros(3, Number_of_Nodes * 2);
            
                for i = 1:Number_of_Nodes
                    Bi(1, 1) = DN(1, i);
                    Bi(2, 2) = DN(2, i);
                    Bi(3, 1) = DN(2, i);
                    Bi(3, 2) = DN(1, i);
            
                    B(:, (i - 1) * 2 + 1 : i * 2) = Bi;
                end                
            elseif dim==3
                alphav=zeros(6,1);
                alphav(1:3,1)=[Dax,Day,Daz];   
               alphavdt=zeros(6,1);
                alphavdt(1:3,1)=[Daxdt,Daydt,Dazdt];
                C = DE./((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                    0 0 0 0 0 (1-2*nu)/2];

                Bi=zeros(6,3);B=zeros(6,Number_of_Nodes*3);
                    for i=1:Number_of_Nodes
                        Bi(1,1) = DN(1,i);
                        Bi(2,2) = DN(2,i); 
                        Bi(3,3) = DN(3,i);
                        Bi(4,:) = [DN(2,i),DN(1,i),0];
                        Bi(5,:) = [0,DN(3,i),DN(2,i)]; 
                        Bi(6,:) = [DN(3,i),0,DN(1,i)];
                        B(:,(i-1)*3+1:(i)*3)=Bi;
        
                    end
            end
                     F_TR=detJ*(B'*C*alphav*N);
                     F_T= detJ*(B'*C*alphav*N + B'*(C*alphavdt)*(N*(Tee-str2double(reader.T0))')*N);%*(Tee-str2double(reader.T0))
                     emp=[];
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [KUU, emp] = stiffnessintegrationFunction(~,natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)
            %% Thermoelectricity Simulation
            % This function calculates thermoelectric properties using finite element analysis.
            % Inputs:
            %
            % Outputs:
            %   - KJ: Jacobian matrix
            %   - R: Residues
            [dim] = mesh.retrieveelementdimension(etype);
            if dim==2
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), 0);
                element_coordinates=element_coordinates(1:2,:);
                N=N';
            else
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));
            end
                Number_of_Nodes = length(element_coordinates(1,:));
                JM = dShape' * element_coordinates';
                %Jacinv = inv(JM);
                DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
                %DN1 = JM \ dShape'; % FIXME and check it is the same!

                % FIXME, calculate from all dofs input
                Th = N * Tee';
            
                % FIXME: Calculate material properties
                %Dep = reader.getmaterialproperty(element_material_index,'ElectricalConductivity');
                %Dkp = reader.getmaterialproperty(element_material_index,'ThermalConductivity');
                %Dap = reader.getmaterialproperty(element_material_index,'Seebeck');
                %Dalphap = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient');
                DEp = reader.getmaterialproperty(element_material_index,'YoungModulus');
                nu = reader.getmaterialproperty(element_material_index,'PoissonRatio');

                %[De,Dde]=CalculateMaterialProperties(1e-6,Dep,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
                %[Da,Dda]=CalculateMaterialProperties(1e-6,Dap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
                %[Dk,Ddk]=CalculateMaterialProperties(reader.kmin,Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_YoungModulus'),reader.getmaterialproperty(element_material_index,'Tmax_YoungModulus')];
                [DE,DdE]=CalculateMaterialProperties(1e-6,DEp,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
                %[Dalpha,Ddalpha]=CalculateMaterialProperties(1e-6,Dalphap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalExpansionCoefficient'));
                
                %alphav=zeros(6,1);
                %alphav(1:3,1)=Dalpha;
                
                % Preallocate memory for B-Operators
 
            if dim==2

                C = DE / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2]; % plane stress
                % C = DE / (1 + nu) / (1 - 2 * nu) * [1 - nu, nu, 0; nu, 1
                % - nu, 0; 0, 0, (1 - 2 * nu) / 2]; % plane strain

                Bi = zeros(3, 2);
                B = zeros(3, Number_of_Nodes * 2);
            
                for i = 1:Number_of_Nodes
                    Bi(1, 1) = DN(1, i);
                    Bi(2, 2) = DN(2, i);
                    Bi(3, 1) = DN(2, i);
                    Bi(3, 2) = DN(1, i);
            
                    B(:, (i - 1) * 2 + 1 : i * 2) = Bi;
                end                
            elseif dim==3 

                C = DE./((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                    0 0 0 0 0 (1-2*nu)/2];

                Bi=zeros(6,3);B=zeros(6,Number_of_Nodes*3);
                    for i=1:Number_of_Nodes
                        Bi(1,1) = DN(1,i);
                        Bi(2,2) = DN(2,i); 
                        Bi(3,3) = DN(3,i);
                        Bi(4,:) = [DN(2,i),DN(1,i),0];
                        Bi(5,:) = [0,DN(3,i),DN(2,i)]; 
                        Bi(6,:) = [DN(3,i),0,DN(1,i)];
                        B(:,(i-1)*3+1:(i)*3)=Bi;
        
                    end
            end

                KUU=det(JM)*(B'*C*B);
                     emp=[];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function SolveLinearSystemInParallel(obj, desiredsolution,bcinit)
            if strcmp(desiredsolution,'decoupledthermoelectromechanical')
                % Extract the necessary variables
                dofs_free = bcinit.dofs_free_mech;
                KT_Distr = distributed(obj.KStiff(dofs_free, dofs_free)); 
                R_Distr=distributed(obj.Residual_mech(dofs_free));
                %[userview, sysview] = memory;
                %fprintf('Memory used by MATLAB before solving: %.2f GB\n', userview.MemUsedMATLAB / 1e9);      
                %spmd
                %    [userview, sysview] = memory;
                %    fprintf('Worker %d memory used before solving: %.2f GB\n', labindex, userview.MemUsedMATLAB / 1e9);
                %end
                dU =( KT_Distr ) \ ( R_Distr );  % Calculation of step
                obj.soldofs_mech(dofs_free)=obj.soldofs_mech(dofs_free)-dU;  % Calculation of step
            else
                % Extract the necessary variables
                dofs_free = bcinit.dofs_free_;
            
                KT_Distr = distributed(obj.KT(dofs_free, dofs_free)); 
                R_Distr=distributed(obj.Residual(dofs_free));
                %[userview, sysview] = memory;
                %fprintf('Memory used by MATLAB before solving (Mech): %.2f GB\n', userview.MemUsedMATLAB / 1e9);   
                %spmd
                %    [userview, sysview] = memory;
                %    fprintf('Worker %d memory used before solving: %.2f GB\n', labindex, userview.MemUsedMATLAB / 1e9);
                %end
                dU =( KT_Distr ) \ ( R_Distr );  % Calculation of step

                obj.soldofs(dofs_free)=obj.soldofs(dofs_free)+dU;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function residual_norm = runNewtonRaphson(obj,reader, mesh,bcinit)
                meshelements=mesh.retrieveElementalSelection(reader.MeshEntityName);
                etype=mesh.data.ElementTypes{meshelements(1)};
                dim = mesh.retrieveelementdimension(etype); 

            total_load = bcinit.loadVector_;
            for loadstep = 1: reader.NR_nloads
                bcinit.loadVector_=total_load/reader.NR_nloads*loadstep;
                % Initialize some parameters and initial guess
                residual = obj.tolerance * 1000;
                threshold=reader.solver_threshold;
                
                for iter = 1:obj.max_iterations
                    % Assembly the system matrix
                    [obj.KT, obj.Residual] = obj.Assembly_Thermoelectricity(reader, mesh);
                    %[userview, sysview] = memory;
                    %fprintf('Memory used by MATLAB at Assembly step: %.2f GB\n', userview.MemUsedMATLAB / 1e9);
                    %if (residual < obj.tolerance && iter > 1)
                        % Converged, return both the solution and the final residual
                    %    fprintf(append('### NR. CONVERGED with tolerance:',num2str(obj.tolerance),'\n'));
                    %    break;
                    %end
                
                    % Solve the system using the tangential matrix and residual: KT*dU=R->dU
                    obj.SolveLinearSystemInParallel('', bcinit)
                
                    % Calculate the residual (error)
                    dofs_free = bcinit.dofs_free_;
                    residual_norm = normest(obj.Residual(dofs_free)); % Calculate the norm
                
                    fprintf('### NR. Iteration %d residual %f\n', iter, residual_norm);
                    Tneg=0;
                    % Check if the change in residual is small
                    if abs(residual_norm) < threshold
                        fprintf('### NR. Convergence with residual change under: %f.\n',threshold);
                        break;
                    elseif (abs(residual_norm)>1000000 || isnan(abs(residual_norm)) || Tneg==1) && dim==3
                        fprintf('### NR. Diverged with residual: %f.\n',residual_norm);
                        odd_numbers = 1:2:length(obj.soldofs);
                        fprintf('### Min Temp: %f.\n',min(obj.soldofs(odd_numbers)));
                        fprintf('### Max Temp: %f.\n',max(obj.soldofs(odd_numbers)));
                        break
                    elseif (abs(residual_norm)>100000000 || isnan(abs(residual_norm)) || Tneg==1) && dim==2
                        fprintf('### NR. Diverged with residual: %f.\n',residual_norm);
                        odd_numbers = 1:2:length(obj.soldofs);
                        fprintf('### Min Temp: %f.\n',min(obj.soldofs(odd_numbers)));
                        fprintf('### Max Temp: %f.\n',max(obj.soldofs(odd_numbers)));
                        break
                    end
                end 
            end

            if strcmp(reader.physics,'decoupledthermoelectromechanical')
                [obj.KStiff,KThermalLoad,KTheta] = obj.Assembly_DecoupledThermoMech(reader, mesh);
                Temperature_solution = obj.soldofs(1:2:end);
                obj.KUT=KThermalLoad;
                obj.loadVector_mech=obj.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                %[userview, sysview] = memory;
                %fprintf('Memory used by MATLAB at Assembly step (Mech): %.2f GB\n', userview.MemUsedMATLAB / 1e9);
                %obj.SolveLinearSystemInParallel(reader.physics,bcinit)
                %%%%
                residual = obj.tolerance * 1000;
                threshold=reader.solver_threshold;
                    %aa=[obj.KStiff*obj.soldofs_mech,obj.loadVector_mech]

if isequal(KThermalLoad, KTheta)
    disp('The matrices are the same.');
else
    disp('The matrices are different.');
end                
                    obj.Residual_mech=(obj.KStiff*obj.soldofs_mech - obj.loadVector_mech);

                for iter = 1:obj.max_iterations
                    % Assembly the system matrix
                    %[userview, sysview] = memory;
                    %fprintf('Memory used by MATLAB at Assembly step: %.2f GB\n', userview.MemUsedMATLAB / 1e9);
                    %if (residual < obj.tolerance && iter > 1)
                        % Converged, return both the solution and the final residual
                    %    fprintf(append('### NR. CONVERGED with tolerance:',num2str(obj.tolerance),'\n'));
                    %    break;
                    %end
                
                    % Solve the system using the tangential matrix and residual: KT*dU=R->dU
                    obj.SolveLinearSystemInParallel(reader.physics,bcinit)
                    obj.Residual_mech=(obj.KStiff*obj.soldofs_mech - obj.loadVector_mech);
                    %aa=[obj.KStiff*obj.soldofs_mech,obj.loadVector_mech]
                    % Calculate the residual (error)
                    dofs_free = bcinit.dofs_free_mech;
                    residual_norm = normest(obj.Residual_mech(dofs_free)); % Calculate the norm
                
                    fprintf('### NR. Iteration %d residual %f\n', iter, residual_norm);
                    % Check if the change in residual is small
                    if abs(residual_norm) < threshold
                        fprintf('### NR. Convergence with residual change under: %f.\n',threshold);
                        break;
                    end
                end 
                
            end
            if iter==obj.max_iterations
            % If we reach here, the Newton-Raphson method did not converge
            warning('Newton-Raphson did not converge: iteration limit');
            residual_norm=1E6;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function residual = runModifiedNewtonRaphsonSolver(obj)
            % Initialize some parameters and initial guess
                % Assembly the system matrix
            [Ksparse,Rsparse] = Assembly_Thermoelectricity();
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
        function [R,Q,K] = funArcLen(obj,reader,bcinit,mesh,t)
             dofs_free = bcinit.dofs_free_;
             obj.soldofs(dofs_free)=t(1:end-1);
             %obj.soldofs=t(1:end-1);
             lambda=t(end);         
             [obj.KT, obj.Residual, R] = obj.Assembly_Thermoelectricity(reader, mesh);
             K=-obj.KT(dofs_free,dofs_free);
             R=R-obj.loadVector_copy*lambda;
             R=R(dofs_free);
             Q=obj.loadVector(dofs_free);
        end
        function [R,K] = funArcLen_LamMorley(obj,reader,bcinit,mesh,t)
             dofs_free = bcinit.dofs_free_;
             obj.soldofs(dofs_free)=t;
             %obj.soldofs=t(1:end-1);
             [K, Ral, R] = obj.Assembly_Thermoelectricity(reader, mesh);
             obj.KT=K;
             obj.Residual=Ral;
             K=-K(dofs_free,dofs_free);
             R=R(dofs_free);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = SingleSolve(obj,reader,mesh,bcinit)
            % Specify the source file and destination folder
            sourceFile = reader.filename;
            destinationFolder = append(reader.rst_folder,'input.txt');
            copyfile(sourceFile, destinationFolder);
            sourceFile = reader.meshFileName;
            destinationFolder = append(reader.rst_folder,'mesh.inp');
            copyfile(sourceFile, destinationFolder);
            if strcmp(reader.solver,'NR')
                    residual_norm=obj.runNewtonRaphson(reader, mesh, bcinit);
                    odd_numbers = 1:2:length(obj.soldofs);
                    Tdofs=obj.soldofs(odd_numbers);
                    if residual_norm>10000  || not(isempty(Tdofs(Tdofs<0)))% divergence in NR catch
                        warning('NR DIVERGED, changing to Arc-len!!!');
                        for i=1:length(bcinit.dofs_free_)
                            df=bcinit.dofs_free_(i);
                            if mod(df, 2)==0
                                obj.soldofs(df)=0.0;
                            else 
                                obj.soldofs(df)=0;
                            end
                        end
                        funAL = @(t) obj.funArcLen(reader,bcinit,mesh,t);
                        [ufree] = arc_length_Crisfield(funAL,obj.soldofs(bcinit.dofs_free_));
                        obj.soldofs(bcinit.dofs_free_)=ufree(:,end);
                        if strcmp(reader.physics,'decoupledthermoelectromechanical')
                            % Extract the necessary variables
                            [obj.KStiff,KThermalLoad] = obj.Assembly_DecoupledThermoMech(reader, mesh);
                            Temperature_solution = obj.soldofs(1:2:end);
                            obj.KUT=KThermalLoad;
                            obj.loadVector_mech=obj.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                            obj.SolveLinearSystemInParallel(reader.physics,bcinit)
                        end
                    end
                elseif strcmp(reader.obj,'Arc-len')
                        for i=1:length(obj.soldofs)/2
                            if obj.soldofs(i*2-1)==0
                                obj.soldofs(i*2-1)=str2double(reader.T0);
                            elseif obj.soldofs(i*2)==0
                                obj.soldofs(i*2)=0.01;
                            end
                        end
                    funAL = @(t) obj.funArcLen(reader,bcinit,mesh,t);
                    [ufree] = arc_length_Crisfield(funAL,obj.soldofs(bcinit.dofs_free_));
                    obj.soldofs(bcinit.dofs_free_)=ufree(:,end);
                    if strcmp(reader.physics,'decoupledthermoelectromechanical')
                        % Extract the necessary variables
                        [obj.KStiff,KThermalLoad] = obj.Assembly_DecoupledThermoMech(reader, mesh);
                        Temperature_solution = obj.soldofs(1:2:end);
                        obj.KUT=KThermalLoad;
                        obj.loadVector_mech=obj.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                        obj.SolveLinearSystemInParallel(reader.physics,bcinit)
                    end
                else
                    warning('No obj recognized, changing to Arc-len!!!');
                        for i=1:length(obj.soldofs)/2
                            if obj.soldofs(i*2-1)==0
                                obj.soldofs(i*2-1)=str2double(reader.T0);
                            end
                            if obj.soldofs(i*2)==0
                                obj.soldofs(i*2)=0.01;
                            end
                        end
                    funAL = @(t) obj.funArcLen(reader,bcinit,mesh,t);
                    [ufree] = arc_length_Crisfield(funAL,obj.soldofs(bcinit.dofs_free_));
                    obj.soldofs(bcinit.dofs_free_)=ufree(:,end);
                    if strcmp(reader.physics,'decoupledthermoelectromechanical')
                        % Extract the necessary variables
                        [obj.KStiff,KThermalLoad] = obj.Assembly_DecoupledThermoMech(reader, mesh);
                        Temperature_solution = obj.soldofs(1:2:end);
                        obj.KUT=KThermalLoad;
                        obj.loadVector_mech=obj.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                        obj.SolveLinearSystemInParallel(reader.physics,bcinit)
                    end
            end
        end
    end
    
    methods (Access = private)
        function result = thermoelectricityintegration(obj, natcoords, coords, dofs, elementTag)
            % Implement your code here
        end
    end
end

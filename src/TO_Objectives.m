classdef TO_Objectives < handle
    %TO_OBJECTIVES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ObjectiveName
        MeshName
        fval
        dfdx
        TOEL
        objective_dofs
        freedofs
        Number_of_dofs
    end
    
    methods
        function obj = TO_Objectives()
            % Initialize mesh TO parameters
            obj.TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            objective_nodes=mesh.retrieveNodalSelection(reader.TopOpt_ObjectiveSelection);
            obj.objective_dofs=(objective_nodes-1)*2+1;
            obj.freedofs=bcinit.dofs_free_;
            obj.Number_of_dofs=length(mesh.data.NODE)*2;
        end
        
        function CalculateObjective(obj,mesh,solver)
            switch obj.ObjectiveName
                case 'AverageTemperature'
                    obj.fval_AverageTemp(mesh,solver)
                    obj.dfdx_AverageTemp(mesh,solver)
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function GaussIntegration_dx(obj,dimension, order, elementTag, mesh, initialdofs,reader,etype)
            
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

            integrationFunction = @(natcoords) obj.integration_AvgTemp_dx(natcoords, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype, mesh.elements_density(elementTag));      
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fval_AverageTemp(obj,reader,mesh,solver)
            objective_nodes=mesh.retrieveNodalSelection(reader.TopOpt_ObjectiveSelection);
            obj.fval=sum(solver.soldofs((objective_nodes-1)*2+1))/length(solver.soldofs((objective_nodes-1)*2+1));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dfdx_AverageTemp(~,reader,mesh,solver,dfdx_index)

            % Initialize solver matrixes
            KJa=-solver.KT;
            AdjT=zeros(obj.Number_of_dofs,1);
            LAdj=zeros(obj.Number_of_dofs,1);
            LAdj(obj.objective_dofs)=1/length(obj.objective_dofs);
            element_sensitivities=zeros(length(obj.TOEL),1);

            % Solve adjoint equation
            A = distributed((-KJa(obj.freedofs,obj.freedofs))'); 
            B=distributed((LAdj(obj.freedofs)));
            AdjT(obj.freedofs)=A\B;

            % Calculate material derivatives
            for ii=1:length(obj.TOEL)
                element_Tag=obj.TOEL(ii);

                [R_dx,element_dofs]=GaussIntegration_dx(3, 14, element_Tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_Tag}) ; 
            
                ADJ_element=AdjT(element_dofs)';%
                element_sensitivities(ii)=ADJ_element*R_dx;
            
            end
                
            % Store elemental sensitivities
            obj.dfdx(1:length(obj.TOEL),dfdx_index)=element_sensitivities;

            % Calculate remaining sensititivies
            % if voltage constraint:
            %doforderV=(Vfnod-1)*2+2;
            %dU1=zeros(length(AdjT),1);
            %dU1(doforderV)=dV;
            %prodF=KJa*dU1;
            %Overall_sensitivities(end,1)=LAdj'*dU1+AdjT'*(+prodF);

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function[Rx]=integration_AvgTemp_dx(~,natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)        
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));

                JM = dShape' * element_coordinates';

                %Jacinv = inv(JM);
                DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
                % FIXME, calculate from all dofs input
                Th = N * Tee';

                Dep = reader.getmaterialproperty(element_material_index,'ElectricalConductivity');
                Dkp = reader.getmaterialproperty(element_material_index,'ThermalConductivity');
                Dap = reader.getmaterialproperty(element_material_index,'Seebeck');
                
                [De]=CalculateMaterialProperties(Dep,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
                [Da]=CalculateMaterialProperties(Dap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
                %[Dk,Ddk]=CalculateMaterialProperties(Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));
                
                % notice that mat(T) and T=f(U) and we only need the
                % partial to respect to x. Furthermore, as we do the
                % derivative delta_R/delat_x, the temperature derivatives
                % do not influence the result.
                [De_dx]=CalculateMaterial_XDerivative(Dep,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
                [Da_dx]=CalculateMaterial_XDerivative(Dap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
                [Dk_dx]=CalculateMaterial_XDerivative(Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));
                
                Vee=Vee';
                Tee=Tee';
                detJ = det(JM);
                
                % Calculate current density and heat flux
                je = -De * DN * Vee - Da * De * DN * Tee; 
                %qe = Da * (N * Tee) * je - Dk * DN * Tee; % Not needed.

                % all matrix
                jx=-De_dx*DN*Vee-Da_dx*De*DN*Tee-Da*Ded*DN*Tee;
                qx=Da_dx*(N'*Tee)*je+Da*(N'*Tee)*jx-Dk_dx*DN*Tee; 
        
                RAx=detJ*(-DN'*qx+N*(jx'*DN*Vee));
                RBx=detJ*(-DN'*jx);

                Rx=[RAx
                   RBx];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end


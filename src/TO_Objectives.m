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
        n
    end
    
    methods
        function obj = TO_Objectives(reader,mesh,bcinit)
            % Initialize mesh TO parameters
            obj.TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            obj.n =length(obj.TOEL)+length(reader.TObcval);
            objective_nodes=mesh.retrieveNodalSelection(reader.TopOpt_ObjectiveSelection);
            obj.objective_dofs=(objective_nodes-1)*2+1;
            obj.freedofs=bcinit.dofs_free_;
            obj.Number_of_dofs=length(mesh.data.NODE)*2;
            obj.dfdx=zeros(obj.n,1);
        end
        
        function CalculateObjective(obj,reader,mesh,solver)
            switch reader.TopOpt_Objective
                case 'AverageTemperature'
                    obj.fval_AverageTemp(reader,mesh,solver)
                    obj.dfdx_AverageTemp(reader,mesh,solver)
                case 'PnormTemperature'
                    obj.fval_PnormTemp(reader,mesh,solver)
                    obj.dfdx_PnormTemp(reader,mesh,solver)
                case 'Heat'
                    obj.fval_Heat(reader,mesh,solver)
                    obj.dfdx_Heat(reader,mesh,solver)
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [K,element_dof_indexes] = GaussIntegration_dx(obj,dimension, order, elementTag, mesh, initialdofs,reader,etype)
                                        %3, 14, element_Tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_Tag}
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

            integrationFunction = @(natcoords) obj.integration_R_dx(natcoords, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype, mesh.elements_density(elementTag));      
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flagGloop=0;
             if dimension == 1
                % 1D integration using a single loop.
                natcoords = zeros(1, 1);
                for i = 1:size(weights, 1)
                    natcoords(1) = gaussPoints(i);
                    % Explicitly use the element-wise multiplication .* for arrays
                    Ke= integrationFunction(natcoords) ;
                    Ke=Ke.* weights(i);
                    K = K + Ke;
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
                        K = K + Ke;
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
                        Ke = integrationFunction(natcoords) ;
                        Ke=Ke .* (weights(i));
                        if flagGloop==1
                                    K = K + Ke;
                        else
                                    flagGloop=1;
                                    K=Ke;
                        end       
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
                                Ke= integrationFunction(natcoords) ;
                                Ke=Ke .* (weights(i) * weights(j) * weights(k));
                                if flagGloop==1
                                            K = K + Ke;
                                else
                                            flagGloop=1;
                                            K=Ke;
                                end       
                            end
                        end
                    end
                end
            else
                fprintf('Invalid dimension for Gauss integration.\n');
                K = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fvalue]=fval_AverageTemp(obj,reader,mesh,solver)

            objective_nodes=mesh.retrieveNodalSelection(reader.TopOpt_ObjectiveSelection);
            fvalue=sum(solver.soldofs((objective_nodes-1)*2+1))/length(solver.soldofs((objective_nodes-1)*2+1));
            obj.fval=fvalue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fvalue]=fval_PnormTemp(obj,reader,mesh,solver)
            objective_nodes=mesh.retrieveNodalSelection(reader.TopOpt_ObjectiveSelection);
            odd_indices = 1:2:length(solver.soldofs); % Get even indices
            Tsol_ref=solver.soldofs(odd_indices)/str2double(reader.T0);
            Tsol_ref=Tsol_ref.^reader.KSUp;
            L=zeros(length(Tsol_ref),1);
            L(objective_nodes)=1/length(objective_nodes);
            fvalue=(L'*Tsol_ref)^(1/reader.KSUp);
            obj.fval=fvalue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fvalue]=fval_Heat(obj,reader,mesh,solver)
            str = reader.TopOpt_Objective;
            parts = split(str, '_');
            direction = str2double(parts(2:4));  
            fvalue=CalculateHeat_direction(reader,mesh,solver,reader.TopOpt_ObjectiveSelection,direction);
            obj.fval=fvalue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dfdx_AverageTemp(obj,reader,mesh,solver)

            % Initialize solver matrixes
            AdjT=zeros(obj.Number_of_dofs,1);
            LAdj=zeros(obj.Number_of_dofs,1);
            LAdj(obj.objective_dofs)=1/length(obj.objective_dofs);
            element_sensitivities=zeros(length(obj.TOEL),1);

            % Solve adjoint equation
            A = distributed((solver.KT(obj.freedofs,obj.freedofs))'); 
            B=distributed((LAdj(obj.freedofs)));
            AdjT(obj.freedofs)=A\B;
            % Calculate material derivatives
            parfor ii=1:length(obj.TOEL)
                element_Tag=obj.TOEL(ii);

                [R_dx,element_dofs]=obj.GaussIntegration_dx(3, 5, element_Tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_Tag}) ; 
            
                ADJ_element=AdjT(element_dofs)';%
                element_sensitivities(ii)=ADJ_element*R_dx;
            
            end
                
            % Store elemental sensitivities
            obj.dfdx(1:length(obj.TOEL))=element_sensitivities;

            % Calculate remaining sensititivies
            bcvariablenames=reader.TObctype;
            for i=1:length(bcvariablenames)
                if strcmp(bcvariablenames(i), 'Voltage')
                    bcvariable_nodes = mesh.retrieveNodalSelection(reader.TObcloc(i));
                    bcvariable_dofs=bcvariable_nodes*2; % Voltage fixed nodes (pairs)
                    dx_bc=zeros(obj.Number_of_dofs,1);
                    dx_bc(bcvariable_dofs)=reader.TObcmaxval(i)-reader.TObcminval(i);
                    prodF=-solver.KT*dx_bc;
                    obj.dfdx(length(obj.TOEL)+i)=LAdj'*dx_bc+AdjT'*(+prodF);
                end
            end


        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dfdx_PnormTemp(obj,reader,mesh,solver)

            % Initialize solver matrixes
            AdjT=zeros(obj.Number_of_dofs,1);
            LAdj=zeros(obj.Number_of_dofs,1);

            objective_nodes=mesh.retrieveNodalSelection(reader.TopOpt_ObjectiveSelection);
            odd_indices = 1:2:length(solver.soldofs); % Get even indices
            Tsol_ref=solver.soldofs(odd_indices)/str2double(reader.T0);

            L=zeros(length(Tsol_ref),1);
            L(objective_nodes)=1/length(objective_nodes);
            LAdj(obj.objective_dofs)=(L'*(Tsol_ref.^reader.KSUp)) ^ (1/reader.KSUp-1) * L'*(Tsol_ref.^(reader.KSUp-1));
            element_sensitivities=zeros(length(obj.TOEL),1);

            % Solve adjoint equation
            A = distributed((solver.KT(obj.freedofs,obj.freedofs))'); 
            B=distributed((LAdj(obj.freedofs)));
            AdjT(obj.freedofs)=A\B;
            % Calculate material derivatives
            parfor ii=1:length(obj.TOEL)
                element_Tag=obj.TOEL(ii);

                [R_dx,element_dofs]=obj.GaussIntegration_dx(3, 5, element_Tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_Tag}) ; 
            
                ADJ_element=AdjT(element_dofs)';%
                element_sensitivities(ii)=ADJ_element*R_dx;
            
            end
                
            % Store elemental sensitivities
            obj.dfdx(1:length(obj.TOEL))=element_sensitivities;

            % Calculate remaining sensititivies
            bcvariablenames=reader.TObctype;
            for i=1:length(bcvariablenames)
                if strcmp(bcvariablenames(i), 'Voltage')
                    bcvariable_nodes = mesh.retrieveNodalSelection(reader.TObcloc(i));
                    bcvariable_dofs=bcvariable_nodes*2; % Voltage fixed nodes (pairs)
                    dx_bc=zeros(obj.Number_of_dofs,1);
                    dx_bc(bcvariable_dofs)=reader.TObcmaxval(i)-reader.TObcminval(i);
                    prodF=-solver.KT*dx_bc;
                    obj.dfdx(length(obj.TOEL)+i)=LAdj'*dx_bc+AdjT'*(+prodF);
                end
            end


        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function[Rx,flag]=integration_R_dx(~,natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)        
                flag=[];
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
                jx=-De_dx*DN*Vee-Da_dx*De*DN*Tee-Da*De_dx*DN*Tee;
                qx=Da_dx*(Th)*je+Da*(Th)*jx-Dk_dx*DN*Tee; 
        
                RAx=detJ*(-DN'*qx+(N*(jx'*DN*Vee))');
                RBx=detJ*(-DN'*jx);

                Rx=[RAx
                   RBx];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dfdx_Heat(obj,reader,mesh,solver,dfdx_index)
            % Initialize solver matrixes
            str = reader.TopOpt_Objective;
            parts = split(str, '_');
            direction = str2double(parts(2:4)); 

            element_sensitivities=zeros(length(obj.TOEL),1);
            ADJP=zeros(obj.Number_of_dofs,1);
            dofs_per_node = 2;
            number_of_TO_elements = length(obj.TOEL);
            nodes_per_element = length(mesh.data.ELEMENTS{obj.TOEL(1)});
            LP_dU_element=zeros(number_of_TO_elements,nodes_per_element*dofs_per_node);
            LP_dU_element_dofs=zeros(number_of_TO_elements,nodes_per_element*dofs_per_node);
            LP_U_element=zeros(number_of_TO_elements,1);
            LP_dU=zeros(obj.Number_of_dofs,1);
            LP_U=zeros(nodes_per_element,1);

            %% Derivatives to Vf
            GaussfunctionTag=@(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) obj.integration_Heat_dx(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
            parfor  ii=1:length(obj.TOEL)
                element_Tag = obj.TOEL(ii);
                [L_dU,L_U,element_dofs]=obj.GaussIntegration_dx(2, 5, element_Tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_Tag},GaussfunctionTag) ;
                % assembly in global residual and jacobian matrix in sparse format
                LP_dU_element_dofs(ii,:)=element_dofs;
                LP_dU_element(ii,:)=L_dU*direction;
                LP_U_element(ii)=L_U*direction;
            end

            % Summation of all element integrations
            for ii=1:length(obj.TOEL)
                LP_dU(LP_dU_element_dofs(ii,:),1)=LP_dU(LP_dU_element_dofs(ii,:),1)+LP_dU_element(ii,:)';
                LP_U(obj.TOEL(ii))=LP_U_element(ii);
            end

            %% Power sensitivity
            ADJP(obj.freedofs)=(solver.KT(obj.freedofs,obj.freedofs))'\LP_dU(obj.freedofs);

            GaussfunctionTag=@(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) obj.integration_R_dx(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
            parfor ii=1:length(obj.TOEL)
                element_Tag=obj.TOEL(ii);
                [Rx,flag,element_dofs]=obj.GaussIntegration_dx(3, 5, element_Tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_Tag},GaussfunctionTag) ;
                element_sensitivities(ii)=LP_U(element_Tag)+ADJP(element_dofs)'*Rx;
            end

            obj.dfdx(dfdx_index,1:length(obj.TOEL))=element_sensitivities;

            % Calculate remaining bc sensititivies
            bcvariablenames=reader.TObctype;
            for i=1:length(bcvariablenames)
                if strcmp(bcvariablenames(i), 'Voltage')
                    bcvariable_nodes = mesh.retrieveNodalSelection(reader.TObcloc(i));
                    bcvariable_dofs=bcvariable_nodes*2;
                    dx_bc=zeros(obj.Number_of_dofs,1);
                    dx_bc(bcvariable_dofs)=reader.TObcmaxval(i)-reader.TObcminval(i);
                    prodF=-solver.KT*dx_bc;
                    obj.dfdx(dfdx_index,length(obj.TOEL)+i)=LP_dU(bcvariable_dofs)'*dx_bc(bcvariable_dofs)...
                        +ADJP(obj.freedofs)'*(+prodF(obj.freedofs));
                end
            end

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function[L_U,L_dU]=integration_Heat_dx(~,natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)
            [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));

            JM = dShape' * element_coordinates';

            %Jacinv = inv(JM);
            DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
            % FIXME, calculate from all dofs input
            Th = N * Tee';

            Dep = reader.getmaterialproperty(element_material_index,'ElectricalConductivity');
            %Dkp = reader.getmaterialproperty(element_material_index,'ThermalConductivity');
            Dap = reader.getmaterialproperty(element_material_index,'Seebeck');

            %[De,De_DT]=CalculateMaterialProperties(Dep,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
            %[Da,Da_DT]=CalculateMaterialProperties(Dap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
            [Dk,Ddk]=CalculateMaterialProperties(Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));

            % notice that mat(T) and T=f(U) and we only need the
            % partial to respect to x. Furthermore, as we do the
            % derivative delta_R/delat_x, the temperature derivatives
            % do not influence the result.
            %[De_dx]=CalculateMaterial_XDerivative(Dep,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
            %[Da_dx]=CalculateMaterial_XDerivative(Dap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
            [Dk_dx]=CalculateMaterial_XDerivative(Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));

            Vee=Vee';
            Tee=Tee';
            detJ = det(JM);

            % Calculate current density and heat flux
            %je = -De * DN * Vee - Da * De * DN * Tee;
            %qe = Da * (N * Tee) * je - Dk * DN * Tee; % Not needed.
            %qeth =  - Dk * DN * Tee; % Not needed.

            dqthdxi_t=- Dk * DN;
            dqthdxi_v=zeros(length(dqdxi_t),1);

            qth_dx=detJ*(-Dk_dx* DN * Tee);

            K11=detJ*(dqthdxi_t);
            K12=detJ*(dqthdxi_v);

            L_U=qth_dx;
            L_dU=[K11 K12];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end
end


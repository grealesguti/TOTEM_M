classdef TO_Constraints < handle
    %TO_OBJECTIVES Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ConstraintNames
        MeshName
        fval
        dfdx
        TOEL
        freedofs
        Number_of_dofs
        Elements_volume
        V_TOT
    end

    methods
        function obj = TO_Constraints(reader,mesh,bcinit)
            % Initialize size of obj. variables
            m=length(reader.TopOpt_ConstraintName);
            obj.fval=zeros(m,1);
            obj.TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            dV=0;
            if dV~=0
                n=length(obj.TOEL)+1;
            else
                n=length(obj.TOEL);
            end
            obj.dfdx=zeros(m,n);
            % Initialize mesh TO parameters
            obj.freedofs=bcinit.dofs_free_;
            obj.Number_of_dofs=length(mesh.data.NODE)*2;
            obj.Elements_volume = obj.CalculateAllElementVolume(mesh); 
            obj.V_TOT=sum(obj.Elements_volume);
        end

        function CalculateConstraint(obj,reader,mesh,solver)
            index_constraint = 1;
            for i=1:length(reader.TopOpt_ConstraintName)
                ConstraintName=reader.TopOpt_ConstraintName{i};
                switch ConstraintName
                    case 'Power'
                        obj.fval_Power(reader,mesh,solver,index_constraint)
                        obj.dfdx_Power(reader,mesh,solver,index_constraint)
                        index_constraint=index_constraint+1;
                    case 'Volume'
                        obj.fval_Volume(reader,mesh,solver,index_constraint)
                        obj.dfdx_Volume(reader,index_constraint)
                        index_constraint=index_constraint+1;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [K,R,element_dof_indexes] = GaussIntegration_dx(~,dimension, order, elementTag, mesh, initialdofs,reader,etype,integrationFunctionhandle)
            %3, 14, element_Tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_Tag}
            if dimension < 1 || order < 1
                fprintf('Invalid dimension or order for Gauss integration.\n');
                K = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
                return;
            end

            [weights, gaussPoints] = getGaussWeightsAndPoints(order);

            if isempty(weights) || isempty(gaussPoints)
                fprintf('Invalid order for Gauss integration.\n');
                K = zeros(1, 1); % Initialize result to a 1x1 matrix with zero value.
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

            integrationFunction = @(natcoords) integrationFunctionhandle(natcoords, element_coordinates, Tee, Vee, element_material_index,reader, mesh, etype, mesh.elements_density(elementTag));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flagGloop=0;
            if dimension == 1
                % 1D integration using a single loop.
                natcoords = zeros(1, 1);
                for i = 1:size(weights, 1)
                    natcoords(1) = gaussPoints(i);
                    % Explicitly use the element-wise multiplication .* for arrays
                    [Ke,Re] = integrationFunction(natcoords) ;
                    Ke=Ke.* weights(i);
                    Re=Re .* (weights(i) );
                    if flagGloop==1
                        K = K + Ke;
                        R = R +Re;
                    else
                        flagGloop=1;
                        K=Ke;
                        R =Re;
                    end
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
                        Re=Re .* (weights(i) * weights(j));
                        if flagGloop==1
                            K = K + Ke;
                            R = R +Re;
                        else
                            flagGloop=1;
                            K=Ke;
                            R =Re;
                        end
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
                        if flagGloop==1
                            K = K + Ke;
                            R = R +Re;
                        else
                            flagGloop=1;
                            K=Ke;
                            R =Re;
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
                                [Ke,Re]= integrationFunction(natcoords) ;
                                Ke=Ke .* (weights(i) * weights(j) * weights(k));
                                Re=Re .* (weights(i) * weights(j) * weights(k));
                                if flagGloop==1
                                    K = K + Ke;
                                    R = R +Re;
                                else
                                    flagGloop=1;
                                    K=Ke;
                                    R =Re;
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
        function fvalue=fval_Power(obj,reader,mesh,solver,index_con)
            fvalue=CalculatePower(reader,mesh,solver)/reader.TopOpt_ConstraintValue(index_con)-1;
            obj.fval(index_con)=fvalue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dfdx_Power(obj,reader,mesh,solver,dfdx_index)
            % Initialize solver matrixes
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
            Pobj = reader.TopOpt_ConstraintValue(dfdx_index);
            
            %% Derivatives to Vf
            GaussfunctionTag=@(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) obj.integration_Power_dx_1(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
            parfor  ii=1:length(obj.TOEL)
                element_Tag = obj.TOEL(ii);
                [LJ,dPdxi_c,element_dofs]=obj.GaussIntegration_dx(3, 14, element_Tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_Tag},GaussfunctionTag) ;
                % assembly in global residual and jacobian matrix in sparse format
                LP_dU_element_dofs(ii,:)=element_dofs;
                LP_dU_element(ii,:)=LJ;
                LP_U_element(ii)=dPdxi_c;
            end
            
            % Summation of all element integrations
            for ii=1:length(obj.TOEL)
                LP_dU(LP_dU_element_dofs(ii,:),1)=LP_dU(LP_dU_element_dofs(ii,:),1)+LP_dU_element(ii,:)';
                LP_U(obj.TOEL(ii))=LP_U_element(ii);
            end

            LP_U=LP_U/Pobj;
            LP_dU=LP_dU/Pobj;

            %% Power sensitivity
            ADJP(obj.freedofs)=(solver.KT(obj.freedofs,obj.freedofs))'\LP_dU(obj.freedofs);

            GaussfunctionTag=@(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) obj.integration_R_dx(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
            parfor ii=1:length(obj.TOEL)
                element_Tag=obj.TOEL(ii);
                [Rx,flag,element_dofs]=obj.GaussIntegration_dx(3, 14, element_Tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_Tag},GaussfunctionTag) ;
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
        function[K_dx,P_dx]=integration_Power_dx_1(~,natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)
            [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));

            JM = dShape' * element_coordinates';

            %Jacinv = inv(JM);
            DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
            % FIXME, calculate from all dofs input
            Th = N * Tee';

            Dep = reader.getmaterialproperty(element_material_index,'ElectricalConductivity');
            %Dkp = reader.getmaterialproperty(element_material_index,'ThermalConductivity');
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
            %[Dk_dx]=CalculateMaterial_XDerivative(Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));

            Vee=Vee';
            Tee=Tee';
            detJ = det(JM);

            % Calculate current density and heat flux
            je = -De * DN * Vee - Da * De * DN * Tee;
            %qe = Da * (N * Tee) * je - Dk * DN * Tee; % Not needed.

            dPdxi_t=Vee'*DN'*Da*De*DN;
            dPdxi_v=Vee'*DN'*De*DN-je'*DN;

            P_dx=detJ*(-Vee'*DN'*(-De_dx*DN*Vee-Da_dx*De*DN*Tee-Da*De_dx*DN*Tee));

            K11=detJ*(dPdxi_t);
            K12=detJ*(dPdxi_v);

            K_dx=[K11 K12];
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
        function [Elements_volume]=CalculateAllElementVolume(obj,mesh)
            Elements_volume=zeros(length(obj.TOEL),1);
            for i=1:length(obj.TOEL)
                element_Tag=obj.TOEL(i);
                element_nodes=mesh.data.ELEMENTS{element_Tag};
                coordinates=zeros(3,8);
                for j=1:8
                    coordinates(:,j)=mesh.data.NODE{element_nodes(j)};
                end
                Elements_volume(i)=CalculateHexVolume(coordinates');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fvalue=fval_Volume(obj,reader,mesh,solver,dfdx_index)
            Vpobj = reader.TopOpt_ConstraintValue(dfdx_index);
            Vx=obj.Elements_volume/(Vpobj*obj.V_TOT);
            fvalue=(Vx'*mesh.elements_density(obj.TOEL)')-1;
            obj.fval(dfdx_index)=fvalue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dfdx_Volume(obj,reader,dfdx_index)
            Vpobj = reader.TopOpt_ConstraintValue(dfdx_index);
            Vx=obj.Elements_volume/(Vpobj*obj.V_TOT);
            obj.dfdx(dfdx_index,1:length(obj.TOEL))=Vx';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end

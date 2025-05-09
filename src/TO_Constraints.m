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
        freedofs_mech
        Number_of_dofs
        Number_of_nodes
        Elements_volume
        V_TOT
        SVM
        dim
    end

    methods
        function obj = TO_Constraints(reader,mesh,bcinit)
            % Initialize size of obj. variables
            meshelements=mesh.retrieveElementalSelection(reader.MeshEntityName);
            etype=mesh.data.ElementTypes{meshelements(1)};
            obj.dim = mesh.retrieveelementdimension(etype); 
            m=length(reader.TopOpt_ConstraintName);
            obj.fval=zeros(m,1);
            if reader.TopOpt_DesignElements==""
                obj.TOEL=[];
            else
                obj.TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            end            
            %obj.TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            n =length(obj.TOEL)+length(reader.TObcval);
            obj.dfdx=zeros(m,n);
            % Initialize mesh TO parameters
            obj.freedofs=bcinit.dofs_free_;
            obj.freedofs_mech=bcinit.dofs_free_mech;
            obj.Number_of_dofs=length(mesh.data.NODE)*2;
            obj.Number_of_nodes=length(mesh.data.NODE);
            obj.Elements_volume = obj.CalculateAllElementVolume(mesh);
            obj.V_TOT=sum(obj.Elements_volume);

        end
    
        function CalculateConstraint(obj, reader, mesh, solver)
            index_constraint = 1;
            for i = 1:length(reader.TopOpt_ConstraintName)
                ConstraintName = reader.TopOpt_ConstraintName{i};
                if startsWith(ConstraintName, 'Power')
                    obj.fval_Power(reader, mesh, solver, index_constraint);
                    obj.dfdx_Power(reader, mesh, solver, index_constraint);
                    index_constraint = index_constraint + 1;
                elseif startsWith(ConstraintName, 'Volume')
                    obj.fval_Volume(reader, mesh, solver, index_constraint);
                    obj.dfdx_Volume(reader, index_constraint);
                    index_constraint = index_constraint + 1;
                elseif startsWith(ConstraintName, 'Stress')
                    obj.fval_Stress(reader, mesh, solver, index_constraint);
                    obj.dfdx_Stress(reader, mesh, solver, index_constraint);
                    index_constraint = index_constraint + 1;
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
            if isempty(obj.TOEL)
                nodes_per_element=20;
            else
                nodes_per_element = length(mesh.data.ELEMENTS{obj.TOEL(1)});
            end

            mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);
            etype=mesh.data.ElementTypes{mesh_elements(1)};
            dim = mesh.retrieveelementdimension(etype); 
            totalelem = length(mesh_elements);

            LP_dU_element=zeros(totalelem,nodes_per_element*dofs_per_node);
            LP_dU_element_dofs=zeros(totalelem,nodes_per_element*dofs_per_node);
            LP_U_element=zeros(totalelem,1);
            LP_dU=zeros(obj.Number_of_dofs,1);
            LP_U=zeros(nodes_per_element,1);
            Pobj = reader.TopOpt_ConstraintValue(dfdx_index);

            %% Derivatives to Vf

            if isempty(obj.TOEL) % FIXME: in cases we only use voltage optim, no densities
               etype="C3D20";
            else
                etype=mesh.data.ElementTypes{obj.TOEL(1)};
            end           % etype=mesh.data.ElementTypes{obj.TOEL(1)};
            dim = mesh.retrieveelementdimension(etype); 
            GaussfunctionTag=@(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) obj.integration_Power_dx_1(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
            %mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);



            parfor  ii=1:totalelem
                element_Tag = mesh_elements(ii);
                [LJ,dPdxi_c,element_dofs]=obj.GaussIntegration_dx(dim,  reader.GI_order, element_Tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_Tag},GaussfunctionTag) ;
                % assembly in global residual and jacobian matrix in sparse format
                LP_dU_element_dofs(ii,:)=element_dofs;
                LP_dU_element(ii,:)=LJ;
                LP_U_element(ii)=dPdxi_c;
            end

            % Summation of all element integrations
            for ii=1:totalelem
                LP_dU(LP_dU_element_dofs(ii,:),1)=LP_dU(LP_dU_element_dofs(ii,:),1)+LP_dU_element(ii,:)';
                LP_U(mesh_elements(ii))=LP_U_element(ii);
            end

            LP_U=LP_U/Pobj;
            LP_dU=LP_dU/Pobj;

            %% Power sensitivity
            ADJP(obj.freedofs)=distributed(solver.KT(obj.freedofs,obj.freedofs))'\distributed(LP_dU(obj.freedofs));
            if isempty(obj.TOEL) % FIXME: in cases we only use voltage optim, no densities
               etype="C3D20";
            else
                etype=mesh.data.ElementTypes{obj.TOEL(1)};
            end           % etype=mesh.data.ElementTypes{obj.TOEL(1)};
            dim = mesh.retrieveelementdimension(etype); 
            GaussfunctionTag=@(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) obj.integration_R_dx(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
            parfor ii=1:length(obj.TOEL)
                element_Tag=obj.TOEL(ii);
                [Rx,flag,element_dofs]=obj.GaussIntegration_dx(dim,  reader.GI_order, element_Tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_Tag},GaussfunctionTag) ;
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
           % [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));
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
            % FIXME, calculate from all dofs input
            Th = N * Tee';

            Dep = reader.getmaterialproperty(element_material_index,'ElectricalConductivity');
            %Dkp = reader.getmaterialproperty(element_material_index,'ThermalConductivity');
            Dap = reader.getmaterialproperty(element_material_index,'Seebeck');
            Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_ElectricalConductivity'),reader.getmaterialproperty(element_material_index,'Tmax_ElectricalConductivity')];
            [De,De_DT]=CalculateMaterialProperties(1e-6,Dep,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
            [De_dx]=CalculateMaterial_XDerivative(1e-6,Dep,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
            Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_Seebeck'),reader.getmaterialproperty(element_material_index,'Tmax_Seebeck')];
            [Da,Da_DT]=CalculateMaterialProperties(1e-6,Dap,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
            %[Dk,Ddk]=CalculateMaterialProperties(reader.kmin,Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));

            % notice that mat(T) and T=f(U) and we only need the
            % partial to respect to x. Furthermore, as we do the
            % derivative delta_R/delat_x, the temperature derivatives
            % do not influence the result.
            [Da_dx]=CalculateMaterial_XDerivative(1e-6,Dap,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
            %[Dk_dx]=CalculateMaterial_XDerivative(reader.kmin,Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));

            Vee=Vee';
            Tee=Tee';
            detJ = det(JM);

            % Calculate current density and heat flux
            je = -De * DN * Vee - Da * De * DN * Tee;
            %qe = Da * (N * Tee) * je - Dk * DN * Tee; % Not needed.

            djdt = -De*(Da*DN+(Da_DT*DN*Tee-je/De^2*De_DT)*N);
            djdv = -De*DN;

            dPdxi_t=-Vee'*DN'*djdt;
            dPdxi_v=-Vee'*DN'*djdv-je'*DN;

            djdx = (-De_dx*DN*Vee-Da_dx*De*DN*Tee-Da*De_dx*DN*Tee);

            P_dx=detJ*(-Vee'*DN'*djdx);

            K11=detJ*(dPdxi_t);
            K12=detJ*(dPdxi_v);

            K_dx=[K11 K12];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function[Rx,flag]=integration_R_dx(~,natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)
            flag=[];
            % [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));
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
            % FIXME, calculate from all dofs input
            Th = N * Tee';

            Dep = reader.getmaterialproperty(element_material_index,'ElectricalConductivity');
            Dkp = reader.getmaterialproperty(element_material_index,'ThermalConductivity');
            Dap = reader.getmaterialproperty(element_material_index,'Seebeck');
            Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_ElectricalConductivity'),reader.getmaterialproperty(element_material_index,'Tmax_ElectricalConductivity')];
            [De]=CalculateMaterialProperties(1e-6,Dep,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
            [De_dx]=CalculateMaterial_XDerivative(1e-6,Dep,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
            Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_Seebeck'),reader.getmaterialproperty(element_material_index,'Tmax_Seebeck')];
            [Da]=CalculateMaterialProperties(1e-6,Dap,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
            [Da_dx]=CalculateMaterial_XDerivative(1e-6,Dap,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));

            %[Dk,Ddk]=CalculateMaterialProperties(reader.kmin,Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));

            % notice that mat(T) and T=f(U) and we only need the
            % partial to respect to x. Furthermore, as we do the
            % derivative delta_R/delat_x, the temperature derivatives
            % do not influence the result.
            Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_ThermalConductivity'),reader.getmaterialproperty(element_material_index,'Tmax_ThermalConductivity')];
            [Dk_dx]=CalculateMaterial_XDerivative(reader.kmin,Dkp,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));

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
                if obj.dim == 2
                    coordinates=zeros(3,4);
                    for j=1:4
                        coordinates(:,j)=mesh.data.NODE{element_nodes(j)};
                    end                    
                    Elements_volume(i)=CalculateQuadArea(coordinates');
                else
                    coordinates=zeros(3,8);
                    for j=1:8
                        coordinates(:,j)=mesh.data.NODE{element_nodes(j)};
                    end
                Elements_volume(i)=CalculateHexVolume(coordinates');
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fvalue=fval_Volume(obj,reader,mesh,solver,dfdx_index)
            Vpobj = reader.TopOpt_ConstraintValue(dfdx_index);
            Vx=obj.Elements_volume/(Vpobj*obj.V_TOT);
            if isempty(obj.TOEL)
                fvalue=-1;
            else
                fvalue=(Vx'*mesh.elements_density(obj.TOEL)')-1;
            end
            obj.fval(dfdx_index)=fvalue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dfdx_Volume(obj,reader,dfdx_index)
            Vpobj = reader.TopOpt_ConstraintValue(dfdx_index);
            Vx=obj.Elements_volume/(Vpobj*obj.V_TOT);
            obj.dfdx(dfdx_index,1:length(obj.TOEL))=Vx';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fvalue=fval_Stress(obj,reader,mesh,solver,index_con)
            if strcmp(reader.TopOpt_ConstraintName{index_con},'Stress_Pnorm')
                obj.SVM=CalculateStressVM_MeshElements(reader,mesh,solver,reader.TopOpt_DesignElements);
                fvalue=(1/length(obj.SVM)*sum((obj.SVM./reader.TopOpt_ConstraintValue(index_con)).^reader.KSUp))^(1/reader.KSUp)-1;
                obj.fval(index_con)=fvalue;
            else
                obj.SVM=CalculateStressVM_MeshElements(reader,mesh,solver,reader.TopOpt_DesignElements);
                fvalue= KSU(obj.SVM/reader.TopOpt_ConstraintValue(index_con)-1,reader.KSUp); % KSU OK
                obj.fval(index_con)=fvalue;
            end

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dfdx_Stress(obj,reader,mesh,solver,index_con)
            
            element_sensitivities=zeros(length(obj.TOEL),1);

            mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);
            etype=mesh.data.ElementTypes{mesh_elements(1)};
            dim = mesh.retrieveelementdimension(etype); 
            totalelem = length(mesh_elements);

            %KJa=-solver.KT;
            %% Derivative of KSU
            Number_Of_M_Dofs = obj.Number_of_nodes*obj.dim;
            Number_Of_TV_Dofs = obj.Number_of_nodes*2;
            AdjU=zeros(Number_Of_M_Dofs,1);
            AdjTV=zeros(Number_Of_TV_Dofs,1);
            if strcmp(reader.TopOpt_ConstraintName{index_con},'Stress_Pnorm')
                [LdU,LdT,Luel] = obj.CalculateStressVM_Derivatives_MeshElements_Pnorm(reader,mesh,solver,index_con);
            else
                [sVM0,LdU,LdT,Luel] = obj.CalculateStressVM_Derivatives_MeshElements(reader,mesh,solver,index_con);
            end
            %% ADJOINT
            %%%% modificar systema a resolver para los adjuntos
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% U Adjoint        % Kuu^T lambda_u=-Ld
            A11 =distributed(solver.KStiff(obj.freedofs_mech,obj.freedofs_mech)');
            B11=distributed(-LdU(obj.freedofs_mech));
            AdjU(obj.freedofs_mech)=A11\B11;

            %% TV Adjoint       % [KJ]^T x lambda_TV = -[LdT,0]-(-K_UT)^T*lambda_u % check
            LdTV=zeros(obj.Number_of_nodes*2,1);
            odd_indices = 1:2:length(LdTV); % Get even indices
            LdTV(odd_indices)=-LdT+solver.KUT'*AdjU;
            A22 = distributed(-(solver.KT(obj.freedofs,obj.freedofs))'); % we store the negative KT
            B22 = distributed(LdTV(obj.freedofs));%% LT is in U T dofs, change into TV and get free ones!!!
            AdjTV(obj.freedofs)=A22\B22;

            %% Elemental wise
            GaussfunctionTag_KUU=@(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) obj.Integration_KuDerivative(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
            GaussfunctionTag_KUT=@(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) obj.Integration_KutDerivative(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
            GaussfunctionTag_Rx=@(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) obj.integration_R_dx(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
            %%% parfor
            parfor ii=1:length(obj.TOEL)
                ElementTag=obj.TOEL(ii);
                
                 element_nodes = mesh.data.ELEMENTS{ElementTag};
                 number_of_nodes=length(element_nodes);
                 orderdof_U=zeros(number_of_nodes*obj.dim,1);
                for nd=1:length(element_nodes)
                    orderdof_U(nd*obj.dim-(obj.dim-1))=element_nodes(nd)*obj.dim-(obj.dim-1);
                    orderdof_U(nd*obj.dim-(obj.dim-2))=element_nodes(nd)*obj.dim-(obj.dim-2);
                    if obj.dim==3
                    orderdof_U(nd*obj.dim)=element_nodes(nd)*obj.dim;
                    end
                end
%                 etype=mesh.data.ElementTypes{obj.TOEL(1)};
                dim = mesh.retrieveelementdimension(etype); 
                %%%% Derivatives
                [KUUd,flag,element_dofs_U]=obj.GaussIntegration_dx(dim,  reader.GI_order, ElementTag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{ElementTag},GaussfunctionTag_KUU) ;
                [KUTd,flag,element_dofs_T]=obj.GaussIntegration_dx(dim,  reader.GI_order, ElementTag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{ElementTag},GaussfunctionTag_KUT) ;
                [Rx,flag,element_dofs_Rx]=obj.GaussIntegration_dx(dim,  reader.GI_order, ElementTag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{ElementTag},GaussfunctionTag_Rx) ;
                Tee=solver.soldofs(element_dofs_T(1:number_of_nodes));
                Uee=solver.soldofs_mech(orderdof_U)';

                %Adji=AdjU(orderdofU(:,ii)');
                %index = mesh_elements == ElementTag;
                element_sensitivities(ii)= sum(Luel(ii,:))...
                    +AdjU(orderdof_U)'*(KUUd*Uee')...
                    -AdjU(orderdof_U)'*(KUTd*(Tee-str2double(reader.T0)))...
                    +AdjTV(element_dofs_Rx)'*Rx;%...

            end

            obj.dfdx(index_con,1:length(obj.TOEL))=element_sensitivities;

            % Calculate remaining bc sensititivies
            bcvariablenames=reader.TObctype;
            for i=1:length(bcvariablenames)
                if strcmp(bcvariablenames(i), 'Voltage')
                    bcvariable_nodes = mesh.retrieveNodalSelection(reader.TObcloc(i));
                    bcvariable_dofs=bcvariable_nodes*2;
                    dx_bc=zeros(obj.Number_of_dofs,1);
                    dx_bc(bcvariable_dofs)=reader.TObcmaxval(i)-reader.TObcminval(i);

                    obj.dfdx(index_con,length(obj.TOEL)+i)=AdjTV'*(-(solver.KT)*dx_bc);
                end
            end


        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [sigmaVM0_ordered,Ld,LdT,Luel] = CalculateStressVM_Derivatives_MeshElements(obj,reader,mesh,solver,index_con)

            PKS=reader.KSUp;
            Sobj=reader.TopOpt_ConstraintValue(index_con);
            mesh_elements = obj.TOEL;
            total_number_of_elements = length(mesh_elements);
            total_number_of_nodes = length(mesh.data.NODE);

            [node_el, etype_element] = mesh.retrievemeshtype(reader);
            sigmaVM0=zeros(total_number_of_nodes,1);

            initialdofs = solver.soldofs;
            daVM=zeros(6,1);
            evTi_summation=0;
            Ld=zeros(total_number_of_nodes*3,1);
            LdT=zeros(total_number_of_nodes*1,1);
            Luel=zeros(total_number_of_elements,2);

            for i = 1:total_number_of_elements

                % Recover each element tag
                elementTag = mesh_elements(i);
                element_nodes = mesh.data.ELEMENTS{elementTag};
                number_of_nodes = length(element_nodes);
                element_coordinates=zeros(3,number_of_nodes);
                dof_per_node=2;
                Tee=zeros(1,number_of_nodes);
                Uee=zeros(1,number_of_nodes*3);
                element_dof_indexes_TV=zeros(number_of_nodes*2,1);
                element_material_index=mesh.elements_material(elementTag);
                xx= mesh.elements_density(elementTag);
                Number_of_Nodes = length(element_coordinates(1,:));
                element_dof_indexes_M=zeros(number_of_nodes*3,1);

                for nd=1:number_of_nodes
                    element_dof_indexes_TV(nd)=element_nodes(nd)*2-1;
                    element_dof_indexes_TV(number_of_nodes+nd)=element_nodes(nd)*2;
                    element_dof_indexes_M((nd-1)*3+1)=(element_nodes(nd)-1)*3+1;
                    element_dof_indexes_M((nd-1)*3+2)=(element_nodes(nd)-1)*3+2;
                    element_dof_indexes_M((nd-1)*3+3)=(element_nodes(nd)-1)*3+3;
                    element_coordinates(:,nd)=mesh.data.NODE{element_nodes(nd)};
                    Tee(nd)=initialdofs(element_nodes(nd)*2-1);
                    Uee((nd-1)*3+1)=solver.soldofs_mech((element_nodes(nd)-1)*3+1);
                    Uee((nd-1)*3+2)=solver.soldofs_mech((element_nodes(nd)-1)*3+2);
                    Uee((nd-1)*3+3)=solver.soldofs_mech((element_nodes(nd)-1)*3+3);
                end

                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype_element, 0, 0, 0);

                JM = dShape' * element_coordinates';
                %Jacinv = inv(JM);
                DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
                % FIXME, calculate from all dofs input
                Th = N * Tee';

                Dalphapx = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_x');
                Dalphapy = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_y');
                Dalphapz = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_z');
                DEp = reader.getmaterialproperty(element_material_index,'YoungModulus');
                nu = reader.getmaterialproperty(element_material_index,'PoissonRatio');%%%check
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_YoungModulus'),reader.getmaterialproperty(element_material_index,'Tmax_YoungModulus')];
                [DE,DdE]=CalculateMaterialProperties(1e-6,DEp,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
                %[Dalpha,Ddalpha]=CalculateMaterialProperties(1e-6,Dalphap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalExpansionCoefficient'));
                [DE_dx]=CalculateMaterial_XDerivative(1e-6,DEp,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
                alphav=zeros(6,1);
                alphav(1:3,1)=[Dalphapx,Dalphapy,Dalphapz];

                C = DE/((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                    0 0 0 0 0 (1-2*nu)/2];

                dCx=DE_dx/((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                    0 0 0 0 0 (1-2*nu)/2];
                % Preallocate memory for B-Operators
                Bi=zeros(6,3);B=zeros(6,Number_of_Nodes*3);
                for bi=1:Number_of_Nodes
                    Bi(1,1) = DN(1,bi);
                    Bi(2,2) = DN(2,bi);
                    Bi(3,3) = DN(3,bi);
                    Bi(4,:) = [DN(2,bi),DN(1,bi),0];
                    Bi(5,:) = [0,DN(3,bi),DN(2,bi)];
                    Bi(6,:) = [DN(3,bi),0,DN(1,bi)];
                    B(:,(bi-1)*3+1:(bi)*3)=Bi;
                end

                eps=B*Uee'+alphav*N*(Tee-str2double(reader.T0))';
                se0=C*B*Uee'-C*alphav*N*(Tee-str2double(reader.T0))';
                %se0=C*eps   -C*alphav*N*(Tee-str2double(reader.T0))';
                sVM0=sqrt(se0(1)^2+se0(2)^2+se0(3)^2-se0(1)*se0(2)-se0(1)*se0(3)-se0(2)*se0(3)+3*se0(4)^2+3*se0(5)^2+3*se0(6)^2);
                %% Derivative of VM stress with each of the stress components
                daVM=[ 1/(2*sVM0)*(2*se0(1)-se0(2)-se0(3));...
                    1/(2*sVM0)*(2*se0(2)-se0(1)-se0(3));...
                    1/(2*sVM0)*(2*se0(3)-se0(1)-se0(2));...
                    3/(sVM0)*se0(4);...
                    3/(sVM0)*se0(5);...
                    3/(sVM0)*se0(6);];
                sigmaVM0(i) = sVM0;

                %%%%%%%%%%%%%%
                evTi=exp(PKS*(sVM0/Sobj-1));
                evTi_summation=evTi_summation+evTi;
                %%%%%%%%%%%%%%
                %% FIXME: for parfor the Li and Li2 is a solution, another way is storing all coefficients and adding after the parfor loop
                Li=zeros(number_of_nodes*3,3*total_number_of_nodes);
                for ll=1:number_of_nodes*3 % transformation matrix for parfor loop
                    Li(ll,element_dof_indexes_M(ll))=1;
                end
                aa=evTi*daVM'*C*B;
                Ld=Ld+(aa*Li)'; % L for mech dofs, which multiplies dU/dx

                Li2=zeros(number_of_nodes,total_number_of_nodes);
                for ll=1:number_of_nodes  % transformation matrix for parfor loop
                    Li2(ll,element_nodes(ll))=1;
                end
                aa=evTi*daVM'*(-C*alphav*N);
                LdT=LdT+(aa*Li2)'; % L for mech dofs, which multiplies d(\Delta \Theta)/dx = d(T-Tref)/dx

                Luel(i,:)=[ evTi*daVM'*dCx*B*Uee';...
                            evTi*daVM'*(-dCx*alphav*N)*(Tee-str2double(reader.T0))'];  % L for mech and T dofs, which multiplies [U;(T-Tref)]
            end

            sigmaVM0_ordered=zeros(total_number_of_elements,1);
            %Luel_ordered = zeros(total_number_of_elements,2);
            for i = 1:total_number_of_elements
                % Recover each element tag
                elementTag = mesh_elements(i);
                sigmaVM0_ordered(elementTag)=sigmaVM0(i);
                %Luel_ordered(elementTag,:)=Luel(i,:);
            end
            Ld=Ld/evTi_summation/Sobj;
            LdT=LdT/evTi_summation/Sobj;
            Luel=Luel/evTi_summation/Sobj;
            %Luel_ordered=Luel_ordered/evTi_summation/Sobj;            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [LdU,LdT,Luel] = CalculateStressVM_Derivatives_MeshElements_Pnorm(obj,reader,mesh,solver,index_con)
            mesh_elements=mesh.retrieveElementalSelection(reader.MeshEntityName);
            etype=mesh.data.ElementTypes{mesh_elements(1)};
            dim = mesh.retrieveelementdimension(etype); 
            
            Sobj=reader.TopOpt_ConstraintValue(index_con);
            mesh_elements = obj.TOEL;
            total_number_of_elements = length(mesh_elements);
            total_number_of_nodes = length(mesh.data.NODE);
            Summation_Pnorm=0;
            [node_el, etype_element] = mesh.retrievemeshtype(reader);

            LdU=zeros(total_number_of_nodes*dim,1);
            LdT=zeros(total_number_of_nodes*1,1);
            Luel=zeros(total_number_of_elements,2);

            %elementTag = mesh_elements(1);
            %element_nodes = mesh.data.ELEMENTS{elementTag};
            %number_of_nodes = length(element_nodes);

            %Ld_U=zeros(total_number_of_elements,number_of_nodes*3*number_of_nodes*3);
            %Ld_U_dofs=zeros(total_number_of_elements,number_of_nodes*3);

            %%% parfor
            parfor i = 1:total_number_of_elements

                % Recover each element tag
                elementTag = mesh_elements(i);
                element_nodes = mesh.data.ELEMENTS{elementTag};
                number_of_nodes = length(element_nodes);
                element_coordinates=zeros(3,number_of_nodes);
                element_dof_indexes_TV=zeros(number_of_nodes*2,1);
                element_dof_indexes_T=zeros(number_of_nodes,1);
                element_material_index=mesh.elements_material(elementTag);
                xx= mesh.elements_density(elementTag);
                Number_of_Nodes = length(element_coordinates(1,:));
                element_dof_indexes_M=zeros(number_of_nodes*dim,1);

                
                for nd=1:number_of_nodes
                    element_dof_indexes_TV(nd)=element_nodes(nd)*2-1;
                    element_dof_indexes_TV(number_of_nodes+nd)=element_nodes(nd)*2;
                    element_dof_indexes_T(nd)=element_nodes(nd);
                    element_dof_indexes_M((nd-1)*dim+1)=(element_nodes(nd)-1)*dim+1;
                    element_dof_indexes_M((nd-1)*dim+2)=(element_nodes(nd)-1)*dim+2;
                    if dim==3
                    element_dof_indexes_M((nd-1)*dim+3)=(element_nodes(nd)-1)*dim+3;
                    end
                    element_coordinates(:,nd)=mesh.data.NODE{element_nodes(nd)};
                end
                Uee=solver.soldofs_mech(element_dof_indexes_M)';
                Tee= solver.soldofs(element_dof_indexes_TV(1:number_of_nodes))';
                if dim==2
                    [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, 0,0,0);
                    element_coordinates=element_coordinates(1:2,:);
                    N=N';
                else
                    [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, 0,0,0);
                end

                JM = dShape' * element_coordinates';
                DN = inv(JM) * dShape'; 
                Th = N * Tee';
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_YoungModulus'),reader.getmaterialproperty(element_material_index,'Tmax_YoungModulus')];
                Dalphapx = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_x');
                [Dax,Daxdt]=CalculateMaterialProperties(1e-6,Dalphapx,Tmat,1,1);
                Dalphapy = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_y');
                [Day,Daydt]=CalculateMaterialProperties(1e-6,Dalphapy,Tmat,1,1);
                Dalphapz = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_z');
                [Daz,Dazdt]=CalculateMaterialProperties(1e-6,Dalphapz,Tmat,1,1);
                %[Dalpha,Ddalpha]=CalculateMaterialProperties(1e-6,Dalphap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalExpansionCoefficient'));
                [DaT_dx]=CalculateMaterial_XDerivative(1e-6,Dalphapx,Tmat,1,1);

                DEp = reader.getmaterialproperty(element_material_index,'YoungModulus');
                nu = reader.getmaterialproperty(element_material_index,'PoissonRatio');

                [DE,DdE]=CalculateMaterialProperties(1e-6,DEp,Tmat,xx,...
                    reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
                %[Dalpha,Ddalpha]=CalculateMaterialProperties(1e-6,Dalphap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalExpansionCoefficient'));
                [DE_dx]=CalculateMaterial_XDerivative(1e-6,DEp,Tmat,xx,...
                    reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));

                
                alphav=zeros(6,1);
                alphav(1:3,1)=[Dax,Day,Daz];
                alphavdt=zeros(6,1);
                alphavdt(1:3,1)=[Daxdt,Daydt,Dazdt];


                if dim==2
               alphav=zeros(3,1);
               alphav(1:2,1)=[Dax,Day];
               alphavdt=zeros(3,1);
                alphavdt(1:2,1)=[Daxdt,Daydt];
                alphav_dx=zeros(3,1);
                alphav_dx(1:2,1)=[DaT_dx,DaT_dx];  


                C = DE / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2]; % plane stress
                % C = DE / (1 + nu) / (1 - 2 * nu) * [1 - nu, nu, 0; nu, 1
                % - nu, 0; 0, 0, (1 - 2 * nu) / 2]; % plane strain
                dCx = DE_dx / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2]; % plane stress

                Bi = zeros(3, 2);
                B = zeros(3, Number_of_Nodes * 2);
            
                for bi = 1:Number_of_Nodes
                    Bi(1, 1) = DN(1, bi);
                    Bi(2, 2) = DN(2, bi);
                    Bi(3, 1) = DN(2, bi);
                    Bi(3, 2) = DN(1, bi);
            
                    B(:, (bi - 1) * 2 + 1 : bi * 2) = Bi;
                end                
            elseif dim==3
                alphav=zeros(6,1);
                alphav(1:3,1)=[Dax,Day,Daz];  
               alphavdt=zeros(6,1);
                alphavdt(1:3,1)=[Daxdt,Daydt,Dazdt];
                alphav_dx=zeros(6,1);
                alphav_dx(1:3,1)=[DaT_dx,DaT_dx,DaT_dx];  

                C = DE/((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                    0 0 0 0 0 (1-2*nu)/2];

                dCx=DE_dx/((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                    0 0 0 0 0 (1-2*nu)/2];

                Bi=zeros(6,3);B=zeros(6,Number_of_Nodes*3);
                    for bi=1:Number_of_Nodes
                        Bi(1,1) = DN(1,bi);
                        Bi(2,2) = DN(2,bi); 
                        Bi(3,3) = DN(3,bi);
                        Bi(4,:) = [DN(2,bi),DN(1,bi),0];
                        Bi(5,:) = [0,DN(3,bi),DN(2,bi)]; 
                        Bi(6,:) = [DN(3,bi),0,DN(1,bi)];
                        B(:,(bi-1)*3+1:(bi)*3)=Bi;
        
                    end
            end
                
            % Compute stress vector (mechanical + thermal)
            se0 = C * B * Uee' - C * alphav * N * (Tee - str2double(reader.T0))';
            
            if dim == 2
                % Voigt notation for 2D: [σ11, σ22, σ12]
                sVM0 = sqrt( ...
                    se0(1)^2 + se0(2)^2 - se0(1)*se0(2) + ...
                    3 * se0(3)^2 ...
                );
            
                Summation_Pnorm = Summation_Pnorm + sVM0^reader.KSUp;
                daVM=[  1/(2*sVM0)*(2*se0(1)-se0(2));...
                        1/(2*sVM0)*(2*se0(2)-se0(1));...
                        3/(sVM0)*se0(3)];            
            elseif dim == 3
                % Voigt notation for 3D: [σ11, σ22, σ33, σ23, σ13, σ12]
                sVM0 = sqrt( ...
                    0.5 * ( ...
                        (se0(1) - se0(2))^2 + ...
                        (se0(2) - se0(3))^2 + ...
                        (se0(3) - se0(1))^2 ) + ...
                    3 * (se0(4)^2 + se0(5)^2 + se0(6)^2) ...
                );
            
                Summation_Pnorm = Summation_Pnorm + sVM0^reader.KSUp;
            
                daVM=[  1/(2*sVM0)*(2*se0(1)-se0(2)-se0(3));...
                        1/(2*sVM0)*(2*se0(2)-se0(1)-se0(3));...
                        1/(2*sVM0)*(2*se0(3)-se0(1)-se0(2));...
                        3/(sVM0)*se0(4);...
                        3/(sVM0)*se0(5);...
                        3/(sVM0)*se0(6);];
            end

% 
%                 se0=C*B*Uee'-C*alphav*N*(Tee-str2double(reader.T0))';
%                 if dim ==2
%                 sVM0=sqrt(se0(1)^2+se0(2)^2 ...
%                     -se0(1)*se0(2) ...
%                     +3*se0(3)^2);
%                 Summation_Pnorm=Summation_Pnorm+(sVM0)^reader.KSUp;
%                 %% Derivative of VM stress with each of the stress components
%                 daVM=[  1/(2*sVM0)*(2*se0(1)-se0(2));...
%                         1/(2*sVM0)*(2*se0(2)-se0(1));...
%                         3/(sVM0)*se0(3)];
% 
%                 elseif dim==3
%                 sVM0=sqrt(se0(1)^2+se0(2)^2+se0(3)^2 ...
%                     -se0(1)*se0(2)-se0(1)*se0(3)-se0(2)*se0(3) ...
%                     +3*se0(4)^2+3*se0(5)^2+3*se0(6)^2);
%                 Summation_Pnorm=Summation_Pnorm+(sVM0)^reader.KSUp;
%                 %% Derivative of VM stress with each of the stress components
%                 daVM=[  1/(2*sVM0)*(2*se0(1)-se0(2)-se0(3));...
%                         1/(2*sVM0)*(2*se0(2)-se0(1)-se0(3));...
%                         1/(2*sVM0)*(2*se0(3)-se0(1)-se0(2));...
%                         3/(sVM0)*se0(4);...
%                         3/(sVM0)*se0(5);...
%                         3/(sVM0)*se0(6);];
%                 end
                %%%%%%%%%%%%%%
                %% FIXME: for parfor the Li and Li2 is a solution, another way is storing all coefficients and adding after the parfor loop
                %Ld_U(i,:)=sVM0^(reader.KSUp-1)*daVM'*C*B; % L for mech dofs, which multiplies dU/dx
                element_multiplier=sVM0^(reader.KSUp-1)*daVM';
                %%
                Li_M=zeros(number_of_nodes*dim,dim*total_number_of_nodes);
                for ll=1:length(element_dof_indexes_M) % transformation matrix for parfor loop
                    Li_M(ll,element_dof_indexes_M(ll))=1;
                end

                Li_T=zeros(number_of_nodes,total_number_of_nodes);
                for ll=1:number_of_nodes  % transformation matrix for parfor loop
                    Li_T(ll,element_nodes(ll))=1;
                end
                %%
                LdU=LdU+(element_multiplier*C*B*Li_M)';
                LdT=LdT-(element_multiplier*(C*alphav*N+(C*alphavdt)*(N*(Tee-str2double(reader.T0))')*N)*Li_T)'; % derivative of CTE with temp %% MODIFICATION!!!
                Luel(i,:)=[element_multiplier*dCx*B*Uee',element_multiplier*((-dCx*alphav*N)*(Tee-str2double(reader.T0))')];% -(C*alphav_dx*N)*(Tee-str2double(reader.T0))'  % L for mech and T dofs, which multiplies [U;(T-Tref)]
            end
            %Pnorm=(Summation_Pnorm/total_number_of_elements)^(1/reader.KSUp);
            dPnorm_common=(Summation_Pnorm/total_number_of_elements)^(1/reader.KSUp-1)*1/total_number_of_elements/Sobj;
            Luel=Luel*dPnorm_common;
            LdU=LdU*dPnorm_common;
            LdT=LdT*dPnorm_common;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function[F_T,flag]=Integration_KutDerivative(~,natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)
                flag=[];
                dim = mesh.retrieveelementdimension(etype); 
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
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_ElectricalConductivity'),reader.getmaterialproperty(element_material_index,'Tmax_ElectricalConductivity')];
                Dalpha_x = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_x');
                [Dax,Daxdt]=CalculateMaterialProperties(1e-6,Dalpha_x,Tmat,1,1);
                Dalpha_y = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_y');
                [Day,Daydt]=CalculateMaterialProperties(1e-6,Dalpha_y,Tmat,1,1);
                Dalpha_z = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_z');
                [Daz,Dazdt]=CalculateMaterialProperties(1e-6,Dalpha_z,Tmat,1,1);
                DEp = reader.getmaterialproperty(element_material_index,'YoungModulus');
                nu = reader.getmaterialproperty(element_material_index,'PoissonRatio');
                %[DE,DdE]=CalculateMaterialProperties(1e-6,DEp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
                %[Dalpha,Ddalpha]=CalculateMaterialProperties(1e-6,Dalphap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalExpansionCoefficient'));
                [DE_dx]=CalculateMaterial_XDerivative(1e-6,DEp,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
                alphav=zeros(6,1);
                alphav(1:3,1)=[Dax,Day,Daz];

                %C = DE./((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                %    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                %    0 0 0 0 0 (1-2*nu)/2];
                dCx = DE_dx./((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                    0 0 0 0 0 (1-2*nu)/2];                
                detJ = det(JM);
            

                if dim==2
               alphav=zeros(3,1);
               alphav(1:2,1)=[Dax,Day];
              alphavdt=zeros(3,1);
                alphavdt(1:2,1)=[Daxdt,Daydt];
                % C = DE / (1 + nu) / (1 - 2 * nu) * [1 - nu, nu, 0; nu, 1
                % - nu, 0; 0, 0, (1 - 2 * nu) / 2]; % plane strain
                dCx = DE_dx / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2]; % plane stress

                Bi = zeros(3, 2);
                B = zeros(3, Number_of_Nodes * 2);
            
                for bi = 1:Number_of_Nodes
                    Bi(1, 1) = DN(1, bi);
                    Bi(2, 2) = DN(2, bi);
                    Bi(3, 1) = DN(2, bi);
                    Bi(3, 2) = DN(1, bi);
            
                    B(:, (bi - 1) * 2 + 1 : bi * 2) = Bi;
                end                
            elseif dim==3
                alphav=zeros(6,1);
                alphav(1:3,1)=[Dax,Day,Daz];  
               alphavdt=zeros(6,1);
                alphavdt(1:3,1)=[Daxdt,Daydt,Dazdt];


                dCx=DE_dx/((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                    0 0 0 0 0 (1-2*nu)/2];

                Bi=zeros(6,3);B=zeros(6,Number_of_Nodes*3);
                    for bi=1:Number_of_Nodes
                        Bi(1,1) = DN(1,bi);
                        Bi(2,2) = DN(2,bi); 
                        Bi(3,3) = DN(3,bi);
                        Bi(4,:) = [DN(2,bi),DN(1,bi),0];
                        Bi(5,:) = [0,DN(3,bi),DN(2,bi)]; 
                        Bi(6,:) = [DN(3,bi),0,DN(1,bi)];
                        B(:,(bi-1)*3+1:(bi)*3)=Bi;
        
                    end
            end
                     %         B'*(dCx*alphav)*N'
                     %F_T=detJ*(B'*(dCx*alphav)*N+B'*(dCx*alphavdt)*N);
                     F_T=detJ*(B'*(dCx*alphav)*N);
        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function[KUU,flag]=Integration_KuDerivative(~,natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)
                flag=[];
                dim = mesh.retrieveelementdimension(etype); 
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

                DEp = reader.getmaterialproperty(element_material_index,'YoungModulus');
                nu = reader.getmaterialproperty(element_material_index,'PoissonRatio');

                %[DE,DdE]=CalculateMaterialProperties(1e-6,DEp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_YoungModulus'),reader.getmaterialproperty(element_material_index,'Tmax_YoungModulus')];
                [DE_dx]=CalculateMaterial_XDerivative(1e-6,DEp,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));

              
                detJ = det(JM);
            

                if dim==2


                % C = DE / (1 + nu) / (1 - 2 * nu) * [1 - nu, nu, 0; nu, 1
                % - nu, 0; 0, 0, (1 - 2 * nu) / 2]; % plane strain
                dCx = DE_dx / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2]; % plane stress

                Bi = zeros(3, 2);
                B = zeros(3, Number_of_Nodes * 2);
            
                for bi = 1:Number_of_Nodes
                    Bi(1, 1) = DN(1, bi);
                    Bi(2, 2) = DN(2, bi);
                    Bi(3, 1) = DN(2, bi);
                    Bi(3, 2) = DN(1, bi);
            
                    B(:, (bi - 1) * 2 + 1 : bi * 2) = Bi;
                end                
            elseif dim==3
 



                dCx=DE_dx/((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                    0 0 0 0 0 (1-2*nu)/2];

                Bi=zeros(6,3);B=zeros(6,Number_of_Nodes*3);
                    for bi=1:Number_of_Nodes
                        Bi(1,1) = DN(1,bi);
                        Bi(2,2) = DN(2,bi); 
                        Bi(3,3) = DN(3,bi);
                        Bi(4,:) = [DN(2,bi),DN(1,bi),0];
                        Bi(5,:) = [0,DN(3,bi),DN(2,bi)]; 
                        Bi(6,:) = [DN(3,bi),0,DN(1,bi)];
                        B(:,(bi-1)*3+1:(bi)*3)=Bi;
        
                    end
            end
                        

                KUU=detJ*(B'*dCx*B);
        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
end

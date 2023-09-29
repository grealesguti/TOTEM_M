classdef TO_Constraints
    %TO_OBJECTIVES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ConstraintNames
        MeshName
        fval
        dfdx
    end
    
    methods
        function obj = TO_Constraints(reader)
            %TO_OBJECTIVES Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function CalculateConstraint(obj,mesh,solver)
            for i=1:length(obj.ConstraintNames)
                ConstraintName=obj.ConstraintNames{i};
                switch ConstraintName
                    case 'Power'
                        obj.fval_AverageTemp(mesh,solver,i)
                        obj.dfdx_AverageTemp(mesh,solver,i)
                    case 'Stress_KSU'
                        obj.fval_AverageTemp(mesh,solver)
                        obj.dfdx_AverageTemp(mesh,solver)
                    case 'Volume'
                        obj.fval_AverageTemp(mesh,solver)
                        obj.dfdx_AverageTemp(mesh,solver)
                    case 'Displacement'
                        obj.fval_AverageTemp(mesh,solver)
                        obj.dfdx_AverageTemp(mesh,solver)
                end
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

            integrationFunction = @(natcoords) obj.integration_AvgTemp(natcoords, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype, mesh.elements_density(elementTag));      
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
        function fval_Power(obj,mesh,solver)
            f0val=sum(U((Objnd-1)*2+1))/length(U((Objnd-1)*2+1));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dfdx_Power(obj,mesh,solver)
            ST=zeros(length(TOel)+1,1);
            KJa=-solver.KT;
            AdjT=zeros(length(LO),1);
            T2v=zeros(length(LO),1);
            T2v((Objnd-1)*2+1)=1/length(T2v((Objnd-1)*2+1));
            LAdj=T2v;
            A = distributed((-KJa(freedof,freedof)-K_conv(freedof,freedof))'); 
            B=distributed((LAdj(freedof)));
            AdjT(freedof)=A\B;
            sf=zeros(length(TOel),1);
            for ii=1:length(TOel)
                 jj=TOel(ii);
                 doforderT=(order(jj,:)-1)*2+1;doforderV=(order(jj,:)-1)*2+2;
                 doforder=[doforderT doforderV];

                Te=U(doforderT);Ve=U(doforderV);
                [Rx]=KderTHEL(coord,order(jj,:),matp,matv,p,xx(jj),1e-9,localsys,sysv,jj,Te,Ve,seebp,rhop);
                ADJi=AdjT(doforder)';%Ui=U(doforder,1);%Ui(1:20)=Ui(1:20);
                sf(ii)=ADJi*Rx;
            
            end

                ST(1:end-1,1)=sf;
                doforderV=(Vfnod-1)*2+2;
                dU1=zeros(length(AdjT),1);
                dU1(doforderV)=dV;
                prodF=KJa*dU1;
                ST(end,1)=LAdj'*dU1+AdjT'*(+prodF);
                ti=toc(ticTxiVf)-tocTxiVf;fprintf("End T vs xi,Vf sensitivity: %f\n",ti)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function[Rx]=integration_Power_dx(~,natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)        
            cooro=coord(order,:);
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));

                JM = dShape' * element_coordinates';

                %Jacinv = inv(JM);
                DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
                % FIXME, calculate from all dofs input
                Th = N * Tee';

                DN = (Jacinv*dShape);

                Dep = reader.getmaterialproperty(element_material_index,'ElectricalConductivity');
                Dkp = reader.getmaterialproperty(element_material_index,'ThermalConductivity');
                Dap = reader.getmaterialproperty(element_material_index,'Seebeck');
                
                [De,Dde]=CalculateMaterialProperties(Dep,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
                [Da,Dda]=CalculateMaterialProperties(Dap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
                [Dk,Ddk]=CalculateMaterialProperties(Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));

                [De_dx,Dde_dx]=CalculateMaterial_XDerivative(Dep,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
                [Da_dx,Dda_dx]=CalculateMaterial_XDerivative(Dap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
                [Dk_dx,Ddk_dx]=CalculateMaterial_XDerivative(Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));
                
                Vee=Vee';
                Tee=Tee';
                detJ = det(JM);
                
                % Calculate current density and heat flux
                je = -De * DN * Vee - Da * De * DN * Tee;
                % all matrix
                jx=-De_dx*DN*Vee-Da_dx*De*DN*Tee-Da*Ded*DN*Tee;
                qx=Da_dx*(N'*Tee)*je+Da*(N'*Tee)*jx-Dk_dx*DN*Tee;
        
                RAx=detJ*(-DN'*qx+N*(jx'*DN*Vee));
                RBx=detJ*(-DN'*jx);

                Rx=[RAx
                   RBx];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fval_Volume(obj,mesh,solver)
            Vx=VTOel/(Vpobj*sum(VTOel));
            fval(i)=(Vx'*xx(TOEL))-1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dfdx_Volume(obj,mesh,solver)
            Vx=VTOel/(Vpobj*sum(VTOel));
            dfdx(i,:)=[Vx'];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end

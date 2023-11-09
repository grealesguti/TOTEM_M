classdef Filtering < handle
    %TO Summary of this class goes here
    %   Detailed explanation goes here

    properties
        nele
        nele_TO
        nd
        rd
        nnv
        nnod
        Newman
        TOdofs
        Kf
        L
        Tt
        rho0
        Hn1
        TOEL
        nele_total
        m
    end

    methods
        function obj = Filtering(reader,mesh)
            obj.TOEL = mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            r=mesh.Element_size*3;
            obj.rd=( r/2/sqrt(3) )^2;
            obj.nele=length( mesh.retrieveElementalSelection(reader.MeshEntityName));
            obj.nele_total =length(mesh.data.ELEMENTS);
            obj.nele_TO=length(obj.TOEL);
            if mesh.dim==2
                obj.nnv=4*4;
                obj.nd=4;
            elseif mesh.dim==3
                obj.nnv=8*8;
                obj.nd=8;
            end
            obj.nnod=length(mesh.data.NODE);

            Hnf=zeros(length(mesh.data.NODE),1);
            for jj=1:length(obj.TOEL)
                nodes_element=mesh.data.ELEMENTS(obj.TOEL(jj));
                nodes=nodes_element{1};
                Hnf(nodes,1)=1;
            end
            obj.TOdofs=find(Hnf==1);
            obj.Hn1=find(Hnf==0);
            obj.Newman=0;
            obj.m=length(reader.TopOpt_ConstraintName);
            obj.Helmholtz_PDE_Initialization(reader,mesh)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Helmholtz_PDE_Initialization(obj,reader,mesh)
            KJvnz=zeros(obj.nnv,obj.nele_TO);
            KJvc=zeros(obj.nnv,obj.nele_TO);
            KJvr=zeros(obj.nnv,obj.nele_TO);
            Tvnz=zeros(obj.nd,obj.nele_TO);
            Tvc=zeros(obj.nd,obj.nele_TO);
            Tvr=zeros(obj.nd,obj.nele_TO);
            
            % 14 point integration for element of any shape
            %GeneralGaussIntegration(dimension, order, elementTag, mesh, initialdofs,reader,etype, integrationFunctionhandle)
            %Integration_InitHelmholtz_1(   natural_coordinates, element_coordinates, mesh, etype,  rd)                                (natural_coordinates, element_coordinates, mesh, etype,  rd) 
            GaussfunctionTag        =@(     natural_coordinates, element_coordinates, mesh, etype,  rd) obj.Integration_InitHelmholtz_1(natural_coordinates, element_coordinates, mesh, etype, obj.rd);
            mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);

            % find elements in mesh not in TOEL
            mesh_element_index=zeros(obj.nele_total,1);
            mesh_element_index(mesh_elements)=1;
            mesh_TOEL_index=zeros(obj.nele_total,1);
            mesh_TOEL_index(obj.TOEL)=1;
            no_TOEL_index=mesh_element_index-mesh_TOEL_index;
            noTOEL=find(no_TOEL_index==1); 

            for jj=1:length(obj.TOEL)
                elementTag=obj.TOEL(jj);
                [Ke,Te,element_dof_indexes]= obj.GaussIntegration_H(mesh.dim,  reader.GI_order, elementTag, mesh, mesh.data.ElementTypes{elementTag},GaussfunctionTag) ;
                KJvnz(:,jj)=reshape(Ke,obj.nnv,1);
                KJvc(:,jj)=reshape(repmat(element_dof_indexes',obj.nd,1),obj.nnv,1);
                KJvr(:,jj)=reshape(repmat(element_dof_indexes,obj.nd,1)',obj.nnv,1);
                Tvnz(:,jj)=Te';
                Tvc(:,jj)=elementTag*ones(1,obj.nd);
                Tvr(:,jj)=element_dof_indexes';
            end
            
            KJvnzr=reshape(KJvnz,obj.nele_TO*obj.nnv,1);
            KJvcr=reshape(KJvc,obj.nele_TO*obj.nnv,1);
            KJvrr=reshape(KJvr,obj.nele_TO*obj.nnv,1); % is this needed??
            Kfs3=sparse(KJvrr,KJvcr,KJvnzr,obj.nnod,obj.nnod);
            Kfs31=Kfs3(obj.TOdofs,obj.TOdofs);
            
            Tvnzr=reshape(Tvnz,obj.nele_TO*obj.nd,1);
            Tvcr=reshape(Tvc,obj.nele_TO*obj.nd,1);
            Tvrr=reshape(Tvr,obj.nele_TO*obj.nd,1); % is this needed??
            T=sparse(Tvrr,Tvcr,Tvnzr,obj.nnod,obj.nele_total);

            obj.Tt=T(obj.TOdofs,obj.TOEL);
            obj.Kf=(Kfs31);
            obj.L = distributed(ichol(Kfs31,struct('michol','on')));
            KEF=Kfs3(obj.TOdofs,obj.Hn1);
            TEt=T(obj.TOdofs,noTOEL);
            obj.rho0=+Kfs31\(TEt*ones(length(noTOEL),1)-KEF*ones(length(obj.Hn1),1));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Ke,Te] = Integration_InitHelmholtz_1(obj,natural_coordinates, element_coordinates, mesh, etype, rd)
                        [dim] = mesh.retrieveelementdimension(etype);
                        if dim==2
                            [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), 0);
                            element_coordinates=element_coordinates(1:2,:);
                            N=N';
                        else
                            [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));
                        end        
                        JM = dShape' * element_coordinates';
                        detJ=det(JM);
                        %Jacinv = inv(JM);
                        DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
                        % FIXME, calculate from all dofs input
                        %[De]=CalculateMaterialProperties(Dep,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
                        %[Da]=CalculateMaterialProperties(Dap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
                        Ke1=detJ*(rd*(DN'*DN));
                        Ke2=detJ*(N'*N);
                        Ke=Ke1+Ke2;
                        Te=detJ*(N);        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [K,R,element_dof_indexes] = GaussIntegration_H(obj,dimension, order, elementTag, mesh, etype,integrationFunctionhandle)
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
            element_dof_indexes=element_nodes;
            for i=1:number_of_nodes
                element_coordinates(:,i)=mesh.data.NODE{element_nodes(i)};
            end
            %Integration_InitHelmholtz_1                                (natural_coordinates, element_coordinates, mesh, etype, rd)
            integrationFunction = @(natcoords) integrationFunctionhandle(natcoords          , element_coordinates, mesh, etype, obj.rd);
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
        function xx=Helmholtz_xx(obj,mesh,xx0)
            B=obj.Tt*xx0';
            xx=zeros(length(xx0),1);
            [rhofnd0,fl,rr,it] = pcg(obj.Kf,B,1e-9,100,obj.L,obj.L');
            rhoall=ones(obj.nnod,1);
            rhoall(obj.TOdofs)=rhofnd0;%+obj.rho0;
            for jj=1:length(obj.TOEL)
                ii=obj.TOEL(jj);
                nodes_element=mesh.data.ELEMENTS{ii};
                rhoallii=rhoall(nodes_element);
                if mesh.dim==2
                    [N, dShape] = mesh.selectShapeFunctionsAndDerivatives("CPS4", 0, 0, 0);
                elseif mesh.dim==3
                    [N, dShape] = mesh.selectShapeFunctionsAndDerivatives("C3D8", 0, 0, 0);
                end
                xx(jj)=N'*rhoallii;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function filter_densities(obj,reader,mesh)
            if reader.Filter_Helmholtz_Density==1 
                        xx=mesh.elements_density(obj.TOEL);
                        mesh.elements_density(obj.TOEL)=obj.Helmholtz_xx(mesh,xx);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function filter_sensitivities(obj,reader,mesh,TOO,TOC)
            if reader.Filter_Helmholtz_Density==1  
                 xx=TOO.dfdx(1:length(obj.TOEL))';
                 TOO.dfdx(1:length(obj.TOEL))=obj.Helmholtz_xx(mesh,xx);
                 for j=1:obj.m
                        xx=TOC.dfdx(j,1:length(obj.TOEL));
                        TOC.dfdx(j,1:length(obj.TOEL))=obj.Helmholtz_xx(mesh,xx);
                 end
            end
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end
end


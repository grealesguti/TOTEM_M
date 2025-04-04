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
        Lin_number_of_nodes
        etype
        mu
        beta
        xe_range
        post
        it
        folderName
        r
        el_multiplier
        Element_size
    end

    methods
        function obj = Filtering(reader,mesh)
            obj.TOEL = mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            obj.Element_size= mesh.Element_size;
            obj.el_multiplier = 3;
            obj.r=max([obj.Element_size*obj.el_multiplier,0.00011]);
            obj.rd=( obj.r/2/sqrt(3) )^2;
            obj.nele=length( mesh.retrieveElementalSelection(reader.MeshEntityName));
            obj.nele_total =length(mesh.data.ELEMENTS);
            obj.nele_TO=length(obj.TOEL);
            if mesh.dim==2
                obj.etype="CPS4";
            elseif mesh.dim==3
                obj.etype="C3D8";
            end
            [obj.nd] = mesh.retrieveelementnumberofnodes(obj.etype);
            obj.nnv=obj.nd^2;
            obj.nnod=length(mesh.data.NODE);

            Hnf=zeros(length(mesh.data.NODE),1);
            obj.Lin_number_of_nodes=obj.retrieveelementnumberofnodes(obj.etype);


            for jj=1:length(obj.TOEL) % if quadratic elements only the first linear ones!!!, depends on element type
                nodes_element=mesh.data.ELEMENTS(obj.TOEL(jj));
                nodes=nodes_element{1};
                nodes=nodes(1:obj.Lin_number_of_nodes);
                Hnf(nodes,1)=1;
            end
            obj.TOdofs=find(Hnf==1);
            obj.Hn1=find(Hnf==0);
            obj.Newman=0;
            obj.m=length(reader.TopOpt_ConstraintName);
            obj.Helmholtz_PDE_Initialization(reader,mesh)
            obj.mu=0.5;
            obj.beta=30;
            obj.xe_range=zeros(length(obj.TOEL)-1,1);
            obj.post = Postprocessing();
            obj.it=0;
            obj.folderName='';


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

            parfor jj=1:length(obj.TOEL)
                elementTag=obj.TOEL(jj);
                [Ke,Te,element_dof_indexes]= obj.GaussIntegration_H(mesh.dim,  reader.GI_order, elementTag, mesh, obj.etype,GaussfunctionTag) ;
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
            obj.Kf=distributed(Kfs31);
            %obj.L = distributed(ichol(Kfs31,struct('michol','on')));
            obj.L = distributed(chol(Kfs31, 'lower'));

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
            number_of_nodes = obj.Lin_number_of_nodes;
            element_coordinates=zeros(3,number_of_nodes);
            element_dof_indexes=element_nodes(1:number_of_nodes);
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
                rhoallii=rhoall(nodes_element(1:obj.Lin_number_of_nodes));
                if mesh.dim==2
                    [N, dShape] = mesh.selectShapeFunctionsAndDerivatives("CPS4", 0, 0, 0);
                elseif mesh.dim==3
                    [N, dShape] = mesh.selectShapeFunctionsAndDerivatives("C3D8", 0, 0, 0);
                    N=N';
                end
                xx(jj)=N'*rhoallii;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function xx=Helmholtz_dxx(obj,mesh,xx0)
            B=obj.Tt*xx0';
            xx=zeros(length(xx0),1);
            [rhofnd0,fl,rr,it] = pcg(obj.Kf,B,1e-9,100,obj.L,obj.L');
            rhoall=ones(obj.nnod,1);
            rhoall(obj.TOdofs)=rhofnd0;%+obj.rho0;
            for jj=1:length(obj.TOEL)
                ii=obj.TOEL(jj);
                nodes_element=mesh.data.ELEMENTS{ii};
                rhoallii=rhoall(nodes_element(1:obj.Lin_number_of_nodes));
                if mesh.dim==2
                    [N, dShape] = mesh.selectShapeFunctionsAndDerivatives("CPS4", 0, 0, 0);
                elseif mesh.dim==3
                    [N, dShape] = mesh.selectShapeFunctionsAndDerivatives("C3D8", 0, 0, 0);
                    N=N';
                end
                xx(jj)=N'*rhoallii;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function filter_densities(obj,reader,mesh)
            if reader.Filter>0
            % Helmholtz
                        xx=mesh.elements_density(obj.TOEL);
                        mesh.elements_density(obj.TOEL)=obj.Helmholtz_xx(mesh,xx);
            end

            obj.post.initVTK(reader,mesh);
            if reader.Filter==2 ||reader.Filter==-2
            %        obj.post.VTK_Mesh_xx(append([obj.folderName,'\','xx_Helmholtz_','_',num2str(obj.it)]),mesh)
            % Heaviside
                        for i =1:length(obj.TOEL)
                            el_index=obj.TOEL(i);
                            xe= mesh.elements_density(el_index);
                            obj.xe_range(i)=xe;
                            new_xe=( tanh(obj.beta*obj.mu)+tanh(obj.beta*(obj.xe_range(i)-obj.mu)) )/( tanh(obj.beta*obj.mu)+tanh(obj.beta*(1-obj.mu)) );
                            mesh.elements_density(el_index)=new_xe;
                        end
                    %obj.post.VTK_Mesh_xx(append([obj.folderName,'\','xx_Heaviside_','_',num2str(obj.it)]),mesh)
            end

            if reader.Filter==3 ||reader.Filter==-3
            %        obj.post.VTK_Mesh_xx(append([obj.folderName,'\','xx_Helmholtz_','_',num2str(obj.it)]),mesh)
            % Heaviside
                        for i =1:length(obj.TOEL)
                            el_index=obj.TOEL(i);
                            xe= mesh.elements_density(el_index);
                            if xe<obj.mu
                                obj.xe_range(i)=xe;
                                mesh.elements_density(el_index)=obj.mu*(exp(-obj.beta*(1-xe/obj.mu))-(1-xe/obj.mu)*exp(-obj.beta));
                            else
                                obj.xe_range(i)=xe;
                                mesh.elements_density(el_index)=(1-obj.mu)*(1-exp(-obj.beta*(xe-obj.mu)/(1*-obj.mu))+(xe-obj.mu)*exp(-obj.beta)/(1-obj.mu))+obj.mu;
                            end
                        end
                    %obj.post.VTK_Mesh_xx(append([obj.folderName,'\','xx_Heaviside_','_',num2str(obj.it)]),mesh)
            end

            if reader.Filter == 4 || reader.Filter == -4
                % Initialize parameters
                N_e = obj.TOEL;   % The indices of the elements in the topology optimization
                
                % Iterate through the topology optimization element list
                for i = 1:length(N_e)
                    % Get the current element index
                    el_index = N_e(i);
                    
                    % Get the current density values for the element
                    rho = mesh.elements_density(N_e);  % Extract density values for all elements in N_e
                    
                    % Apply the new expression to compute rho_e_tilde
                    numerator = sum(exp(obj.beta * rho));         % Summation of exponential terms
                    denominator = sum(ones(size(N_e)));      % Summation of ones (count of elements)
                    rho_e_tilde = log(numerator / denominator) / obj.beta;  % Compute the final value
                    
                    % Update the element density
                    mesh.elements_density(el_index) = rho_e_tilde;
                end
                
                % Optional: Write VTK mesh data for visualization (uncomment if needed)
                % obj.post.VTK_Mesh_xx(append([obj.folderName, '\', 'xx_Filter4_', '_', num2str(obj.it)]), mesh);
            end

            if reader.Filter == 5 || reader.Filter == -5
                % Heav dilation
                        for i =1:length(obj.TOEL)
                            el_index=obj.TOEL(i);
                            xe= mesh.elements_density(el_index);
                            obj.xe_range(i)=xe;
                            new_xe=exp(-obj.beta * (1 - obj.xe_range(i))) - (1 - obj.xe_range(i)) * exp(-obj.beta);
                            mesh.elements_density(el_index)=new_xe;
                        end                

            end

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function filter_sensitivities(obj,reader,mesh,TOO,TOC)

            if reader.Filter==2||reader.Filter==-2
                
                     for i =1:length(obj.TOEL)
                            %dxe=-(obj.beta * (tanh(obj.beta * (obj.mu - obj.xe_range(i)))^2 - 1)) / (2 * tanh(obj.beta * obj.mu));
                            dxe=(obj.beta*(tanh(obj.beta*(obj.mu - obj.xe_range(i)))^2 - 1))/(tanh(obj.beta*(obj.mu - 1)) - tanh(obj.beta*obj.mu));
                            TOO.dfdx(i)=TOO.dfdx(i)*dxe;
                     end
                 % Constraints Heaviside
                 for j=1:obj.m
                     for i =1:length(obj.TOEL)
                            %dxe=( tanh(obj.beta*obj.mu)+tanh(obj.beta*(obj.xe_range(i)-obj.mu)) )/( tanh(obj.beta*obj.mu)+tanh(obj.beta*(1-obj.mu)) );
                            dxe=(obj.beta*(tanh(obj.beta*(obj.mu - obj.xe_range(i)))^2 - 1))/(tanh(obj.beta*(obj.mu - 1)) - tanh(obj.beta*obj.mu));
                            TOC.dfdx(j,i)=TOC.dfdx(j,i)*dxe;
                     end
                 end       
            end

            if reader.Filter==3||reader.Filter==-3
            % Heaviside
                     for i =1:length(obj.TOEL)
                            if obj.xe_range(i)<obj.mu
                                xde=obj.beta*exp(-obj.beta*(1-obj.xe_range(i)/obj.mu))+exp(-obj.beta);
                            else
                                xde=obj.beta*exp(-obj.beta*(1-obj.xe_range(i)/obj.mu)/(1-obj.mu))+exp(-obj.beta);
                            end
                            TOO.dfdx(i)=TOO.dfdx(i)*xde;
                     end
                 % Constraints Heaviside
                 for j=1:obj.m
                     for i =1:length(obj.TOEL)
                            if obj.xe_range(i)<obj.mu
                                xde=obj.beta*exp(-obj.beta*(1-obj.xe_range(i)/obj.mu))+exp(-obj.beta);
                            else
                                xde=obj.beta*exp(-obj.beta*(1-obj.xe_range(i)/obj.mu)/(1-obj.mu))+exp(-obj.beta);
                            end
                            TOC.dfdx(j,i)=TOC.dfdx(j,i)*xde;
                     end
                 end       
            end

            if reader.Filter == 5 || reader.Filter == -5
                % Heav dilation
                     for i =1:length(obj.TOEL)
                                xde=obj.beta * exp(-obj.beta * (1 - obj.xe_range(i))) + exp(-obj.beta);
                            TOO.dfdx(i)=TOO.dfdx(i)*xde;
                     end
                 % Constraints Heaviside
                 for j=1:obj.m
                     for i =1:length(obj.TOEL)
                                xde=obj.beta * exp(-obj.beta * (1 - obj.xe_range(i))) + exp(-obj.beta);
                            TOC.dfdx(j,i)=TOC.dfdx(j,i)*xde;
                     end
                 end             

            end

        if reader.Filter == 4 || reader.Filter == -4
            % Initialize parameters
            %beta = obj.beta;          % Beta value from the object
            N_e = obj.TOEL;           % Element indices in the topology optimization
            rho = mesh.elements_density(N_e);  % Current density values for elements in N_e
            
            % Precompute terms for the derivative
            Z = sum(exp(obj.beta * rho));  % Sum of exponentials in the numerator
            C = length(N_e);          % Total number of elements (constant denominator)

                     for i =1:length(obj.TOEL)
                        xde = exp(obj.beta * rho(i)) / Z;  % Derivative term
                        TOO.dfdx(i)=TOO.dfdx(i)*xde;
                     end
                 % Constraints Heaviside
                 for j=1:obj.m
                     for i =1:length(obj.TOEL)
                        % Get the current element index
                        % Compute the derivative of rho_e_tilde with respect to rho_j
                        xde = exp(obj.beta * rho(i)) / Z;  % Derivative term
                        TOC.dfdx(j,i)=TOC.dfdx(j,i)*xde;
                     end
                 end      
            
        end


            if reader.Filter>0
                 xx=TOO.dfdx(1:length(obj.TOEL))';
                 TOO.dfdx(1:length(obj.TOEL))=obj.Helmholtz_xx(mesh,xx);

                 % Constraints Helmholtz
                 for j=1:obj.m
                        xx=TOC.dfdx(j,1:length(obj.TOEL));
                        TOC.dfdx(j,1:length(obj.TOEL))=obj.Helmholtz_xx(mesh,xx);
                 end
            end

        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [nodes] = retrieveelementnumberofnodes(~,etype_element)


                if etype_element == "CPS4"      % DIM-2 4-node quadrangle
                    nodes=4;      
                elseif etype_element == "T3D2"  % DIM-1 2-node line
                    nodes=2;
                elseif etype_element == "T3D3"  % DIM-1 3-node line
                    nodes=2;
                elseif etype_element == "CPS8"  % DIM-2 8-node second-order quadrangle
                    nodes=4;
                elseif etype_element == "C3D8"  % DIM-3 8-node Hexahedral 8 node element
                    nodes=8;
                elseif etype_element == "C3D20" % DIM-3 20-node Hexahedral 20 node element
                    nodes=8;
                else
                    nodes=-1;
                end
                          
        end       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = updateElMultiplier(obj, new_multiplier)
            % Function to update el_multiplier and recalculate r and rd
            obj.el_multiplier = new_multiplier;
            obj.r = max([obj.Element_size * obj.el_multiplier, 0.00011]);
            obj.rd = (obj.r / 2 / sqrt(3))^2;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end


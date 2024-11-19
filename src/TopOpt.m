classdef TopOpt
    %TO Summary of this class goes here
    %   Detailed explanation goes here

    properties
        m
        n
        xmin
        xmax
        dfdx
        low
        upp
        c
        d
        a0
        a
        outeriter
        maxiter
        kkttol
        xold1
        xold2
        fval
        xval
        TOEL
        f0val
        df0dx
        fval_iter
        f0val_iter
        xbc_iter
        Voltage_value
        Voltage_initial
        onlyvol
        mma_move
        mma_incr
        mma_decr
        mma_init
        Hev_update
        Hev_max
        Hev_init
    end

    methods
        function obj = TopOpt(reader,mesh)
            obj.m=length(reader.TopOpt_ConstraintName);
            obj.outeriter = 0;
            obj.maxiter = 100;
            if reader.TopOpt_DesignElements==""
                obj.TOEL=[];
            else
                obj.TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            end
            obj.n =length(obj.TOEL)+length(reader.TObcval);
            obj.xval=zeros(obj.n,1);
            if isempty(reader.TopOpt_Initial_x)
                reader.TopOpt_Initial_x=1;
                obj.xval(1:length(obj.TOEL))=mesh.elements_density(obj.TOEL);
            elseif (isnumeric(reader.TopOpt_Initial_x))
                mesh.elements_density(obj.TOEL)=ones(length(mesh.elements_density(obj.TOEL)),1)*reader.TopOpt_Initial_x;
                obj.xval(1:length(obj.TOEL))=mesh.elements_density(obj.TOEL);
            elseif reader.TopOpt_DesignElements==""
                xx_init=[];
            elseif isstring(reader.TopOpt_Initial_x) || ischar(reader.TopOpt_Initial_x)
                if contains(reader.TopOpt_Initial_x, '_2filter')
                    disp('The string contains "_2filter".');
                    % Perform actions for strings containing '_2filter'
                    xin=read_vtk_cell_data(reader.TopOpt_Initial_x);
                    xin(xin>0.95)=1;
                    xin(xin<0.95)=0.0001;
                    mesh.elements_density = xin;
                else
                    disp('The string does not contain "_2filter".');
                    % Perform other actions
                    mesh.elements_density = read_vtk_cell_data(reader.TopOpt_Initial_x);
                end
            else
                %xx_init = obj.ReadVTKxx(reader.TopOpt_Initial_x);
                %mesh.elements_density(obj.TOEL)=xx_init;
                obj.xval(1:length(obj.TOEL))=mesh.elements_density(obj.TOEL);
            end

            for i=1:length(reader.TObcval)
                obj.xval(length(obj.TOEL)+i)=(reader.TObcval(i)-reader.TObcminval(i))/(reader.TObcmaxval(1)-reader.TObcminval(1));
            end
            obj.xmin=zeros(obj.n,1);
            obj.xmax=ones(obj.n,1);
            obj.xold1=obj.xval;
            obj.xold2=obj.xval;
            %obj.xmax(end)=reader.TObcmaxval(1);
            %obj.xmin(end)=reader.TObcminval(1);
            obj.dfdx=zeros(obj.m,obj.n);
            obj.fval=zeros(obj.m,1);
            obj.low     = 0.3;
            obj.upp     = 0.7;
            obj.c       = ones(obj.m,1)*20000;
            obj.d       = ones(obj.m,1);
            obj.a0      = 1;
            obj.a       = zeros(obj.m,1);
            obj.kkttol = 1e-8;
            obj.f0val_iter=zeros(obj.maxiter+1,1);
            obj.fval_iter=zeros(obj.maxiter+1,obj.m);
            obj.xbc_iter=zeros(obj.maxiter+1,length(reader.TObcval));

            dofs_TO=zeros(length(mesh.data.NODE)*2,1);
            for i=1:length(obj.TOEL)
                element_nodes=mesh.data.ELEMENTS{i};
                for j=1:length(element_nodes)
                    dofs_TO(j*2-1)=1;
                    dofs_TO(j*2)=1;
                end
            end
            obj.mma_move=0.5;
            obj.mma_incr=1.3;
            obj.mma_decr=0.7;
            obj.mma_init=0.5;
            obj.Hev_update = 8;
            obj.Hev_max=100;
            obj.Hev_init=1;

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function  obj=runMMA(obj,reader,mesh)
            close all
            bcinit = BCInit(reader, mesh);
                if reader.Filter>0
                    filtering = Filtering(reader,mesh);
                    filtering.beta=obj.Hev_init;
                end
            solver = Solver(mesh, bcinit);

            % Specify the source file and destination folder
            sourceFile = reader.filename;
            destinationFolder = append(reader.rst_folder,'input.txt');
            copyfile(sourceFile, destinationFolder);
            sourceFile = reader.meshFileName;
            destinationFolder = append(reader.rst_folder,'mesh.inp');
            copyfile(sourceFile, destinationFolder);

            if strcmp(reader.solver,'NR')
                    residual_norm=solver.runNewtonRaphson(reader, mesh, bcinit);
                    odd_numbers = 1:2:length(solver.soldofs);
                    Tdofs=solver.soldofs(odd_numbers);
                    if residual_norm>10000  || not(isempty(Tdofs(Tdofs<0)))% divergence in NR catch
                        warning('NR DIVERGED, changing to Arc-len!!!');
                        for i=1:length(bcinit.dofs_free_)
                            df=bcinit.dofs_free_(i);
                            if mod(df, 2)==0
                                solver.soldofs(df)=0.0;
                            else 
                                solver.soldofs(df)=0;
                            end
                        end
                        funAL = @(t) solver.funArcLen(reader,bcinit,mesh,t);
                        [ufree] = arc_length_Crisfield(funAL,solver.soldofs(bcinit.dofs_free_));
                        solver.soldofs(bcinit.dofs_free_)=ufree(:,end);
                        if strcmp(reader.physics,'decoupledthermoelectromechanical')
                            % Extract the necessary variables
                            [solver.KStiff,KThermalLoad] = solver.Assembly_DecoupledThermoMech(reader, mesh);
                            Temperature_solution = solver.soldofs(1:2:end);
                            solver.KUT=KThermalLoad;
                            solver.loadVector_mech=solver.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                            solver.SolveLinearSystemInParallel(reader.physics,bcinit)
                        end
                    end
                elseif strcmp(reader.solver,'Arc-len')
                        for i=1:length(solver.soldofs)/2
                            if solver.soldofs(i*2-1)==0
                                solver.soldofs(i*2-1)=str2double(reader.T0);
                            elseif solver.soldofs(i*2)==0
                                solver.soldofs(i*2)=0.01;
                            end
                        end
                    funAL = @(t) solver.funArcLen(reader,bcinit,mesh,t);
                    [ufree] = arc_length_Crisfield(funAL,solver.soldofs(bcinit.dofs_free_));
                    solver.soldofs(bcinit.dofs_free_)=ufree(:,end);
                    if strcmp(reader.physics,'decoupledthermoelectromechanical')
                        % Extract the necessary variables
                        [solver.KStiff,KThermalLoad] = solver.Assembly_DecoupledThermoMech(reader, mesh);
                        Temperature_solution = solver.soldofs(1:2:end);
                        solver.KUT=KThermalLoad;
                        solver.loadVector_mech=solver.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                        solver.SolveLinearSystemInParallel(reader.physics,bcinit)
                    end
                else
                    warning('No solver recognized, changing to Arc-len!!!');
                        for i=1:length(solver.soldofs)/2
                            if solver.soldofs(i*2-1)==0
                                solver.soldofs(i*2-1)=str2double(reader.T0);
                            end
                            if solver.soldofs(i*2)==0
                                solver.soldofs(i*2)=0.01;
                            end
                        end
                    funAL = @(t) solver.funArcLen(reader,bcinit,mesh,t);
                    [ufree] = arc_length_Crisfield(funAL,solver.soldofs(bcinit.dofs_free_));
                    solver.soldofs(bcinit.dofs_free_)=ufree(:,end);
                    if strcmp(reader.physics,'decoupledthermoelectromechanical')
                        % Extract the necessary variables
                        [solver.KStiff,KThermalLoad] = solver.Assembly_DecoupledThermoMech(reader, mesh);
                        Temperature_solution = solver.soldofs(1:2:end);
                        solver.KUT=KThermalLoad;
                        solver.loadVector_mech=solver.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                        solver.SolveLinearSystemInParallel(reader.physics,bcinit)
                    end
                end


            for i = 1:length(reader.TObcval)
                if(strcmp(reader.TObctype,'Voltage'))
                    obj.Voltage_initial=reader.TObcval(i);
                end
            end
            TOO = TO_Objectives(reader,mesh,bcinit);
            TOO.CalculateObjective(reader,mesh,solver)
            TOC = TO_Constraints(reader,mesh,bcinit);
            TOC.CalculateConstraint(reader,mesh,solver);
            post = Postprocessing();
            post.initVTK(reader,mesh);
            currentDate = datestr(now, 'yyyy_mm_dd_HH_MM');
            folderName = fullfile(reader.rst_folder, append(reader.Rst_name,'_', currentDate));
            filtering.folderName=folderName;
            mkdir(folderName);
            post.VTK_x_TV(mesh,solver,append([folderName,'/',reader.Rst_name, '_TVMMA_',currentDate,'_',num2str(1000+obj.outeriter),'.vtk']))
            if strcmp(reader.physics,'decoupledthermoelectromechanical')
                post.VTK_x_U(mesh,solver,append([folderName,'/',reader.Rst_name, '_UMMA_',currentDate,'_',num2str(1000+obj.outeriter),'.vtk']))
            end
            %% New derivatives
            TOO.CalculateObjective(reader,mesh,solver)
            TOC.CalculateConstraint(reader,mesh,solver);
            obj.f0val = TOO.fval;
            obj.df0dx = TOO.dfdx;
            obj.fval = TOC.fval;
            obj.dfdx = TOC.dfdx;
            obj.f0val_iter(obj.outeriter+1)= TOO.fval;
            obj.fval_iter(obj.outeriter+1,:)= TOC.fval;

            if not(isempty(reader.TObcval))
                obj.xbc_iter(obj.outeriter+1,length(reader.TObcval))= max(obj.xold1(length(obj.TOEL)+1:length(obj.xold1)),[0]);
            end
            %obj.initMMA() % including dfdx!!! running sensitivities should return dfdx, and filters
            kktnorm = 1000;
            lowv=obj.low;
            uppv=obj.upp;

            while kktnorm > obj.kkttol && obj.outeriter < obj.maxiter 
                obj.outeriter = obj.outeriter+1;
                %postprocess.save()
                if obj.onlyvol==1
                    obj.df0dx(1:end-1)=0;
                    obj.dfdx(:,1:end-1)=0;
                end
                [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,lowv,uppv] = ...
                    mmasub(obj.m,obj.n,obj.outeriter,obj.xval,obj.xmin,obj.xmax,obj.xold1,obj.xold2, ...
                    obj.f0val,obj.df0dx,obj.fval,obj.dfdx,lowv,uppv,obj.a0,obj.a,obj.c,obj.d,...
                    obj.mma_move,obj.mma_incr,obj.mma_decr,obj.mma_init);

                if obj.onlyvol==1
                    xmma(1:end-1)=obj.xold1(1:end-1);
                end                
                
                %% Filter densities
                for i = 1:length(reader.TObcval)
                    if(strcmp(reader.TObctype,'Voltage'))
                        obj.Voltage_value=reader.TObcminval(i)+xmma(length(obj.TOEL)+i)*(reader.TObcmaxval(i)-reader.TObcminval(i));
                        reader.TObcval(i)=obj.Voltage_value;
                        reader.bcval(length(reader.bcval)-length(reader.TObcval)+i)=obj.Voltage_value;
                    end
                end
                %mesh_1 = Mesh(reader);
                for i=1:length(obj.TOEL)
                    filtering.it=obj.outeriter;
                    mesh.elements_density(obj.TOEL(i))=xmma(i);
                end
                if reader.Filter>0
                    if mod(obj.outeriter, obj.Hev_update) == 0
                        filtering.beta = filtering.beta * 2;
                        if filtering.beta> obj.Hev_max
                            filtering.beta=obj.Hev_max;
                        end
                    end
                    filtering.filter_densities(reader,mesh)
                end

                %% New NR starting point
                % New Voltage drop
                %obj.modifyNRStartingpoint()
                odd_numbers = 1:2:length(solver.soldofs);
                even_numbers = 2:2:length(solver.soldofs);
                prevdofs_odd=solver.soldofs(odd_numbers);
                prevdofs_even=solver.soldofs(even_numbers);

                %bcinit1 = BCInit(reader, mesh);
                solver = Solver(mesh, bcinit);

                for i=1:length(reader.TObcval)
                    nodes=mesh.retrieveNodalSelection(reader.TObcloc(i));
                    if(strcmp(reader.TObctype,'Voltage'))
                        %solver.soldofs(odd_numbers)=prevdofs_odd;
                        %solver.soldofs(even_numbers)=prevdofs_even*obj.Voltage_value/obj.Voltage_initial;
                        solver.soldofs(nodes*2)=obj.Voltage_value;
                    end
                end


                %% New Solve
                if strcmp(reader.solver,'NR')
                    residual_norm=solver.runNewtonRaphson(reader, mesh, bcinit);
                    Tdofs=solver.soldofs(odd_numbers);
                    if residual_norm>10000  || not(isempty(Tdofs(Tdofs<0)))% divergence in NR catch
                        warning('NR DIVERGED, changing to Arc-len!!!');
                        solver = Solver(mesh, bcinit);
        
                        for i=1:length(reader.TObcval)
                            nodes=mesh.retrieveNodalSelection(reader.TObcloc(i));
                            if(strcmp(reader.TObctype,'Voltage'))
                                solver.soldofs(odd_numbers)=prevdofs_odd;
                                solver.soldofs(even_numbers)=prevdofs_even*obj.Voltage_value/obj.Voltage_initial;
                                solver.soldofs(nodes*2)=obj.Voltage_value;
                            end
                        end
                        %[ufree] = arc_length_Lam_Morley(funALM,bcinit.loadVector_(solver.soldofs(bcinit.dofs_free_)),solver.soldofs(bcinit.dofs_free_));
                        funAL = @(t) solver.funArcLen(reader,bcinit,mesh,t);
                        [ufree] = arc_length_Crisfield(funAL,solver.soldofs(bcinit.dofs_free_));
                        solver.soldofs(bcinit.dofs_free_)=ufree(:,end);
                        if strcmp(reader.physics,'decoupledthermoelectromechanical')
                            % Extract the necessary variables
                            [solver.KStiff,KThermalLoad] = solver.Assembly_DecoupledThermoMech(reader, mesh);
                            Temperature_solution = solver.soldofs(1:2:end);
                            solver.KUT=KThermalLoad;
                            solver.loadVector_mech=solver.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                            solver.SolveLinearSystemInParallel(reader.physics,bcinit)
                        end
                    end
                elseif strcmp(reader.solver,'Arc-len')
                    funAL = @(t) solver.funArcLen(reader,bcinit,mesh,t);
                    [ufree] = arc_length_Crisfield(funAL,solver.soldofs(bcinit.dofs_free_));
                    solver.soldofs(bcinit.dofs_free_)=ufree(:,end);
                    if strcmp(reader.physics,'decoupledthermoelectromechanical')
                        % Extract the necessary variables
                        [solver.KStiff,KThermalLoad] = solver.Assembly_DecoupledThermoMech(reader, mesh);
                        Temperature_solution = solver.soldofs(1:2:end);
                        solver.KUT=KThermalLoad;
                        solver.loadVector_mech=solver.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                        solver.SolveLinearSystemInParallel(reader.physics,bcinit)
                    end
                else
                    warning('No solver recognized, changing to Arc-len!!!');
                    funAL = @(t) solver.funArcLen(reader,bcinit,mesh,t);
                    [ufree] = arc_length_Crisfield(funAL,solver.soldofs(bcinit.dofs_free_));
                    solver.soldofs(bcinit.dofs_free_)=ufree(:,end);
                    if strcmp(reader.physics,'decoupledthermoelectromechanical')
                        % Extract the necessary variables
                        [solver.KStiff,KThermalLoad] = solver.Assembly_DecoupledThermoMech(reader, mesh);
                        Temperature_solution = solver.soldofs(1:2:end);
                        solver.KUT=KThermalLoad;
                        solver.loadVector_mech=solver.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                        solver.SolveLinearSystemInParallel(reader.physics,bcinit)
                    end
                end

                post.VTK_x_TV(mesh,solver,append([folderName,'/',reader.Rst_name, '_TVMMA_',currentDate,'_',num2str(1000+obj.outeriter),'.vtk']))
                if strcmp(reader.physics,'decoupledthermoelectromechanical')
                    post.VTK_x_U(mesh,solver,append([folderName,'/',reader.Rst_name, '_UMMA_',currentDate,'_',num2str(1000+obj.outeriter),'.vtk']))
                end
            %% New objective, constraints and derivatives
                TOC = TO_Constraints(reader,mesh,bcinit);
                TOO = TO_Objectives(reader,mesh,bcinit);
                TOO.CalculateObjective(reader,mesh,solver)
                TOC.CalculateConstraint(reader,mesh,solver);
                if reader.Filter>0
                    filtering.filter_sensitivities(reader,mesh,TOO,TOC)
                end
                %f0valold=obj.f0val;
                obj.f0val = TOO.fval;
                obj.df0dx = TOO.dfdx;
                obj.fval = TOC.fval;
                obj.dfdx = TOC.dfdx;

                obj.f0val_iter(obj.outeriter+1)= TOO.fval;
                obj.fval_iter(obj.outeriter+1,:)= TOC.fval;
                if not(isempty(reader.TObcval))
                    for bc=1:length(reader.TObcval)
                        obj.xbc_iter(obj.outeriter+1,bc)= xmma(length(obj.TOEL)+bc);
                    end
                end
                %obj.FilteringSensitivities()

                %% MMA parameters update
                obj.xold2=obj.xold1;
                obj.xold1=obj.xval;
                obj.xval=xmma;

                %% write results
                %postprocesing.save()


                kktnorm=norm((obj.xold2-obj.xval)./obj.xval)/length(obj.xval);
                post.PlotIter(1,reader,obj.outeriter+1,obj.f0val_iter,obj.fval_iter,obj.xbc_iter)
                %saveas(1, append([reader.rst_folder,reader.Rst_name,'_',currentDate,'.png']), 'png')
                saveas(1, append([folderName,'/',reader.Rst_name, 'MMA_',currentDate,'_',num2str(1000+obj.outeriter),'.png']), 'png')
                %post.SaveIterCSV(append([reader.rst_folder,reader.Rst_name,'_',currentDate,'.csv']),reader,obj.outeriter+1,obj.f0val_iter,obj.fval_iter,obj.xbc_iter)
                %saveas(1, append([folderName,'/',reader.Rst_name, 'MMA_',currentDate,'_',num2str(1000+obj.outeriter),'.vtk']), 'png')
                post.SaveIterCSV(append([folderName,'/',reader.Rst_name, 'MMA_',currentDate,'_','.csv']),reader,obj.outeriter+1,obj.f0val_iter,obj.fval_iter,obj.xbc_iter)
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function FiniteDifferences(obj,reader,mesh,bcinit,solver,fun)
            mesh_copy=mesh;
            tol=1e-6;
            FDmaxiter=100;
            density=0.9;density_old=1;
            mesh_copy.element_densities(obj.TOEL(1))=density;
            f_old=fun(reader,mesh,solver);
            solver=Solver(reader,mesh,bcinit);
            err=1;
            iter=1;
            while err>tol && iter<FDmaxiter
                solver.Assembly(reader,mesh_copy,bcinit)
                solver.runNewtonRaphson();
                f=fun(reader,mesh,solver);
                df=(f-f_old)/(density-density_old);
                if iter>1
                    err=(df-df_old)/df;
                end
                df_old=df;
                density_old=density;
                density=(1-density)/2+density;
                mesh_copy.element_densities(obj.TOEL(1))=density;
                iter=iter+1;

            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function  xx=ReadVTKxx(~,filename)
            A=regexp(fileread((filename)),'\n','split');
            An=find(startsWith(A,'CELL_DATA')==1);
            B=strtrim(strsplit(char(A(startsWith(A,'CELL_DATA')==1)),' '));
            Nxxdata=str2double(char(B(2)));
            xx=zeros(Nxxdata,1);
            for i=1:Nxxdata
                xx(i)=str2double(char(A(An+2+i)));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end


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
    end
    
    methods
        function obj = TopOpt(reader,mesh)
            obj.m=length(reader.TopOpt_ConstraintName);
            obj.outeriter = 0;
            obj.maxiter = 50;
            obj.TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            obj.n =length(obj.TOEL)+length(reader.TObcval);
            obj.xval=zeros(obj.n,1);
            obj.xval(1:length(obj.TOEL))=mesh.elements_density(obj.TOEL);
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
            obj.c       = ones(obj.m,1)*1000;
            obj.d       = ones(obj.m,1);
            obj.a0      = 1;
            obj.a       = zeros(obj.m,1);
            obj.kkttol = 1e-6;
            obj.f0val_iter=zeros(obj.maxiter+1,1);
            obj.fval_iter=zeros(obj.maxiter+1,obj.n);

            dofs_TO=zeros(length(mesh.data.NODE)*2,1);
            for i=1:length(obj.TOEL)
                element_nodes=mesh.data.ELEMENTS{i};
                for j=1:length(element_nodes)
                    dofs_TO(j*2-1)=1;
                    dofs_TO(j*2)=1;
                end
            end

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function  runMMA(obj,reader,mesh)
            bcinit = BCInit(reader, mesh);
            solver = Solver(mesh, bcinit);
            solver.runNewtonRaphson(reader, mesh, bcinit);
            TOO = TO_Objectives(reader,mesh,bcinit);
            TOO.CalculateObjective(reader,mesh,solver)
            TOC = TO_Constraints(reader,mesh,bcinit);
            TOC.CalculateConstraint(reader,mesh,solver);
            post = Postprocessing();
            post.initVTK(reader,mesh);
            currentDate = datestr(now, 'yyyy_mm_dd_HH_MM');
            post.VTK_x_TV(mesh,solver,append([reader.rst_folder, 'MMA_date_',currentDate,'_iter_',num2str(1000+obj.outeriter),'.vtk']))
            
            %% New derivatives
            TOO.CalculateObjective(reader,mesh,solver)
            TOC.CalculateConstraint(reader,mesh,solver);
            obj.f0val = TOO.fval;
            obj.df0dx = TOO.dfdx;
            obj.fval = TOC.fval;
            obj.dfdx = TOC.dfdx;
            obj.f0val_iter(obj.outeriter+1)= TOO.fval;
            obj.fval_iter(obj.outeriter+1,:)= TOC.fval;
            %obj.initMMA() % including dfdx!!! running sensitivities should return dfdx, and filters
            kktnorm = 1000;
            while kktnorm > obj.kkttol && obj.outeriter < obj.maxiter 
                obj.outeriter = obj.outeriter+1;
                %postprocess.save()
                [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,obj.low,obj.upp] = ...
                    mmasub(obj.m,obj.n,obj.outeriter,obj.xval,obj.xmin,obj.xmax,obj.xold1,obj.xold2, ...
                    obj.f0val,obj.df0dx,obj.fval,obj.dfdx,obj.low,obj.upp,obj.a0,obj.a,obj.c,obj.d);
                %% Filter densities
                mesh.elements_density(obj.TOEL)=xmma(1:length(obj.TOEL));
                %obj.FilteringDensitites()
    
                %% New NR starting point
                % New Voltage drop
                %obj.modifyNRStartingpoint()
                solver.soldofs=bcinit.initialdofs_;
                for i=1:length(reader.TObcval)
                    nodes=mesh.retrieveNodalSelection(reader.TObcloc(i));
                    if(strcmp(reader.TObctype,'Voltage'))
                        Voltage_value=reader.TObcminval(i)+xmma(length(obj.TOEL)+i)*(reader.TObcmaxval(i)-reader.TObcminval(i));
                        solver.soldofs(nodes*2)=Voltage_value;
                    end
                end
    
                %% New Solve
                solver.runNewtonRaphson(reader, mesh, bcinit);
                post.VTK_x_TV(mesh,solver,append([reader.rst_folder, 'MMA_date_',currentDate,'_iter_',num2str(1000+obj.outeriter),'.vtk']))

                %% New objective, constraints and derivatives
                TOO.CalculateObjective(reader,mesh,solver)
                TOC.CalculateConstraint(reader,mesh,solver);
                obj.f0val = TOO.fval;
                obj.df0dx = TOO.dfdx;
                obj.fval = TOC.fval;
                obj.dfdx = TOC.dfdx;
                obj.f0val_iter(obj.outeriter+1)= TOO.fval;
                obj.fval_iter(obj.outeriter+1,:)= TOC.fval;
                %obj.FilteringSensitivities()
    
                %% MMA parameters update
                obj.xold2=obj.xold1;
                obj.xold1=obj.xval;
                obj.xval=xmma;
    
                %% write results
                %postprocesing.save()
    
                %% Convergence
                if obj.outeriter>5
                    kktnorm=norm((obj.xval-obj.xold1)./obj.xval);
                end
                post.PlotIter(1,reader,obj.outeriter+1,obj.f0val_iter,obj.fval_iter)
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
            solver1=Solver(reader,mesh,bcinit);
            err=1;
            iter=1;
            while err>tol && iter<FDmaxiter
                solver1.Assembly(reader,mesh_copy,bcinit)
                solver1.runNewtonRaphson();
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
        function TestSensitivities(obj,reader,mesh,solver)
            objectives = Objectives();
            list_of_functions={};
            for i=1:length(list_of_functions)
                list_of_functions(i);
                obj.FiniteDifferences()
            end

            constraints = Constraints();
            list_of_functions={};
            for i=1:length(list_of_functions)
                list_of_functions(i);
                obj.FiniteDifferences()
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    end
end


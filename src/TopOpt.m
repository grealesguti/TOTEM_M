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
        dV
        vx0
        TOEL
    end
    
    methods
        function obj = TopOpt(reader,mesh)
            m=length(reader.TopOpt_ConstraintName);
            obj.outeriter=0;
            obj.TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            obj.n =length(TOEL)+length(reader.TObcval);
            xval=zeros(n,1);
            xval(1:length(TOEL))=mesh.elements_density(obj.TOEL);
            for i=1:length(reader.TObcval)
                xval(length(TOEL)+i)=reader.TObcval(i);
            end
            obj.xmin=zeros(n,1);
            obj.xmax=ones(n,1);
            obj.xold1=xval;
            obj.xold2=xval;
            obj.xmax(end)=scalev;
            obj.dfdx=zeros(m,n);
            obj.fval=zeros(m,1);
            obj.low     = 0.3;
            obj.upp     = 0.7;
            obj.c       = ones(m,1)*1000;
            obj.d       = ones(m,1);
            obj.a0      = 1;
            obj.a       = zeros(m,1);
            obj.kkttol = 1e-6;

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
            fprintf('Initialized InputReader with filename: %s\n', inputfilename);
            fprintf('Initialized Mesh\n');
            bcinit = BCInit(reader, mesh);
            solver = Solver(reader,mesh,bcinit);
            solver.runNewtonRaphson();
            constraints=TO_Constraints();
            objective = TO_Objectives();
            %% New derivatives
            obj.f0val = objective.CalculateObjective(reader,mesh,solver);
            obj.df0dx = objective.Calculate_dfdx(reader,mesh,solver);
            obj.fval = constraints.CalculateConstraint(reader,mesh,solver);
            obj.dfdx = constraints.calculate_dfdx(reader,mesh,solver);
            %obj.initMMA() % including dfdx!!! running sensitivities should return dfdx, and filters
    
            while kktnorm > obj.kkttol && obj.outeriter < obj.maxouteriter 
                obj.outeriter = obj.outeriter+1;
                postprocess.save()
                [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,obj.low,obj.upp] = ...
                    mmasub(obj.m,obj.n,obj.outeriter,obj.xval,obj.xmin,obj.xmax,obj.xold1,obj.xold2, ...
                    f0val,df0dx,obj.fval,obj.dfdx,obj.low,obj.upp,obj.a0,obj.a,obj.c,obj.d);
                %% Filter densities
                xx=xmma;
                mesh.elements_density(obj.TOEL)=xmma;
                %obj.FilteringDensitites()
    
                %% New NR starting point
                % New Voltage drop
                %obj.modifyNRStartingpoint()
                solver.soldofs=bcinit.initialdofs_;
    
                %% New Solve
                solver.Assembly()
                solver.runNewtonRaphson();
    
                %% New objective, constraints and derivatives
                obj.f0val = objective.CalculateObjective(reader,mesh,solver);
                obj.df0dx = objective.Calculate_dfdx(reader,mesh,solver);
                obj.fval = constraints.CalculateConstraint(reader,mesh,solver);
                obj.dfdx = constraints.calculate_dfdx(reader,mesh,solver);
                %obj.FilteringSensitivities()
    
                %% MMA parameters update
                obj.xold2=obj.xold1;obj.xold1=obj.xval;
                obj.xval=xmma;
    
                %% write results
                postprocesing.save()
    
                %% Convergence
                kktnorm=changexval(obj.outeriter);
                postprocesing.plot()
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


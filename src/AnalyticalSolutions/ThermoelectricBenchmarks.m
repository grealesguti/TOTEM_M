classdef ThermoelectricBenchmarks < handle
    %THERMOELECTRICBENCHMARKS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Benchmarks
        BenchmarksFunctions
        diffFEM_ctemat
        FD_vals_ctemat
        diffFEM_nonlinmat
        FD_vals_nonlinmat
        diffFEM_ctemat_TEC
        FD_vals_ctemat_TEC
        diffFEM_nonlinmat_TEC
        FD_vals_nonlinmat_TEC
    end
    
    methods
        function obj = ThermoelectricBenchmarks()
            close all;
            %THERMOELECTRICBENCHMARKS Construct an instance of this class
            Benchmark_Perez_Aparicio_LinUncoupSEffect_HexLinear = "Benchmarks/Elements/Benchmark_HexLinear/input_LinearUncoupledSeebeck.txt";
            Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexLinear = "Benchmarks/Elements/Benchmark_HexLinear/input_NonLinCouplSEffect.txt";
            Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity = "Benchmarks/Elements/Benchmark_HexSerendipity/input_LinearUncoupledSeebeck.txt";
            Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexSerendipity = "Benchmarks/Elements/Benchmark_HexSerendipity/input_NonLinCouplSEffect.txt";

            obj.Benchmarks={
            Benchmark_Perez_Aparicio_LinUncoupSEffect_HexLinear;
            Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexLinear;
            Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity;
            Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexSerendipity};

            obj.BenchmarksFunctions=[1,2,1,2];

            %Benchmark_sensitivities ="Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect.txt";
            %[obj.diffFEM_ctemat, obj.FD_vals_ctemat] = obj.run_SingleFEM_diff(Benchmark_sensitivities);
            %[diffFEM_ctemat_TEC, FD_vals_ctemat_TEC] = thb.run_SingleFEM_diff("Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect.txt"); 

            %Benchmark_sensitivities ="Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect_nonlinmat.txt";
            %[obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat] = obj.run_SingleFEM_diff(Benchmark_sensitivities);

            %Benchmark_sensitivities ="TECTO/input_TECTO_StressConstrained_cte.txt";
            %[obj.diffFEM_ctemat_TEC, obj.FD_vals_ctemat_TEC] = obj.run_SingleFEM_diff_TEC(Benchmark_sensitivities); 
            %[diffFEM_ctemat_TEC, FD_vals_ctemat_TEC] = thb.run_SingleFEM_diff_TEC("TECTO/input_TECTO_StressConstrained_cte.txt"); 


            %Benchmark_sensitivities ="TECTO/input_TECTO_StressConstrained_noncte.txt";
            %[obj.diffFEM_nonlinmat_TEC, obj.diffFEM_nonlinmat_TEC] = obj.run_SingleFEM_diff_TEC(Benchmark_sensitivities); 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function []=runPhysicsBenchmarks(obj)
            for i=1:length(obj.BenchmarksFunctions)
                if(obj.BenchmarksFunctions(i)==1)
                    Perez_Aparicio_LinearUncoupledSeebeck(obj,obj.Benchmarks{i},i)
                elseif(obj.BenchmarksFunctions(i)==2)
                    Perez_Aparicio_CoupledSeebeck(obj,obj.Benchmarks{i},i)
                elseif(obj.BenchmarksFunctions(i)==3)
                    Perez_Aparicio_CoupledPeltier(obj,obj.Benchmarks{i},i)
                end
            end        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % https://link.springer.com/article/10.1007/s00466-006-0080-7
        function [xv,Tx,Vx,Ux] = Perez_Aparicio_LinearUncoupledSeebeck(~,filename,index)

                reader = InputReader(filename);
                fprintf('Initialized InputReader with filename: %s\n', filename);
                mesh = Mesh(reader);
                fprintf('Initialized Mesh\n');
                bcinit = BCInit(reader, mesh);
                fprintf('Initialized Loads\n');
                solver = Solver(mesh, bcinit);
                solver.runNewtonRaphson(reader, mesh, bcinit);
                fprintf('Postprocessing\n');
                postprocessing = Postprocessing();
                postprocessing.initVTK(reader,mesh);

                L=1.524e-3;
                steps=100;
                x_step=L/steps;
                xv=zeros(steps+1,1);
                Tx=zeros(steps+1,1);
                Vx=zeros(steps+1,1);
                Ux=zeros(steps+1,1);% int((T(x)-298.15)*alpha dx) (symbolic calculation)
                alpha=1e-3;
                Th=273.15;
                Tc=298.15;
                a=-4.26E-3;
                b=1.754;
                c=3.598E-7;
                d=1.804;
            
                kh=-0.0043*Th+2.9176;
                kc=-0.0043*Tc+2.9176;
                c1=(kh^2-kc^2)/L;
                               
                for i=1:101
                    xv(i)=x_step*(i-1);
                    Tx(i)=412-sqrt(169587-1.31e7*xv(i))+273.15;
                    Vx(i)= sqrt(0.018 - 1.41*xv(i)) + 2.36*xv(i) - 0.13;
                    Ux(i)=alpha*(412*xv(i) - (56529*169587^(1/2))/6550000 + (169587 - 13100000*xv(i))^(3/2)/19650000);
                end

                figure(index)
                hold on
                subplot(1, 3,1);
                hold on
                TFEM = postprocessing.Benchmark_T_PLOT_axis(index,solver,2);
                plot(xv,Tx)
                subplot(1, 3,2);
                hold on
                VFEM = postprocessing.Benchmark_V_PLOT_axis(index,solver,2);
                plot(xv,Vx)
                subplot(1, 3,3);
                hold on
                UFEM = postprocessing.Benchmark_U_PLOT_axis(index,solver,2);
                plot(xv,Ux)

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % https://link.springer.com/article/10.1007/s00466-006-0080-7
        function [xv,Tx,Vx,Ux,Power_Bench] = Perez_Aparicio_CoupledSeebeck(~,filename,index)

                reader = InputReader(filename);
                fprintf('Initialized InputReader with filename: %s\n', filename);
                mesh = Mesh(reader);
                fprintf('Initialized Mesh\n');
                bcinit = BCInit(reader, mesh);
                fprintf('Initialized Loads\n');
                solver = Solver(mesh, bcinit);
                solver.runNewtonRaphson(reader, mesh, bcinit);
                fprintf('Postprocessing\n');
                postprocessing = Postprocessing();
                postprocessing.initVTK(reader,mesh);


                L=1.524e-3;
                steps=100;
                x_step=L/steps;
                xv=zeros(steps,1);
                Tx=zeros(steps,1);
                Vx=zeros(steps,1);
                Ux=zeros(steps,1); % int((T(x)-298.15)*alpha dx)
                alpha=1e-3;

                for i=1:100
                    xv(i)=x_step*(i-1);
                    Tx(i)=3.794e7*xv(i)*(1.524e-3-xv(i))+298.15;
                    Vx(i)= 5.788E-2 - 4.913E1*xv(i) + 7.315E3*xv(i)^2;
                    Ux(i)=(-xv(i)^2*((37940000*xv(i))/3 - 8332820879051308801875/288230376151711744))*alpha;
                end
                j=3.199e6;
                Power = 5.788e-2*(j*0.0014^2);
                Power_FEM = CalculatePower(reader,mesh,solver);
                Power_Bench = abs(Power-Power_FEM);
                fprintf('Power benchmark, Analytical: %s FEM: %s\n', [Power,Power_FEM]);

                figure(index)
                hold on
                subplot(1, 3,1);
                hold on
                postprocessing.Benchmark_T_PLOT_axis(index,solver,2)
                plot(xv,Tx)
                subplot(1, 3,2);
                hold on
                postprocessing.Benchmark_V_PLOT_axis(index,solver,2)
                plot(xv,Vx)
                subplot(1, 3,3);
                hold on
                postprocessing.Benchmark_U_PLOT_axis(index,solver,2)
                plot(xv,Ux)

        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [xv,Tx,Vx] = Perez_Aparicio_CoupledPeltier(~,filename,index)
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [xv,Tv,Vv] = Benchmark1_2_3(~,reader,mesh,solver,postprocessing)
            %METHOD1 Summary of this method goes here
                g1 = reader.getmaterialproperty(1,'ElectricalConductivity');
                k1 = reader.getmaterialproperty(1,'ThermalConductivity');
                a1 = reader.getmaterialproperty(1,'Seebeck');    
                coordinates=zeros(3,length(mesh.data.NODE));
                for i=1:length(mesh.data.NODE)
                    coordinates(:,i)=mesh.data.NODE{i};
                end
                L=abs(max(coordinates(2,:))-min(coordinates(2,:)));
                l1=abs(max(coordinates(1,:))-min(coordinates(1,:)));
                l2=abs(max(coordinates(3,:))-min(coordinates(3,:)));
                A=l1*l2;
                boundaryConditions = reader.bctype;
                % Iterate through the cell array to find the index location of the string
                for i = 1:numel(boundaryConditions)
                    if strcmp(boundaryConditions{i}, 'heat_n')
                        qin = reader.bcval(i);
                    elseif strcmp(boundaryConditions{i}, 'Temperature')
                        Th = reader.bcval(i);
                    elseif strcmp(boundaryConditions{i}, 'Voltage')
                        if(abs(reader.bcval(i))>0)
                            Vc=reader.bcval(i);
                        end
                    end
                end
                V0=0;
                syms C31 Tc j1 
                eqTc1 = -j1^2 / (k1*g1) * L^2/2 + C31*L + Th == Tc;
                eqVc1 = -j1/g1*L + a1*j1^2/(g1*k1)*L^2/2 - a1*C31*L == Vc;              
                eqq = (a1*j1*Tc-k1*(-j1^2/(k1*g1)*L+C31))*A == -qin;
            
                Ss=[eqTc1, eqVc1, eqq];
                % SoLve the equations
                S = solve(Ss, ...
                         [Tc, C31, j1]); 
                TCs=eval(S.Tc)
                C31s=eval(S.C31)
                Vcs=eval(S.Vc)
                steps=100;
                x_step=L/steps;
                xv=zeros(steps,1);
                Tv=zeros(steps,1);
                Vv=zeros(steps,1);
                for i=1:100
                    xv(i)=x_step*(i-1);
                    Tv(i)=-j1^2 / (k1*g1) * xv(i)^2/2 + C31*xv(i) + Th;
                    Vv(i)=-j1/g1*xv(i) + a1*j1^2/(g1*k1)*xv(i)^2/2 - a1*C31*xv(i) + V0;
                end

                figure(1)
                plot(xv,Tv)

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [diff, err] = Finite_Differences_DensityElement(~, reader, mesh, bcinit, solver, eval_fun, idx)
            % Get the test element for density evaluation
            TOEL = mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            test_element = TOEL(idx);
            
            % Initialize density and iterate variables
            xx = mesh.elements_density(test_element);
            xx_iter = xx * 0.97;
            mesh_1 = Mesh(reader);
                if isempty(reader.TopOpt_Initial_x)
                    reader.TopOpt_Initial_x=1;
                else
                    mesh_1.elements_density(TOEL)=ones(length(mesh_1.elements_density(TOEL)),1)*reader.TopOpt_Initial_x;
                end            
            % Calculate the initial objective function value
            value_0 = eval_fun(reader, mesh, solver);
            % Initialize history arrays to store err and diff values
                err_history = [];
                diff_history = [];
    
            % Set convergence tolerance and counters
            Tol = 1e-6;
            cc = 1;
            err = 100;
            diff_old = 0; % Initialize diff_old
            max_iterations = 15;
            
            while (err > Tol && cc <= max_iterations) || cc<3
                % Calculate finite difference and update density
                xx_iter = (xx + xx_iter) / 2;
                mesh_1.elements_density(test_element) = xx_iter;
                
                % Create a new solver and run the Newton-Raphson method
                solver_1 = Solver(mesh_1,bcinit);
                solver_1.runNewtonRaphson(reader, mesh_1, bcinit);
                
                % Calculate the objective function value with the updated density
                value = eval_fun(reader, mesh_1, solver_1);
                
                % Calculate finite difference and check for convergence
                diff = (value - value_0) / (xx_iter - xx);
                if cc > 1
                    err = abs((diff_old - diff) / diff_old);
                end
                diff_old = diff;
                % Print information to the screen
                fprintf('FD Iteration: %d, Error: %e, value: %e, value_0: %e, diff: %e, xx: %e\n', cc, err, value, value_0, diff, xx_iter);     

                % Store err and diff values in history arrays
                err_history(cc) = err;
                diff_history(cc) = diff;

                % Increment the counter
                cc = cc + 1;
            end

                % Plot err and diff history at the end of iterations
                figure;
                subplot(2, 1, 1);
                plot(1:length(err_history), err_history);
                xlabel('Iteration');
                ylabel('Error');
                title('Convergence Plot - Error');
                
                subplot(2, 1, 2);
                plot(1:length(diff_history), diff_history);
                xlabel('Iteration');
                ylabel('Difference');
                title('Convergence Plot - Difference');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [diff, err, cc] = Finite_Differences_bc(~, filename, reader, mesh, solver, eval_fun, index,index_TO)
            % Calculate the initial objective function value
            value_0 = eval_fun(reader, mesh, solver);
            TOEL = mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);

            % Get the test element for density evaluation
            bcval = reader.bcval(index);
            bcval_iter=bcval*0.95;
            
            reader_1 = InputReader(filename);
            bcmax =  reader.TObcmaxval(index_TO);
            bcmin = reader.TObcminval(index_TO);
            dbc=bcmax-bcmin;
            xbcval0=(bcval-bcmin)/dbc;
            xbcval_iter=xbcval0*0.95;


            % Set convergence tolerance and counters
            Tol = 1e-3;
            cc = 1;
            err = 1;
            diff_old = 0; % Initialize diff_old

                % Store err and diff values in history arrays
                err_history = [];
                diff_history= [];

            max_iterations = 15;
            
            while err > Tol && cc <= max_iterations
                % Calculate finite difference and update density
                xbcval_iter = (xbcval0 + xbcval_iter) / 2;
                bcvalue=bcmin+dbc*xbcval_iter;
                reader_1.bcval(index)=bcvalue;
                mesh_1 = Mesh(reader_1);
                mesh_1.elements_density=mesh.elements_density; 

                bcinit_1 = BCInit(reader_1, mesh_1);
                
                % Create a new solver and run the Newton-Raphson method
                solver_1 = Solver(mesh_1,bcinit_1);
                solver_1.runNewtonRaphson(reader_1, mesh_1, bcinit_1);
                
                % Calculate the objective function value with the updated density
                value = eval_fun(reader_1, mesh_1, solver_1);
                
                % Calculate finite difference and check for conve-rgence
                diff = (value - value_0) / (xbcval_iter - xbcval0);
                if cc > 2
                    err = abs((diff_old - diff) / diff_old);
                end
                diff_old = diff;
                fprintf('FD Iteration: %d, Error: %e, value: %e, value_0: %e, diff: %e, bcval_iter: %e\n', cc, err, value, value_0, diff, bcval_iter);     
                
                % Store err and diff values in history arrays
                err_history(cc) = err;
                diff_history(cc) = diff;

                % Increment the counter
                cc = cc + 1;
            end

                % Plot err and diff history at the end of iterations
                figure;
                subplot(2, 1, 1);
                plot(1:length(err_history), err_history);
                xlabel('Iteration');
                ylabel('Error');
                title('Convergence Plot - Error');
                
                subplot(2, 1, 2);
                plot(1:length(diff_history), diff_history);
                xlabel('Iteration');
                ylabel('Difference');
                title('Convergence Plot - Difference');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [diff, err, cc] = Finite_Differences_bcTO(~, filename, reader, mesh, solver, eval_fun, index,index_TO)
            % Calculate the initial objective function value
            value_0 = eval_fun(reader, mesh, solver);
            TOEL = mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);

            % Get the test element for density evaluation
            bcval = xmma(length(obj.TOEL)+i);
            bcval_iter=bcval*0.95;
            
            reader_1 = InputReader(filename);
            bcmax =  reader.TObcmaxval(index_TO);
            bcmin = reader.TObcminval(index_TO);
            dbc=bcmax-bcmin;
            xbcval0=(bcval-bcmin)/dbc;
            xbcval_iter=xbcval0*0.95;

            
            % Set convergence tolerance and counters
            Tol = 1e-3;
            cc = 1;
            err = 1;
            diff_old = 0; % Initialize diff_old

                % Store err and diff values in history arrays
                err_history = [];
                diff_history= [];

            max_iterations = 15;
            
            while err > Tol && cc <= max_iterations
                % Calculate finite difference and update density
                xbcval_iter = (xbcval0 + xbcval_iter) / 2;
                bcvalue=bcmin+dbc*xbcval_iter;
                reader_1.bcval(index)=bcvalue;
                mesh_1 = Mesh(reader_1);
                mesh_1.elements_density=mesh.elements_density;
                if isempty(reader_1.TopOpt_Initial_x)
                    reader.TopOpt_Initial_x=1;
                else
                    mesh_1.elements_density(TOEL)=ones(length(mesh_1.elements_density(TOEL)),1)*reader_1.TopOpt_Initial_x;
                end      

                bcinit_1 = BCInit(reader_1, mesh_1);
                
                % Create a new solver and run the Newton-Raphson method
                solver_1 = Solver(mesh_1,bcinit_1);
                solver_1.runNewtonRaphson(reader_1, mesh_1, bcinit_1);
                
                % Calculate the objective function value with the updated density
                value = eval_fun(reader_1, mesh_1, solver_1);
                
                % Calculate finite difference and check for conve-rgence
                diff = (value - value_0) / (xbcval_iter - xbcval0);
                if cc > 2
                    err = abs((diff_old - diff) / diff_old);
                end
                diff_old = diff;
                fprintf('FD Iteration: %d, Error: %e, value: %e, value_0: %e, diff: %e, bcval_iter: %e\n', cc, err, value, value_0, diff, bcval_iter);     
                
                % Store err and diff values in history arrays
                err_history(cc) = err;
                diff_history(cc) = diff;

                % Increment the counter
                cc = cc + 1;
            end

                % Plot err and diff history at the end of iterations
                figure;
                subplot(2, 1, 1);
                plot(1:length(err_history), err_history);
                xlabel('Iteration');
                ylabel('Error');
                title('Convergence Plot - Error');
                
                subplot(2, 1, 2);
                plot(1:length(diff_history), diff_history);
                xlabel('Iteration');
                ylabel('Difference');
                title('Convergence Plot - Difference');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [diffFEM, FD_vals] = run_FEM_diff(~, filename)
            reader = InputReader(filename);
            mesh = Mesh(reader);
            bcinit = BCInit(reader, mesh);
            solver = Solver(mesh, bcinit);
            solver.runNewtonRaphson(reader, mesh, bcinit);
            TOO = TO_Objectives(reader,mesh,bcinit);
            TOO.CalculateObjective(reader,mesh,solver);
            TOC = TO_Constraints(reader,mesh,bcinit);
            TOC.CalculateConstraint(reader,mesh,solver)

            cc=1; n_con = length(reader.TopOpt_ConstraintValue); n_bc = length(reader.TObctype);
            diffFEM=zeros(1+length(n_con),length(TOO.TOEL));
            diffFEM(1,:)=TOO.dfdx;
            diffFEM(2,:,:)=TOC.dfdx;
            idx=1;

            TOO_1 = TO_Objectives(reader,mesh,bcinit);
            if (reader.TopOpt_Objective=="AverageTemperature")
                        eval_fun=@(reader,mesh,solver) TOO_1.fval_AverageTemp(reader,mesh,solver);
                        [FD_vals(cc), err, cc_FD] = Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,idx);
            end
            cc=cc+1;

            FD_vals=zeros(ncon+n_bc,1);
            for i=1:n_con
            TOC_1 = TO_Constraints(reader,mesh,bcinit);

                if (reader.TopOpt_ConstraintName=="Power")
                        eval_fun=@(reader,mesh,solver) CalculatePower(reader,mesh,solver);
                        [FD_vals(cc), err, cc_FD] = Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,idx);
                elseif (reader.TopOpt_ConstraintName=="Volume")
                        eval_fun=@(reader,mesh,solver) TOC_1.fval_Volume(reader,mesh,solver,1);
                        [FD_vals(cc), err, cc_FD] = Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,idx);
                end

                cc=cc+1;
            end

            for i=1:n_bc
                    if (reader.TObctype=="Voltage")
                        eval_fun=@(reader,mesh,solver) TOO_1.fval_AverageTemp(reader,mesh,solver);
                        [FD_vals(cc), err, cc_FD] = Finite_Differences_bc( reader, mesh, bcinit, solver, eval_fun);
                    end
                    cc=cc+1;
            end

            TOO_1 = TO_Objectives(reader,mesh,bcinit);
            for i=1:n_con
                if (reader.TopOpt_ConstraintName=="Power")
                        eval_fun=@(reader,mesh,solver) CalculatePower(reader,mesh,solver);
                        [FD_vals(cc), err, cc_FD] = Finite_Differences_bc( reader, mesh, bcinit, solver, eval_fun);
                elseif (reader.TopOpt_ConstraintName=="Volume")
                        FD_vals(cc)=0;
                end
                    cc=cc+1;
            end

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [diffFEM, FD_vals] = run_SingleFEM_diff(obj, filepath)
                % filepath="Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect.txt";
                %reader = InputReader("Benchmarks/Elements/Benchmark1_HexLinear/input_Benchmark1_LINEARHEX_PARAM.txt");
                %reader = InputReader("Benchmarks/Elements/Benchmark_HexSerendipity/input_NonLinCouplSEffect.txt");
                reader = InputReader(filepath);
                fprintf('Initialized InputReader with filename: %s\n', filepath);
                mesh = Mesh(reader);
                TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
                if isempty(reader.TopOpt_Initial_x)
                    reader.TopOpt_Initial_x=1;
                else
                    mesh.elements_density(TOEL)=ones(length(mesh.elements_density(TOEL)),1)*reader.TopOpt_Initial_x;
                end
                fprintf('Initialized Mesh\n');
                bcinit = BCInit(reader, mesh);
                fprintf('Initialized Loads\n');
                solver = Solver(mesh, bcinit);
                solver.runNewtonRaphson(reader, mesh, bcinit);
                fprintf('Postprocessing\n');
                TOO = TO_Objectives(reader,mesh,bcinit);
                TOO.CalculateObjective(reader,mesh,solver)
                TOC = TO_Constraints(reader,mesh,bcinit);
                TOC.CalculateConstraint(reader,mesh,solver);
                ncon = length(reader.TopOpt_ConstraintValue);
                %bench = ThermoelectricBenchmarks();
        
                %TOO_1 = TO_Objectives(reader,mesh,bcinit);
                %eval_fun=@(reader,mesh,solver) TOO_1.fval_AverageTemp(reader,mesh,solver);
                %[FD_vals(1), err] = obj.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); % OK
        
                %TOO_1 = TO_Objectives(reader,mesh,bcinit);
                %eval_fun=@(reader,mesh,solver) TOO_1.fval_AverageTemp(reader,mesh,solver);
                %[FD_vals(2), err] = obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, 6,1); %  OK. Needs to adjust the value of the constraint (6=default) in the overall bc values vector in reader
        
                %TOC_1 = TO_Constraints(reader,mesh,bcinit);
                %eval_fun=@(reader,mesh,solver) TOC_1.fval_Volume(reader,mesh,solver,2); % matters which index is given!!!
                %[FD_vals(3), err] = obj.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); % OK
        
                %TOC_1 = TO_Constraints(reader,mesh,bcinit);
                %eval_fun=@(reader,mesh,solver) TOC_1.fval_Power(reader,mesh,solver,1);
                %[FD_vals(4), err] = obj.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); %  OK

                %TOC_1 = TO_Constraints(reader,mesh,bcinit);
                %eval_fun=@(reader,mesh,solver) TOC_1.fval_Power(reader,mesh,solver,1);
                %[FD_vals(5), err] = obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, 6,1); % OK

                TOC_1 = TO_Constraints(reader,mesh,bcinit);
                eval_fun=@(reader,mesh,solver) TOC_1.fval_Stress(reader,mesh,solver,3);
                % OK for all Penaly 1 and xx = 1, FD tolerance to 1e-6
                [FD_vals(6), err] = obj.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); % NOT OK 23.10.19

                TOC_1 = TO_Constraints(reader,mesh,bcinit);
                eval_fun=@(reader,mesh,solver) TOC_1.fval_Stress(reader,mesh,solver,3);
                [FD_vals(7), err] = obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, 6,1); % OK

                diffFEM=zeros(1+ncon,length(TOO.TOEL)+1);
                diffFEM(1,:)=TOO.dfdx;
                diffFEM(2:end,:)=TOC.dfdx;
        end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [diffFEM, FD_vals] = run_SingleFEM_diff_TEC(obj, filepath)
                % filepath="Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect.txt";
                %reader = InputReader("Benchmarks/Elements/Benchmark1_HexLinear/input_Benchmark1_LINEARHEX_PARAM.txt");
                %reader = InputReader("Benchmarks/Elements/Benchmark_HexSerendipity/input_NonLinCouplSEffect.txt");
                reader = InputReader(filepath);
                fprintf('Initialized InputReader with filename: %s\n', filepath);
                mesh = Mesh(reader);
                TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
                if isempty(reader.TopOpt_Initial_x)
                    reader.TopOpt_Initial_x=1;
                else
                    mesh.elements_density(TOEL)=ones(length(mesh.elements_density(TOEL)),1)*reader.TopOpt_Initial_x;
                end
                fprintf('Initialized Mesh\n');
                bcinit = BCInit(reader, mesh);
                fprintf('Initialized Loads\n');
                solver = Solver(mesh, bcinit);
                solver.runNewtonRaphson(reader, mesh, bcinit);
                fprintf('Postprocessing\n');
                TOO = TO_Objectives(reader,mesh,bcinit);
                TOO.CalculateObjective(reader,mesh,solver)
                TOC = TO_Constraints(reader,mesh,bcinit);
                TOC.CalculateConstraint(reader,mesh,solver);
                ncon = length(reader.TopOpt_ConstraintValue);
                %bench = ThermoelectricBenchmarks();
        
                TOO_1 = TO_Objectives(reader,mesh,bcinit);
                eval_fun=@(reader,mesh,solver) TOO_1.fval_AverageTemp(reader,mesh,solver);
                [FD_vals(1), err] = obj.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); % OK
        
                TOO_1 = TO_Objectives(reader,mesh,bcinit);
                eval_fun=@(reader,mesh,solver) TOO_1.fval_AverageTemp(reader,mesh,solver);
                [FD_vals(2), err] = obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, 8,1); %  OK. Needs to adjust the value of the constraint (6=default) in the overall bc values vector in reader
        
                TOC_1 = TO_Constraints(reader,mesh,bcinit);
                eval_fun=@(reader,mesh,solver) TOC_1.fval_Volume(reader,mesh,solver,2); % matters which index is given!!!
                [FD_vals(3), err] = obj.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); % OK
        
                TOC_1 = TO_Constraints(reader,mesh,bcinit);
                eval_fun=@(reader,mesh,solver) TOC_1.fval_Power(reader,mesh,solver,1);
                [FD_vals(4), err] = obj.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,126); %  OK

                TOC_1 = TO_Constraints(reader,mesh,bcinit);
                eval_fun=@(reader,mesh,solver) TOC_1.fval_Power(reader,mesh,solver,1);
                [FD_vals(5), err] = obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, 8,1); % OK

                diffFEM=zeros(1+ncon,length(TOO.TOEL)+1);
                diffFEM(1,:)=TOO.dfdx;
                diffFEM(2:end,:)=TOC.dfdx;
        end
    end
end


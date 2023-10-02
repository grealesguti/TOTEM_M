classdef ThermoelectricBenchmarks < handle
    %THERMOELECTRICBENCHMARKS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Benchmarks
        BenchmarksFunctions
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

            for i=1:length(obj.BenchmarksFunctions)
                if(obj.BenchmarksFunctions(i)==1)
                    Perez_Aparicio_LinearUncoupledSeebeck(obj,obj.Benchmarks{i},i)
                elseif(obj.BenchmarksFunctions(i)==2)
                    Perez_Aparicio_CoupledSeebeck(obj,obj.Benchmarks{i},i)
                elseif(obj.BenchmarksFunctions(i)==3)
                    Perez_Aparicio_CoupledPeltier(obj,obj.Benchmarks{i},i)
                end
            end

            Benchmark_sensitivities ="Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect.txt";
            %[diffFEM, FD_vals] = obj.run_FEM_diff( Benchmark_sensitivities)

            

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % https://link.springer.com/article/10.1007/s00466-006-0080-7
        function [xv,Tx,Vx] = Perez_Aparicio_LinearUncoupledSeebeck(~,filename,index)

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

                figure(index)
                hold on
                postprocessing.Benchmark_T_PLOT_axis(index,solver,2)
                L=1.524e-3;
                steps=100;
                x_step=L/steps;
                xv=zeros(steps,1);
                Tx=zeros(steps,1);
                Vx=zeros(steps,1);
                for i=1:100
                    xv(i)=x_step*(i-1);
                    Tx(i)=412-sqrt(169587-1.31e7*xv(i))+273.15;
                    Vx(i)=sqrt(0.018-1.31*xv(i))+2.36*xv(i)-0.13;
                end
                plot(xv,Tx)

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % https://link.springer.com/article/10.1007/s00466-006-0080-7
        function [xv,Tx,Vx] = Perez_Aparicio_CoupledSeebeck(~,filename,index)

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

                figure(index)
                hold on
                postprocessing.Benchmark_T_PLOT_axis(index,solver,2)
                L=1.524e-3;
                steps=100;
                x_step=L/steps;
                xv=zeros(steps,1);
                Tx=zeros(steps,1);
                Vx=zeros(steps,1);
                for i=1:100
                    xv(i)=x_step*(i-1);
                    Tx(i)=3.794e7*xv(i)*(1.524e-3-xv(i))+298.15;
                    Vx(i)=5.788e-2-4.913*10*xv(i)+7.315e3*xv(i)^2;
                end
                j=3.199e6;
                Power = 5.788e-2*(j*0.0014^2);
                Power_FEM = CalculatePower(reader,mesh,solver);
                Power_Bench = abs(Power-Power_FEM);
                fprintf('Power benchmark, Analytical: %s FEM: %s\n', [Power,Power_FEM]);
                plot(xv,Tx)

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
            xx_iter = xx * 0.9;
            mesh_1 = Mesh(reader);
            
            % Calculate the initial objective function value
            value_0 = eval_fun(reader, mesh, solver);
            % Initialize history arrays to store err and diff values
                err_history = [];
                diff_history = [];
    
            % Set convergence tolerance and counters
            Tol = 1e-3;
            cc = 1;
            err = 1;
            diff_old = 0; % Initialize diff_old
            max_iterations = 15;
            
            while err > Tol && cc <= max_iterations
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
        function [diff, err, cc] = Finite_Differences_bc(~, reader, mesh, bcinit, solver, eval_fun, index)
            % Calculate the initial objective function value
            value_0 = eval_fun(reader, mesh, solver);

            % Get the test element for density evaluation
            bcval = reader.bcval(index);
            bcval_iter=bcval_0*0.95;

            % Initialize density and iterate variables

            

            % Set convergence tolerance and counters
            Tol = 1e-6;
            cc = 1;
            err = 1;
            diff_old = 0; % Initialize diff_old
            
            while err > Tol
                % Calculate finite difference and update density
                bcval_iter = (bcval - bcval_iter) / 2;
                reader_1.bcval(index)=bcval_iter;
                bcinit_1 = BCInit(reader, mesh);

                
                % Create a new solver and run the Newton-Raphson method
                solver_1 = Solver(mesh,bcinit_1);
                solver_1.runNewtonRaphson(reader_1, mesh, bcinit_1);
                
                % Calculate the objective function value with the updated density
                value = eval_fun(reader_1, mesh, solver_1);
                
                % Calculate finite difference and check for convergence
                diff = (value - value_0) / (bcval_iter - bcval);
                if cc > 2
                    err = abs((diff_old - diff) / diff_old);
                end
                diff_old = diff;
                
                % Increment the counter
                cc = cc + 1;
            end
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

            TOO_1 = TO_Objectives(reader,mesh,bcinit);
            if (reader.TopOpt_Objective=="AverageTemperature")
                        eval_fun=@(reader,mesh,solver) TOO_1.fval_AverageTemp(reader,mesh,solver);
                        [FD_vals(cc), err, cc_FD] = Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun);
            end
            cc=cc+1;

            FD_vals=zeros(ncon+n_bc,1);
            for i=1:n_con
            TOC_1 = TO_Constraints(reader,mesh,bcinit);

                if (reader.TopOpt_ConstraintName=="Power")
                        eval_fun=@(reader,mesh,solver) CalculatePower(reader,mesh,solver);
                        [FD_vals(cc), err, cc_FD] = Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun);
                elseif (reader.TopOpt_ConstraintName=="Volume")
                        eval_fun=@(reader,mesh,solver) TOC_1.fval_Volume(reader,mesh,solver,1);
                        [FD_vals(cc), err, cc_FD] = Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun);
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

    
    end
end


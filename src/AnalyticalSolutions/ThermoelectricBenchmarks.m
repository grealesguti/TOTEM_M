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
        bench
        errhist
    end

    methods
        function obj = ThermoelectricBenchmarks()

            close all;%clear all; clc;
            srcPath = 'src';
            utilsPath = 'utils';
            obj.bench = 2;
            % Add the generated paths to the MATLAB path
            addpath(srcPath);
            addpath(utilsPath);
            if obj.bench == 0
                % THERMOELECTRICBENCHMARKS Construct an instance of this class

                Benchmark_Perez_Aparicio_LinUncoupSEffect_HexLinear = ...
                    "Benchmarks/Elements/Benchmark_HexLinear/input_LinearUncoupledSeebeck.txt"; % Hex DIM-3 nodes-8
                Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexLinear = ...
                    "Benchmarks/Elements/Benchmark_HexLinear/input_NonLinCouplSEffect.txt";
                Benchmark_Perez_Aparicio_LinUncoupSEffect_HexLinear_Q = ...
                    "Benchmarks/Elements/Benchmark_HexLinear/input_LinearUncoupledSeebeck_Q.txt";
                Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity = ...
                    "Benchmarks/Elements/Benchmark_HexSerendipity/input_LinearUncoupledSeebeck.txt"; % Hex DIM-3 nodes-20
                Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexSerendipity = ...
                    "Benchmarks/Elements/Benchmark_HexSerendipity/input_NonLinCouplSEffect.txt";
                Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity_Q = ...
                    "Benchmarks/Elements/Benchmark_HexSerendipity/input_LinearUncoupledSeebeck_Q.txt";

                obj.Benchmarks = {
                    Benchmark_Perez_Aparicio_LinUncoupSEffect_HexLinear;
                    Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexLinear;
                    Benchmark_Perez_Aparicio_LinUncoupSEffect_HexLinear_Q;
                    Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity;
                    Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexSerendipity;
                    Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity_Q
                    };

                obj.BenchmarksFunctions = [1, 2, 3, 1, 2, 3];

            elseif obj.bench == 1
                % THERMOELECTRICBENCHMARKS Construct an instance of this class

                Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity = ...
                    "Benchmarks/Elements/Benchmark_HexSerendipity/input_LinearUncoupledSeebeck.txt";
                Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexSerendipity = ...
                    "Benchmarks/Elements/Benchmark_HexSerendipity/input_NonLinCouplSEffect.txt";
                Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity_Q = ...
                    "Benchmarks/Elements/Benchmark_HexSerendipity/input_LinearUncoupledSeebeck_Q.txt";

                obj.Benchmarks = {
                    Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity;
                    Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexSerendipity;
                    Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity_Q
                    };

                obj.BenchmarksFunctions = [1, 2, 3];

            elseif obj.bench == 2
                % THERMOELECTRICBENCHMARKS Construct an instance of this class

                Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexSerendipity = ...
                    "Benchmarks/Elements/Benchmark_HexSerendipity/input_NonLinCouplSEffect.txt";

                obj.Benchmarks = {
                    Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexSerendipity
                    };

                obj.BenchmarksFunctions = [2];

            end


            %obj.runPhysicsBenchmarks()

            %Benchmark_sensitivities ="Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect.txt";
            %[obj.diffFEM_ctemat, obj.FD_vals_ctemat] = obj.run_SingleFEM_diff(Benchmark_sensitivities);
            %[diffFEM_ctemat_TEC, FD_vals_ctemat_TEC] = thb.run_SingleFEM_diff("Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect.txt");
            
            %%%% next one is the basic!!! best to run
            %Benchmark_sensitivities ="Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect_nonlinmat.txt";
            %idxvoltage=7;idxelement=1;
            %[obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);
            %%[diffFEM_ctemat_TEC, FD_vals_ctemat_TEC] = thb.run_SingleFEM_diff("Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect_nonlinmat.txt");

           % Benchmark_sensitivities ="TECTO/input_TECTO_StressConstrained_cte.txt";
            %            idxvoltage=9;
            %[obj.diffFEM_ctemat_TEC, obj.FD_vals_ctemat_TEC] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage);
            %[diffFEM_ctemat_TEC, FD_vals_ctemat_TEC] = thb.run_SingleFEM_diff_TEC("TECTO/input_TECTO_StressConstrained_cte.txt");

            %Benchmark_sensitivities ="TECTO/input_TECTO_StressConstrained_noncte.txt";
            %idxvoltage=7;
            %[obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage);
             %Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D_reviewalphabench.txt";
%              Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D_reviewctebench_volumepenalty.txt";
%                         idxvoltage=10;idxelement=5;
%              %Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal3D_reviewAir.txt";
%             %idxvoltage=12;
%             [obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);
%              Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D_reviewctebench_seebeck.txt";
%                         idxvoltage=10;idxelement=5;
%              %Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal3D_reviewAir.txt";
%             %idxvoltage=12;
%             [obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);
%              Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D_reviewctebench_seebeckpenalty.txt";
%                         idxvoltage=10;idxelement=5;
%              %Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal3D_reviewAir.txt";
%             %idxvoltage=12;
%             [obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);
%              Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D_reviewctebench_nonlin.txt";
%                         idxvoltage=10;idxelement=5;
             %Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal3D_reviewAir.txt";
            %idxvoltage=12;
%            [obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);

            %[obj.diffFEM_nonlinmat_TEC, obj.diffFEM_nonlinmat_TEC] = obj.run_SingleFEM_diff(Benchmark_sensitivities);
            %[obj.diffFEM_nonlinmat_TEC, obj.diffFEM_nonlinmat_TEC] = obj.run_SingleFEM_diff_TEC(Benchmark_sensitivities);
           
             Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D_reviewctebench_contactsmat.txt";
                        idxvoltage=10;idxelement=5;
             %Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal3D_reviewAir.txt";
            %idxvoltage=12;
            [obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);

             Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D_reviewctebench_penaltyfinal.txt";
                        idxvoltage=10;idxelement=5;
             %Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal3D_reviewAir.txt";
            %idxvoltage=12;
            [obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);
             Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D_reviewctebench_nonlinallpenalty3.txt";
                        idxvoltage=10;idxelement=5;
             %Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal3D_reviewAir.txt";
            %idxvoltage=12;
            [obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);

            [obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);
             Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D_reviewctebench_nonlinallpenalty5.txt";
                        idxvoltage=10;idxelement=5;
             %Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal3D_reviewAir.txt";
            %idxvoltage=12;
            [obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);


            [obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);
             Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D_reviewctebench_seebeckmintemp.txt";
                        idxvoltage=10;idxelement=5;
             %Benchmark_sensitivities ="TECTO/input_TECTO_Thermoel_Serend_240124_Journal3D_reviewAir.txt";
            %idxvoltage=12;
            [obj.diffFEM_nonlinmat, obj.FD_vals_nonlinmat, obj.errhist] = obj.run_SingleFEM_diff(Benchmark_sensitivities,idxvoltage,idxelement);



        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function []=runPhysicsBenchmarks(obj)
            for i=1:length(obj.BenchmarksFunctions)
                if(obj.BenchmarksFunctions(i)==1)
                    Perez_Aparicio_LinearUncoupledSeebeck(obj,obj.Benchmarks{i},i)
                elseif(obj.BenchmarksFunctions(i)==2)
                    Perez_Aparicio_CoupledSeebeck(obj,obj.Benchmarks{i},i)
                elseif(obj.BenchmarksFunctions(i)==3)
                    Perez_Aparicio_LinearUncoupledSeebeck_Q(obj,obj.Benchmarks{i},i)
                    %Perez_Aparicio_CoupledPeltier(obj,obj.Benchmarks{i},i)
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
        function [xv,Tx,Vx,Ux] = Perez_Aparicio_LinearUncoupledSeebeck_Q(~,filename,index)

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
            subplot(1, 3,1);
            hold on
            TFEM = postprocessing.Benchmark_T_PLOT_axis(index,solver,2);
            subplot(1, 3,2);
            hold on
            VFEM = postprocessing.Benchmark_V_PLOT_axis(index,solver,2);
            subplot(1, 3,3);
            hold on
            UFEM = postprocessing.Benchmark_U_PLOT_axis(index,solver,2);

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
        % https://link.springer.com/article/10.1007/s00466-006-0080-7
        function [xv,Tx,Vx,Ux,Power_Bench] = GR_HeatInjection(~,filename,index)

            %reader = InputReader(filename);
            %fprintf('Initialized InputReader with filename: %s\n', filename);
            %mesh = Mesh(reader);
            %fprintf('Initialized Mesh\n');
            %bcinit = BCInit(reader, mesh);
            %fprintf('Initialized Loads\n');
            %solver = Solver(mesh, bcinit);
            %solver.runNewtonRaphson(reader, mesh, bcinit);
            %fprintf('Postprocessing\n');
            %postprocessing = Postprocessing();
            %postprocessing.initVTK(reader,mesh);


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
function [diff, err, cc, diff_history, err_history, epsilon_history,epsilon_vector] = Finite_Differences_DensityElement(obj, reader, mesh, bcinit, solver, eval_fun, idx, do_plot, stop_at_eps, min_perturbation, max_iterations, Tol)
    % Optional input defaults
    if nargin < 8, do_plot = true; end
    if nargin < 9, stop_at_eps = false; end
    if nargin < 10, min_perturbation = 1e-4; end
    if nargin < 11, max_iterations = 15; end
    if nargin < 12, Tol = 1e-4; end

    % Get the test element for density evaluation
    TOEL = mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
    test_element = TOEL(idx);

    % Initial density and perturbed density
    xx = mesh.elements_density(test_element);

    % Automatically determine the order of xx
    order_of_xx = floor(log10(xx));  % Find the order of xx by using log10 and floor it to get the integer part
    
    % Set initial epsilon based on the order of xx
    if order_of_xx == -1  % For values like 0.5
        epsilon = 0.01;
    elseif order_of_xx == 0  % For values like 1, or close to it
        epsilon = 0.1;
    elseif order_of_xx > 0  % For values larger than 1 (e.g., 10, 100)
        epsilon = 0.1 * 10^(-order_of_xx);  % Dynamically scale epsilon based on order of xx
    else  % For values smaller than 1 but greater than 0 (e.g., 0.01, 0.1)
        epsilon = 0.01 * 10^(-order_of_xx);  % Adjust epsilon based on order of xx
    end

    % Pre-calculate all epsilon values based on a constant step (e.g., halving epsilon)
    epsilon_vector = zeros(1, max_iterations);
    epsilon_vector(1) = epsilon;  % Initialize the first epsilon value
    for i = 2:max_iterations
        epsilon_vector(i) = epsilon_vector(i-1) / (2);  % Constant step reduction
    end


    % Initial objective function value
    value_0 = eval_fun(reader, mesh, solver);

    % History arrays
    err_history = zeros(1, max_iterations);
    diff_history = zeros(1, max_iterations);
    epsilon_history = zeros(1, max_iterations);

    % Convergence setup
    cc = 1;
    err = 100;
    diff_old = 0;

    while cc <= max_iterations
        % Stopping criteria
        if stop_at_eps && epsilon < min_perturbation
            fprintf('Stopping: Relative Epsilon = %.2e < %.2e\n', epsilon, min_perturbation);
            break;
        end
        if ~stop_at_eps && err <= Tol && cc > 2
            break;
        end

        % Perturbation update using the set epsilon
        xx_iter = xx + epsilon;  % Apply the perturbation to xx

        mesh_1 = Mesh(reader);
        if isempty(reader.TopOpt_Initial_x)
            reader.TopOpt_Initial_x = 1;
        else
            mesh_1.elements_density(TOEL) = ones(length(mesh_1.elements_density(TOEL)), 1) * reader.TopOpt_Initial_x;
        end
        mesh_1.elements_density(test_element) = xx_iter;

        bcinit_1 = BCInit(reader, mesh_1);
        solver_1 = Solver(mesh_1, bcinit_1);
        solver_1.runNewtonRaphson(reader, mesh_1, bcinit_1);

        value = eval_fun(reader, mesh_1, solver_1);
        diff = (value - value_0) / (xx_iter - xx);

        if cc > 1
            err = abs((diff_old - diff) / diff_old);
        end
        diff_old = diff;

        % Print info
        fprintf('FD Iteration: %d, Epsilon: %.2e, Error: %e, value: %e, diff: %e, xx_iter: %e\n', ...
            cc, epsilon, err, value, diff, xx_iter);

        % Store history
        err_history(cc) = err;
        diff_history(cc) = diff;
        epsilon_history(cc) = epsilon;

        % Reduce epsilon by half an order (divide by √10)
        epsilon = epsilon / (2);

        cc = cc + 1;
    end

    % Trim history arrays
    err_history = err_history(1:cc-1);
    diff_history = diff_history(1:cc-1);
    epsilon_history = epsilon_history(1:cc-1);

    % Plot if requested
    if do_plot
        figure;
        subplot(3, 1, 1);
        plot(1:length(err_history), err_history);
        xlabel('Iteration');
        ylabel('Error');
        title('Convergence Plot - Error');

        subplot(3, 1, 2);
        plot(1:length(diff_history), diff_history);
        xlabel('Iteration');
        ylabel('Difference');
        title('Convergence Plot - Difference');

        subplot(3, 1, 3);
        plot(1:length(epsilon_history), epsilon_history);
        xlabel('Iteration');
        ylabel('Relative Perturbation');
        title('Relative Perturbation History');
    end
end
function [diff, err, cc, diff_history, err_history, epsilon_history, epsilon_vector] = Finite_Differences_bc(obj, filename, reader, mesh, solver, eval_fun, index, index_TO, do_plot, stop_at_eps, min_perturbation, max_iterations, Tol)
    % Optional input defaults
    if nargin < 9, do_plot = true; end
    if nargin < 10, stop_at_eps = false; end
    if nargin < 11, min_perturbation = 1e-4; end
    if nargin < 12, max_iterations = 15; end
    if nargin < 13, Tol = 1e-4; end

    % Initial objective function value
    value_0 = eval_fun(reader, mesh, solver);
    TOEL = mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);

    % Setup
    bcval = reader.bcval(index);
    bcval_iter = bcval * 0.95;
    reader_1 = InputReader(filename);
    bcmax = reader.TObcmaxval(index_TO);
    bcmin = reader.TObcminval(index_TO);
    dbc = bcmax - bcmin;
    xbcval0 = (bcval - bcmin) / dbc;
    xbcval_iter = xbcval0 * 0.95;

    % Calculate the order of xbcval0 (log scale)
    order_of_xbc = floor(log10(xbcval0));  % Find the order of xbcval0 by using log10

    % Set initial epsilon based on the order of xbcval0
    if order_of_xbc == -1  % For values like 0.5
        epsilon = 0.01;
    elseif order_of_xbc == 0  % For values like 1, or close to it
        epsilon = 0.1;
    elseif order_of_xbc > 0  % For values larger than 1 (e.g., 10, 100)
        epsilon = 0.1 * 10^(-order_of_xbc);  % Dynamically scale epsilon based on order of xbcval0
    else  % For values smaller than 1 but greater than 0 (e.g., 0.01, 0.1)
        epsilon = 0.01 * 10^(-order_of_xbc);  % Adjust epsilon based on order of xbcval0
    end

    % Pre-calculate all epsilon values based on a constant step (e.g., halving epsilon)
    epsilon_vector = zeros(1, max_iterations);
    epsilon_vector(1) = epsilon;  % Initialize the first epsilon value
    for i = 2:max_iterations
        epsilon_vector(i) = epsilon_vector(i-1) / (2);  % Constant step reduction
    end

    % Convergence setup
    cc = 1;
    err = 1;
    diff_old = 0;

    % History arrays
    err_history = zeros(1, max_iterations);
    diff_history = zeros(1, max_iterations);
    epsilon_history = zeros(1, max_iterations);

    while cc <= max_iterations
        % Get the epsilon value for this iteration
        epsilon = epsilon_vector(cc);

        % Calculate epsilon (relative perturbation)
        xbcval_iter = xbcval0 + epsilon;  % Apply the perturbation to xbcval0

        % Exit condition: epsilon-based stopping
        if stop_at_eps && epsilon < min_perturbation
            fprintf('Stopping: Relative Epsilon = %.2e is below threshold %.2e\n', epsilon, min_perturbation);
            break;
        end

        % Exit condition: convergence-based stopping
        if ~stop_at_eps && err <= Tol && cc > 2
            break;
        end

        % Midpoint update of perturbation
        xbcval_iter = (xbcval0 + xbcval_iter) / 2;
        bcvalue = bcmin + dbc * xbcval_iter;
        reader_1.bcval(index) = bcvalue;

        mesh_1 = Mesh(reader_1);
        mesh_1.elements_density = mesh.elements_density;
        bcinit_1 = BCInit(reader_1, mesh_1);

        solver_1 = Solver(mesh_1, bcinit_1);
        solver_1.runNewtonRaphson(reader_1, mesh_1, bcinit_1);

        value = eval_fun(reader_1, mesh_1, solver_1);
        diff = (value - value_0) / (xbcval_iter - xbcval0);

        if cc > 2
            err = abs((diff_old - diff) / diff_old);
        end
        diff_old = diff;

        fprintf('FD Iteration: %d, Rel. Eps: %.2e, Error: %e, value: %e, diff: %e\n', cc, epsilon, err, value, diff);

        % Store history
        err_history(cc) = err;
        diff_history(cc) = diff;
        epsilon_history(cc) = epsilon;

        cc = cc + 1;
    end

    % Trim arrays
    err_history = err_history(1:cc-1);
    diff_history = diff_history(1:cc-1);
    epsilon_history = epsilon_history(1:cc-1);
    epsilon_vector = epsilon_vector(1:cc-1);  % Trim epsilon_vector to match the number of iterations

    % Plot if requested
    if do_plot
        figure;
        subplot(3, 1, 1);
        plot(1:length(err_history), err_history);
        xlabel('Iteration');
        ylabel('Error');
        title('Convergence Plot - Error');

        subplot(3, 1, 2);
        plot(1:length(diff_history), diff_history);
        xlabel('Iteration');
        ylabel('Difference');
        title('Convergence Plot - Difference');

        subplot(3, 1, 3);
        plot(1:length(epsilon_history), epsilon_history);
        xlabel('Iteration');
        ylabel('Epsilon');
        title('Perturbation (Epsilon) History');
    end
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
        function [diffFEM, FD_vals] = run_SingleFEM_diff_2D(obj, filepath)
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
            [FD_vals(2), err] = obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, 6,1); %  OK. Needs to adjust the value of the constraint (6=default) in the overall bc values vector in reader

            TOC_1 = TO_Constraints(reader,mesh,bcinit);
            eval_fun=@(reader,mesh,solver) TOC_1.fval_Volume(reader,mesh,solver,2); % matters which index is given!!!
            [FD_vals(3), err] = obj.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); % OK

            TOC_1 = TO_Constraints(reader,mesh,bcinit);
            eval_fun=@(reader,mesh,solver) TOC_1.fval_Power(reader,mesh,solver,1);
            [FD_vals(4), err] = obj.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); %  OK

            TOC_1 = TO_Constraints(reader,mesh,bcinit);
            eval_fun=@(reader,mesh,solver) TOC_1.fval_Power(reader,mesh,solver,1);
            [FD_vals(5), err] = obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, 6,1); % OK

            TOC_1 = TO_Constraints(reader,mesh,bcinit);
            eval_fun=@(reader,mesh,solver) TOC_1.fval_Stress(reader,mesh,solver,3);
            [FD_vals(6), err] = obj.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); % OK 23.10.24

            TOC_1 = TO_Constraints(reader,mesh,bcinit);
            eval_fun=@(reader,mesh,solver) TOC_1.fval_Stress(reader,mesh,solver,3);
            [FD_vals(7), err] = obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, 6,1); % OK

            diffFEM=zeros(1+ncon,length(TOO.TOEL)+1);
            diffFEM(1,:)=TOO.dfdx;
            diffFEM(2:end,:)=TOC.dfdx;
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
        function [diffFEM, FD_vals, errorhisttofem] = run_SingleFEM_diff(obj, filepath, idxvoltage,idxelement, do_plot)
            if nargin < 5
                do_plot = true;
            end
            reader = InputReader(filepath);
            fprintf('Initialized InputReader with filename: %s\n', filepath);
            mesh = Mesh(reader);
            TOEL = mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            if isempty(reader.TopOpt_Initial_x)
                reader.TopOpt_Initial_x = 1;
            else
                mesh.elements_density(TOEL) = ones(length(mesh.elements_density(TOEL)), 1) * reader.TopOpt_Initial_x;
            end
            fprintf('Initialized Mesh\n');
            bcinit = BCInit(reader, mesh);
            fprintf('Initialized Loads\n');
            solver = Solver(mesh, bcinit);
            solver.runNewtonRaphson(reader, mesh, bcinit);
            fprintf('Postprocessing\n');

            TOO = TO_Objectives(reader, mesh, bcinit);
            TOO.CalculateObjective(reader, mesh, solver);
            TOC = TO_Constraints(reader, mesh, bcinit);
            TOC.CalculateConstraint(reader, mesh, solver);
            ncon = length(reader.TopOpt_ConstraintValue);

            FD_vals = struct();
            errorhisttofem = struct();
% Set parameters
FD_Tol = 1e-8;
FD_max_iter = 20;
epsmin = 1e-8;
cc=0;
% 
% 
% %             % Objective: AverageTemp (Density Element)
%              TOO_1 = TO_Objectives(reader, mesh, bcinit); cc = 1;
%              eval_fun = @(reader, mesh, solver) TOO_1.fval_AverageTemp(reader, mesh, solver);
%             [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, FD_vals(cc).diff_history, FD_vals(cc).err_history, FD_vals(cc).eps_history,FD_vals(cc).eps_vect] = ...
%                 obj.Finite_Differences_DensityElement(reader, mesh, bcinit, solver, eval_fun, 1, false, true, epsmin, FD_max_iter, FD_Tol);
%             errorhisttofem.Trho = abs((TOO.dfdx(1) - FD_vals(cc).diff_history) ./ TOO.dfdx(1));
%             if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Trho, 'Trho', cc == 1); end
% 
%             %Objective: AverageTemp (BC)
%             TOO_1 = TO_Objectives(reader, mesh, bcinit); cc = cc + 1;
%             eval_fun = @(reader, mesh, solver) TOO_1.fval_AverageTemp(reader, mesh, solver);
%             [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, FD_vals(cc).diff_history, FD_vals(cc).err_history, FD_vals(cc).eps_history,FD_vals(cc).eps_vect] = ...
%                 obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, idxvoltage, idxelement, false, true, epsmin, FD_max_iter, FD_Tol);
%             errorhisttofem.Tbc = abs((TOO.dfdx(end) - FD_vals(cc).diff_history) ./ TOO.dfdx(end));
% if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Tbc, 'Tbc', cc == 1); end
% 
% %             % Constraint: Volume (Density Element)
              TOC_1 = TO_Constraints(reader, mesh, bcinit); cc = cc + 1;
%              eval_fun = @(reader, mesh, solver) TOC_1.fval_Volume(reader, mesh, solver, 2);
%             [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, FD_vals(cc).diff_history, FD_vals(cc).err_history, FD_vals(cc).eps_history,FD_vals(cc).eps_vect] = ...
%                 obj.Finite_Differences_DensityElement(reader, mesh, bcinit, solver, eval_fun, 1, false, true, epsmin, FD_max_iter, FD_Tol);
%             errorhisttofem.Vrho = abs((TOC.dfdx(2,1) - FD_vals(cc).diff_history) ./ TOC.dfdx(2,1));
% if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Vrho, 'Vrho', cc == 1); end

%             %Constraint: Power (Density Element)
%             eval_fun = @(reader, mesh, solver) TOC_1.fval_Power(reader, mesh, solver, 1); cc = cc + 1;
%             [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, FD_vals(cc).diff_history, FD_vals(cc).err_history, FD_vals(cc).eps_history,FD_vals(cc).eps_vect] = ...
%                 obj.Finite_Differences_DensityElement(reader, mesh, bcinit, solver, eval_fun, idxelement, false, true, epsmin, FD_max_iter, FD_Tol);
%             errorhisttofem.Prho = abs((TOC.dfdx(1,idxelement) - FD_vals(cc).diff_history) ./ TOC.dfdx(1,idxelement));
% if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Prho, 'Prho', cc == 1); end
% 
%             % Constraint: Power (BC)
%             eval_fun = @(reader, mesh, solver) TOC_1.fval_Power(reader, mesh, solver, 1); cc = cc + 1;
%             [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, FD_vals(cc).diff_history, FD_vals(cc).err_history, FD_vals(cc).eps_history,FD_vals(cc).eps_vect] = ...
%                 obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, idxvoltage, 1, false, true, epsmin, FD_max_iter, FD_Tol);
%             errorhisttofem.Pbc = abs((TOC.dfdx(1,end) - FD_vals(cc).diff_history) ./ TOC.dfdx(1,end));
% if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Pbc, 'Pbc', cc == 1); end

            % Constraint: Stress (Density Element)
            eval_fun = @(reader, mesh, solver) TOC_1.fval_Stress(reader, mesh, solver, 3); cc = cc + 1;
            [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, FD_vals(cc).diff_history, FD_vals(cc).err_history, FD_vals(cc).eps_history,FD_vals(cc).eps_vect] = ...
                obj.Finite_Differences_DensityElement(reader, mesh, bcinit, solver, eval_fun, idxelement, false, true, epsmin, FD_max_iter, FD_Tol);
            errorhisttofem.Srho = abs((TOC.dfdx(3,idxelement) - FD_vals(cc).diff_history) ./ TOC.dfdx(3,idxelement));
if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Srho, 'Srho', cc == 1); end

            % Constraint: Stress (BC)
            eval_fun = @(reader, mesh, solver) TOC_1.fval_Stress(reader, mesh, solver, 3); cc = cc + 1;
            [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, FD_vals(cc).diff_history, FD_vals(cc).err_history, FD_vals(cc).eps_history,FD_vals(cc).eps_vect] = ...
                obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, idxvoltage, 1, false, true, epsmin, FD_max_iter, FD_Tol);
            errorhisttofem.Sbc = abs((TOC.dfdx(3,end) - FD_vals(cc).diff_history) ./ TOC.dfdx(3,end));
if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Sbc, 'Sbc', cc == 1); end

            % FEM sensitivities
            diffFEM = zeros(1 + ncon, length(TOO.TOEL) + 1);
            diffFEM(1, :) = TOO.dfdx;
            diffFEM(2:end, :) = TOC.dfdx;

            error_fields = fieldnames(errorhisttofem);

% Optional plotting
if do_plot
   
    % Save the plot as an image
    saveas(gcf, 'relative_error_vs_epsilon_log_plot.png');
end



            % Prepare filename and delete if exists
        csv_filename = fullfile("C:\Archive\Programming\MATLAB\TOTEM_M\Benchmarks", ...
            append("error_vs_epsilon_simple_", datestr(now, 'yyyymmdd'), ".csv"));
            if isfile(csv_filename)
                delete(csv_filename);
            end
            % Step 1: Determine maximum length across all fields
            maxlen = 0;
            for i = 1:length(error_fields)
                field = error_fields{i};
                errors = errorhisttofem.(field)(:);
                epsilons = FD_vals(i).eps_history(:);
                maxlen = max(maxlen, max(length(errors), length(epsilons)));
            end
            
            % Step 2: Initialize table
            T = table();
            
            % Step 3: Fill the table with padded columns
            for i = 1:length(error_fields)
                field = error_fields{i};
                errors = errorhisttofem.(field)(:);
                epsilons = FD_vals(i).eps_history(:);
            
                % Pad with NaNs to the common maxlen
                errors(end+1:maxlen) = NaN;
                epsilons(end+1:maxlen) = NaN;
            
                % Add to table
                T.([field '_error']) = errors;
                T.([field '_epsilon']) = epsilons;
            end
            
            % Step 4: Write to CSV
            writetable(T, csv_filename);
            fprintf('Saved error history and epsilons to %s\n', csv_filename);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [diffFEM, FD_vals, errorhisttofem] = run_SingleFEM_diff_TEC(obj, filepath, do_plot)
    if nargin < 3
        do_plot = true;
    end

    reader = InputReader(filepath);
    fprintf('Initialized InputReader with filename: %s\n', filepath);
    mesh = Mesh(reader);
    TOEL = mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
    
    if isempty(reader.TopOpt_Initial_x)
        reader.TopOpt_Initial_x = 1;
    else
        mesh.elements_density(TOEL) = ones(length(mesh.elements_density(TOEL)), 1) * reader.TopOpt_Initial_x;
    end
    fprintf('Initialized Mesh\n');
    
    bcinit = BCInit(reader, mesh);
    fprintf('Initialized Loads\n');
    
    solver = Solver(mesh, bcinit);
    solver.runNewtonRaphson(reader, mesh, bcinit);
    fprintf('Postprocessing\n');
    
    TOO = TO_Objectives(reader, mesh, bcinit);
    TOO.CalculateObjective(reader, mesh, solver);
    
    TOC = TO_Constraints(reader, mesh, bcinit);
    TOC.CalculateConstraint(reader, mesh, solver);
    
    ncon = length(reader.TopOpt_ConstraintValue);

    % FD parameters
    FD_Tol = 1e-20;
    FD_max_iter = 25;
    epsmin = 1e-9;

    % Testing
    FD_Tol = 1e-4;
    FD_max_iter = 5;
    epsmin = 1e-5;

    FD_vals = struct();
    errorhisttofem = struct();
    cc = 0;

            TOO_1 = TO_Objectives(reader, mesh, bcinit);% cc = 1;
            TOC_1 = TO_Constraints(reader, mesh, bcinit);% cc = cc + 1;

    % ========== Objective: AverageTemp (Density) ==========
%      eval_fun = @(r, m, s) TOO_1.fval_AverageTemp(r, m, s); cc = cc + 1;
%     [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, ...
%      FD_vals(cc).diff_history, FD_vals(cc).err_history, ...
%      FD_vals(cc).eps_history, FD_vals(cc).eps_vect] = ...
%         obj.Finite_Differences_DensityElement(reader, mesh, bcinit, solver, ...
%             eval_fun, 1, false, true, epsmin, FD_max_iter, FD_Tol);
%     errorhisttofem.Trho = abs((TOO.dfdx(1) - FD_vals(cc).diff_history) ./ TOO.dfdx(1));
%     if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Trho, 'Trho', cc == 1); end
% 
%     % ========== Objective: AverageTemp (BC) ==========
%             TOO_1 = TO_Objectives(reader, mesh, bcinit); %cc = 1;
%             TOC_1 = TO_Constraints(reader, mesh, bcinit); %cc = cc + 1;
%             cc = cc + 1;
%     [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, ...
%      FD_vals(cc).diff_history, FD_vals(cc).err_history, ...
%      FD_vals(cc).eps_history, FD_vals(cc).eps_vect] = ...
%         obj.Finite_Differences_bc(filepath, reader, mesh, solver, ...
%             eval_fun, 12, 1, false, true, epsmin, FD_max_iter, FD_Tol);
%     errorhisttofem.Tbc = abs((TOO.dfdx(end) - FD_vals(cc).diff_history) ./ TOO.dfdx(end));
%     if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Tbc, 'Tbc', false); end
% 
%     % ========== Constraint: Volume (Density) ==========
%             TOO_1 = TO_Objectives(reader, mesh, bcinit);% cc = 1;
%             TOC_1 = TO_Constraints(reader, mesh, bcinit);% cc = cc + 1;
%             eval_fun = @(r, m, s) TOC_1.fval_Volume(r, m, s, 2); cc = cc + 1;
%     [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, ...
%      FD_vals(cc).diff_history, FD_vals(cc).err_history, ...
%      FD_vals(cc).eps_history, FD_vals(cc).eps_vect] = ...
%         obj.Finite_Differences_DensityElement(reader, mesh, bcinit, solver, ...
%             eval_fun, 1, false, true, epsmin, FD_max_iter, FD_Tol);
%     errorhisttofem.Vrho = abs((TOC.dfdx(2,1) - FD_vals(cc).diff_history) ./ TOC.dfdx(2,1));
%     if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Vrho, 'Vrho', false); end

    % ========== Constraint: Power (Density) ==========
            TOO_1 = TO_Objectives(reader, mesh, bcinit);% cc = 1;
            TOC_1 = TO_Constraints(reader, mesh, bcinit);% cc = cc + 1;

    eval_fun = @(r, m, s) TOC_1.fval_Power(r, m, s, 1); cc = cc + 1;
    [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, ...
     FD_vals(cc).diff_history, FD_vals(cc).err_history, ...
     FD_vals(cc).eps_history, FD_vals(cc).eps_vect] = ...
        obj.Finite_Differences_DensityElement(reader, mesh, bcinit, solver, ...
            eval_fun, 1, false, true, epsmin, FD_max_iter, FD_Tol);
    errorhisttofem.Prho = abs((TOC.dfdx(1,1) - FD_vals(cc).diff_history) ./ TOC.dfdx(1,1));
    if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Prho, 'Prho', false); end

    % ========== Constraint: Power (BC) ==========
            TOO_1 = TO_Objectives(reader, mesh, bcinit);% cc = 1;
            TOC_1 = TO_Constraints(reader, mesh, bcinit);% cc = cc + 1;

    cc = cc + 1;
    [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, ...
     FD_vals(cc).diff_history, FD_vals(cc).err_history, ...
     FD_vals(cc).eps_history, FD_vals(cc).eps_vect] = ...
        obj.Finite_Differences_bc(filepath, reader, mesh, solver, ...
            eval_fun, 12, 1, false, true, epsmin, FD_max_iter, FD_Tol);
    errorhisttofem.Pbc = abs((TOC.dfdx(1,end) - FD_vals(cc).diff_history) ./ TOC.dfdx(1,end));
    if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Pbc, 'Pbc', false); end


            % Constraint: Stress (Density Element)
            TOO_1 = TO_Objectives(reader, mesh, bcinit);% cc = 1;
            TOC_1 = TO_Constraints(reader, mesh, bcinit);% cc = cc + 1;

             cc = cc + 1;
            eval_fun = @(r, m, s) TOC_1.fval_Stress(reader, mesh, solver, 3);
            [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, FD_vals(cc).diff_history, FD_vals(cc).err_history, FD_vals(cc).eps_history,FD_vals(cc).eps_vect] = ...
                obj.Finite_Differences_DensityElement(reader, mesh, bcinit, solver, eval_fun, 1, false, true, epsmin, FD_max_iter, FD_Tol);
            errorhisttofem.Srho = abs((TOC.dfdx(3,1) - FD_vals(cc).diff_history) ./ TOC.dfdx(3,1));
if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Srho, 'Srho', cc == 1); end

            % Constraint: Stress (BC)
            TOO_1 = TO_Objectives(reader, mesh, bcinit);% cc = 1;
            TOC_1 = TO_Constraints(reader, mesh, bcinit);% cc = cc + 1;

             cc = cc + 1;
            eval_fun = @(reader, mesh, solver) TOC_1.fval_Stress(reader, mesh, solver, 3); %cc = cc + 1;
            [FD_vals(cc).diff, FD_vals(cc).err, FD_vals(cc).cc, FD_vals(cc).diff_history, FD_vals(cc).err_history, FD_vals(cc).eps_history,FD_vals(cc).eps_vect] = ...
                obj.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, 12, 1, false, true, epsmin, FD_max_iter, FD_Tol);
            errorhisttofem.Sbc = abs((TOC.dfdx(3,end) - FD_vals(cc).diff_history) ./ TOC.dfdx(3,end));
if do_plot, obj.plot_relative_error(FD_vals(cc).eps_history, errorhisttofem.Sbc, 'Sbc', cc == 1); end


    % ========== FEM sensitivities ==========
    diffFEM = zeros(1 + ncon, length(TOO.TOEL) + 1);
    diffFEM(1,:) = TOO.dfdx;
    diffFEM(2:end,:) = TOC.dfdx;

    % ========== Export CSV ==========
    error_fields = fieldnames(errorhisttofem);
    csv_filename = fullfile(filepath, 'error_vs_epsilon_TEC.csv');
    if isfile(csv_filename), delete(csv_filename); end
    
    T = table;
    for i = 1:length(error_fields)
        field = error_fields{i};
        errors = errorhisttofem.(field)(:);
        epsilons = FD_vals(i).eps_history(:);
        maxlen = max(length(errors), length(epsilons));
        errors(end+1:maxlen) = NaN;
        epsilons(end+1:maxlen) = NaN;
        T.([field '_error']) = errors;
        T.([field '_epsilon']) = epsilons;
    end
    writetable(T, csv_filename);
    fprintf('Saved error history and epsilons to %s\n', csv_filename);
end


function []=plot_relative_error(obj,epsilons, errors, label_name, is_first_plot)
    persistent hFig
    if isempty(hFig) || ~isvalid(hFig)
        hFig = figure('Name', 'Relative Error to FEM vs Epsilon');
        set(hFig, 'NumberTitle', 'off');
        hold on;
        grid on;
        title('Relative Error to FEM vs Epsilon');
        xlabel('\epsilon');
        ylabel('Relative Error');
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');
        set(gca, 'XDir', 'reverse');

    end

    figure(hFig); % Bring figure to front
    loglog(epsilons, errors, '-o', 'LineWidth', 1.5, 'DisplayName', label_name);
    legend('-DynamicLegend', 'Interpreter', 'none', 'Location', 'best');
    drawnow; % Force update
end


    end
end


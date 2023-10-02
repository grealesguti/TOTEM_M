%TOTEM_M Script
%Inputs
close all; clear all; clc;
inputfilename = "Benchmarks/Elements/Benchmark1_HexLinear/input_Benchmark1_LINEARHEX_PARAM.txt";
    
    % Define the path to the source folder (change this as needed)
    srcFolder = 'src';
    addpath(srcFolder);
    utilsFolder = 'utils';
    addpath(utilsFolder);
    
    % Create an instance of the InputReader class
        filepath="Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect.txt";
        %reader = InputReader("Benchmarks/Elements/Benchmark1_HexLinear/input_Benchmark1_LINEARHEX_PARAM.txt");
        %reader = InputReader("Benchmarks/Elements/Benchmark_HexSerendipity/input_NonLinCouplSEffect.txt");
        reader = InputReader(filepath);
        fprintf('Initialized InputReader with filename: %s\n', inputfilename);
        mesh = Mesh(reader);
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

        bench = ThermoelectricBenchmarks();

        TOO_1 = TO_Objectives(reader,mesh,bcinit);
        eval_fun=@(reader,mesh,solver) TOO_1.fval_AverageTemp(reader,mesh,solver);
        [FD_vals(1), err] = bench.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); % OK

        TOO_1 = TO_Objectives(reader,mesh,bcinit);
        eval_fun=@(reader,mesh,solver) TOO_1.fval_AverageTemp(reader,mesh,solver);
        [FD_vals(2), err] = bench.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, 4,1); %  OK

        TOC_1 = TO_Constraints(reader,mesh,bcinit);
        eval_fun=@(reader,mesh,solver) TOC_1.fval_Volume(reader,mesh,solver,2); % matters which index is given!!!
        [FD_vals(3), err] = bench.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); % OK

        TOC_1 = TO_Constraints(reader,mesh,bcinit);
        eval_fun=@(reader,mesh,solver) TOC_1.fval_Power(reader,mesh,solver,1);
        [FD_vals(4), err] = bench.Finite_Differences_DensityElement( reader, mesh, bcinit, solver, eval_fun,1); % Not OK

        TOC_1 = TO_Constraints(reader,mesh,bcinit);
        eval_fun=@(reader,mesh,solver) TOC_1.fval_Power(reader,mesh,solver,1);
        [FD_vals(5), err] = bench.Finite_Differences_bc(filepath, reader, mesh, solver, eval_fun, 4,1); % Not OK

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %postprocessing = Postprocessing();
        %postprocessing.initVTK(reader,mesh);
        %postprocessing.Benchmark_T_PLOT_axis(1,solver,2)
        %postprocessing.VTK_TV(solver,"Benchmarks/Elements/input_Benchmark2_QUADSERENDIPITYHEX_PARAM")
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TO = TopOpt(reader,mesh);

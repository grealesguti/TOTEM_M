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

        TO = TopOpt(reader,mesh);

        TO.runMMA(reader,mesh)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %postprocessing = Postprocessing();
        %postprocessing.initVTK(reader,mesh);
        %postprocessing.Benchmark_T_PLOT_axis(1,solver,2)
        %postprocessing.VTK_TV(solver,"Benchmarks/Elements/input_Benchmark2_QUADSERENDIPITYHEX_PARAM")
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TO = TopOpt(reader,mesh);

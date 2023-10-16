%TOTEM_M Script
%Inputs
close all; clear all; clc;
inputfilename = "Benchmarks/Elements/Benchmark1_HexLinear/input_Benchmark1_LINEARHEX_PARAM.txt";
    
    % Define the folder paths
    srcFolder = 'src';
    utilsFolder = 'utils';
    
    % Generate paths including subfolders
    srcPath = genpath(srcFolder);
    utilsPath = genpath(utilsFolder);
    
    % Add the generated paths to the MATLAB path
    addpath(srcPath);
    addpath(utilsPath);
    
    % Create an instance of the InputReader class
        %filepath="Benchmarks/Elements/Benchmark_TO/input_NonLinCouplSEffect.txt";
        filepath="TECTO/input_TECTO_StressConstrained_noncte.txt";
        %filepath="Benchmarks/Elements/Benchmark_DecoupledThermoElectroMech/input_NonLinCouplSEffect.txt";
        %reader = InputReader("Benchmarks/Elements/Benchmark1_HexLinear/input_Benchmark1_LINEARHEX_PARAM.txt");
        %reader = InputReader("Benchmarks/Elements/Benchmark_HexSerendipity/input_NonLinCouplSEffect.txt");
        reader = InputReader(filepath);
        fprintf('Initialized InputReader with filename: %s\n', inputfilename);
        mesh = Mesh(reader);
        fprintf('Initialized Mesh\n');

        postprocessing = Postprocessing();
        postprocessing.initVTK(reader,mesh);
        postprocessing.VTK_Mesh("TECTO/Results/Mesh")

        bcinit = BCInit(reader, mesh);
        fprintf('Initialized Loads\n');
        %solver = Solver(mesh, bcinit);
        %solver.runNewtonRaphson(reader, mesh, bcinit);
        %fprintf('Postprocessing\n');
        
        %postprocessing.VTK_TV(solver,"TECTO/Results/VTK_TV")

        %TOO = TO_Objectives(reader,mesh,bcinit);
        %TOO.CalculateObjective(reader,mesh,solver)
        %TOC = TO_Constraints(reader,mesh,bcinit);
        %TOC.CalculateConstraint(reader,mesh,solver);

        TO = TopOpt(reader,mesh);

        TO.runMMA(reader,mesh)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %postprocessing.Benchmark_T_PLOT_axis(1,solver,2)
        postprocessing.VTK_TV(solver,"Benchmarks/Elements/input_Benchmark2_QUADSERENDIPITYHEX_PARAM")
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TO = TopOpt(reader,mesh);

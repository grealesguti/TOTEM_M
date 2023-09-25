%TOTEM_M Script
%Inputs
close all; clear all; clc;
inputfilename = "Benchmarks/Elements/Benchmark2_SERENDIPITYHEX_PARAM_M.txt";
    
    % Define the path to the source folder (change this as needed)
    srcFolder = 'src';
    addpath(srcFolder);
    utilsFolder = 'utils';
    addpath(utilsFolder);
    
    % Create an instance of the InputReader class
        reader = InputReader("Benchmarks/Elements/input_Benchmark2_QUADSERENDIPITYHEX_PARAM.txt");
        fprintf('Initialized InputReader with filename: %s\n', inputfilename);
        mesh = Mesh(reader);
        fprintf('Initialized Mesh\n');
        bcinit = BCInit(reader, mesh);
        fprintf('Initialized Loads\n');
        solver = Solver(reader, mesh, bcinit);

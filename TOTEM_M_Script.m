%TOTEM_M Script
%Inputs
close all; clear all; clc;
inputfilename = "Benchmarks/Elements/input_Benchmark1_LINEARHEX_PARAM.txt";
    
    % Define the path to the source folder (change this as needed)
    srcFolder = 'src';
    addpath(srcFolder);

    % Create an instance of the InputReader class
        reader = InputReader("Benchmarks/Elements/input_Benchmark1_LINEARHEX_PARAM.txt");
        fprintf('Initialized InputReader with filename: %s\n', inputfilename);
        mesh = Mesh(reader);
        fprintf('Initialized Mesh');
        bcinit = BCInit(reader, mesh);

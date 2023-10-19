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
        filepath="TECTO/input_TECTO_StressConstrained_cte.txt";
        %filepath="Benchmarks/Elements/Benchmark_DecoupledThermoElectroMech/input_NonLinCouplSEffect.txt";
        %reader = InputReader("Benchmarks/Elements/Benchmark1_HexLinear/input_Benchmark1_LINEARHEX_PARAM.txt");
        %reader = InputReader("Benchmarks/Elements/Benchmark_HexSerendipity/input_NonLinCouplSEffect.txt");
        reader = InputReader(filepath);
        fprintf('Initialized InputReader with filename: %s\n', inputfilename);
        mesh = Mesh(reader);
        fprintf('Initialized Mesh\n');

        postprocessing = Postprocessing();
        postprocessing.initVTK(reader,mesh);
        %postprocessing.VTK_Mesh("TECTO/Results/Mesh")
        %bcinit = BCInit(reader, mesh);
        %solver = Solver(mesh, bcinit);
        %solver.runNewtonRaphson(reader, mesh, bcinit);
        %fprintf('Postprocessing\n');
                TO = TopOpt(reader,mesh);
        
                TO.runMMA(reader,mesh)
                %postprocessing.VTK_TV(solver,"TECTO/Results/VTK_TV")

        %TOO = TO_Objectives(reader,mesh,bcinit);
        %TOO.CalculateObjective(reader,mesh,solver)
        %TOC = TO_Constraints(reader,mesh,bcinit);
        %TOC.CalculateConstraint(reader,mesh,solver);
        
        qmin=6000;qinval=qmin;qn=5;qinstep=(30000-qmin)/qn;
        pmin=0.011;powval=pmin;pn=4;powstep=(0.05-pmin)/pn;
        for i=1:qn
            for j=1:pn
                reader.bcval(6)=qinval;
                reader.TopOpt_ConstraintValue(1)=powval;

                %bcinit = BCInit(reader, mesh);

                close all
                reader.rst_folder=append(reader.rst_folder,'Q',num2str(qinval),'_P',num2str(powval));
    


                qinval=qinval+qinstep;
                powval=powval+powstep;
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %postprocessing.Benchmark_T_PLOT_axis(1,solver,2)
        postprocessing.VTK_TV(solver,"Benchmarks/Elements/input_Benchmark2_QUADSERENDIPITYHEX_PARAM")
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TO = TopOpt(reader,mesh);

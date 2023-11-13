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
        %filepath="TECTO/input_TECTO_StressConstrained_cte.txt";
        %filepath="TECTO/input_TECTO_StressConstrained_cte_Al2O3v2.txt";
        %filepath="Benchmarks/Elements/Benchmark_HexLinear/input_LinearUncoupledSeebeck.txt"
        %filepath="Benchmarks/Elements/Benchmark_HexSerendipity/input_LinearUncoupledSeebeck.txt"
        %reader = InputReader("Benchmarks/Elements/Benchmark1_HexLinear/input_Benchmark1_LINEARHEX_PARAM.txt");
        %reader = InputReader("Benchmarks/Elements/Benchmark_HexSerendipity/input_LinearUncoupledSeebeck.txt");
        powval=15;
        qinval=5000;
        filepath="TECTO/input_TECTO_Thermoel_Quad_cte_th.txt";      
        reader = InputReader(filepath);
        reader.Rst_name=append('Al2O3_nonlin_Q_',num2str(qinval),'_P',num2str(powval));
        %reader.bcval(9)=qinval;
        %reader.TopOpt_ConstraintValue(1)=powval;
        mesh = Mesh(reader);
        TO = TopOpt(reader,mesh);
        TO.runMMA(reader,mesh)

        filepath="TECTO/input_TECTO_Thermoel_Quad_cte.txt";      
        reader = InputReader(filepath);
        reader.Rst_name=append('Al2O3_nonlin_0topU_Q_',num2str(qinval),'_P',num2str(powval));
        %reader.bcval(9)=qinval;
        %reader.TopOpt_ConstraintValue(1)=powval;
        mesh = Mesh(reader);
        TO = TopOpt(reader,mesh);
        TO.runMMA(reader,mesh)

        filepath="TECTO/input_TECTO_Thermoel_Serend_cte.txt";   
        reader = InputReader(filepath);
        reader.Rst_name=append('Al2O3_nonlin_0topU_HF_3D_Q_',num2str(qinval),'_P',num2str(powval));
        %reader.bcval(10)=qinval;
        %reader.TopOpt_ConstraintValue(1)=powval;
        mesh = Mesh(reader);
        TO = TopOpt(reader,mesh);
        TO.runMMA(reader,mesh)

        filepath="TECTO/input_TECTO_Thermoel_Serend_cte.txt";   
        reader = InputReader(filepath);
        reader.Rst_name=append('Al2O3_nonlin_0topU_3D_Q_',num2str(qinval),'_P',num2str(powval));
        %reader.bcval(10)=qinval;
        %reader.TopOpt_ConstraintValue(1)=powval;
        mesh = Mesh(reader);
        %bcinit = BCInit(reader, mesh);
        %solver = Solver(mesh, bcinit);
        %residual_norm=solver.runNewtonRaphson(reader, mesh, bcinit);
        TO = TopOpt(reader,mesh);
        TO.runMMA(reader,mesh)


        %funAL = @(t) solver.funArcLen(reader,bcinit,mesh,t); 
        %funAL2 = @(t) solver.funArcLen_LamMorley(reader,bcinit,mesh,t); 
        %u0=bcinit.initialdofs_;
        %odd_numbers = 1:2:length(u0);
        %even_numbers = 2:2:length(u0);
        %u0(odd_numbers)=str2double(reader.T0);
        %u0(even_numbers)=max(u0(even_numbers));
        %[ufree] = arc_length_Crisfield(funAL,u0(bcinit.dofs_free_));
        %[ufree] = arc_length_Lam_Morley_modified(funAL,u0(bcinit.dofs_free_));
        %solver.runNewtonRaphson(reader, mesh, bcinit);
        %filter=Filtering(reader,mesh);
        %filter.filter_densities(reader,mesh)
        %fprintf('Postprocessing\n');
        %    post = Postprocessing();
        %    post.initVTK(reader,mesh);
        %   post.VTK_TV(solver,"TECTO/Results/VTK_TV")
        %    post.VTK_U(solver,"TECTO/Results/VTK_U")
        %   post.VTK_matidx_TV(mesh,solver,"TECTO/Results/VTK_matTV")
        %    post.VTK_x_TV(mesh,solver,"TECTO/Results/xxfilter.vtk")
        %post.VTK_qth(reader,mesh,solver,"TECTO/Results/VTK_qth")
        %post.VTK_q(reader,mesh,solver,"TECTO/Results/VTK_q")
        %post.VTK_j(reader,mesh,solver,"TECTO/Results/VTK_j")
        %TOO = TO_Objectives(reader,mesh,bcinit);
        %TOO.CalculateObjective(reader,mesh,solver)
        %TOC = TO_Constraints(reader,mesh,bcinit);
        %TOC.CalculateConstraint(reader,mesh,solver);
        %postprocessing.VTK_TV(solver,'test.vtk')
        %filter.filter_sensitivities(reader,mesh,TOO,TOC)

        qmin=6000;qinval=qmin;qn=3;qinstep=(40000-qmin)/qn;
        scale=1000;
        pmin=0.015*scale;powval=pmin;pn=3;powstep=(0.1*scale-pmin)/pn;
        initial_folder=reader.rst_folder;
        for i=1:qn
            powval=pmin;
            for j=1:pn
                
                reader.bcval(9)=qinval;
                reader.TopOpt_ConstraintValue(1)=powval;
                mesh = Mesh(reader);
                bcinit = BCInit(reader, mesh);

                close all
                reader.Rst_name=append('Al2O3_nonlin_0topdispl_Q_',num2str(qinval),'_P',num2str(powval));
                %reader.Rst_name=append('Q',num2str(qinval),'_P',num2str(powval));
                TO = TopOpt(reader,mesh);
        
                TO.runMMA(reader,mesh)
                powval=powval+powstep;
            end
            qinval=qinval+qinstep;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %postprocessing.Benchmark_T_PLOT_axis(1,solver,2)
        postprocessing.VTK_TV(solver,"Benchmarks/Elements/input_Benchmark2_QUADSERENDIPITYHEX_PARAM")
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TO = TopOpt(reader,mesh);

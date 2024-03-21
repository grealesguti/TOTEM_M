%TOTEM_M Script
%Inputs
close all; clear all; clc;    
    % Define the folder paths
    srcFolder = 'src';
    utilsFolder = 'utils';
    
    % Generate paths including subfolders
    srcPath = genpath(srcFolder);
    utilsPath = genpath(utilsFolder);
    
    % Add the generated paths to the MATLAB path
    addpath(srcPath);
    addpath(utilsPath);
    
    % 2D results
    % case study lists, Replace with your desired values
    volumes =[1,1,1,1,1,1,0.5,0.5];
    qinvals = [3000, 3000, 5000, 5000, 3000, 3000,3000,5000]; % Replace with your desired values
    powvals = [15,15,15,15,30,30,15,15]; % Replace with your desired values
    filterval = [0,1,0,1,0,1,1,1]; % Replace with your desired values

    volumes =[0.6];
    qinvals = [7500]; % Replace with your desired values
    powvals = [15]; % Replace with your desired values
    filterval = [1]; % Replace with your desired values
    stresses = [10E6]; % Replace with your desired values

    filepath="TECTO/AnsysInp/input_AnsysSolutionRun.txt";   
    reader = InputReader(filepath);
    pcon_loc = 1;
    vcon_loc = 2;
    qinloc = 9;

    mesh = Mesh(reader);
    post = Postprocessing();
    post.initVTK(reader,mesh);
    %post.VTK_Mesh("TECTO/AnsysInp/test")
    %VTKapplyDensity(reader,mesh,1)
    %post.VTK_Mesh_xx("TECTO/AnsysInp/testxx",mesh)
    bcinit = BCInit(reader, mesh);
    post.VTK_freedofs(bcinit,"TECTO/AnsysInp/bcinit.vtk")
    %post.VTK_x_TV(mesh,solver,"TECTO/AnsysInp/inittestTV.vtk")
    solver = Solver(mesh, bcinit);
    post.VTK_matidx_TV(mesh,solver,"TECTO/AnsysInp/matinit.vtk")
    residual_norm=solver.runNewtonRaphson(reader, mesh, bcinit);
    post.VTK_x_TV(mesh,solver,"TECTO/AnsysInp/testTV.vtk")
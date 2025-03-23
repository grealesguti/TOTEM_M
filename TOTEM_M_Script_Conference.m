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
    volumes =[1,1,1];
    qinvals = [3000,3000,3000]; % Replace with your desired values
    powvals = [15,15,15]; % Replace with your desired values
    filterval = [0,0,1]; % Replace with your desired values
    stresses = [1e100,10E6,10E6]; % Replace with your desired values
    volumes =[1];
    qinvals = [3000]; % Replace with your desired values
    powvals = [15]; % Replace with your desired values
    filterval = [0]; % Replace with your desired values
    stresses = [1E100]; % Replace with your desired values
    %volumes =[1];
    %qinvals = [3000]; % Replace with your desired values
    %powvals = [15]; % Replace with your desired values
    %filterval = [0]; % Replace with your desired values
    %stresses = [50e7]; % Replace with your desired values

    filepath="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D.txt";   
    filepath="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2DConf.txt";   

    reader = InputReader(filepath);
    pcon_loc = 1;
    vcon_loc = 2;
    qinloc = 9;
    for i = 1:length(qinvals)
            reader.bcval(qinloc) = qinvals(i);
            reader.TopOpt_ConstraintValue(pcon_loc) = powvals(i);
            reader.TopOpt_ConstraintValue(vcon_loc) = volumes(i);
            reader.Filter = filterval(i);
            reader.TopOpt_ConstraintValue(3) = stresses(i);
            mesh = Mesh(reader);
            bcinit = BCInit(reader, mesh);
            close all
            reader.Rst_name = append('Journal_Al2O3_nonlin_0topdispl_Q_', num2str(qinvals(i)), '_P', num2str(powvals(i)),'_F',num2str(filterval(i)),'_v',num2str(volumes(i)));
            TO = TopOpt(reader, mesh);
            TO.runMMA(reader, mesh)
    end


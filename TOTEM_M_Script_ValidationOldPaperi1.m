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

    folderPath = "C:/Archive/Programming/MATLAB/TOTEM_M/TECTO/AnsysInp/xxvalues";  % Replace 'path_to_your_folder' with the actual folder path
    files = dir(fullfile(folderPath, '*.vtk'));  % Change '*.txt' to match the file extension of the files you want to read
  %% Step 1
    filepath="TECTO/AnsysInp/inputs/input_AnsysSolutionRun.txt";   
    fvals_old=zeros(length(files),1);
i=15;
reader = InputReader(filepath);
    
       reader.TopOpt_Initial_x=append(files(i).folder,"\",files(i).name); % to modify in loop
    
        modifiedFileName0= reader.TopOpt_Initial_x;
        % Substitute '.vtk' by '.csv' and 'inputs' by 'mesh'
        modifiedFileName = strrep(modifiedFileName0, '.vtk', '.csv');
        modifiedFileName = strrep(modifiedFileName, 'xxvalues', 'csv');

        %filename = 'folder/HeatStudy__heat_20000_pcon_15_.vtk';
        
        % Define the regular expression pattern to match the expected format
        pattern = '^.*heat_(\d+)_pcon_(\d+)_.*\.vtk$';
        
        % Check if the filename matches the pattern
        if regexp(modifiedFileName0, pattern)
            % If the filename matches the pattern, extract heat and pcon values
            match = regexp(modifiedFileName0, pattern, 'tokens');
            heat = str2double(match{1}{1});  % Convert to numeric value
            pcon = str2double(match{1}{2});  % Convert to numeric value
            
            % Display the extracted values
            disp(['Heat: ', num2str(heat)]);
            disp(['Pcon: ', num2str(pcon)]);

            reader.bcval(4)=heat; % keep 1 constant
        else
            disp('No heat change');
        end
        
        % Read the subsequent file
        data = readtable(modifiedFileName); % Assuming the file is a CSV file, adjust this function accordingly if the file format is different
        voltages = data{:, 4};
        outeriter = max(data{:,1});
        endvoltage = voltages(outeriter);
        reader.TObcval(1)=endvoltage;
        reader.bcval(1)=endvoltage; % keep 1 constant
    
        mesh = Mesh(reader);
        post = Postprocessing();
        post.initVTK(reader,mesh);
        VTKapplyDensity(reader,mesh,1)
        bcinit = BCInit(reader, mesh);
        solver = Solver(mesh, bcinit);
        residual_norm=solver.runNewtonRaphson(reader, mesh, bcinit);
    
        outputname = strrep(modifiedFileName0, 'xxvalues', 'output');
        post.VTK_x_TV(mesh,solver,outputname)
        % calculate objective!!! and store!!!
        TOO = TO_Objectives(reader,mesh,bcinit);
        TOO.CalculateObjective(reader,mesh,solver)    
        fvals_old(i)=TOO.fval;
        Pvalue_old=CalculatePower(reader,mesh,solver);

  %%% Step 2
%     % nonlinear properties!!!
    filepath="TECTO/AnsysInp/inputs/input_AnsysSolutionRun_nonlin.txt";   
    fvals_nonlin=zeros(length(files),1);
        reader = InputReader(filepath);
    
        reader.TopOpt_Initial_x=append(files(i).folder,"\",files(i).name); % to modify in loop
    
        modifiedFileName0= reader.TopOpt_Initial_x;
        % Substitute '.vtk' by '.csv' and 'inputs' by 'mesh'
        modifiedFileName = strrep(modifiedFileName0, '.vtk', '.csv');
        modifiedFileName = strrep(modifiedFileName, 'xxvalues', 'csv');

        %filename = 'folder/HeatStudy__heat_20000_pcon_15_.vtk';
        
        % Define the regular expression pattern to match the expected format
        pattern = '^.*heat_(\d+)_pcon_(\d+)_.*\.vtk$';
        
        % Check if the filename matches the pattern
        if regexp(modifiedFileName0, pattern)
            % If the filename matches the pattern, extract heat and pcon values
            match = regexp(modifiedFileName0, pattern, 'tokens');
            heat = str2double(match{1}{1});  % Convert to numeric value
            pcon = str2double(match{1}{2});  % Convert to numeric value
            
            % Display the extracted values
            disp(['Heat: ', num2str(heat)]);
            disp(['Pcon: ', num2str(pcon)]);

            reader.bcval(4)=heat; % keep 1 constant
        else
            disp('No heat change');
        end
        
        % Read the subsequent file
        data = readtable(modifiedFileName); % Assuming the file is a CSV file, adjust this function accordingly if the file format is different
        voltages = data{:, 5};
        endvoltage = voltages(end);
    
        reader.bcval(1)=endvoltage; % keep 1 constant
    
        mesh = Mesh(reader);
        post = Postprocessing();
        post.initVTK(reader,mesh);
        VTKapplyDensity(reader,mesh,1)
        bcinit = BCInit(reader, mesh);
        solver = Solver(mesh, bcinit);
        residual_norm=solver.runNewtonRaphson(reader, mesh, bcinit);
    
        outputname = strrep(modifiedFileName0, 'xxvalues', 'output');
        outputname = strrep(outputname, '.vtk', '_nonlin.vtk');
        post.VTK_x_TV(mesh,solver,outputname)
        % calculate objective!!! and store!!!
        TOO = TO_Objectives(reader,mesh,bcinit);
        TOO.CalculateObjective(reader,mesh,solver)    
        fvals_nonlin(i)=TOO.fval;
        Pvalue_nonlin(i)=CalculatePower(reader,mesh,solver);

  %%% Step 3
%non-linear & linear optimal results full volume, only for heat tests!!!
    filepath="TECTO/AnsysInp/inputs/input_AnsysSolutionRun.txt"; 
    fvals_TOCte=zeros(length(files),1);
i=1;
reader = InputReader(filepath);
    
        reader.TopOpt_Initial_x=append(files(i).folder,"\",files(i).name); % to modify in loop
        modifiedFileName0= reader.TopOpt_Initial_x;
        reader.TopOpt_Initial_x=1;
        
        % Define the regular expression pattern to match the expected format
        pattern = '^.*heat_(\d+)_pcon_(\d+)_.*\.vtk$';
        
        % Check if the filename matches the pattern
        if regexp(modifiedFileName0, pattern)
            % If the filename matches the pattern, extract heat and pcon values
            match = regexp(modifiedFileName0, pattern, 'tokens');
            heat = str2double(match{1}{1});  % Convert to numeric value
            pcon = str2double(match{1}{2});  % Convert to numeric value
            
            % Display the extracted values
            disp(['Heat: ', num2str(heat)]);
            disp(['Pcon: ', num2str(pcon)]);

            reader.bcval(4)=heat; % keep 1 constant
            reader.TopOpt_ConstraintValue(1)=pcon/1000;
        else
            disp('No heat change');
        end
        
        mesh = Mesh(reader);
        TO = TopOpt(reader, mesh);
        test=TO.runMMA(reader, mesh);
        fvals_TOCte(i)=test.f0val;

    %%% Step 4

        filepath="TECTO/AnsysInp/inputs/input_AnsysSolutionRun_nonlin.txt";   
        fvals_TONonlin=zeros(length(files),1);

        reader = InputReader(filepath);
    
        reader.TopOpt_Initial_x=append(files(i).folder,"\",files(i).name); % to modify in loop
        modifiedFileName0= reader.TopOpt_Initial_x;
        reader.TopOpt_Initial_x=1;
        
        % Define the regular expression pattern to match the expected format
        pattern = '^.*heat_(\d+)_pcon_(\d+)_.*\.vtk$';
        
        % Check if the filename matches the pattern
        if regexp(modifiedFileName0, pattern)
            % If the filename matches the pattern, extract heat and pcon values
            match = regexp(modifiedFileName0, pattern, 'tokens');
            heat = str2double(match{1}{1});  % Convert to numeric value
            pcon = str2double(match{1}{2});  % Convert to numeric value
            
            % Display the extracted values
            disp(['Heat: ', num2str(heat)]);
            disp(['Pcon: ', num2str(pcon)]);

            reader.bcval(4)=heat; % keep 1 constant
            reader.TopOpt_ConstraintValue(1)=pcon/1000;
        else
            disp('No heat change');
        end
        
        mesh = Mesh(reader);
        TO = TopOpt(reader, mesh);
        test=TO.runMMA(reader, mesh);
        fvals_TONonlin(i)=test.f0val;

        Pvals_TONonlin(i)=CalculatePower(reader,mesh,solver);




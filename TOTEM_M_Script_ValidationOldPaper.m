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
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 0
    filepath="TECTO/AnsysInp/inputs/input_AnsysSolutionRun.txt";   
    fvals_old=zeros(length(files),1);
    for i=1:length(files)
        reader = InputReader(filepath);
    
       reader.TopOpt_Initial_x=append(files(i).folder,"\",files(i).name); % to modify in loop
    
        modifiedFileName0= reader.TopOpt_Initial_x;
        % Substitute '.vtk' by '.csv' and 'inputs' by 'mesh'
        modifiedFileName = strrep(modifiedFileName0, '.vtk', '.csv');
        modifiedFileName = strrep(modifiedFileName, 'xxvalues', 'csv');

        %filename = 'folder/HeatStudy__heat_20000_pcon_15_.vtk';
        
        % Define the regular expression pattern to match the expected format
        pattern = 'heat_(\d+)';
        % Check if the filename matches the pattern
        if regexp(files(i).name, 'heat_(\d+)')
            % If the filename matches the pattern, extract heat and pcon values
            h=regexp(files(i).name, 'heat_(\d+)', 'tokens', 'once');
                        v=regexp(files(i).name, 'Vol_0.(\d+)', 'tokens', 'once');
        if isempty(v)
                    p=regexp(files(i).name, 'pcon_0.0(\d+)', 'tokens', 'once');
                    v=1;
        else
                    p=regexp(files(i).name, 'pcon_0.0(\d+).vtk', 'tokens', 'once');
                    v=str2double(v);
        end
            heat = str2double(h);  % Convert to numeric value
            pcon = str2double(p);  % Convert to numeric value
        if isempty(p)
            p=regexp(files(i).name, 'pcon_(\d+)', 'tokens', 'once');
                        pcon = str2double(p);  % Convert to numeric value

        end
        
        if pcon<5
                pcon=pcon*10;
            end

            % Display the extracted values
            disp(['Heat: ', num2str(heat)]);
            disp(['Pcon: ', num2str(pcon)]);

            reader.bcval(4)=heat; % keep 1 constant
        else
            heat=0;
            pcon=0;
            v=0;
            disp('No heat change');
        end
        if pcon==0
            disp('pcon = 0')
        end
        
        % Read the subsequent file
        data = readtable(modifiedFileName); % Assuming the file is a CSV file, adjust this function accordingly if the file format is different
        voltages = data{:, 4};
        outeriter = max(data{:,1});
        endvoltage = voltages(outeriter);
        reader.TObcval(1)=endvoltage;
        reader.bcval(1)=endvoltage; % keep 1 constant
        reader.bcval(5)=endvoltage; % keep 1 constant
        reader.TObcminval(1)=endvoltage;

        data_M(i,:)=[outeriter,endvoltage,heat,pcon,v,data{outeriter, 3}];
        file_M{i}=files(i).name;
    end
        file_M=file_M';

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Step 1
%     n=length(data_M(1,:));
% 
%     filepath="TECTO/AnsysInp/inputs/input_AnsysSolutionRun.txt";   
%     fvals_old=zeros(length(files),1);
%     for i=1:length(files)
%         fprintf('### STEP 1: %i out of %i\n', i, length(files));
%         reader = InputReader(filepath);
%     
%         reader.TopOpt_Initial_x=append(files(i).folder,"\",files(i).name); % to modify in loop
%     
%         modifiedFileName0= reader.TopOpt_Initial_x;
%         % Substitute '.vtk' by '.csv' and 'inputs' by 'mesh'
%         modifiedFileName = strrep(modifiedFileName0, '.vtk', '.csv');
%         modifiedFileName = strrep(modifiedFileName, 'xxvalues', 'csv');
% 
%         %filename = 'folder/HeatStudy__heat_20000_pcon_15_.vtk';
%         
%         % Define the regular expression pattern to match the expected format
%         pattern = '^.*heat_(\d+)_pcon_(\d+)_.*\.vtk$';
%         
%         % Check if the filename matches the pattern
%          if regexp(files(i).name, 'heat_(\d+)')
%             % If the filename matches the pattern, extract heat and pcon values
%             h=regexp(files(i).name, 'heat_(\d+)', 'tokens', 'once');
%                         v=regexp(files(i).name, 'Vol_0.(\d+)', 'tokens', 'once');
%         if isempty(v)
%                     p=regexp(files(i).name, 'pcon_(\d+)', 'tokens', 'once');
%                     v=0;
%         
%         else
%                     p=regexp(files(i).name, 'pcon_0.0(\d+).vtk', 'tokens', 'once');
%                     v=str2double(v);
%         end
%             heat = str2double(h);  % Convert to numeric value
%             pcon = str2double(p);  % Convert to numeric value
%             if pcon<5
%                 pcon=pcon*10;
%             end
% 
%             % Display the extracted values
%             disp(['Heat: ', num2str(heat)]);
%             disp(['Pcon: ', num2str(pcon)]);
% 
%             reader.bcval(4)=heat; % keep 1 constant
%         else
%             heat=0;
%             pcon=0;
%             v=0;
%             disp('No heat change');
%         end
%         
%         
%         % Read the subsequent file
%         data = readtable(modifiedFileName); % Assuming the file is a CSV file, adjust this function accordingly if the file format is different
%         voltages = data{:, 4};
%         outeriter = max(data{:,1});
%         endvoltage = voltages(outeriter);
%         reader.TObcval(1)=endvoltage;
%         reader.bcval(1)=endvoltage; % keep 1 constant
%         reader.bcval(5)=endvoltage; % keep 1 constant
%         reader.TObcminval(1)=endvoltage;
% 
%         mesh = Mesh(reader);
%         post = Postprocessing();
%         post.initVTK(reader,mesh);
%         VTKapplyDensity(reader,mesh,1)
%         bcinit = BCInit(reader, mesh);
%         solver = Solver(mesh, bcinit);
%         residual_norm=solver.runNewtonRaphson(reader, mesh, bcinit);
%     
%         outputname = strrep(modifiedFileName0, 'xxvalues', 'output');
%         post.VTK_x_TV(mesh,solver,outputname)
%         % calculate objective!!! and store!!!
%         TOO = TO_Objectives(reader,mesh,bcinit);
%         TOO.CalculateObjective(reader,mesh,solver)    
%         fvals_old(i)=TOO.fval;
%         Pvalue_old(i)=CalculatePower(reader,mesh,solver);
%         data_M(i,n+1)=fvals_old(i);
%         data_M(i,n+2)=Pvalue_old(i);
%     end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %% Step 2
%   load('C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\AnsysInp\output\Step1_data.mat')
%       n=length(data_M(1,:));
% 
% %     % nonlinear properties!!!
%     filepath="TECTO/AnsysInp/inputs/input_AnsysSolutionRun_nonlin.txt";   
%     fvals_nonlin=zeros(length(files),1);
%     for i=1:length(files)
%         fprintf('### STEP 2: %i out of %i\n', i, length(files));
%         reader = InputReader(filepath);
%     
%         reader.TopOpt_Initial_x=append(files(i).folder,"\",files(i).name); % to modify in loop
%     
%         modifiedFileName0= reader.TopOpt_Initial_x;
%         % Substitute '.vtk' by '.csv' and 'inputs' by 'mesh'
%         modifiedFileName = strrep(modifiedFileName0, '.vtk', '.csv');
%         modifiedFileName = strrep(modifiedFileName, 'xxvalues', 'csv');
% 
%         %filename = 'folder/HeatStudy__heat_20000_pcon_15_.vtk';
%         
%         % Define the regular expression pattern to match the expected format
%         pattern = '^.*heat_(\d+)_pcon_(\d+)_.*\.vtk$';
%         
%         % Check if the filename matches the pattern
%          if regexp(files(i).name, 'heat_(\d+)')
%             % If the filename matches the pattern, extract heat and pcon values
%             h=regexp(files(i).name, 'heat_(\d+)', 'tokens', 'once');
%                         v=regexp(files(i).name, 'Vol_0.(\d+)', 'tokens', 'once');
%         if isempty(v)
%                     p=regexp(files(i).name, 'pcon_(\d+)', 'tokens', 'once');
%                     v=0;
%         
%         else
%                     p=regexp(files(i).name, 'pcon_0.0(\d+).vtk', 'tokens', 'once');
%                     v=str2double(v);
%         end
%             heat = str2double(h);  % Convert to numeric value
%             pcon = str2double(p);  % Convert to numeric value
%             if pcon<5
%                 pcon=pcon*10;
%             end
% 
%             % Display the extracted values
%             disp(['Heat: ', num2str(heat)]);
%             disp(['Pcon: ', num2str(pcon)]);
% 
%             reader.bcval(4)=heat; % keep 1 constant
%         else
%             heat=0;
%             pcon=0;
%             v=0;
%             disp('No heat change');
%         end
%        
%         
%         % Read the subsequent file
%         data = readtable(modifiedFileName); % Assuming the file is a CSV file, adjust this function accordingly if the file format is different
%         voltages = data{:, 4};
%         outeriter = max(data{:,1});
%         endvoltage = voltages(outeriter);
%         reader.TObcval(1)=endvoltage;
%         reader.bcval(5)=endvoltage; % keep 1 constant
%         reader.bcval(1)=endvoltage; % keep 1 constant
%     
%         mesh = Mesh(reader);
%         post = Postprocessing();
%         post.initVTK(reader,mesh);
%         VTKapplyDensity(reader,mesh,1)
%         bcinit = BCInit(reader, mesh);
%         solver = Solver(mesh, bcinit);
%         residual_norm=solver.runNewtonRaphson(reader, mesh, bcinit);
%     
%         outputname = strrep(modifiedFileName0, 'xxvalues', 'output');
%         outputname = strrep(outputname, '.vtk', '_nonlin.vtk');
%         post.VTK_x_TV(mesh,solver,outputname)
%         % calculate objective!!! and store!!!
%         TOO = TO_Objectives(reader,mesh,bcinit);
%         TOO.CalculateObjective(reader,mesh,solver)    
%         fvals_nonlin(i)=TOO.fval;
%         Pvalue_nonlin(i)=CalculatePower(reader,mesh,solver);
%         data_M(i,n+1)=fvals_nonlin(i);
%         data_M(i,n+2)=Pvalue_nonlin(i);        
%     end
% 
%         fvalname = strrep(modifiedFileName0, 'xxvalues', 'output');
%         fvalname = strrep(fvalname, '.vtk', '_datavol.mat');
%         save(fvalname, 'data_');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Step 2.5
  %load('C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\AnsysInp\output\Step2_data.mat')
      n=length(data_M(1,:));
%non-linear & linear optimal results full volume, only for heat tests!!!
    filepath="TECTO/AnsysInp/inputs/input_AnsysSolutionRun_nonlinGoodPenal.txt";   
    fvals_TOCte=zeros(length(files),1);
    for i=7:length(files)
        fprintf('### STEP 2.5: %i out of %i\n', i, length(files));
        reader = InputReader(filepath);
            
        reader.TopOpt_Initial_x=append(files(i).folder,"\",files(i).name); % to modify in loop
        modifiedFileName0= reader.TopOpt_Initial_x;
        % Substitute '.vtk' by '.csv' and 'inputs' by 'mesh'
        modifiedFileName = strrep(modifiedFileName0, '.vtk', '.csv');
        modifiedFileName = strrep(modifiedFileName, 'xxvalues', 'csv');        %.TopOpt_Initial_x=1;
        
        % Define the regular expression pattern to match the expected format
        pattern = '^.*heat_(\d+)_pcon_(\d+)_.*\.vtk$';
        
        % Check if the filename matches the pattern
         if regexp(files(i).name, 'heat_(\d+)')
            % If the filename matches the pattern, extract heat and pcon values
            h=regexp(files(i).name, 'heat_(\d+)', 'tokens', 'once');
                        v=regexp(files(i).name, 'Vol_0.(\d+)', 'tokens', 'once');
        if isempty(v)
                    p=regexp(files(i).name, 'pcon_(\d+)', 'tokens', 'once');
                    v=0;
        
        else
                    p=regexp(files(i).name, 'pcon_0.0(\d+).vtk', 'tokens', 'once');
                    v=str2double(v);
        end
            heat = str2double(h);  % Convert to numeric value
            pcon = str2double(p);  % Convert to numeric value
            if pcon<5
                pcon=pcon*10;
            end

            % Display the extracted values
            disp(['Heat: ', num2str(heat)]);
            disp(['Pcon: ', num2str(pcon)]);

            reader.bcval(4)=heat; % keep 1 constant
            reader.TopOpt_ConstraintValue(1)=pcon/1000;
        else
            heat=0;
            pcon=0;
            v=0;
            disp('No heat change');
         end

         % Read the subsequent file
         data = readtable(modifiedFileName); % Assuming the file is a CSV file, adjust this function accordingly if the file format is different
         voltages = data{:, 4};
         outeriter = max(data{:,1});
         endvoltage = voltages(outeriter);
         reader.TObcval(1)=0.01;
        mesh = Mesh(reader);
        VTKapplyDensity(reader,mesh,1)
        %reader.TopOpt_DesignElements="";
        TO = TopOpt(reader, mesh);
        TO.onlyvol=1;
        test=TO.runMMA(reader, mesh);
        fvals_TOCte(i)=test.f0val;
        Pvalue_TOnonlin(i)=(test.fval(1)+1)*reader.TopOpt_ConstraintValue(1);
        data_M(i,n+1)=test.f0val;
        data_M(i,n+2)=Pvalue_TOnonlin(i);
        data_M(i,n+3)=test.Voltage_value;
    end
        fvalname = strrep(modifiedFileName0, 'xxvalues', 'output');
        fvalname = strrep(fvalname, '.vtk', '_datavolonly.mat');
        save(fvalname, 'data_M');
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Step 3
  load('C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\AnsysInp\output\Step25_data.mat')
  n=length(data_M(1,:));
%non-linear & linear optimal results full volume, only for heat tests!!!
%     filepath="TECTO/AnsysInp/inputs/input_AnsysSolutionRun.txt"; 
%     fvals_TOCte=zeros(length(files),1);
%     for i=1:length(files)
%         reader = InputReader(filepath);
%             fprintf('### STEP 3: %i out of %i\n', i, length(files));
% 
%         reader.TopOpt_Initial_x=append(files(i).folder,"\",files(i).name); % to modify in loop
%         modifiedFileName0= reader.TopOpt_Initial_x;
%         reader.TopOpt_Initial_x=1;
%         
%         % Define the regular expression pattern to match the expected format
%         pattern = '^.*heat_(\d+)_pcon_(\d+)_.*\.vtk$';
%         
%          % Check if the filename matches the pattern
%          if regexp(files(i).name, 'heat_(\d+)')
%             % If the filename matches the pattern, extract heat and pcon values
%             h=regexp(files(i).name, 'heat_(\d+)', 'tokens', 'once');
%                         v=regexp(files(i).name, 'Vol_0.(\d+)', 'tokens', 'once');
%         if isempty(v)
%                     p=regexp(files(i).name, 'pcon_(\d+)', 'tokens', 'once');
%                     v=0;
%         
%         else
%                     p=regexp(files(i).name, 'pcon_0.0(\d+).vtk', 'tokens', 'once');
%                     v=str2double(v);
%         end
%             heat = str2double(h);  % Convert to numeric value
%             pcon = str2double(p);  % Convert to numeric value
%             if pcon<5
%                 pcon=pcon*10;
%             end
% 
%             % Display the extracted values
%             disp(['Heat: ', num2str(heat)]);
%             disp(['Pcon: ', num2str(pcon)]);
% 
%             reader.bcval(4)=heat; % keep 1 constant
%             reader.TopOpt_ConstraintValue(1)=pcon/1000;
%         else
%             heat=0;
%             pcon=0;
%             v=0;
%             disp('No heat change');
%          end       
%         mesh = Mesh(reader);
%         TO = TopOpt(reader, mesh);
%         TO.onlyvol=1;
% 
%         test=TO.runMMA(reader, mesh);
%         fvals_TOCte(i)=test.f0val;
% 
%         Pvalue_TOnonlin(i)=(test.fval(1)+1)*reader.TopOpt_ConstraintValue(1);
%         data_M(i,n+1)=test.f0val;
%         data_M(i,n+2)=Pvalue_TOnonlin(i);
%         data_M(i,n+3)=test.Voltage_value;
%     end
%         fvalname = strrep(modifiedFileName0, 'xxvalues', 'output');
%         fvalname = strrep(fvalname, '.vtk', '_3data.mat');
%         save(fvalname, 'data_M');
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          % Step 4
%         n=length(data_M(1,:));
%         filepath="TECTO/AnsysInp/inputs/input_AnsysSolutionRun_nonlin.txt";   
%         fvals_TONonlin=zeros(length(files),1);
% 
%     for i=1:length(files)
%         reader = InputReader(filepath);
%             fprintf('### STEP 4: %i out of %i\n', i, length(files));
% 
%         reader.TopOpt_Initial_x=append(files(i).folder,"\",files(i).name); % to modify in loop
%         modifiedFileName0= reader.TopOpt_Initial_x;
%         reader.TopOpt_Initial_x=1;
%         
%        % Define the regular expression pattern to match the expected format
%         pattern = '^.*heat_(\d+)_pcon_(\d+)_.*\.vtk$';
%               % Check if the filename matches the pattern
%          if regexp(files(i).name, 'heat_(\d+)')
%             % If the filename matches the pattern, extract heat and pcon values
%             h=regexp(files(i).name, 'heat_(\d+)', 'tokens', 'once');
%                         v=regexp(files(i).name, 'Vol_0.(\d+)', 'tokens', 'once');
%         if isempty(v)
%                     p=regexp(files(i).name, 'pcon_(\d+)', 'tokens', 'once');
%                     v=0;
%         
%         else
%                     p=regexp(files(i).name, 'pcon_0.0(\d+).vtk', 'tokens', 'once');
%                     v=str2double(v);
%         end
%             heat = str2double(h);  % Convert to numeric value
%             pcon = str2double(p);  % Convert to numeric value
%             if pcon<5
%                 pcon=pcon*10;
%             end
% 
%             % Display the extracted values
%             disp(['Heat: ', num2str(heat)]);
%             disp(['Pcon: ', num2str(pcon)]);
% 
%             reader.bcval(4)=heat; % keep 1 constant
%             reader.TopOpt_ConstraintValue(1)=pcon/1000;
%         else
%             heat=0;
%             pcon=0;
%             v=0;
%             disp('No heat change');
%          end    
%         mesh = Mesh(reader);
%         TO = TopOpt(reader, mesh);
%                 TO.onlyvol=1;
% 
%         test=TO.runMMA(reader, mesh);
%         fvals_TONonlin(i)=test.f0val;
% 
%        Pvalue_TOnonlin(i)=(test.fval(1)+1)*reader.TopOpt_ConstraintValue(1);
%         data_M(i,n+1)=test.f0val;
%         data_M(i,n+2)=Pvalue_TOnonlin(i);
%         data_M(i,n+3)=test.Voltage_value;
%     end
% 
%         fvalname = strrep(modifiedFileName0, 'xxvalues', 'output');
%         fvalname = strrep(fvalname, '.vtk', '_4data.mat');
%         save(fvalname, 'data_M');


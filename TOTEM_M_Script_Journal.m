%TOTEM_M Script
%Inputs
close all; clear all; clc;
disp('Script is starting execution...'); % Print message to screen
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

qinvals = [7500,7500,7500,7500,7500]; % Replace with your desired values
powvals = [15,15,15,15,15]; % Replace with your desired values
filterval = [0,0,1,2,2]; % Replace with your desired values

initvolt = [0.5,0.5,0.5,0.5,0.5];
volumes =[1,1,1,1,0.5];
stresses = [1000E6,10E6,10E6,10E6,10E6]; % Replace with your desired values
filepath="TECTO/input_TECTO_Thermoel_Serend_240124_Journal2D_s.txt";
reader = InputReader(filepath);
pcon_loc = 1;
vcon_loc = 2;
qinloc = 9;
xinit = [0.5,0.5,0.5,0.4,0.4,0.4];
moveval = [0.5,0.5,0.5,0.02,0.02,0.02,0.02]; % Replace with your desired values
incrval = [1.3,1.3,1.3,1.02,1.02,1.02,1.02];
decreval = [0.7,0.7,0.7,0.98,0.98];
initval = [0.5,0.5,0.5,0.02,0.02];
Hevmu=0.4;
for i = 1:length(qinvals)
    reader.bcval(qinloc) = qinvals(i);
    reader.TopOpt_ConstraintValue(pcon_loc) = powvals(i);
    reader.TopOpt_ConstraintValue(vcon_loc) = volumes(i);
    reader.Filter = filterval(i);
    %reader.TopOpt_Objective = 'AverageTemperature';
    reader.KSUp = 3;
    reader.TopOpt_Initial_x = xinit(i);
    reader.TopOpt_ConstraintValue(3) = stresses(i);
    mesh = Mesh(reader);
    bcinit = BCInit(reader, mesh);
    close all
    Hmult = 3;
    Hev_init = 64;
    Hev_upt = 25;
    reader.Rst_name = append('BB_Journal2D_Q_', num2str(qinvals(i)), ...
        '_P', num2str(powvals(i)), ...
        '_F', num2str(filterval(i)), ...
        '_v', num2str(volumes(i)), ...
        '_x', num2str(xinit(i)), ...
        '_mv', num2str(moveval(i)),...
        '_Hs_', num2str(Hmult),...
        '_Hu_',num2str(Hev_upt),...
        '_H_',num2str(Hev_init),...
        '_V_',num2str(initvolt(i)),...
        '_Hm_',num2str(Hevmu));
    TO = TopOpt(reader, mesh);
    TO.mma_incr = incrval(i);
    TO.mma_decr = decreval(i);
    TO.mma_init = initval(i);
    TO.Helm_mult =Hmult ;
    TO.maxiter = 150;
    TO.Hev_update = Hev_upt;
    TO.Hev_init=Hev_init;
    TO.Hev_update = Hev_upt;
    TO.Hev_max=64;
    TO.Hev_mu = Hevmu;
    TO.xval(end)=initvolt(i);

    TO.xold1=TO.xval;
    TO.xold2=TO.xval;
    TO.runMMA(reader, mesh)
end


filepath="TECTO/input_TECTO_Thermoel_Serend_240124_Journal3D_s.txt";
reader = InputReader(filepath);
pcon_loc = 1;
vcon_loc = 2;
qinloc = 9;

xinit = [0.4,0.4];
volumes =[1,0.4];
qinvals = [7500,7500]; % Replace with your desired values
powvals = [0.0075,0.015 ]; % Replace with your desired values
filterval = [2,2]; % Repilace with your desired values
stresses = [10E6,10E6]; % Replace with your desired values
moveval = [0.02,0.02]; % Replace with your desired values
incrval = [1.02,1.02];
decreval = [0.98,0.98];
initval = [0.02,0.02];
initvolt = [0.5,0.5];


for i = 1:length(qinvals)
    reader.bcval(qinloc) = qinvals(i);
    reader.TopOpt_ConstraintValue(pcon_loc) = powvals(i);
    reader.TopOpt_ConstraintValue(vcon_loc) = volumes(i);
    reader.Filter = filterval(i);
    reader.KSUp = 3;
    reader.TopOpt_Initial_x = xinit(i);
    mesh = Mesh(reader);
    bcinit = BCInit(reader, mesh);
    close all
    Hmult = 3;
    Hevmu = 0.4;
    Hev_init = 64;
    Hev_upt = 300;
    reader.Rst_name = append('CCC_Journal3D_Q_', num2str(qinvals(i)), ...
        '_P', num2str(powvals(i)), ...
        '_F', num2str(filterval(i)), ...
        '_v', num2str(volumes(i)), ...
        '_xinit', num2str(xinit(i)), ...
        '_mv', num2str(moveval(i)),...
        '_mult_', num2str(Hmult),...
        '_Hu_',num2str(Hev_upt),...
        '_Hi_',num2str(Hev_init), ...
        '_Vi_',num2str(initvolt(i)));

    TO = TopOpt(reader, mesh);
    TO.mma_incr = incrval(i);
    TO.mma_decr = decreval(i);
    TO.mma_init = initval(i);
    TO.Helm_mult =Hmult ;
    TO.maxiter = 150;
    TO.Hev_update = Hev_upt;
    TO.Hev_init=Hev_init;
    TO.Hev_update = Hev_upt;
    TO.Hev_max=64;
    TO.Hev_mu = Hevmu;
    TO.xval(end)=initvolt(i);

    TO.xold1=TO.xval;
    TO.xold2=TO.xval;
    TO.runMMA(reader, mesh)
end

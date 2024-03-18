% Visualizing the Loss Landscape of neural nets
% Visualizing high-dimensional loss landscapes with Hessian directions
%custom_f = @(x) sum(x.^2); % replace this with your actual loss function
% Assume theta is your vector of parameters
    load('opt_test_mat.mat');

    % Define the folder paths
    srcFolder = 'src';
    utilsFolder = 'utils';
    
    % Generate paths including subfolders
    srcPath = genpath(srcFolder);
    utilsPath = genpath(utilsFolder);
    
    % Add the generated paths to the MATLAB path
    addpath(srcPath);
    addpath(utilsPath);
    

filename="TECTO/input_TECTO_Thermoel_Quad_cteF_Rnd.txt";
[N] = getTOEL(filename);
N=N+1;
theta = ones(1, N)*0.95;
xx1 = ones(1, N-1)*0.95;
xx(1:N-1)=0.2;

% Choose the direction vectors
delta = (0.5-(rand(1, N)))*0.05 ;
eta = (0.5-(rand(1, N)))*0.05;
%for i = 1:N
    % Convert random values to the desired range
%    delta(i) = (-1 + 2 * delta(i)) * xx(i) / 2;
%    eta(i) = (-1 + 2 * eta(i)) * xx(i) / 2;
%end
disp([delta(1),eta(1)])
% Compute the alpha and beta values
Np=7;
[alpha,beta] = meshgrid(linspace(0,1,Np),linspace(0,1,Np));

% Preallocate the values array
values = zeros(size(alpha));
values1 = zeros(size(alpha));
values2 = zeros(size(alpha));

% Compute the loss values along these directions
%values = arrayfun(@(a, b) custom_f(theta + a * delta + b * eta), alpha, beta);
% Compute the loss values along these directions in parallel
total_time_start = tic; % Start the timer for the entire loop

for i = 1:Np
    % Start the timer for the current outer loop iteration
    loop_start = tic;
    % Start the timer    
    for j = 1:Np
        a = alpha(i, j);
        b = beta(i, j);
        
        % Start the timer for the specific calculation
        t_start = tic;
        
        [values(i, j), fval] = custom_f(xx + a * delta + b * eta, filename);
        values1(i, j)=fval(1);
        values2(i, j)=fval(2);
        % Stop the timer for the specific calculation
        t_elapsed = toc(t_start);
        %save('DesignSpace_2312.mat', 'values');
        fprintf('Time for iteration %d, calculation %d: %f seconds\n', i, j, t_elapsed);
    end

    % Stop the timer for the current outer loop iteration
    loop_elapsed = toc(loop_start);
    fprintf('Time for iteration %d: %f seconds\n', i, loop_elapsed);
end


% Stop the timer for the entire loop
total_time_elapsed = toc(total_time_start);

% Calculate and display the average time per iteration
average_time_per_iteration = total_time_elapsed / Np^2;
fprintf('Average time per iteration: %f seconds\n', average_time_per_iteration);

% Get the current date and time
current_datetime = datestr(now, 'yyyymmdd_HHMMSS');

    % Save alpha, beta, and values
    % Append date and direction information to the data_filename
    data_filename = ['TECTO\Results\LossFun\loss_landscape_data_',current_datetime, ...
        '_delta_', num2str(delta(1)), '_eta_', num2str(eta(1)), '.mat'];

    % Save alpha, beta, and values
    save(data_filename, 'alpha', 'beta', 'values', 'values1', 'values2');

plotAndSaveLossLandscape(alpha, beta, values, delta, eta, [current_datetime,'_Temp_'])
plotAndSaveLossLandscape(alpha, beta, values1, delta, eta, [current_datetime,'_val1_'])
plotAndSaveLossLandscape(alpha, beta, values2, delta, eta, [current_datetime,'_val2_'])

plotAndSaveCustomLossLandscape(alpha, beta, values, delta, eta, [current_datetime,'_Temp_'])
plotAndSaveCustomLossLandscape(alpha, beta, values1, delta, eta, [current_datetime,'_val1_'])
plotAndSaveCustomLossLandscape(alpha, beta, values2, delta, eta, [current_datetime,'_val2_'])


function [loss,cons] = custom_f(x,filename)

        %filepath="TECTO/input_TECTO_Thermoel_Quad_cteF.txt";   
        filepath=filename;   
        reader = InputReader(filepath);

        mesh = Mesh(reader);
        TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            bcinit = BCInit(reader, mesh);
                if reader.Filter>0
                    filtering = Filtering(reader,mesh);
                end

                %% Filter densities
                for i = 1:length(reader.TObcval)
                    if(strcmp(reader.TObctype,'Voltage'))
                        Voltage_value=reader.TObcminval(i)+x(length(TOEL)+i)*(reader.TObcmaxval(i)-reader.TObcminval(i));
                        reader.TObcval(i)=Voltage_value;
                        reader.bcval(length(reader.bcval)-length(reader.TObcval)+i)=Voltage_value;
                    end
                end
                %mesh_1 = Mesh(reader);
                for i=1:length(TOEL)
                    mesh.elements_density(TOEL(i))=x(i);
                end
                if reader.Filter>0
                    filtering.filter_densities(reader,mesh)
                end

                %bcinit1 = BCInit(reader, mesh);
                solver = Solver(mesh, bcinit);

                for i=1:length(reader.TObcval)
                    nodes=mesh.retrieveNodalSelection(reader.TObcloc(i));
                    if(strcmp(reader.TObctype,'Voltage'))
                        solver.soldofs(nodes*2)=Voltage_value;
                    end
                end

            if strcmp(reader.solver,'NR')
                    residual_norm=solver.runNewtonRaphson(reader, mesh, bcinit);
                    odd_numbers = 1:2:length(solver.soldofs);
                    Tdofs=solver.soldofs(odd_numbers);
                    if residual_norm>10000  || not(isempty(Tdofs(Tdofs<0)))% divergence in NR catch
                        warning('NR DIVERGED, changing to Arc-len!!!');
                        for i=1:length(bcinit.dofs_free_)
                            df=bcinit.dofs_free_(i);
                            if mod(df, 2)==0
                                solver.soldofs(df)=0.0;
                            else 
                                solver.soldofs(df)=273.15;
                            end
                        end
                        funAL = @(t) solver.funArcLen(reader,bcinit,mesh,t);
                        [ufree] = arc_length_Crisfield(funAL,solver.soldofs(bcinit.dofs_free_));
                        solver.soldofs(bcinit.dofs_free_)=ufree(:,end);
                        if strcmp(reader.physics,'decoupledthermoelectromechanical')
                            % Extract the necessary variables
                            [solver.KStiff,KThermalLoad] = solver.Assembly_DecoupledThermoMech(reader, mesh);
                            Temperature_solution = solver.soldofs(1:2:end);
                            solver.KUT=KThermalLoad;
                            solver.loadVector_mech=solver.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                            solver.SolveLinearSystemInParallel(reader.physics,bcinit)
                        end
                    end
                end

            TOO = TO_Objectives(reader,mesh,bcinit);
            TOO.CalculateObjective(reader,mesh,solver)
            TOC = TO_Constraints(reader,mesh,bcinit);
            TOC.CalculateConstraint(reader,mesh,solver);

            %% New derivatives
            f0val = TOO.fval;
            fval = TOC.fval;

            loss = f0val;
            cons = fval;
end

function plotAndSaveLossLandscape(alpha, beta, values, delta, eta,name)
    % Get the current date and time
    % Plot the loss values
    figure;
    contourf(alpha, beta, values);
    colorbar;

    % Add title with random direction information
    title(sprintf('Loss Landscape\nRandom Direction: delta=%f, eta=%f', delta(1), eta(1)));

    % Save the figure
    fig_filename = ['TECTO\Results\LossFun\loss_landscape_plot', name, ...
        '_delta_', num2str(delta(1)), '_eta_', num2str(eta(1)), '.png'];
    saveas(gcf, fig_filename);


    % Display a message
    disp(['Figure saved as: ' fig_filename]);

    % Close the current figure
    close(gcf);
end

function plotAndSaveCustomLossLandscape(alpha, beta, values, delta, eta, name)

    % Plot another figure with specific size and axis locations
    fig = figure('Position', [100, 100, 800, 600]); % Adjust the size as needed

    % Create axes
    ax = axes(fig);

    % Specify the position of the axes
    set(ax, 'Position', [0 0 1 1]);

    % Plot contour plot
    contourf(ax, alpha, beta, values);
    colormap(inferno)

    % Save the second figure
    fig_filename_custom = ['TECTO\Results\LossFun\loss_landscape_custom_plot',name, ...
        '_delta_', num2str(delta(1)), '_eta_', num2str(eta(1)), '.png'];
    saveas(gcf, fig_filename_custom);

    % Display a message for the second figure
    disp(['Customized Figure saved as: ' fig_filename_custom]);
    %lose the current figure
    close(gcf);
end

% SCRIPT RESULTS POSTPROCESSING

% Define the filename of the CSV file
filename = 'C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\Results\Journal_Al2O3_nonlin_0topdispl_Q_3000_P15_F0_v1_2024_07_18_14_06\Journal_Al2O3_nonlin_0topdispl_Q_3000_P15_F0_v1MMA_2024_07_18_14_06_.csv';
imgfolder = 'C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\Results\animationtests';

% Read the CSV file into a table
data = readtable(filename);

% Display the contents of the table
disp(data);

% Define the video file name
videoFilename = 'animation.mp4';
% Check if the video file already exists and delete it if it does
if isfile(videoFilename)
    delete(videoFilename);
end

% Initialize video writer with specified frame rate (frames per second)
frameRate = 10; % Adjust frame rate as needed
writerObj = VideoWriter(videoFilename, 'MPEG-4');
writerObj.FrameRate = frameRate; % Set the frame rate
open(writerObj);

% Get a list of all JPG files in the folder and sort them alphabetically
imgFiles = dir(fullfile(imgfolder, '*.jpg'));
[~, idx] = sort({imgFiles.name});
imgFiles = imgFiles(idx);

minfact = 0.99;
maxfact = 1.01;
alpha_start = 0.1; % Starting alpha value
alpha_end = 1.0;   % Ending alpha value
alpha_steps = 10;  % Number of steps to increase alpha

% Iterate over increasing values of Iteration
for i = 1:length(data.Iteration)
    % Load the current image from the sorted list
    imgFile = fullfile(imgfolder, imgFiles(i).name);
    img = imread(imgFile);
    
    % Internal loop for smooth transparency transition
    for alpha_i = linspace(alpha_start, alpha_end, alpha_steps)
        % Create subplots for indexes [1:i]
        %createSubplots(data.Iteration(1:i), data.Power_right_value(1:i), data.Volume_right_value(1:i), data.Voltage_right_value(1:i), data.PnormTemperature(1:i), length(data.Iteration),...
        %    [min(data.Power_right_value) * minfact, max(data.Power_right_value) * maxfact], [min(data.Volume_right_value) * minfact, max(data.Volume_right_value) * maxfact], [min(data.Voltage_right_value) * minfact, max(data.Voltage_right_value) * maxfact], [min(data.PnormTemperature) * minfact, max(data.PnormTemperature) * maxfact],...
        %    alpha_i);
        
        % Create subplots and add the image
        createSubplotsWithImage(data.Iteration(1:i), data.Power_right_value(1:i), data.Volume_right_value(1:i), data.Voltage_right_value(1:i), data.PnormTemperature(1:i), length(data.Iteration),...
            [min(data.Power_right_value) * minfact, max(data.Power_right_value) * maxfact], [min(data.Volume_right_value) * minfact, max(data.Volume_right_value) * maxfact], [min(data.Voltage_right_value) * minfact, max(data.Voltage_right_value) * maxfact], [min(data.PnormTemperature) * minfact, max(data.PnormTemperature) * maxfact],...
            alpha_i, img, img);

        % Capture the current figure
        frame = getframe(gcf);
        
        % Write the frame to the video file
        writeVideo(writerObj, frame);
        
        % Close the current figure to clear for the next iteration
        close(gcf);
    end
end

% Close the video writer
close(writerObj);

close all;


function createSubplotsWithImage(x, y1, y2, y3, y4, xlen, y1len, y2len, y3len, y4len, alpha, img, imgprev)
    % Create a figure
    totalImageSize=[1000,600];
    figure('Color', 'k', 'Position', [100, 100, totalImageSize(1), totalImageSize(2)]); % Set figure background color to black

    %figure('Color', 'k'); % Set figure background color to black
    lineWidth = 1.5; % Default line width if not provided

    % Create subplots for the data
    subplot('Position', [0.7 0.55 0.28 0.28]); % Adjust position and size as needed
    hold on;
    xlim([1, xlen]); % Set x-axis limits from 0 to length of x
    ylim(y1len); % Set y-axis limits
    if length(x) > 1
        plot(x(1:end-1), y1(1:end-1), 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line
        plot(x, y1, 'w', 'Color', [1 1 1 alpha], 'LineWidth', lineWidth); % Transparent line
        scatter(x, y1, 'w', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    else
        plot(x, y1, 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line for single point
    end
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w'); % Set axes color
    xlabel('$\mathrm{Iteration}$', 'Color', 'w', 'Interpreter', 'latex'); % Set x-label color
    ylabel('$\mathrm{Power}$', 'Color', 'w', 'Interpreter', 'latex'); % Set y-label color

    % Subplot 2: y2 vs x
    subplot('Position', [0.7 0.1 0.28 0.28]); % Adjust position and size as needed
    hold on;
    xlim([1, xlen]); % Set x-axis limits from 0 to length of x
    ylim(y2len); % Set y-axis limits
    if length(x) > 1
        plot(x(1:end-1), y2(1:end-1), 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line
        plot(x, y2, 'w', 'Color', [1 1 1 alpha], 'LineWidth', lineWidth); % Transparent line
        scatter(x, y2, 'w', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    else
        plot(x, y2, 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line for single point
    end
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w'); % Set axes color
    xlabel('$\mathrm{Iteration}$', 'Color', 'w', 'Interpreter', 'latex'); % Set x-label color
    ylabel('$\mathrm{Volume}$', 'Color', 'w', 'Interpreter', 'latex'); % Set y-label color

    % Subplot 3: y3 vs x
    subplot('Position', [0.35 0.55 0.28 0.28]); % Adjust position and size as needed
    hold on;
    xlim([1, xlen]); % Set x-axis limits from 0 to length of x
    ylim(y3len); % Set y-axis limits
    if length(x) > 1
        plot(x(1:end-1), y3(1:end-1), 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line
        plot(x, y3, 'w', 'Color', [1 1 1 alpha], 'LineWidth', lineWidth); % Transparent line
        scatter(x, y3, 'w', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    else
        plot(x, y3, 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line for single point
    end
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w'); % Set axes color
    xlabel('$\mathrm{Iteration}$', 'Color', 'w', 'Interpreter', 'latex'); % Set x-label color
    ylabel('$\mathrm{Voltage}$', 'Color', 'w', 'Interpreter', 'latex'); % Set y-label color

    % Subplot 4: y4 vs x
    subplot('Position', [0.35 0.1 0.28 0.28]); % Adjust position and size as needed
    hold on;
    xlim([1, xlen]); % Set x-axis limits from 0 to length of x
    ylim(y4len); % Set y-axis limits
    if length(x) > 1
        plot(x(1:end-1), y4(1:end-1), 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line
        plot(x, y4, 'w', 'Color', [1 1 1 alpha], 'LineWidth', lineWidth); % Transparent line
        scatter(x, y4, 'w', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    else
        plot(x, y4, 'w-o', 'Color', [1 1 1 alpha], 'LineWidth', lineWidth, 'MarkerFaceColor', 'w'); % Transparent line
    end
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w'); % Set axes color
    xlabel('$\mathrm{Iteration}$', 'Color', 'w', 'Interpreter', 'latex'); % Set x-label color
    ylabel('$\mathrm{PnormTemperature}$', 'Color', 'w', 'Interpreter', 'latex'); % Set y-label color

    % Create an axes for the image
    hImg = axes('Position', [-0.05 0.35 0.35 0.35]); % Adjust position and size as needed
    imshow(img);
    set(gca, 'Visible', 'off'); % Hide the axes
end
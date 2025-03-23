% SCRIPT RESULTS POSTPROCESSING

% Define the filename of the CSV file
% Only Thermoel
filename = 'C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\Results\Conference\Journal_Al2O3_nonlin_0topdispl_Q_3000_P15_F0_v1_2024_07_19_15_20\thermoel.csv';
imgfolder = 'C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\Results\animationtests\thermoel';

% Filter
%filename = 'C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\Results\Conference\Journal_Al2O3_nonlin_0topdispl_Q_3000_P15_F1_v1_2024_07_21_12_34\thermoel.csv';
%imgfolder = 'C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\Results\animationtests\thermoelectromech_F';

% No Filter
%filename = 'C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\Results\Conference\Journal_Al2O3_nonlin_0topdispl_Q_3000_P15_F0_v1_2024_07_20_23_32\thermoel.csv';
%imgfolder = 'C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\Results\animationtests\thermoelectromech_noF';

% Bad penalization
%filename = 'C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\Results\Conference\Journal_Al2O3_nonlin_0topdispl_Q_3000_P15_F0_v1_2024_07_22_20_17\thermoel.csv';
%imgfolder = 'C:\Archive\Programming\MATLAB\TOTEM_M\TECTO\Results\animationtests\badcoeff';


% Read the CSV file into a table
data = readtable(filename);

% Display the contents of the table
disp(data);

% Define the video file name
videoFilename = 'animation0.mp4';
% Check if the video file already exists and delete it if it does
if isfile(videoFilename)
    delete(videoFilename);
end

% Initialize video writer with specified frame rate (frames per second)
frameRate = 30; % Adjust frame rate as needed
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
for i = 1:length(data.Iteration)-1
    % Load the current image from the sorted list
    imgFile = fullfile(imgfolder, imgFiles(i).name);
    img1 = imread(imgFile);
    
    imgFile = fullfile(imgfolder, imgFiles(i+1).name);
    img2 = imread(imgFile);
    % Internal loop for smooth transparency transition
    for alpha_i = linspace(alpha_start, alpha_end, alpha_steps)
        % Create subplots for indexes [1:i]
        %createSubplots(data.Iteration(1:i), data.Power_right_value(1:i), data.Volume_right_value(1:i), data.Voltage_right_value(1:i), data.PnormTemperature(1:i), length(data.Iteration),...
        %    [min(data.Power_right_value) * minfact, max(data.Power_right_value) * maxfact], [min(data.Volume_right_value) * minfact, max(data.Volume_right_value) * maxfact], [min(data.Voltage_right_value) * minfact, max(data.Voltage_right_value) * maxfact], [min(data.PnormTemperature) * minfact, max(data.PnormTemperature) * maxfact],...
        %    alpha_i);
        
        % Create subplots and add the image
        createSubplotsWithImage(data.Iteration(1:i), data.Power_right_value(1:i), data.Volume_right_value(1:i), data.Voltage_right_value(1:i), data.PnormTemperature(1:i), length(data.Iteration),...
            [min(data.Power_right_value) * minfact, max(data.Power_right_value) * maxfact], [min(data.Volume_right_value) * minfact, max(data.Volume_right_value) * maxfact], [min(data.Voltage_right_value) * minfact, max(data.Voltage_right_value) * maxfact], [min(data.PnormTemperature) * minfact, max(data.PnormTemperature) * maxfact],...
            alpha_i, img1, img2, ...
            [],[])
            %data.Stress_Pnorm_right_value(1:i),[min(data.Stress_Pnorm_right_value) * minfact, max(data.Stress_Pnorm_right_value) * maxfact]);

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


function createSubplotsWithImage(x, y1, y2, y3, y4, xlen, y1len, y2len, y3len, y4len, alpha, img1, img2, y5,y5len)
    % Create a figure
    totalImageSize = [1000, 600];
    figure('Color', 'k', 'Position', [100, 100, totalImageSize(1), totalImageSize(2)]); % Set figure background color to black

    % Default line width
    lineWidth = 1.5;

    % Subplot 1: y1 vs x
    subplot('Position', [0.7 0.55 0.27 0.27]); % Adjust position and size as needed
    hold on;
    xlim([1, xlen]);
    ylim(y1len);
    if length(x) > 1
        plot(x(1:end-1), y1(1:end-1), 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line
        plot(x, y1, 'w', 'Color', [1 1 1 alpha], 'LineWidth', lineWidth); % Transparent line
        scatter(x, y1, 'w', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    else
        plot(x, y1, 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line for single point
    end
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    xlabel('$\mathrm{Iteration}$', 'Color', 'w', 'Interpreter', 'latex');
    ylabel('$\mathrm{Power}$', 'Color', 'w', 'Interpreter', 'latex');

    % Subplot 2: y2 vs x
    subplot('Position', [0.7 0.1 0.27 0.27]);
    hold on;
    xlim([1, xlen]);
    ylim(y2len);
    if length(x) > 1
        plot(x(1:end-1), y2(1:end-1), 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line
        plot(x, y2, 'w', 'Color', [1 1 1 alpha], 'LineWidth', lineWidth); % Transparent line
        scatter(x, y2, 'w', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    else
        plot(x, y2, 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line for single point
    end
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    xlabel('$\mathrm{Iteration}$', 'Color', 'w', 'Interpreter', 'latex');
    ylabel('$\mathrm{Volume}$', 'Color', 'w', 'Interpreter', 'latex');

    % Subplot 3: y3 vs x
    subplot('Position', [0.37 0.55 0.27 0.27]);
    hold on;
    xlim([1, xlen]);
    ylim(y3len);
    if length(x) > 1
        plot(x(1:end-1), y3(1:end-1), 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line
        plot(x, y3, 'w', 'Color', [1 1 1 alpha], 'LineWidth', lineWidth); % Transparent line
        scatter(x, y3, 'w', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    else
        plot(x, y3, 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line for single point
    end
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    xlabel('$\mathrm{Iteration}$', 'Color', 'w', 'Interpreter', 'latex');
    ylabel('$\mathrm{Voltage}$', 'Color', 'w', 'Interpreter', 'latex');

    % Subplot 4: y4 vs x
    subplot('Position', [0.37 0.1 0.27 0.27]);
    hold on;
    xlim([1, xlen]);
    ylim(y4len);
    if length(x) > 1
        plot(x(1:end-1), y4(1:end-1), 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line
        plot(x, y4, 'w', 'Color', [1 1 1 alpha], 'LineWidth', lineWidth); % Transparent line
        scatter(x, y4, 'w', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    else
        plot(x, y4, 'w-o', 'Color', [1 1 1 alpha], 'LineWidth', lineWidth, 'MarkerFaceColor', 'w'); % Transparent line
    end
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    xlabel('$\mathrm{Iteration}$', 'Color', 'w', 'Interpreter', 'latex');
    ylabel('$\mathrm{PnormTemperature}$', 'Color', 'w', 'Interpreter', 'latex');

    % Conditional Subplot 5: y5 vs x (if y5 is not empty)
    if ~isempty(y5)
        subplot('Position', [0.05 0.55 0.27 0.27]);
        hold on;
        xlim([1, xlen]);
        ylim(y5len); % Set y5 limits
        if length(x) > 1
            plot(x(1:end-1), y5(1:end-1), 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line
            plot(x, y5, 'w', 'Color', [1 1 1 alpha], 'LineWidth', lineWidth); % Transparent line
            scatter(x, y5, 'w', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
        else
            plot(x, y5, 'w-o', 'MarkerFaceColor', 'w', 'LineWidth', lineWidth); % Solid line for single point
        end
        set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
        xlabel('$\mathrm{Iteration}$', 'Color', 'w', 'Interpreter', 'latex');
        ylabel('$\mathrm{Stress Pnorm}$', 'Color', 'w', 'Interpreter', 'latex');
    end

    % Ensure both images are the same size
    assert(all(size(img1) == size(img2)), 'Images must be the same size.');

    % Interpolate between the images
    img1_double = double(img1);
    img2_double = double(img2);
    interpolated_img_double = (1 - alpha) * img1_double + alpha * img2_double;
    interpolated_img = uint8(interpolated_img_double);

    % Create an axes for the image
    hImg = axes('Position', [0.0 0.1 0.325 0.325]); % Adjust position and size as needed
    imshow(interpolated_img, 'Parent', hImg);
    set(hImg, 'Visible', 'off'); % Hide the axes
end

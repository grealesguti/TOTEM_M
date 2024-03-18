
data_k=[304.24978602936335, 3.0080321285140563;
        322.3130225821318, 2.859437751004016;
        348.9778787280269, 2.734939759036145;
        373.1796036605438, 2.6305220883534135;
        397.76812166699585, 2.63855421686747;
        422.357462637435, 2.642570281124498;
        448.15902956086643, 2.7309236947791167;
        473.13845546118904, 2.8313253012048194;
        497.288333662519, 2.9799196787148596;
        523.480808479821, 3.1606425702811247];


%% Extracting data
x = data_k(:, 1);
y = data_k(:, 2);
% Example usage:
degree = 6;
xmin = 200;
xmax = 650;
ymin = 3.05;
ymax = 3.2;
minDerivative = 0;

[p, error] = fitPolynomialWithConstraint(x, y, degree, xmin, xmax, ymin, ymax, minDerivative);


function [p, error] = fitPolynomialWithConstraint(x, y, degree, xmin, xmax, ymin, ymax, minDerivative)
    % Objective function for polynomial fitting with constraint
    
    % Objective function: sum of squared errors
    objective = @(p) sum((polyval(p, x) - y).^2);
    
    % Constraint function: minimum absolute derivative at xmax
    constraint = @(p) minDerivative - abs(polyder(p, 1, length(p) - 1) * (xmax - x(end)));
    
    % Initial guess for polynomial coefficients
    initialGuess = zeros(1, degree + 1);
    
    % Polynomial coefficients are subject to the constraint
    options = optimset('fmincon');
    options.Display = 'off'; % Suppress display of optimization information
    
    % Perform constrained optimization
    [p, error] = fmincon(objective, initialGuess, [], [], [], [], [], [], constraint, options);
    
    % Display the results
    disp('Optimal Coefficients:');
    disp(p);
    disp(['Error: ', num2str(error)]);
    
    % Plotting the results
    figure;
    plot(x, y, 'o-', 'LineWidth', 1.5, 'DisplayName', 'Data');
    hold on;
    x_fit = linspace(min(x), max(x), 1000);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Fit');
    legend('show');
    title('Polynomial Fit with Constraint');
    xlabel('X');
    ylabel('Y');
    grid on;
end
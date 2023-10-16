clear;clc;
% Define the original function
%rho = @(T) -7.36e-13 * T^3 + 6.14e-10 * T^2 - 6.35e-8 * T - 1.78e-6;
rho = @(T) -6.83e-13*T^3 + 6.66e-10*T^2 - 1.55e-7*T + 1.81e-5;

% Define the specific temperature point around which you want to approximate
T0 = 300;  % You can choose a different value

% Define margins for the temperature range
lower_margin = T0 - 10;
upper_margin = T0 + 50;

% Create a set of temperature values within the specified range
T = lower_margin:0.1:upper_margin;

% Calculate the corresponding values of 1/rho(T)
sigma = 1 ./ arrayfun(rho, T);

% Fit a cubic polynomial to the data within the specified range
coefficients = polyfit(T - T0, sigma, 3);  % 3 for a cubic polynomial

% Print the coefficients of the fit
disp('Cubic Polynomial Coefficients:');
disp(coefficients);

% Define a new function handle for the cubic polynomial
cubic_poly = @(T) polyval(coefficients, T - T0);

% Print the final formula for the cubic polynomial
disp('Cubic Polynomial Formula:');
fprintf('sigma(T) = %.5e(T - %d)^3 + %.5e(T - %d)^2 + %.5e(T - %d) + %.5e\n', ...
    coefficients(1), T0, coefficients(2), T0, coefficients(3), T0, coefficients(4));

% Evaluate the cubic polynomial at all temperature points
approximated_values = cubic_poly(T);
% Plot the original data and the cubic fit
figure;
plot(T, sigma, 'b', 'LineWidth', 2); % Original data in blue
hold on;
plot(T, approximated_values, 'r--', 'LineWidth', 2); % Cubic fit in red dashed line
xlabel('Temperature');
ylabel('1/rho(T)');
legend('Original Data', 'Cubic Fit');
title('Cubic Fit for 1/rho(T)');
grid on;

syms a b c d Tx

% Define the cubic polynomial formula symbolically
cubic_poly = a * (Tx-T0)^3 + b * (Tx-T0)^2 + c * (Tx-T0) + d;
% Substitute the coefficients and the temperature point into the symbolic formula
cubic_poly = subs(cubic_poly, [a, b, c, d], coefficients);

% Display the simplified formula
disp(['Cubic Polynomial Formula:']);
eqn=eval(simplify(cubic_poly));
% Expand and simplify the equation
expanded_eqn = (expand(eqn))


% Collect terms in aTx^2 + bTx + c form
coeffs_eqn = coeffs(expanded_eqn, Tx);

% Extract coefficients
a = coeffs_eqn(1);  % Coefficient of Tx^3
b = coeffs_eqn(2);  % Coefficient of Tx^2
c = coeffs_eqn(3);  % Coefficient of Tx^1
d = coeffs_eqn(4);  % Constant term

ae = eval(coeffs_eqn(1));  % Coefficient of Tx^3
be = eval(coeffs_eqn(2));  % Coefficient of Tx^2
ce = eval(coeffs_eqn(3));  % Coefficient of Tx^1
de = eval(coeffs_eqn(4));  % Constant term

% Display the simplified equation
simplified_eqn = de * Tx^3 + ce * Tx^2 + be*Tx+ae

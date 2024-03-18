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
data_rhop =  [   197.6047904191617, 0.18379482334634734;
                223.35329341317365, 0.22615557747613835;
                247.9041916167665, 0.2754491017964078;
                277.245508982036, 0.3352946079997776;
                301.79640718562877, 0.4054698028536896;
                327.54491017964074, 0.4756727843616715;
                352.0958083832336, 0.5736902065937735;
                374.2514970059881, 0.6472901065618186;
                398.2035928143713, 0.7313725217777904];

data_rhon =  [  300.5988023952096, 0.45068563568917885;
                322.7544910179641, 0.48948275143448616;
                350.2994011976048, 0.542326020812204;
                373.0538922155689, 0.5950581435736417;
                398.8023952095808, 0.658300568237076;
                423.3532934131737, 0.7110743709796188;
                449.70059880239523, 0.7778109673923619;
                474.2514970059881, 0.8445058838239994;
                500, 0.9147088653319813;
                522.1556886227545, 0.9674270947663839];
S_p=[   198.79032258064518, 74.83443708609272;
        227.82258064516128, 85.43046357615896;
        253.22580645161293, 96.02649006622514;
        272.58064516129036, 101.98675496688742;
        299.7983870967742, 109.93377483443709;
        326.41129032258067, 120.5298013245033;
        355.44354838709677, 131.12582781456953;
        372.98387096774195, 139.0728476821192;
        401.41129032258067, 149.66887417218544;
        428.02419354838713, 152.9801324503311;
        446.1693548387097, 149.66887417218544];

S_n=[   305.0967741935484, -126.66666666666669;
        312.9758064516129, -129.33333333333337;
        329.9354838709678, -133.33333333333337;
        352.3467741935484, -138.66666666666669;
        371.72983870967744, -143.33333333333337;
        398.9717741935484, -147.33333333333337;
        418.3387096774194, -149.33333333333337;
        446.16532258064524, -150;
        466.72580645161304, -149.33333333333337;
        493.91935483870975, -145.33333333333337;
        513.8548387096776, -141.33333333333337];

%% modifying data K
% New data to append
data_start=[200,3.1];
data_end=[650,3.17];
data_k = [data_start; data_k; data_end];

data_start=[150,0.2];
data_end=[475,0.8];
data_rhop = [data_start; data_rhop; data_end];

data_start=[275,0.46];
data_end=[550,0.95];
data_rhon = [data_start; data_rhon; data_end];

data_start=[150,85];
data_end=[500,153];
S_p = [data_start; S_p; data_end];

data_start=[275,-130];
data_end=[550,-145];
S_n = [data_start; S_n; data_end];

%% Extracting data
temperature_k = data_k(:, 1);
thermal_conductivity_k = data_k(:, 2);

temperature_rhop = data_rhop(:, 1);
electrical_resistivity_rhop = 1./data_rhop(:, 2);

temperature_rhon = data_rhon(:, 1);
electrical_resistivity_rhon = 1./data_rhon(:, 2);

temperature_sp = S_p(:, 1);
seebeck_coefficient_sp = S_p(:, 2);

temperature_sn = S_n(:, 1);
seebeck_coefficient_sn = S_n(:, 2);

%% Plotting K
figure;
% Set x-axis limits
xmin = 150;
xmax = 700;

% Degree of the polynomial fit
degree = 6; % You can adjust this as needed

% Subplot 1: Thermal Conductivity
subplot(3, 1, 1);
plot(temperature_k, thermal_conductivity_k, 'o-', 'LineWidth', 2.5);
title('Thermal Conductivity', 'FontSize', 14);
xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Thermal Conductivity', 'FontSize', 14);
grid on;
xlim([xmin, xmax]); % Set x-axis limits
ylim([0, 10]); % Set x-axis limits


% Fit polynomial and plot
p1 = polyfit(temperature_k, thermal_conductivity_k, degree);
x_fit = linspace(xmin, xmax, 100);
y_fit = polyval(p1, x_fit);
hold on;
plot(x_fit, y_fit, 'r--', 'LineWidth', 2.5);
legend('Data', 'Polyfit');

% Using findExtrema function to find minima and maxima
poly_order = degree;
% Find extrema locations
poly_coefficients = p1;
extrema_locations = findExtrema(poly_order, poly_coefficients, xmin, xmax);
for i=1:length(extrema_locations)
    xline(extrema_locations(i), '--', 'Color', 'k','LineWidth', 2.5);
end
extrema_locations_k=extrema_locations;

k_value_minx = polyval(p1, 235.4822);
k_value_maxx = polyval(p1, 605.3485);

hold off;
%% Plotting Rho
% Subplot 2: Electrical Resistivity
subplot(3, 1, 2);
hold on
plot(temperature_rhop, electrical_resistivity_rhop, 'o-', 'LineWidth', 2.5);
plot(temperature_rhon, electrical_resistivity_rhon, 'o-', 'LineWidth', 2.5);
title('Electrical Resistivity', 'FontSize', 14);
xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Electrical Conductivity', 'FontSize', 14);
grid on;
xlim([xmin, xmax]); % Set x-axis limits
ylim([1, 6]); % Set x-axis limits


% Fit polynomials and plot
p2_rhop = polyfit(temperature_rhop, electrical_resistivity_rhop, degree);
p2_rhon = polyfit(temperature_rhon, electrical_resistivity_rhon, degree);
y_fit_rhop = polyval(p2_rhop, x_fit);
y_fit_rhon = polyval(p2_rhon, x_fit);
hold on;
plot(x_fit, y_fit_rhop, 'r--', 'LineWidth', 1.5);
plot(x_fit, y_fit_rhon, 'g--', 'LineWidth', 1.5);
legend('Data (rhop)', 'Data (rhon)', 'Polyfit (rhop)', 'Polyfit (rhon)');
poly_coefficients = p2_rhop;
extrema_locations = findExtrema(poly_order, poly_coefficients, xmin, xmax);
for i=1:length(extrema_locations)
    xline(extrema_locations(i), '--', 'Color', 'k','LineWidth', 2.5);
end
extrema_locations_rhop=extrema_locations;
poly_coefficients = p2_rhon;
extrema_locations = findExtrema(poly_order, poly_coefficients, xmin, xmax);
for i=1:length(extrema_locations)
    xline(extrema_locations(i), '--', 'Color', 'r','LineWidth', 2.5);
end
extrema_locations_rhon=extrema_locations;

hold off;

RHOP_value_minx = polyval(p2_rhop, 173.46);
RHOP_value_maxx = polyval(p2_rhop, 451.56);


RHON_value_minx = polyval(p2_rhon, 289.1785);
RHON_value_maxx = polyval(p2_rhon, 618.6265);

%% Subplot 3: Seebeck Coefficient
subplot(3, 1, 3);
hold on
plot(temperature_sp, seebeck_coefficient_sp, 'o-', 'LineWidth', 2.5);
plot(temperature_sn, seebeck_coefficient_sn, 'o-', 'LineWidth', 2.5);
title('Seebeck Coefficient', 'FontSize', 14);
xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Seebeck Coefficient', 'FontSize', 14);
ylim([-200, 200]); % Set x-axis limits

grid on;
xlim([xmin, xmax]); % Set x-axis limits

% Fit polynomials and plot
p3_sp = polyfit(temperature_sp, seebeck_coefficient_sp, degree);
p3_sn = polyfit(temperature_sn, seebeck_coefficient_sn, degree);
y_fit_sp = polyval(p3_sp, x_fit);
y_fit_sn = polyval(p3_sn, x_fit);
hold on;
plot(x_fit, y_fit_sp, 'r--', 'LineWidth', 2.5);
plot(x_fit, y_fit_sn, 'g--', 'LineWidth', 2.5);
legend('Data (sp)', 'Data (sn)', 'Polyfit (sp)', 'Polyfit (sn)');
poly_coefficients = p3_sp;
extrema_locations = findExtrema(poly_order, poly_coefficients, xmin, xmax);
extrema_locations_sp=extrema_locations;
for i=1:length(extrema_locations)
    xline(extrema_locations(i), '--', 'Color', 'k','LineWidth', 2.5);
end
poly_coefficients = p3_sn;
extrema_locations = findExtrema(poly_order, poly_coefficients, xmin, xmax);
for i=1:length(extrema_locations)
    xline(extrema_locations(i), '--', 'Color', 'r','LineWidth', 2.5);
end
extrema_locations_sn=extrema_locations;

hold off;

SN_value_minx = polyval(p3_sn, 289.1785);
SN_value_maxx = polyval(p3_sn, 529.4214);

SP_value_minx = polyval(p3_sp, 177.1991);
SP_value_maxx = polyval(p3_sp, 478.2327);

%% Adjusting layout
sgtitle('Material Properties vs. Temperature');

disp('Extrema Locations for k:');
disp(sort(extrema_locations_k));
p1=fliplr(p1); 
disp(p1)
disp('Extrema Locations for rhop:');
disp(sort(extrema_locations_rhop));
p2_rhop=fliplr(p2_rhop./0.00001); 
disp(p2_rhop)
disp('Extrema Locations for rhon:');
disp(sort(extrema_locations_rhon));
p2_rhon=fliplr(p2_rhon./0.00001); 
disp(p2_rhon)
disp('Extrema Locations for sp:');
disp(sort(extrema_locations_sp));
p3_sp=fliplr(p3_sp*1e-6); 
disp(p3_sp)
disp('Extrema Locations for snn:');
disp(sort(extrema_locations_sn));
p3_sn=fliplr(p3_sn*1e-6); 
disp(p3_sn)




function extrema_locations = findExtrema(poly_order, poly_coefficients, xmin, xmax)
    % Input:
    %   - poly_order: The order of the polynomial
    %   - poly_coefficients: Coefficients of the polynomial, starting from the highest order
    %   - xmin, xmax: Range to search for extrema

    % Create the polynomial function
    poly_func = @(x) polyval(poly_coefficients, x);

    % Derivative of the polynomial function
    poly_derivative = polyder(poly_coefficients);

    % Find all roots (extrema) using the roots function
    extrema_roots = roots(poly_derivative);

    % Filter out roots outside the specified range
    extrema_locations = extrema_roots(extrema_roots >= xmin & extrema_roots <= xmax);

    % Filter out complex roots with a complex component higher than 0.001
    is_real = abs(imag(extrema_locations)) < 0.001;
    extrema_locations = extrema_locations(is_real);
end



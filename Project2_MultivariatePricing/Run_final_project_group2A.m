% Final FE PRoject
% Group 2A
% 17/05/2024
%
% Description:

% Marco Maspes
% Matteo Torba
% Andrea Tarditi

clear all; close all; clc;

%% load fol
addpath('data');
addpath('forward price');
addpath('calibration');

%% Loading of the matrices
% Loading of the matrices necessary for the projects

data_USA = load("OptionData.mat").mkt;
data_EU = load("OptionData.mat").mkt_EU;

SP500_EUR500 = load("SPXSX5Ereturns.mat").Returns;

%% 

%% POINT 5: Forward Prices computation
% Choice of the flags: flag = 0 [EUROPEAN], flag = 1 [AMERICAN]
flag = 0;

if flag == 1
    data = data_USA;
else
    data = data_EU;
end

% for i = 1:length(data.datesExpiry)

%     date = data.datesExpiry(i);

%     [F_vector, G_vector] = forward_prices(data, date);
% end

%% POINT 6: Model calibration

rho_mkt = zeros(length(data_EU.datesExpiry), 1);

for i = 1:length(data_EU.datesExpiry)

    date = data_EU.datesExpiry(i);

    % compute correlation coefficient between the two series
    rho_mkt(i) = compute_corr_coeff(data_EU,data_USA,date);
end

%% Creation of the constraints for the simulations

% General function handles for the writing

sigma2_i = @(g2_i, a_i, gz) g2_i + a_i^2 * gz;
theta_i = @(beta_i, a_i, bz) beta_i + a_i * bz;
k_i = @(nu_i, nz) (nu_i * nz)/(nu_i + nz);

% Calibrationà

fun = @(a1, a2, b1, n1, g1, b2, n2, g2, bz, nz, gz) ...
    a1*a2*(gz + bz^2 * nz)/ (sqrt(sigma2_i(g1, a1, gz) + k_i(n1, nz) * (theta_i(b1, a1, bz)^2)) * ...
     sqrt(sigma2_i(g2, a2, gz) + k_i(n2, nz) * (theta_i(b2, a2, bz)^2)));

x0 = zeros(11, 1);
lb = [-Inf; -Inf; -Inf; 0; 0; -Inf; 0; 0; -Inf; 0; 0];
ub = [];
options = optimset('Display', 'Off');

params = fmincon(@(a1, a2, b1, n1, g1, b2, n2, g2, bz, nz, gz) fun(a1, a2, b1, n1, g1, b2, n2, g2, bz, nz, gz) - rho_mkt(1), ...
    x0, lb, ub, [], [], [], [], @(a1, a2, b1, n1, g1, b2, n2, g2, bz, nz, gz) nonlinconstr(a1, a2, b1, n1, g1, b2, n2, g2, bz, nz, gz), options);




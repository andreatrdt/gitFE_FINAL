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

sigma2_i = @(gamma2_i, a_i, gamma_z) gamma2_i + a_i^2 * gamma_z;
theta_i = @(beta_i, a_i, beta_z) beta_i + a_i * beta_z;
k_i = @(nu_i, nu_z) (nu_i * nu_z)/(nu_i + nu_z);

% Calibrationà

fun = @(a_1, a_2, beta_1, nu_1, gamma_1, beta_2, nu_2, gamma_2, ...
    beta_z, nu_z, gamma_z) a_1*a_2*(gamma_z + beta_z^2 * nu_z)/ ...
    (sqrt(sigma2_i(gamma_1, a_1, gamma_z) + k_i(nu_1, nu_z) * (theta_i(beta_1, a_1, beta_z)^2)) * ...
     sqrt(sigma2_i(gamma_2, a_2, gamma_z) + k_i(nu_2, nu_z) * (theta_i(beta_2, a_2, beta_z)^2)));

x0 = zeros(11, 1);
lb = [-Inf; -Inf; -Inf; 0; 0; -Inf; 0; 0; -Inf; 0; 0];
ub = [];
options = optimset('Display', 'Off');

params = lsqnonlin(@(a_1, a_2, beta_1, nu_1, gamma_1, beta_2, nu_2, gamma_2, beta_z, nu_z, gamma_z) fun(a_1, a_2, beta_1, nu_1, gamma_1, beta_2, nu_2, gamma_2, ...
    beta_z, nu_z, gamma_z) - rho_mkt(1), x0, lb, ub, [], [], [], [], ...
    @(a_1, a_2, beta_1, nu_1, gamma_1, beta_2, nu_2, gamma_2, beta_z, nu_z, gamma_z) nonlinconstr(a_1, a_2, beta_1, nu_1, gamma_1, beta_2, nu_2, gamma_2, ...
    beta_z, nu_z, gamma_z), options);

%



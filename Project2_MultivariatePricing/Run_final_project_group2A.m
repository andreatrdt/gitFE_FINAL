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

%     [F_vector, G_vector , ~] = forward_prices(data, date);
% end

%% POINT 6: Model calibration

% rho_mkt = zeros(length(data_EU.datesExpiry), 1);

% for i = 1:length(data_EU.datesExpiry)

%     date = data_EU.datesExpiry(i);

%     % compute correlation coefficient between the two series
%     rho_mkt(i) = compute_corr_coeff(data_EU,data_USA,date);
% end

%% Creation of the constraints for the simulations

%% calibration 
alpha = 1/2; %(NIG model)

M = 15;
dz = 0.0025;
t = 1;

data = data_EU;

idx=1;

date = data.datesExpiry(idx);

[F_vector, G_vector , discount_1y ] = forward_prices(data, date);

F_0 = F_vector(1,:);

% compute the log moneyess from the strikes
log_moneyness = log( F_0 ./ data.strikes(idx).value);

% create the distance function to minimize
dist = @(x) abs(sum((E3_callPriceLewis(discount_1y, F_0, alpha, log_moneyness, x(1), x(2), x(3), t, 2) - data.callAsk(idx).prices)));

% create the constraint
const = @(x) constraint(x, alpha);

% calibrate the model using fmincon
% initial guess
x0 = [0.1, 1, 3];

A = [
    -1, 0, 0;
    0, -1, 0;
];
b = [
    0;
    0;
];

% calibration
options = optimoptions('fmincon', 'Display', 'iter');

[x, fval] = fmincon(dist, x0, A, b, [], [], [], [], const, options)



% % General function handles for the writing

% sigma2_i = @(g2_i, a_i, gz) g2_i + a_i^2 * gz;
% theta_i = @(beta_i, a_i, bz) beta_i + a_i * bz;
% k_i = @(nu_i, nz) (nu_i * nz)/(nu_i + nz);

% % Calibration

% fun = @(a1, a2, b1, n1, g1, b2, n2, g2, bz, nz, gz) ...
%     a1*a2*(gz + bz^2 * nz)/ (sqrt(sigma2_i(g1, a1, gz) + k_i(n1, nz) * (theta_i(b1, a1, bz)^2)) * ...
%      sqrt(sigma2_i(g2, a2, gz) + k_i(n2, nz) * (theta_i(b2, a2, bz)^2)));

% x0 = zeros(11, 1);
% lb = [-Inf; -Inf; -Inf; 0; 0; -Inf; 0; 0; -Inf; 0; 0];
% ub = [];
% options = optimset('Display', 'Off');

% params = fmincon(@(a1, a2, b1, n1, g1, b2, n2, g2, bz, nz, gz) fun(a1, a2, b1, n1, g1, b2, n2, g2, bz, nz, gz) - rho_mkt(1), ...
%     x0, lb, ub, [], [], [], [], @(a1, a2, b1, n1, g1, b2, n2, g2, bz, nz, gz) nonlinconstr(a1, a2, b1, n1, g1, b2, n2, g2, bz, nz, gz), options);




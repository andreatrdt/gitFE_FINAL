% Final FE Project
% Group 2A
%
% Description:

% Marco Maspes
% Matteo Torba
% Andrea Tarditi

% start run time
tic;

%% Clearing the workspace
clear all; close all; clc;

%% Load folders

addpath('data');
addpath('forward price');
addpath('calibration');
addpath('general');

%% Loading of the matrices
% Loading of the matrices necessary for the projects

data_USA = load("OptionData.mat").mkt;
data_EU = load("OptionData.mat").mkt_EU;

SP500_EUR500 = load("SPXSX5Ereturns.mat").Returns;

%% Initialization of the base parameters

% Settlement date:
date_settlement = datenum('2023-07-09');

% Dates Vector:
dates_EU = datenum(data_EU.datesExpiry);
dates_USA = datenum(data_USA.datesExpiry);

%% Year frac conventions

conv_ACT360 = 2; conv_ACT365 = 3; conv_30360_EU = 6;

%% Plot of the returns
 
% plot_returns(SP500_EUR500,date_settlement)

%% POINT 5: Forward Prices 

[F0_EU, B_bar_EU] = forward_prices(data_EU, 0);
[F0_USA, B_bar_USA] = forward_prices(data_USA, 0);

%% POINT 6: Calibration

%% Options selection

% Choice of only OTM options for the further calibration
data_EU_OTM = OTM_preprocessing(data_EU, B_bar_EU, F0_EU);
data_USA_OTM = OTM_preprocessing(data_USA, B_bar_USA, F0_USA);

% Computing the Delta of Black & Scholes over the OTM Call/Put in order to
% clean dataset from too far from the ATM point prices

data_calib_EU = dataset_preprocessing(data_EU_OTM, F0_EU, B_bar_EU, date_settlement, 0);
data_calib_USA = dataset_preprocessing(data_USA_OTM, F0_USA, B_bar_USA, date_settlement, 0);

%% Plot of the surface of the implied volatilities

% surface_vols(data_calib_EU,F0_EU);
% surface_vols(data_calib_USA,F0_USA);

%% Calibration of the model parameters

% Quantities of interest
% x0 = [10 2 0.5 10 2 0.5];
x0 = [5 2 0.5 5 2 0.5];
% x0 = [32 0.04 0.36 11.8 0.09 0.37];
% x0 = 0.5 * ones(1, 6);

% Linear inequality constraints 
A = [-1 0 0 0 0 0;
    0 0 -1 0 0 0;
    0 0 0 -1 0 0;
    0 0 0 0 0 -1];

b = [0; 0; 0; 0];

% Linear equality constraints
Aeq = []; beq = [];

% Lower and upper bounds 
lb = [0; -Inf; 0; 0; -Inf; 0];
ub = [];

% Options for the visualization
options = optimset('MaxFunEvals', 3e3, 'Display', 'iter');

% Calibration
params_marginals = fmincon(@(params) new_calibration(params, data_calib_EU, data_calib_USA, ...
    F0_EU, B_bar_EU, F0_USA, B_bar_USA, date_settlement), x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr(params), options);

% Display of the parameters on the console
params_USA = params_marginals(1:3)
params_EU = params_marginals(4:6)

%% Plots of the prices with calibrated values

% Plot over the entire curve
% plot_calls_puts_total(data_EU, F0_EU, B_bar_EU, params_EU, date_settlement);
% plot_calls_puts_total(data_USA, F0_USA, B_bar_USA, params_USA, date_settlement);

% Plot over the filtered options
% plot_calls_puts(data_calib_EU, F0_EU, B_bar_EU, params_EU, date_settlement);
% plot_calls_puts(data_calib_USA, F0_USA, B_bar_USA, params_USA, date_settlement);

% Plot the implied volatilities over the Calls
plot_volatility_smiles(data_calib_EU, F0_EU, B_bar_EU, params_EU, date_settlement)
plot_volatility_smiles(data_calib_USA, F0_USA, B_bar_USA, params_USA, date_settlement)

%% 2nd Calibration over the rho

% Computation of the historical correlation
rho = hist_corr(SP500_EUR500);

% Initialization of the parameters
A = []; b = []; Aeq = []; beq = [];
lb = 0; ub = [];
x0 = 1;

% Recall the parameters
k1 = params_USA(1); k2 = params_EU(1);

% Calibration of the nuZ parameter
nu_z_single = fmincon(@(nu_z) (sqrt(k1 * k2) / nu_z - rho)^2, ...
    x0, A, b, Aeq, beq, lb, ub, [], options)

% Initialization of the parameters
A = []; b = []; Aeq = []; beq = [];
lb = zeros(1, 3); ub = [];
x0 = ones(1, 3);

% Recall the parameters
k1 = params_USA(1); k2 = params_EU(1);

% Calibration of the nuZ parameter
params_singles = fmincon(@(params) (sqrt(params(1) * params(2) / ((params(1) + nu_z_single)*(params(2) + nu_z_single))) - rho)^2, ...
    x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr_corr([params, nu_z_single], k1, k2), options);

nu_1 = params_singles(1)
nu_2 = params_singles(2)

%% Alternative calibration over the rho

% Computation of the historical correlation
rho = hist_corr(SP500_EUR500);

% Initialization of the parameters
A = []; b = []; Aeq = []; beq = [];
lb = zeros(1, 3); ub = [];
x0 = ones(1, 3);

% Recall the parameters
k1 = params_USA(1); k2 = params_EU(1);

% Calibration of the nuZ parameter
params = fmincon(@(params) (sqrt(params(1) * params(2) / ((params(1) + params(3))*(params(2) + params(3)))) - rho)^2, ...
    x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr_corr(params, k1, k2), options);

nu_1 = params(1)
nu_2 = params(2)
nu_z = params(3)

elapsed_time = toc;
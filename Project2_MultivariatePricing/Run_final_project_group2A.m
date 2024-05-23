% Final FE Project
% Group 2A
%
% Description:

% Marco Maspes
% Matteo Torba
% Andrea Tarditi

clear all; close all; clc;

%% Load folders

addpath('data');
addpath('forward price');
addpath('calibration');
addpath('general')

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

%% plot
 
% plot_returns(SP500_EUR500,date_settlement)

%%
%% POINT 5: Forward Prices 

[F0_EU, B_bar_EU] = forward_prices(data_EU, 0);
[F0_USA, B_bar_USA] = forward_prices(data_USA, 0);

%% POINT 6: Calibration

%% Options selection

% Choice of only OTM options for the further calibration
data_EU = OTM_preprocessing(data_EU, B_bar_EU);
data_USA = OTM_preprocessing(data_USA, B_bar_USA);

% Computing the Delta of Black & Scholes over the OTM Call/Put in order to
% clean dataset from too far from the ATM point prices

data_EU = dataset_preprocessing(data_EU, F0_EU, B_bar_EU, date_settlement, 0);
data_USA = dataset_preprocessing(data_USA, F0_USA, B_bar_USA, date_settlement, 0);

%% Calibration of the model parameters

x0 = [20 0.2 0.1 30 0.5 0.3];

% % Quantities of interest
% x0 = [10 2 0.5 10 2 0.5];

% Linear inequality constraints 
A = []; b = [];

% Linear equality constraints
Aeq = []; beq = [];

% Lower and upper bounds 
lb = [0; -Inf; 0; 0; -Inf; 0];
ub = [];

% Options for the visualization
options = optimset('Display', 'iter');

% Calibration
params = fmincon(@(params) new_calibration(params, data_EU, data_USA, ...
    F0_EU, B_bar_EU, F0_USA, B_bar_USA, date_settlement), x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr(params), options)

% params(T=1) =   0.1300   19.3334    0.1200    0.0257   42.8258    0.1182
% params(T=2) =    1.1832   10.0741    0.1473    0.3111   17.4346    0.1307
% CAMBIO DELLA FUNZIONE OBJ (prima avevo dimenticato la radice, quindi era un MSE pesato)
% params(T=3) =   1.0e+03 *    0.0000   -0.6944    0.0002    4.0417    0.0010    0.0042

% load('params_II.mat');
% 
% log_moneyness = log(F0_EU(1) ./ data_EU.strikes(1).value);
% TTM = yearfrac(date_settlement, datenum(data_EU.datesExpiry(1)), conv_ACT365);
% prices = E3_callPriceLewis(B_bar_EU(1), F0_EU(1), 1/2, log_moneyness, params(6), params(4), params(5), TTM, 1, 19, 0.0048)
% mean_call_price = (data_EU.callAsk(1).prices + data_EU.callBid(1).prices)/2

%% 2nd Calibration over the rho

% rho = hist_corr(SP500_EUR500);
% 
% % Initialization of the parameters
% A = []; b = []; Aeq = []; beq = [];
% lb = 0; ub = [];
% x0 = 1;
% 
% % Calibration of the nuZ parameter
% nu_Z = fmincon(@(nu_Z) abs(sqrt(params(1) * params(4))/nu_Z - rho), x0, A, b, Aeq, beq, lb, ub, [], options)

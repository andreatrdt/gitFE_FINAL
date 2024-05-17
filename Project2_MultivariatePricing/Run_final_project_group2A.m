% Final FE PRoject
% Group 2A
% 17/05/2024
%
% Description:

% Marco Maspes
% Matteo Torba
% Andrea Tarditi

clear all; close all; clc;

%% load folders
addpath('data');
addpath('forward price');

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

date = data.datesExpiry(1);

% compute correlation coefficient between the two series
rho_mkt = compute_corr_coeff(data_EU,data_USA,date);

% parameters for the calibration
% sigma , kappa , vega

rho = @(vega_1, vega_2 , vega_z) sqrt(vega_1*vega_2/((vega_z+vega_1)*(vega_z+vega_2)));


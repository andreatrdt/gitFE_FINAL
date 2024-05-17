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

if flag
    data = data_USA;
else
    data = data_EU;
end

for i = 1:length(data.datesExpiry)

    date = data.datesExpiry(i);

    [F_vector, G_vector] = forward_prices(data, date);
end

%% POINT 6: Model calibration








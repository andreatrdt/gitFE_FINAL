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

%%
%% POINT 5: Forward Prices 

[F0_EU, B_bar_EU] = forward_prices(data_EU, 0);
[F0_USA, B_bar_USA] = forward_prices(data_USA, 0);

%% POINT 6: Calibration

%% Options selection

% Choice of only OTM options for the further calibration
data_EU = OTM_preprocessing(data_EU, B0_EU);
data_USA = OTM_preprocessing(data_USA, B0_USA);

% Computing the Delta of Black & Scholes over the OTM Call/Put in order to
% clean dataset from too far from the ATM point prices

data_EU = dataset_preprocessing(data_EU, F0_EU, B_bar_EU, date_settlement, 0);
data_USA = dataset_preprocessing(data_USA, F0_USA, B_bar_USA, date_settlement, 0);


%% Calibration of the model parameters

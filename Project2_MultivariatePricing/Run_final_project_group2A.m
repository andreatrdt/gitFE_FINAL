% Final FE PRoject
% Group 2A
% 17/05/2024
%
% Description:

% Marco Maspes
% Matteo Torba
% Andrea Tarditi

clear all; close all; clc;

%% Loading of the matrices
% Loading of the matrices necessary for the projects

data_USA = load("OptionData.mat").mkt;
data_EU = load("OptionData.mat").mkt_EU;

SP500_EUR500 = load("SPXSX5Ereturns.mat").Returns;

%% Prova punto 5

[Gi, Gi_ask, Gi_bid] = Synthethic_Forward(data_EU.callBid, data_EU.callAsk, data_EU.putBid, data_EU.putAsk, 1);
Ki = data_EU.strikes(1).value;

B_bar = estimation_discount_factor(Gi, Ki);

F = Gi/B_bar + Ki;
F_ask = Gi_ask/B_bar + Ki;
F_bid = Gi_bid/B_bar + Ki;

figure();
plot(Ki, F_ask, '*-'); hold on;
plot(Ki, F, '*-'); hold on;
plot(Ki, F_bid, '*-'); 






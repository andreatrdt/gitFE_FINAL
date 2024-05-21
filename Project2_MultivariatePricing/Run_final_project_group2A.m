% Final FE PRoject
% Group 2A
% 17/05/2024
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

% Settlement date:
formatData = 'YYYY-MM-DD';
date_settlement = datenum('2023-07-09');

% Dates Vector:
dates_EU = datenum(data_EU.datesExpiry);
dates_USA = datenum(data_USA.datesExpiry);

% Act365 convenction for yearfrac:
Act365 = 3;

%% POINT 5: Forward Prices computation
% Choice of the flags: flag = 0 [EUROPEAN], flag = 1 [AMERICAN]
flag = 0;

if flag == 1
    data = data_USA;
else
    data = data_EU;
end

for i = 1:length(data.datesExpiry)

    date = data.datesExpiry(i);

    [F_vector, G_vector , ~] = forward_prices(data, date, 1);
end

%% POINT 6: Model calibration

% Computation of the correlation of the market to impose the non linear
% constraint in the calibration

rho_mkt = zeros(length(data_EU.datesExpiry), 1);

for i = 1:length(data_EU.datesExpiry)

    date = data_EU.datesExpiry(i);
    rho_mkt(i) = compute_corr_coeff(SP500_EUR500, date_settlement, date); 
end

%% Joint calibration
alpha = 1/2; % (NIG model)
idx = 1;

% EU:
data = data_EU;
date = data.datesExpiry(idx);

% compute the forward in 0:
[F_0_EU, ~ , discount_at_expiry_EU] = forward_prices(data, date, 1);

% compute the log moneyess from the strikes
log_moneyness = log(F_0_EU(1,:) ./ data.strikes(idx).value);

% time to maturity
t = yearfrac(date_settlement,dates_EU(idx),Act365);

% create a function that the prices of the call options given the strikes
prices_EU = @(p) callIntegral(discount_at_expiry_EU, F_0_EU(1,:), p(1), p(2), p(3), t, log_moneyness); % CORRECT and use integral

% compute the implied volatilities:
volatility_EU = @(p) blkimpv(F_0_EU(1,:), data.strikes(idx).value, -log(discount_at_expiry_EU)/t, t, prices_EU(p));


% USA:
data = data_USA;
date = data.datesExpiry(idx);

% compute the forward in 0:
[F_0_USA, ~ , discount_at_expiry_USA] = forward_prices(data, date, 1);

% compute the log moneyess from the strikes
log_moneyness = log(F_0_USA(1,:) ./ data.strikes(idx).value);

% time to maturity
t = yearfrac(date_settlement,dates_USA(idx),Act365);

% create a function that the prices of the call options given the strikes
prices_USA = @(p) callIntegral(discount_at_expiry_USA, F_0_USA(1,:), p(1), p(2), p(3), t, log_moneyness);

% compute the implied volatilities:
volatility_USA = @(p) blkimpv(F_0_USA(1,:), data.strikes(idx).value, -log(discount_at_expiry_USA)/t, t, prices_USA(p));


% create the distance function to minimize
% dist = @(p_EU,p_USA) 1/length(data_EU.callAsk(idx).impvol)*sum((volatility_EU(p_EU)' - data_EU.callAsk(idx).impvol).^2) + 1/length(data_USA.callAsk(idx).impvol)*sum((volatility_USA(p_USA)' - data_USA.callAsk(idx).impvol).^2);
mean_call_price_EU = (data_EU.callAsk(idx).prices+data_EU.callBid(idx).prices)/2;
mean_call_price_USA = (data_USA.callAsk(idx).prices+data_USA.callBid(idx).prices)/2;

weights_EU = data_EU.Volume_call(idx).volume./sum(data_EU.Volume_call(idx).volume);
weights_USA = data_USA.Volume_call(idx).volume./sum(data_USA.Volume_call(idx).volume);

dist = @(p_EU,p_USA) 1/length(data_EU.callAsk(idx).prices)*sum(weights_EU' .* (prices_EU(p_EU) - mean_call_price_EU).^2) + ...
    1/length(data_USA.callAsk(idx).prices)*sum(weights_USA' .* (prices_USA(p_USA) - mean_call_price_USA).^2);


% calibrate the model using fmincon
% initial guess
% Quantities of interest
esp_thr = 1e-1;

% Initial values for the initialization
x0 = 1 * ones(11, 1);
% x0 = ones(6, 1);
% x0 = 0.01*ones(6,1);

% Linear inequality constraints on the theta_i
A = [];

% Plain term for the previous matrix
% b = zeros(4, 1);
b = [];

% Unused inequality matrixies
Aeq = []; beq = [];

% Bounds for the single parameter, no ub required
lb = [0; -Inf; 0; 0; -Inf; 0; 0; -Inf; 0; -Inf; -Inf];
% lb = [0; -Inf; 0; 0; -Inf; 0];
ub = [];

% Options for the visualization
options = optimset('Display', 'iter');

x = fmincon(@(x) dist([x(1) x(2) x(3)],[x(4) x(5) x(6)]), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconstr(x, rho_mkt(1), esp_thr), options);

% 
%     0.0094
%    83.6219
%     0.1058
%     0.0590
%    10.8153
%     0.1121
%    16.7909
%    56.5838
%     3.0280
%     0.0008
%     0.3834

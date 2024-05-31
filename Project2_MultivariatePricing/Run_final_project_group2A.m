% Final FE Project
% Group 2A
%
% Description:

% Marco Maspes
% Matteo Torba
% Andrea Tarditi


%% Clearing the workspace
clear all; close all; clc;

% start run time
tic;

% fix random seed
rng(2);

%% Flag

% flag = 1 for plotting enabled
% flag = 0 for plotting disabled

flag = 0;

% save results

save_results = 1;

%% Load folders

addpath('data');
addpath('forward price');
addpath('calibration');
addpath('general');
addpath('pricing_certificate');
addpath('plotting');

%% Loading of the matrices
% Loading of the matrices necessary for the projects

data_USA = load("OptionData.mat").mkt;
data_EU = load("OptionData.mat").mkt_EU;

SP500_EUR50 = load("SPXSX5Ereturns.mat").Returns;

%% Initialization of the base parameters

% Settlement date:
date_settlement = datenum('2023-07-09');

% Dates Vector:
dates_EU = datenum(data_EU.datesExpiry);
dates_USA = datenum(data_USA.datesExpiry);

% Computation of the historical correlation
rho_historical = hist_corr(SP500_EUR50);

%% Year frac conventions

conv_ACT360 = 2; conv_ACT365 = 3; conv_30360_EU = 6;

%% Plot of the penny and initial options

if flag == 1
    dataset_exploration(data_EU, data_USA, date_settlement)
end

%% Plot of the returns

if flag == 1
    plot_returns(SP500_EUR50,date_settlement)
end

%% POINT 5: Forward Prices

[F0_EU, B_EU] = forward_prices(data_EU, flag);
[F0_USA, B_USA] = forward_prices(data_USA, flag);

if flag
    figure;
    plot(dates_EU,B_EU,'-*','Color','b'); grid on;
    title('Discount factor for the EU market');
    ylabel('Discounts');
    datetick('x','dd-mmm-yyyy','keepticks')


    figure;
    plot(dates_USA,B_USA,'-*','Color','r'); grid on;
    title('Discount factor for the USA market');
    ylabel('Discounts');
    datetick('x','dd-mmm-yyyy','keepticks')
end

if flag
    figure;
    plot(dates_EU, F0_EU,'-*', 'Color', 'b'); grid on;
    title('Forwards value EU');
    ylabel('Forwards');
    datetick('x','dd-mmm-yyyy','keepticks');

    yf = yearfrac(date_settlement, data_EU.datesExpiry, conv_ACT365);
    rates = -log(B_EU)./yf;
    F0_sim_EU = data_EU.spot .* exp(rates .* yf);
    hold on;  plot(dates_EU, F0_sim_EU, '-o', 'Color', 'r');

    figure;
    plot(dates_USA, F0_USA, '-*', 'Color', 'b'); grid on;
    title('Forwards value USA');
    ylabel('Forwards');
    datetick('x','dd-mmm-yyyy','keepticks');

    yf = yearfrac(date_settlement, data_USA.datesExpiry, conv_ACT365);
    rates = -log(B_USA)./yf;
    F0_sim_USA = data_USA.spot .* exp(rates .* yf);
    hold on;  plot(dates_USA, F0_sim_USA, '-o', 'Color', 'r');
  
end

%% POINT 6: Calibration

%% Options selection

% Choice of only OTM options for the further calibration
data_EU_OTM = OTM_preprocessing(data_EU, F0_EU);
data_USA_OTM = OTM_preprocessing(data_USA, F0_USA);

% Computing the Delta of Black & Scholes over the OTM Call/Put in order to
% clean dataset from too far from the ATM point prices

data_calib_EU = dataset_preprocessing(data_EU_OTM, F0_EU, B_EU, date_settlement, flag);
data_calib_USA = dataset_preprocessing(data_USA_OTM, F0_USA, B_USA, date_settlement, flag);

%% Plot of the surface of the implied volatilities

if flag == 1
    surface_vols(data_calib_EU,F0_EU);
    surface_vols(data_calib_USA,F0_USA);
end

%% Calibration of the model parameters

% Quantities of interest
% x0 = [10 2 0.5 10 2 0.5];
% x0 = 0.09 * [1 1 0.5 1 1 0.5];
% x0 = [32 0.04 0.36 11.8 0.09 0.37];
% x0 = 0.1 * ones(1, 6);
x0 = [0.3 -0.5 0.15 0.3 -0.5 0.15];

initial_cond = x0;

% Linear inequality constraints
A = [-1 0 0 0 0 0;
    0 0 -1 0 0 0;
    0 0 0 -1 0 0;
    0 0 0 0 0 -1];

b = [0; 0; 0; 0];

% Linear equality constraints
Aeq = []; beq = [];

% Lower and upper bounds
lb = [0.01; -Inf; 0; 0.01; -Inf; 0];
ub = [];

% Options for the visualization
options = optimset('MaxFunEvals', 3e3, 'Display', 'iter');

% Calibration
params_marginals = fmincon(@(params) new_calibration(params, data_calib_EU, data_calib_USA, ...
    F0_EU, B_EU, F0_USA, B_USA, date_settlement), x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr(params), options);

% Display of the parameters on the console
params_USA = params_marginals(1:3);
params_EU = params_marginals(4:6);

% Explicit useful params
k1 = params_USA(1); k2 = params_EU(1);

%% parameters group 2B

% % sigma_EU = 0.11991
% % kappa_EU = 0.0024632
% % theta_EU = 0.021422
% % sigma_US = 0.10851
% % kappa_US = 0.0019843
% % theta_US = 0.021599

% params_USA = [0.0019843 , 0.021599 , 0.10851];
% params_EU = [0.0024632 , 0.021422 , 0.11991];


%% Plots of the prices with calibrated values

if flag == 1

    % Plot over the entire curve
    plot_calls_puts_total(data_EU, F0_EU, B_EU, params_EU, date_settlement);
    plot_calls_puts_total(data_USA, F0_USA, B_USA, params_USA, date_settlement);

    % Plot over the filtered options
    plot_calls_puts(data_calib_EU, F0_EU, B_EU, params_EU, date_settlement);
    plot_calls_puts(data_calib_USA, F0_USA, B_USA, params_USA, date_settlement);

    % Plot the implied volatilities over the Calls
    plot_volatility_smiles(data_calib_EU, F0_EU, B_EU, params_EU, date_settlement)
    plot_volatility_smiles(data_calib_USA, F0_USA, B_USA, params_USA, date_settlement)

end

%%

%% POINT 7: comparison with the historical correlation

%% Calibration over the rho - OLD

% Initialization of the parameters
A = []; b = []; Aeq = []; beq = [];
lb = 0; ub = [];
x0 = 1;

% % Calibration of the nuZ parameter
params = fmincon(@(params) abs(sqrt( k1 * k2)/ params - rho_historical), ...
    x0, A, b, Aeq, beq, lb, ub, [], options);

nu_z = params;

nu_1 = k1*nu_z/(nu_z-k1);
nu_2 = k2*nu_z/(nu_z-k2);

%% Calibration over the rho - MISMATCH CORR

% Initialization of the parameters
A = [-1 0 0; 0 -1 0; 0 0 -1]; b = [0; 0; 0];
Aeq = []; beq = [];
lb = zeros(1, 3); ub = [];

x0 = ones(1, 3);

% Calibration of the nu parameters
% params = fmincon(@(params) abs(sqrt(params(1) * params(2) / ((params(1) + params(3))*(params(2) + params(3)))) - rho_historical), ...
%    x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr_corr(params, k1, k2), options);

params = fmincon(@(params) (sqrt(params(1) * params(2) / ((params(1) + params(3))*(params(2) + params(3)))) - rho_historical)^2, ...
    x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr_corr(params, k1, k2), options);

nu_1 = params(1);
nu_2 = params(2);
nu_z = params(3);

%% disp the calibrated parameters

disp_params(params_marginals, nu_1 ,nu_2 ,nu_z, initial_cond, save_results);

%% Common and idiosynchratic parameters

% kappa_1 = params_USA(1); theta_1 = params_USA(2); sigma_1 = params_USA(3);
% kappa_2 = params_EU(1); theta_2 = params_EU(2); sigma_2 = params_EU(3);

% Compute the rho obtained from the model
rho_model_Levy = sqrt(params(1) * params(2) / ((params(1) + params(3))*(params(2) + params(3))));

% Compute the solution of the system
sol =  marginal_param(params_USA,params_EU , nu_z , rho_model_Levy);

% Explicitate the parameters
a_USA = sol.x(1); a_EU = sol.x(2);
beta_z = sol.x(3); gamma_z = sol.x(4);

% Remaining parameters computation
beta_USA = params_USA(2) - a_USA * beta_z;                       % From condition (10.1)
gamma_USA = sqrt((params_USA(3)^2) - (a_USA^2) * (gamma_z^2));   % From condition (10.2)
% nu_USA = (nu_z * params_USA(1))/(nu_z - params_USA(1)); 
nu_USA = nu_1;

beta_EU = params_EU(2) - a_EU * beta_z;
gamma_EU = sqrt(params_EU(3)^2 - a_EU ^ 2 * gamma_z^2);
% nu_EU = (nu_z * params_EU(1))/(nu_z - params_EU(1));
nu_EU = nu_2;

% Creation of zipped vectors

% Idiosyncratic parameters USA
% sol_USA = [a_USA, beta_USA, gamma_USA];
idiosync_USA = [nu_USA, beta_USA, gamma_USA, a_USA];

% Idiosyncratic parameters EU
% sol_EU = [a_EU , beta_EU, gamma_EU];
idiosync_EU = [nu_EU, beta_EU, gamma_EU, a_EU];

% Systematic parameters
syst_Z = [nu_z , beta_z, gamma_z];

%% Display of the parameters over the command window

% disp_marginal_params(idiosync_USA , idiosync_EU , beta_z, gamma_z,save_results);

%%

%% Point 8: Black model calibration

%% Calibration of the model parameters

% Initialization of the parameters
x0 = 0.5;
x0 = 1e-4;
% Linear inequality constraints
A = []; b = [];
% Linear equality constraints
Aeq = []; beq = [];
% Lower and upper bounds
lb = 0; ub = [];

% Options for the visualization
options = optimset('MaxFunEvals', 3e3, 'Display', 'iter');

% Calibration of sigma EU
sigma_EU = fmincon(@(sigma) blk_calibration(sigma, data_calib_EU, F0_EU, B_EU, date_settlement), ...
    x0, A, b, Aeq, beq, lb, ub, [], options);

% Calibration of sigma USA
sigma_USA = fmincon(@(sigma) blk_calibration(sigma, data_calib_USA, F0_USA, B_USA, date_settlement), ...
    x0, A, b, Aeq, beq, lb, ub, [], options);

%% Plots of the prices with calibrated values

if flag == 1
    % Plot over the entire curve
    blk_plot_calls_puts_total(data_EU, F0_EU, B_EU, sigma_EU, date_settlement);
    blk_plot_calls_puts_total(data_USA, F0_USA, B_USA, sigma_USA, date_settlement);

    % Plot over the filtered options
    blk_plot_calls_puts(data_calib_EU, F0_EU, B_EU, sigma_EU, date_settlement);
    blk_plot_calls_puts(data_calib_USA, F0_USA, B_USA, sigma_USA, date_settlement);

    % Plot the implied volatilities over the Calls
    blk_plot_volatility_smiles(data_calib_EU, F0_EU, B_EU, sigma_EU, date_settlement)
    blk_plot_volatility_smiles(data_calib_USA, F0_USA, B_USA, sigma_USA, date_settlement)

end

%%

%% Point 9: Pricing of the certificate - LEVY

year_to_maturity = 1;

[rate_USA, TTM] = interp_pricing_params(datenum(data_calib_USA.datesExpiry), B_USA, date_settlement, year_to_maturity);
[rate_EU, ~] = interp_pricing_params(datenum(data_calib_EU.datesExpiry), B_EU, date_settlement, year_to_maturity);

yf = yearfrac(date_settlement, data_USA.datesExpiry, conv_ACT365);
interp_date = datenum(busdate(datetime(date_settlement, 'ConvertFrom', 'datenum') - caldays(1) + calyears(year_to_maturity), 1, eurCalendar));

if flag
    figure;
    plot(dates_USA, -log(B_USA)./yf, '-*', 'Color', 'b'); hold on;
    plot(interp_date, rate_USA, 'o', 'Color', 'r'); grid on;
    title('Rates USA');
    ylabel('Rates');
    datetick('x','dd-mmm-yyyy','keepticks');
end

    
S0_USA = data_USA.spot; S0_EU = data_EU.spot;

% % Simulation of theunderlying stock prices
% St_USA = stock_simulation_Levy(params_USA, rate_USA , TTM , S0_USA);
% St_EU = stock_simulation_Levy(params_EU, rate_EU , TTM , S0_EU);

% Computation of the discount at 1y
B0_Levy = exp(-rate_USA * TTM);

S0_Levy = [S0_USA S0_EU];
rates_Levy = [rate_USA rate_EU];

St_Levy = stock_simulation_Levy(idiosync_USA, idiosync_EU, syst_Z, params_USA, params_EU, S0_Levy, rates_Levy, TTM);

% Unpacking the results
St_USA_Levy = St_Levy(:, 1);
St_EU_Levy = St_Levy(:, 2);


% Computation of the pricing certificate payoff 
indicator_Levy = St_EU_Levy < (0.95 * S0_EU);
certificate_payoff_Levy = max(St_USA_Levy - S0_USA, 0) .* indicator_Levy;

% Mean price and confidence interval
[mean_price_Levy, ~, IC_Levy] = normfit(B0_Levy * certificate_payoff_Levy);

%% Point 9: Pricing of the certificate - Brownian Motion

year_to_maturity = 1;

% Computation of the rates and the time to maturity
[rate_USA, TTM] = interp_pricing_params(datenum(data_calib_USA.datesExpiry), B_USA, date_settlement, year_to_maturity);
[rate_EU, ~] = interp_pricing_params(datenum(data_calib_EU.datesExpiry), B_EU, date_settlement, year_to_maturity);

% Stock prices
S0_USA = data_USA.spot; S0_EU = data_EU.spot;

% Forward prices
F01_USA = S0_USA*exp(rate_USA*TTM);
F01_EU = S0_EU*exp(rate_EU*TTM);

% % Simulation of the underlying stock prices
[St_Black, St_Black_AV] = stock_simulation_Black([sigma_USA; sigma_EU], [F01_USA; F01_EU], ...
    [rate_USA; rate_EU], rho_historical, TTM);

% Computation of the discount at 1y
B0_black = exp(-rate_USA * TTM);

%% NORMAL

% Unpacking the results
St_USA_Black = St_Black(:, 1); St_EU_Black = St_Black(:, 2);

% Computation of the pricing certificate payoff 
indicator_Black = St_EU_Black < (0.95 * S0_EU);
certificate_payoff_Black = max(St_USA_Black - S0_USA, 0) .* indicator_Black;

% Mean price and confidence interval
[mean_price_Black, ~, IC_Black] = normfit(B0_black * certificate_payoff_Black);

%% ANTITHETIC

% Unpacking the results
St_USA_Black_AV = St_Black_AV(:, 1); St_EU_Black_AV = St_Black_AV(:, 2);

% Computation of the pricing certificate payoff 
indicator_Black_AV = St_EU_Black_AV < (0.95 * S0_EU);
certificate_payoff_Black_AV = max(St_USA_Black_AV - S0_USA, 0) .* indicator_Black_AV;

certificate_payoff_Black_AV = (certificate_payoff_Black_AV + certificate_payoff_Black)/2;

% Mean price and confidence interval
[mean_price_Black_AV, ~, IC_Black_AV] = normfit(B0_black * certificate_payoff_Black_AV);

%% Closed formula

% Price computed via the closed formula
price_semiclosed = blk_semiclosed(data_USA.spot, rate_USA, rate_EU, sigma_USA, sigma_EU, rho_historical, TTM);

%% Display of the prices:

% mean_price_Levy = 0;
% IC_Levy = [0; 0];
disp_contract_prices(mean_price_Levy,IC_Levy,mean_price_Black,IC_Black,mean_price_Black_AV,IC_Black_AV,price_semiclosed)
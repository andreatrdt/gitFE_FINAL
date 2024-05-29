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
rng(42);

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

if flag == 1
    plot_returns(SP500_EUR500,date_settlement)
end

%% POINT 5: Forward Prices

[F0_EU, B_bar_EU] = forward_prices(data_EU, flag);
[F0_USA, B_bar_USA] = forward_prices(data_USA, flag);

B_EU = B_bar_EU;
B_USA = B_bar_USA;

if flag == 1
    figure;
    plot(dates_EU,B_EU,'-*','Color','b');
    title('Discount factor for the EU market');
    datetick('x','dd-mmm-yyyy','keepticks')


    figure;
    plot(dates_USA,B_USA,'-*','Color','r');
    title('Discount factor for the USA market');
    datetick('x','dd-mmm-yyyy','keepticks')
end



%% POINT 6: Calibration

%% Options selection

% Choice of only OTM options for the further calibration
data_EU_OTM = OTM_preprocessing(data_EU, B_EU, F0_EU);
data_USA_OTM = OTM_preprocessing(data_USA, B_USA, F0_USA);

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
% x0 = [1 1 0.5 1 1 0.5];
% x0 = [32 0.04 0.36 11.8 0.09 0.37];
% x0 = ones(1, 6);
x0 = [0.1 0.02 0.1 0.1 0.02 0.1];

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
lb = [0; -Inf; 0; 0; -Inf; 0];
ub = [];

% Options for the visualization
options = optimset('MaxFunEvals', 50e3, 'Display', 'iter');

% Calibration
params_marginals = fmincon(@(params) new_calibration(params, data_calib_EU, data_calib_USA, ...
    F0_EU, B_EU, F0_USA, B_USA, date_settlement), x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr(params), options);

% Display of the parameters on the console
params_USA = params_marginals(1:3);
params_EU = params_marginals(4:6);


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

%% 2nd Calibration over the rho

% Computation of the historical correlation
rho = hist_corr(SP500_EUR500);

% % Initialization of the parameters
% A = []; b = []; Aeq = []; beq = [];
% lb = zeros(1, 3); ub = [];
% x0 = ones(1, 3);

% Recall the parameters
k1 = params_USA(1); k2 = params_EU(1);

% % Calibration of the nuZ parameter
% params = fmincon(@(params) (sqrt(params(1) * params(2) / ((params(1) + params(3))*(params(2) + params(3)))) - rho)^2, ...
%     x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr_corr(params, k1, k2), options);

% nu_1 = params(1);
% nu_2 = params(2);
% nu_z = params(3);

nu_z = sqrt(k1*k2)/rho;
nu_1 = k1*nu_z/(nu_z-k1);
nu_2 = k2*nu_z/(nu_z-k2);

%% disp the calibrated parameters

disp('Calibrated parameters for the USA market: ')
disp(params_USA)
disp('Calibrated parameters for the EU market: ')
disp(params_EU)


disp('calibrated nu_1: ')
disp(nu_1)
disp('calibrated nu_2: ')
disp(nu_2)
disp('calibrated nu_z: ')
disp(nu_z)

if save_results == 1

    % save results in a txt file

    fileID = fopen('results.txt','w');

    fprintf(fileID,'-----------------------');
    fprintf(fileID,'X0 used :\n');
    fprintf(fileID,'%f \n',initial_cond);
    fprintf(fileID,'-----------------------');
    fprintf(fileID,'Calibrated parameters for the USA market: \n');
    fprintf(fileID,'%f \n',params_USA);
    fprintf(fileID,'-----------------------');
    fprintf(fileID,'Calibrated parameters for the EU market: \n');
    fprintf(fileID,'%f \n',params_EU);
    fprintf(fileID,'-----------------------');
    fprintf(fileID,'calibrated nu_1: \n');
    fprintf(fileID,'%f \n',nu_1);
    fprintf(fileID,'-----------------------');
    fprintf(fileID,'calibrated nu_2: \n');
    fprintf(fileID,'%f \n',nu_2);
    fprintf(fileID,'-----------------------');
    fprintf(fileID,'calibrated nu_z: \n');
    fprintf(fileID,'%f \n',nu_z);
    fprintf(fileID,'-----------------------');


    fclose(fileID);

end


%% Common and idiosynchratic parameters


sol_USA =  marginal_param(params_USA,nu_z)

a_USA = sol_USA.x(1)
Beta_USA = sol_USA.x(4)
gamma_USA = sol_USA.x(5)


sol_EU = marginal_param(params_EU,nu_z)

a_EU = sol_EU.x(1)
Beta_z = sol_EU.x(2)
gamma_z = sol_EU.x(3)
Beta_EU = sol_EU.x(4)
gamma_EU = sol_EU.x(5)


%% Point 8: Black model calibration

%% Calibration of the model parameters

% Initialization of the parameters
x0 = 0.5;
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

    %Plot over the filtered options (TO BE IMPLEMENTED)
    plot_calls_puts(data_calib_EU, F0_EU, B_EU, params_EU, date_settlement);
    plot_calls_puts(data_calib_USA, F0_USA, B_USA, params_USA, date_settlement);

    %Plot the implied volatilities over the Calls (TO BE IMPLEMENTED)
    blk_plot_volatility_smiles(data_calib_EU, F0_EU, B_EU, params_EU, date_settlement)
    blk_plot_volatility_smiles(data_calib_USA, F0_USA, B_USA, params_USA, date_settlement)

end

%% Point 9: Pricing of the certificate - Levy

[rate_USA, interp_F0_USA] = interp_pricing_params(datenum(data_calib_USA.datesExpiry), F0_USA, B_USA, date_settlement);
[rate_EU, interp_F0_EU] = interp_pricing_params(datenum(data_calib_EU.datesExpiry), F0_EU, B_EU, date_settlement);

[St_Levy, S0_Levy] = stock_simulation_Levy(sol_USA, sol_EU, nu_1 , nu_2 , nu_z , [interp_F0_USA; interp_F0_EU] ...
                , [B_USA ; B_EU] , [rate_USA; rate_EU] , date_settlement);

% Unpacking the results
St_USA_Levy = St_Levy(:,1);
St_EU_Levy = St_Levy(:,2);

S0_USA_Levy = S0_Levy(1);
S0_EU_Levy = S0_Levy(2);

% Computation of the pricing certificate payoff


indx_Levy = St_EU_Levy > (0.95 * S0_EU_Levy);

certificate_payoff_Levy = max(St_USA_Levy - S0_USA_Levy, 0) .* indx_Levy;

% Mean price and confidence interval
[mean_price_levy, ~, IC_Levy] = normfit(certificate_payoff_Levy);

% Plot of the histogram of positive payoffs
certificate_reduced_Levy = certificate_payoff_Levy(find(certificate_payoff_Levy));
%histogram(certificate_reduced_Levy);

[mean, ~, IC] = normfit(certificate_reduced_Levy)

%% Point 9: Pricing of the certificate - Brownian Motion

% Computation of the rates and initial forwards
[rate_USA, interp_F0_USA] = interp_pricing_params(datenum(data_calib_USA.datesExpiry), F0_USA, B_USA, date_settlement);
[rate_EU, interp_F0_EU] = interp_pricing_params(datenum(data_calib_EU.datesExpiry), F0_EU, B_EU, date_settlement);

% Simulation of theunderlying stock prices
[St_Black, S0_Black] = stock_simulation_Black([sigma_USA; sigma_EU], [interp_F0_USA; interp_F0_EU], ...
    [rate_USA; rate_EU], rho, date_settlement);

% Unpacking the results
St_USA_Black = St_Black(:, 1); St_EU_Black = St_Black(:, 2);
S0_USA_Black = S0_Black(1); S0_EU_Black = S0_Black(2);

% Computation of the pricing certificate payoff
indicator_Black = St_EU_Black < (0.95 * S0_EU_Black);
certificate_payoff_Black = max(St_USA_Black - S0_USA_Black, 0) .* indicator_Black;

% Mean price and confidence interval
[mean_price_Black, ~, IC_Black] = normfit(certificate_payoff_Black);

% Plot of the histogram of positive payoffs
idx = find(certificate_payoff_Black > 0);
certificate_reduced_Black = certificate_payoff_Black(idx);
histogram(certificate_reduced_Black);

[mean, ~, IC] = normfit(certificate_reduced_Black)
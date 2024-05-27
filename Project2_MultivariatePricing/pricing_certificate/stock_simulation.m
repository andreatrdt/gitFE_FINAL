function [prices, S0] = stock_simulation(dataset, params, F0, B0, date_settlement)
% Pricing of the underlying process Si(t)
% 
% INPUT:
% dataset:              initial dataset
% params:               calibration NIG parameters [k, theta, sigma]
% F0:                   initial forward value
% B0:                   initial discount value
% date_settlement:      initial date of the certificate
% 
% OUTPUT:
% prices:               underlying stock to be simulated
% 
% USES:
% function rate_interpolation()

    %% Conventions

    conv_ACT365 = 3;
    
    %% Unpacking of the parameters

    k = params(1); 
    theta = params(2);
    sigma = params(3);

    %% Computation of the support params

    nSim = 1e6;

    drift_compensator = - 1/k * (1 - sqrt(1 - 2*k*theta - k*sigma^2));

    % Dates for each expiry
    dates = datenum(dataset.datesExpiry);
    
    % Computation of the TTM
    interp_date = datenum(busdate(datetime(dates(1), 'ConvertFrom', 'datenum') - caldays(1) + calyears(1), 1, eurCalendar));
    TTM = yearfrac(date_settlement, interp_date, conv_ACT365);

    %% Computation of interest rate

    rate = rate_interpolation(dates, B0, date_settlement);

    %% Computation of the interpolated F(0, 1y)

    interp_F0 = interp1(dates, F0, interp_date);

    %% Computation of the initial stock price S(0)

    S0 = interp_F0 / exp(rate * TTM);
    
    %% Simulation of the NIG process

    % Stochastic parts

    g = randn(nSim, 1);
    G = random('InverseGaussian', 1, TTM/k, [nSim, 1]);
    
    % Creation of Xt dynamic

    Xt =theta.*G + sigma .* sqrt(TTM .* G) .* g;

    %% Computation of the initial stock

    prices = S0 .* exp(rate - drift_compensator * TTM + Xt);

end % function stock_simulation
function [prices , S0] = stock_simulation_Levy_prova( params, rate , TTM , S0)
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

    %% conventions

    conv_ACT365 = 3;

    %% unpack parameters

    k = params(1);

    theta = params(2);

    sigma = params(3);

    % Computation of the support params

    nSim = 1e6;

    drift_compensator = - 1/k * (1 - sqrt(1 - 2*k*theta - k*sigma^2));
        
    %% Simulation of the NIG process

    % Stochastic parts

    g = randn(nSim, 1);
    G = random('InverseGaussian', 1, TTM/k, [nSim, 1]);
    
    % Creation of Xt dynamic

    Xt =theta.*G + sigma .* sqrt(TTM .* G) .* g;

    %% Computation of the initial stock

    prices = S0 .* exp(rate - drift_compensator * TTM + Xt);

end % function stock_simulation
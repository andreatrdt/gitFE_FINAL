function [prices, pricesAV] = stock_simulation_Black(sigmas, F0, rho, TTM)
% Pricing of the underlying process Si(t)
% 
% INPUT:
% sigmas:               [VECTOR] volatilities of the Black process
% S0:                   [VECTOR] initial spot value
% rates:                [VECTOR] interest rates
% rho:                  [SCALAR] historical correlation
% TTM:                  [SCALAR] time to maturity
% 
% OUTPUT:
% prices:               [MATRIX] underlying stock to be simulated
% S0:                   [VECTOR] initial value of the stock
%
% USES:      none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Initialization

    nSim = 1e7;

    %% Simulation of the NIG process

    % Stochastic part
    covarianceMatrix = [1 rho; rho 1];
    meanVector = [0; 0];

    g = mvnrnd(meanVector, covarianceMatrix, nSim);
    
    % Creation of Xt dynamic

    Xt = -0.5 .* sigmas'.^2 .* TTM + sigmas' .* sqrt(TTM) .* g;
    % Xt = -0.5 .* sigmas'.^2 .* TTM + sigmas' .* A;

    %% Computation of the initial stock

    prices = F0' .* exp(Xt);
    
    %% Computation for the antithetic behaviour

    Xt_AV = [-0.5 .* sigmas'.^2 .* TTM + sigmas' .* sqrt(TTM) .* g(1:nSim/2, :); ...
        -0.5 .* sigmas'.^2 .* TTM - sigmas' .* sqrt(TTM) .* g(nSim/2+1:end, :)] ;

    pricesAV = F0' .* exp(Xt_AV);

end % function stock_simulation_Black
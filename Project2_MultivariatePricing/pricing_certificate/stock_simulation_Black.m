function [prices, pricesAV] = stock_simulation_Black(sigmas, F0, rates, rho, TTM)
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


    %% Initialization

    nSim = 5 * 1e7;

    %% Simulation of the NIG process

    % Stochastic part
    covarianceMatrix = [TTM rho*TTM; rho*TTM TTM];
    cholesky_matrix = chol(covarianceMatrix)';

    % A = (cholesky_matrix * randn(2, nSim))';

    meanVector = [0; 0];

    rng(2);

    g = mvnrnd(meanVector, covarianceMatrix, nSim);
    
    % Creation of Xt dynamic

    Xt = -0.5 .* sigmas'.^2 .* TTM + sigmas' .* sqrt(TTM) .* g;
    % Xt = -0.5 .* sigmas'.^2 .* TTM + sigmas' .* A;

    %% Computation of the initial stock

    % prices = S0' .* exp(rates'*TTM + Xt);

    prices = F0' .* exp(Xt);
    
    %% Computation for the antithetic behaviour

    Xt_AV = -0.5 .* sigmas'.^2 .* TTM - sigmas' .* sqrt(TTM) .* g;
    % Xt_AV = -0.5 .* sigmas'.^2 .* TTM - sigmas' .* A;
    % pricesAV = S0' .* exp(rates'*TTM + Xt_AV);
    pricesAV = F0' .* exp(Xt_AV);
    
    % pricesAV = (prices + pricesAV)/2;

end % function stock_simulation_Black
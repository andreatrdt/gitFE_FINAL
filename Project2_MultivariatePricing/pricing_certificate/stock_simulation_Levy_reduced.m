function [prices , S0] = stock_simulation_Levy_reduced(params_USA,params_EU, rates , TTM , S0, rho)
% 
% INPUT:
% params:               calibration NIG parameters [k, theta, sigma]
% rate                  rate
% TTM                   time to maturity
% SO                    initial value of the stocks
% 
% OUTPUT:
% prices:               underlying stock to be simulated
% 
% USES:     none

% Authors:
% M.Maspes, A.Tarditi, M.Torba

    %% unpack parameters

    k_1 = params_USA(1);

    theta_1 = params_USA(2);

    sigma_1 = params_USA(3);

    k_2 = params_EU(1);

    theta_2 = params_EU(2);

    sigma_2 = params_EU(3);

    % Computation of the support params

    nSim = 1e7;

    drift_compensator = - [1/k_1 * (1 - sqrt(1 - 2*k_1*theta_1 - k_1*sigma_1^2)) 1/k_2 * (1 - sqrt(1 - 2*k_2*theta_2 - k_2*sigma_2^2))];
        
    %% Simulation of the NIG process

    % Stochastic part
    covarianceMatrix = [1 rho; rho 1];
    meanVector = [0; 0];
    
    g = mvnrnd(meanVector, covarianceMatrix, nSim);
        
    G_1 = random('InverseGaussian', 1, TTM/k_1, [nSim, 1]);
    G_2 = random('InverseGaussian', 1, TTM/k_2, [nSim, 1]);

    G=[G_1 G_2];

    % Creation of Xt dynamic

    X_1 =-(1/2 + theta_1).*G_1.*sigma_1^2.*TTM+ sigma_1 .* sqrt(TTM .* G_1) .* g(:,1);

    X_2 =-(1/2 + theta_2).*G_2.*sigma_2^2.*TTM+ sigma_2 .* sqrt(TTM .* G_2) .* g(:,2);

    Xt = [X_1 X_2];
    
    %% Computation of the initial stock

    prices = real(S0 .* exp((rates + drift_compensator) * TTM + Xt));

end % function stock_simulation_Levy_reduced
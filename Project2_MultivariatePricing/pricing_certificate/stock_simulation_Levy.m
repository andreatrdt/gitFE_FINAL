function [prices , S0] = stock_simulation_Levy(sol_USA, sol_EU, nu_1 , nu_2 , nu_z , params_USA , params_EU , S0 , rates , TTM)
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
;
    
    %% Unpacking of the parameters


    kappa_USA = params_USA(1);
    theta_USA = params_USA(2);
    sigma_USA = params_USA(3);


    kappa_EU = params_EU(1);
    theta_EU = params_EU(2);
    sigma_EU = params_EU(3);

  
    a_EU = sol_EU.x(1);
    a_USA = sol_USA.x(1);

    Beta_z = sol_EU.x(2);
    Beta_USA = sol_USA.x(4);
    Beta_EU = sol_EU.x(4);

    gamma_z = sol_EU.x(3);
    gamma_USA = sol_USA.x(5);  
    gamma_EU = sol_EU.x(5);


    %% Computation of the support params

    nSim = 1e6;

    drift_compensator_USA = - 1/kappa_USA * (1 - sqrt(1 - 2*kappa_USA*theta_USA - kappa_USA*sigma_USA^2));

    drift_compensator_EU = - 1/kappa_EU * (1 - sqrt(1 - 2*kappa_EU*theta_EU - kappa_EU*sigma_EU^2));

    drift_compensator = [drift_compensator_USA, drift_compensator_EU];
    
    %% Simulation of the NIG process

    % Stochastic parts

    g = randn(nSim, 1);
    G_1 = random('InverseGaussian', 1, TTM/nu_1, [nSim, 1]);
    G_2 = random('InverseGaussian', 1, TTM/nu_2, [nSim, 1]);
    G_z = random('InverseGaussian', 1, TTM/nu_z, [nSim, 1]);


    Y_1 = Beta_USA.*G_1 + gamma_USA .* sqrt(TTM .* G_1) .* g;
    Y_2 = Beta_EU.*G_2 + gamma_EU .* sqrt(TTM .* G_2) .* g;
    Z = Beta_z.*G_z + gamma_z .* sqrt(TTM .* G_z) .* g;

    %% Computation of the initial stock

    X_1 = Y_1 + a_USA * Z;
    X_2 = Y_2 + a_EU * Z;

    Xt = [X_1 , X_2];

    prices = S0' .* exp(rates' - drift_compensator * TTM + Xt);

    prices = prices';

end % function stock_simulation
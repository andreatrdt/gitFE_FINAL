function [stock , S0] = stock_simulation_Levy(idiosync_USA, idiosync_EU, syst_Z, params_USA, params_EU, S0, rates, TTM)
% Pricing of the underlying process Si(t)
% 
% INPUT:
% idiosync_USA:         idiosyncratic parameters for the USA
% idiosync_EU:          idiosyncratic parameters for the EU
% syst_Z:               systematic parameters for the Z process
% params_USA:           parameters for the USA
% params_EU:            parameters for the EU
% S0:                   initial stock value
% rates:                interest rates
% TTM:                  time to maturity
% 
% OUTPUT:
% stock:                underlying stock to be simulated
    
    %% Unpacking of the parameters

    kappa_USA = params_USA(1);
    theta_USA = params_USA(2);
    sigma_USA = params_USA(3);

    kappa_EU = params_EU(1);
    theta_EU = params_EU(2);
    sigma_EU = params_EU(3);

    nu_z = syst_Z(1);
    nu_USA = idiosync_USA(1);
    nu_EU = idiosync_EU(1);
    
    Beta_z = syst_Z(2);
    Beta_USA = idiosync_USA(2);
    Beta_EU = idiosync_EU(2);

    gamma_z = syst_Z(3);
    gamma_USA = idiosync_USA(3);  
    gamma_EU = idiosync_EU(3);

    a_USA = idiosync_USA(4);
    a_EU = idiosync_EU(4);

    %% Computation of the support params

    nSim = 1e6;

    drift_compensator_USA = - 1/kappa_USA * (1 - sqrt(1 - 2*kappa_USA*theta_USA - kappa_USA*sigma_USA^2));
    drift_compensator_EU = - 1/kappa_EU * (1 - sqrt(1 - 2*kappa_EU*theta_EU - kappa_EU*sigma_EU^2));
    drift_compensator = [drift_compensator_USA drift_compensator_EU];
    
    %% Simulation of the NIG process

    % Stochastic parts
    g_1 = randn(nSim, 1);
    g_2 = randn(nSim, 1);
    g_z = randn(nSim, 1);

    G_1 = random('InverseGaussian', 1, TTM/nu_USA, [nSim, 1]);
    G_2 = random('InverseGaussian', 1, TTM/nu_EU, [nSim, 1]);
    G_z = random('InverseGaussian', 1, TTM/nu_z, [nSim, 1]);

    Y_1 = -Beta_USA .* gamma_USA^2 .*G_1 * TTM + gamma_USA .* sqrt(TTM .* G_1) .* g_1;
    Y_2 = -Beta_EU .* gamma_EU^2.* G_2 * TTM + gamma_EU .* sqrt(TTM .* G_2) .* g_2;
    Z = -Beta_z .* gamma_z^2.* G_z * TTM + gamma_z .* sqrt(TTM .* G_z) .* g_z;

    %% Conjunction of the processes

    % Marginal processes
    X_1 = Y_1 + a_USA * Z;
    X_2 = Y_2 + a_EU * Z;

    % General vector
    Xt = [X_1 X_2];

    stock = S0 .* exp((rates + drift_compensator) * TTM + Xt);

end % function stock_simulation
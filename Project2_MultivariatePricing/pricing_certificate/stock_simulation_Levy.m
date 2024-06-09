function [stock , stock_AV] = stock_simulation_Levy(idiosync_USA, idiosync_EU, syst_Z, params_USA, params_EU, S0, rates, TTM)
% Pricing of the underlying process Si(t)
% 
% INPUT:
% idiosync_USA:         [VECTOR]idiosyncratic parameters for the USA
% idiosync_EU:          [VECTOR]idiosyncratic parameters for the EU
% syst_Z:               [VECTOR]systematic parameters for the Z process
% params_USA:           [VECTOR]parameters for the USA
% params_EU:            [VECTOR]parameters for the EU
% S0:                   [SCALAR]initial stock value
% rates:                [VECTOR]interest rates
% TTM:                  [SCALAR]time to maturity
% 
% OUTPUT:
% stock:                underlying stock to be simulated
% stock_AV:             antithetic underlying
%
% USES:            none

% Authors:
% M.Maspes, A.Tarditi, M.Torba

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

    nSim = 1e7;

    drift_compensator_USA = - 1/kappa_USA * (1 - sqrt(1 - 2*kappa_USA*theta_USA - kappa_USA*sigma_USA^2));
    drift_compensator_EU = - 1/kappa_EU * (1 - sqrt(1 - 2*kappa_EU*theta_EU - kappa_EU*sigma_EU^2));
    drift_compensator = [drift_compensator_USA drift_compensator_EU];

    % drift_compensator_Y_USA = - 1/nu_USA * (1 - sqrt(1 - 2*nu_USA*Beta_USA - nu_USA*gamma_USA^2));
    % drift_compensator_Y_EU = - 1/nu_EU * (1 - sqrt(1 - 2*nu_EU*Beta_EU - nu_EU*gamma_EU^2));
    % drift_compensator_Z = - 1/nu_z * (1 - sqrt(1 - 2*nu_z*Beta_z - nu_z*gamma_z^2));
    
    
    %% Simulation of the NIG process

    % Stochastic parts
    g_1 = randn(nSim, 1);
    g_2 = randn(nSim, 1);
    g_z = randn(nSim, 1);

    G_1 = random('InverseGaussian', 1, TTM/nu_USA, [nSim, 1]);
    G_2 = random('InverseGaussian', 1, TTM/nu_EU, [nSim, 1]);
    G_z = random('InverseGaussian', 1, TTM/nu_z, [nSim, 1]);

    % Y_1 = -(0.5+Beta_USA) * gamma_USA^2 .*G_1 * TTM + gamma_USA .* sqrt(TTM .* G_1) .* g_1;
    % Y_2 = -(0.5+Beta_EU) * gamma_EU^2 .* G_2 * TTM + gamma_EU .* sqrt(TTM .* G_2) .* g_2;
    % Z = -(0.5+Beta_z) * gamma_z^2 .* G_z * TTM + gamma_z .* sqrt(TTM .* G_z) .* g_z;

    Y_1 = Beta_USA .*G_1 * TTM + gamma_USA .* sqrt(TTM .* G_1) .* g_1;
    Y_2 = Beta_EU .* G_2 * TTM + gamma_EU .* sqrt(TTM .* G_2) .* g_2;
    Z = Beta_z .* G_z * TTM + gamma_z .* sqrt(TTM .* G_z) .* g_z;

    %% Conjunction of the processes

    % Marginal processes
    X_1 = Y_1 + a_USA * Z;
    X_2 = Y_2 + a_EU * Z;

    % General vector
    Xt = [X_1 X_2];

    stock = real(S0 .* exp((rates + drift_compensator) * TTM + Xt));

    %% ANTITHETIC FRAMEWORK

    % Y_1_AV = [Y_1(1:nSim/2, :); ...
    %     -(0.5+Beta_USA) * gamma_USA^2 .*G_1(nSim/2 + 1:end, :) * TTM - gamma_USA .* sqrt(TTM .* G_1(nSim/2 + 1:end, :)) .* g_1(nSim/2 + 1:end, :)];
    % Y_2_AV = [Y_2(1:nSim/2, :); ...
    %     -(0.5+Beta_EU) * gamma_EU^2 .* G_2(nSim/2 + 1:end, :) * TTM - gamma_EU .* sqrt(TTM .* G_2(nSim/2 + 1:end, :)) .* g_2(nSim/2 + 1:end, :)];
    % Z_AV = [Z(1:nSim/2, :); ...
    %     -(0.5+Beta_z) * gamma_z^2 .* G_z(nSim/2 + 1:end, :) * TTM - gamma_z .* sqrt(TTM .* G_z(nSim/2 + 1:end, :)) .* g_z(nSim/2 + 1:end, :)];


    Y_1_AV = [Y_1(1:nSim/2, :); ...
        Beta_USA .*G_1(nSim/2 + 1:end, :) * TTM - gamma_USA .* sqrt(TTM .* G_1(nSim/2 + 1:end, :)) .* g_1(nSim/2 + 1:end, :)];
    Y_2_AV = [Y_2(1:nSim/2, :); ...
        Beta_EU .* G_2(nSim/2 + 1:end, :) * TTM - gamma_EU .* sqrt(TTM .* G_2(nSim/2 + 1:end, :)) .* g_2(nSim/2 + 1:end, :)];
    Z_AV = [Z(1:nSim/2, :); ...
        Beta_z .* G_z(nSim/2 + 1:end, :) * TTM - gamma_z .* sqrt(TTM .* G_z(nSim/2 + 1:end, :)) .* g_z(nSim/2 + 1:end, :)];

    % Marginal processes
    X_1_AV = Y_1_AV + a_USA * Z_AV;
    X_2_AV = Y_2_AV + a_EU * Z_AV;

    % General vector
    Xt_AV = [X_1_AV X_2_AV];

    stock_AV = real(S0 .* exp((rates + drift_compensator) * TTM + Xt_AV));


end % function stock_simulation
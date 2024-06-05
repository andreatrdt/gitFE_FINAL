function branch_out_procedure(data_calib_EU, data_calib_USA, F0_EU, B_EU, F0_USA, B_USA, date_settlement)
% Function to show that the calibration with too many dates is erroneous
% and worsen the results approximated
% 
% INPUT:
% 
% OUTPUT:
% 

    %% Removal unused dates
    % In order to simplify the calibration we tried to remove the dates
    % after 1y than, theoretically are not useful to price our derivative:
    % EU tille the 11th date, USA till the 13th date

    [data_calib_EU, F0_EU, B_EU] = removal_expiry(data_calib_EU, F0_EU, B_EU, [12 13]);
    [data_calib_USA, F0_USA, B_USA] = removal_expiry(data_calib_USA, F0_USA, B_USA, [14:1:20]);

    x0 = [0.3 -0.5 0.15 0.3 -0.5 0.15];

    %% First calibration   
    
    % Linear inequality constraints
    A = [-1 0 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 0 -1];
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
    params_USA = params_marginals(1:3); params_EU = params_marginals(4:6);

    % Explicit useful params
    k1 = params_USA(1); k2 = params_EU(1);

    %%

    %% Calibration over the rho 
    
    rho_historical = 0.801;

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
    
    nu_USA = params(1);
    nu_EU = params(2);
    nu_z = params(3);

    %% Common and idiosynchratic parameters
    
    % Compute the rho obtained from the model
    rho_model_Levy = sqrt(params(1) * params(2) / ((params(1) + params(3))*(params(2) + params(3))));
    
    error_rho = abs(rho_model_Levy - rho_historical);
    disp('error for rho calibration: \n ')
    disp(error_rho);
    
    % Compute the solution of the system
    sol =  marginal_param(params_USA,params_EU , nu_z , rho_model_Levy);
    
    % Explicitate the parameters
    a_USA = sol.x(1); a_EU = sol.x(2);
    beta_z = sol.x(3); gamma_z = sol.x(4);
    
    % Remaining parameters computation
    beta_USA = params_USA(2) - a_USA * beta_z;                       % From condition (10.1)
    gamma_USA = sqrt((params_USA(3)^2) - (a_USA^2) * (gamma_z^2));   % From condition (10.2)
    % nu_USA = (nu_z * params_USA(1))/(nu_z - params_USA(1)); 
    
    beta_EU = params_EU(2) - a_EU * beta_z;
    gamma_EU = sqrt(params_EU(3)^2 - a_EU ^ 2 * gamma_z^2);
    % nu_EU = (nu_z * params_EU(1))/(nu_z - params_EU(1));
    
    % Creation of zipped vectors
    
    % Idiosyncratic parameters USA
    idiosync_USA = [nu_USA, beta_USA, gamma_USA, a_USA];
    
    % Idiosyncratic parameters EU
    idiosync_EU = [nu_EU, beta_EU, gamma_EU, a_EU];
    
    % Systematic parameters
    syst_Z = [nu_z , beta_z, gamma_z];

    %%

    %% Point 9: Pricing of the certificate - LEVY
    
    year_to_maturity = 1;
    
    [rate_USA, TTM] = interp_pricing_params(datenum(data_calib_USA.datesExpiry), B_USA, date_settlement, year_to_maturity);
    [rate_EU, ~] = interp_pricing_params(datenum(data_calib_EU.datesExpiry), B_EU, date_settlement, year_to_maturity);
    
    S0_USA = data_calib_USA.spot; S0_EU = data_calib_EU.spot;
    
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

end
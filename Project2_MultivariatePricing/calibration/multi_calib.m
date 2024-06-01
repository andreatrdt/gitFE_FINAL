function multi_calib(data_EU, data_USA, F_EU, B0_EU, F_USA, B0_USA, date_settlement)
% Multi-calibration to study the impact of each single maturity on the
% overall RMSE
% 
% INPUT:
% data_EU:            [STRUCT] dataset of the European mkt
% data_USA:           [STRUCT] dataset of the American mkt
% F0_EU:              [VECTOR] initial forward F(0, T) EU
% B0_EU:              [VECTOR] initial discount B(0, T) EU
% F_USA:              [VECTOR] initial forward F(0, T) USA
% B0_USA:             [VECTOR] initial discount B(0, T) USA 
% date_settlement:    [DATENUM] settlement date of the computation
% 
% OUTPUT:
% 

    %% Initialization of the calib parameters

    x0 = [0.3 -0.5 0.15 0.3 -0.5 0.15];
    
    % Linear inequality constraints
    A = [-1 0 0 0 0 0;
        0 0 -1 0 0 0;
        0 0 0 -1 0 0;
        0 0 0 0 0 -1];
    
    b = [0; 0; 0; 0];
    
    % Linear equality constraints
    Aeq = []; beq = [];
    
    % Lower and upper bounds
    lb = [0.01; -Inf; 0; 0.01; -Inf; 0];
    ub = [];
    
    % Options for the visualization
    options = optimset('MaxFunEvals', 3e3, 'Display', 'iter');
    
    %% Rotation of the dates
    removal_expiry(data_calib_EU, F0_EU, B_EU, 12)

    %% Calibration
    params_marginals = fmincon(@(params) new_calibration(params, data_calib_EU, data_calib_USA, ...
        F0_EU, B_EU, F0_USA, B_USA, date_settlement), x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr(params), options);

    

end % function multi_calib
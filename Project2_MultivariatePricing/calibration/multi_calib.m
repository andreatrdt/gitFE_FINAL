function [params_removal, fun_eval, exit_condition] = multi_calib(data_EU, data_USA, F_EU, B0_EU, F_USA, B0_USA, date_settlement)
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
% params_removal:     [MATRIX] parameters for each calibration
% fun_eval:           [VECTOR] RMSE for each calibration
% exit_condition:     [VECTOR] exit condition of the solver
%
% USES:           calibration() , removal_expiry()

% Authors:
% M.Maspes, A.Tarditi, M.Torba


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
    
    %% Initialization of the structures to get all the calibrations

    parameters_removal_EU = zeros(length(data_EU.datesExpiry), 6);
    function_eval_EU = zeros(length(data_EU.datesExpiry), 1);

    parameters_removal_USA = zeros(length(data_USA.datesExpiry), 6);
    function_eval_USA = zeros(length(data_USA.datesExpiry), 1);

    exit_condition_EU = zeros(length(data_EU.datesExpiry), 1);
    exit_condition_USA = zeros(length(data_USA.datesExpiry), 1);

    %% Calibration

    for ii = 1:length(data_EU.datesExpiry)

        [data, F0, B0] = removal_expiry(data_EU, F_EU, B0_EU, ii);
        
        [parameters_removal_EU(ii, :), function_eval_EU(ii), exit_condition_EU(ii), ~]  = fmincon(@(params) calibration(params, data, data_USA, ...
            F0, B0, F_USA, B0_USA, date_settlement), x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr(params), options);
    end

    for ii = 1:length(data_USA.datesExpiry)

        [data, F0, B0] = removal_expiry(data_USA, F_USA, B0_USA, ii);
        
        [parameters_removal_USA(ii, :), function_eval_USA(ii), exit_condition_USA(ii), ~]  = fmincon(@(params) calibration(params, data_EU, data, ...
            F_EU, B0_EU, F0, B0, date_settlement), x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr(params), options);
    end

    %% Creation of the return structures

    params_removal = [parameters_removal_EU; parameters_removal_USA];
    fun_eval = [function_eval_EU; function_eval_USA];
    exit_condition = [exit_condition_EU; exit_condition_USA];

end % function multi_calib
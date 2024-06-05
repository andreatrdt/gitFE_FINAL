function find_x0(data_calib_EU, data_calib_USA, F0_EU, B_EU, F0_USA, B_USA, date_settlement)
% Find the initial parameters for the calibration
%
% INPUT:
% data_calib_EU:        [STRUCT] dataset of the options used for the
%                       calibration for the EU market   
% data_calib_USA:       [STRUCT] dataset of the options used for the
%                       calibration for the USA market
% F0_EU:                [VECTOR] initial forward value F(0, T) for the EU
%                       market
% B_EU:                 [VECTOR] initial discounts B(0, T) for the EU
%                       market
% F0_USA:               [VECTOR] initial forward value F(0, T) for the USA
%                       market
% B_USA:                [VECTOR] initial discounts B(0, T) for the USA
%                       market
% date_settlement:      [DATENUM] initial date
% OUTPUT:
% None
% USES:  calibration(), nonlinconstr(), disp_params(), nonlinconstr_corr(),
%        marginal_param(), disp_marginal_params()

% Authors:
% M.Maspes, A.Tarditi, M.Torba




vect = [0.1:0.1:1]';
matrix = repmat(vect,1,6);



matrix(:,2) = -matrix(:,2);
matrix(:,5) = -matrix(:,5);


for i = 1:length(vect)
    x0 = matrix(i,:)
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

    % Calibration
    params_marginals = fmincon(@(params) calibration(params, data_calib_EU, data_calib_USA, ...
        F0_EU, B_EU, F0_USA, B_USA, date_settlement), x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr(params), options);

    % Display of the parameters on the console
    params_USA = params_marginals(1:3);
    params_EU = params_marginals(4:6);

    disp_params(params_marginals, x0, 1);

    % Initialization of the parameters
    A = [-1 0 0; 0 -1 0; 0 0 -1]; b = [0; 0; 0];
    Aeq = []; beq = [];
    lb = zeros(1, 3); ub = [];
    
    x0 = ones(1, 3);
    
    k1 = params_USA(1); k2 = params_EU(1);
    rho_historical = 0.801;

    params = fmincon(@(params) (sqrt(params(1) * params(2) / ((params(1) + params(3))*(params(2) + params(3)))) - rho_historical)^2, ...
        x0, A, b, Aeq, beq, lb, ub, @(params) nonlinconstr_corr(params, k1, k2), options);
    
    nu_USA = params(1);
    nu_EU = params(2);
    nu_z = params(3);


    % Compute the rho obtained from the model
    rho_model_Levy = sqrt(params(1) * params(2) / ((params(1) + params(3))*(params(2) + params(3))));

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

    % Creation of vectors

    % Idiosyncratic parameters USA
    idiosync_USA = [nu_USA, beta_USA, gamma_USA, a_USA];

    % Idiosyncratic parameters EU
    idiosync_EU = [nu_EU, beta_EU, gamma_EU, a_EU];

    % Systematic parameters
    syst_Z = [nu_z , beta_z, gamma_z];

    disp_marginal_params(idiosync_USA , idiosync_EU , beta_z, gamma_z, nu_z, 1);

end
function [sol] = marginal_param(params_USA,params_EU,nu_z)
% This function computes the common and idiosyncratic parameters for the
% USA and EU markets
%
% INPUTS
% params_USA:      [VECTOR]the calibrated parameters for USA mkt
% params_EU:       [VECTOR]calibrated parameters for EU mkt
% nu_z:            [SCALAR]the calibrated nu_z
%
% OUTPUTS
% sol: the solution of the optimization problem
%
% USES: none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Unpacking of the parameters
    kappa_1 = params_USA(1);
    theta_1 = params_USA(2);
    sigma_1 = params_USA(3);

    kappa_2 = params_EU(1);
    theta_2 = params_EU(2);
    sigma_2 = params_EU(3);

    %% Initialization
    % The parameters supposed are
    % a_1 = x(1) - - - - a_2 = x(2)
    % Beta_z = x(3) - - - - gamma_z = x(4)

    x = optimvar('x',4);

    %% Composition of the equations

    % Conditions (9.1) BB
    eq1 = x(1) * x(3) - (kappa_1 * theta_1 / nu_z) == 0;
    eq2 = x(2) * x(3) - (kappa_2 * theta_2 / nu_z) == 0;

    % Conditions (9.2) BB
    eq3 = kappa_1 * sigma_1^2 - nu_z * x(1)^2 * x(4) ^2  == 0;
    eq4 = kappa_2 * sigma_2^2 - nu_z * x(2)^2 * x(4) ^2  == 0;

    % eq4 = x(1) * x(2) *( x(4)^2 + x(3)^2 * nu_z) / (sqrt(sigma_1^2 + theta_1 ^ 2 * kappa_1) * sqrt(sigma_2^2 + theta_2 ^ 2 * kappa_2)) - rho == 0;

    %% Solving

    prob = eqnproblem;
    prob.Equations.eq1 = eq1;
    prob.Equations.eq2 = eq2;
    prob.Equations.eq3 = eq3;
    prob.Equations.eq4 = eq4;

    % Initial value
    x0.x = ones(4,1);

    % Solution
    options = optimset('MaxFunEval', 5e3);
    
    sol = solve(prob,x0);

end % function marginal_param()
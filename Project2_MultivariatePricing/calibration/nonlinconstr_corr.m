function [c, ceq] = nonlinconstr_corr(params, params_USA, params_EU)
% Computation of the non linear constraints for the nu_Z, both the
% equality and the inequality ones
% 
% INPUT:
% params:       [VECTOR] [nu_1, nu_2, nu_z]
% k1:           [SCALAR] value calibrated for the USA mkt
% k2:           [SCALAR] value calibrated for the EU mkt
% 
% OUTPUT:
% c:            inequality constraints
% ceq:          equality constraints
%
% USES:        none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Unpacking the constraints
    nu_1 = params(1);
    nu_2 = params(2);
    nu_z = params(3);

    %% Unpacking of the other parameters
    
    k1 = params_USA(1);
    k2 = params_EU(1);

    theta1 = params_USA(2);
    tehta2 = params_EU(2);
    
    sigma1 = params_USA(3);
    sigma2 = params_EU(3);

    %% Computation of the other parameters

    sol = marginal_param(params_USA,params_EU,nu_z);

    % Explicitate the parameters
    a_USA = sol.x(1); a_EU = sol.x(2);
    beta_z = sol.x(3); gamma_z = sol.x(4);

    % Remaining parameters computation
    beta_USA = params_USA(2) - a_USA * beta_z;                       % From condition (10.1)
    gamma_USA = sqrt((params_USA(3)^2) - (a_USA^2) * (gamma_z^2));   % From condition (10.2)

    beta_EU = params_EU(2) - a_EU * beta_z;
    gamma_EU = sqrt(params_EU(3)^2 - a_EU ^ 2 * gamma_z^2);

    %% Constraints on the equalities

    ceq = [% nu_1 * nu_z /(nu_1 + nu_z) - k1;
        % nu_2 * nu_z /(nu_2 + nu_z) - k2;
        sqrt((nu_1 * nu_2) / ((nu_1 + nu_z) * (nu_2 + nu_z))) - (sqrt(k1 * k2)/ nu_z);
        sigma1^2/(theta1^2*k1)-gamma_z^2/(beta_z^2*nu_z)];
           
    %% Constraints on the inequalities

    c = [-gamma_USA;
        -gamma_EU;
        -gamma_z;
        -a_EU*a_USA;
        -((sigma1)^2 - (a_USA^2) * (gamma_z^2));
        -(sigma2^2 - a_EU ^ 2 * gamma_z^2)];

end % function nonlinconstr_corr
function [c, ceq] = nonlinconstr_corr(params, k1, k2)
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

    %% Constraints on the equalities

    ceq = [nu_1 * nu_z /(nu_1 + nu_z) - k1;
        nu_2 * nu_z /(nu_2 + nu_z) - k2;
        sqrt((nu_1 * nu_2) / ((nu_1 + nu_z) * (nu_2 + nu_z))) - (sqrt(k1 * k2)/ nu_z)];
           
    %% Constraints on the inequalities

    c = [];

end % function nonlinconstr_corr
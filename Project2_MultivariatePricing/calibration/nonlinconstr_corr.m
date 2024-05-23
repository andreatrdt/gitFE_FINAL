function [c, ceq] = nonlinconstr_corr(params, k1, k2)
% Computation of the non linear constraints for the nu_Z, both the
% equality and the inequality ones
% 
% INPUT:
% params:       [VECTOR] [nu_1, T1, S1, K2, T2, S2]
% k1:           [SCALAR] value calibrated for the USA mkt
% k2:           [SCALAR] value calibrated for the EU mkt
% 
% OUTPUT:
% c:            inequality constraints
% ceq:          equality constraints

    %% Constraints on the equalities

    % Constraint to create the entire equality given on the final
    % parameters
    ceq = x(3)^2/(x(1)*x(2)^2) - x(6)^2/(x(4)*x(5)^2);

    %% Constraints on the inequalities

    c = [];
end
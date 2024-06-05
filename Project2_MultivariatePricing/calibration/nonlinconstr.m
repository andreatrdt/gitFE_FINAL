function [c, ceq] = nonlinconstr(x)
% Computation of the non linear constraints for the solver, both the
% equality and the inequality ones
% 
% INPUT:
% x:            [VECTOR] [K1, T1, S1, K2, T2, S2]
% 
% OUTPUT:
% c:            inequality constraints
% ceq:          equality constraints
%
% USES:           none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Unpacking the parameters
    k1 = x(1); theta1 = x(2); sigma1 = x(3);
    k2 = x(4); theta2 = x(5); sigma2 = x(6);

    %% Constraints on the equalities

    % Constraint to create the entire equality given on the final
    % parameters
    ceq = (sigma1^2/(k1*theta1^2)) - (sigma2^2/(k2*theta2^2));

    %% Constraints on the inequalities

    c = [];
    
end % function nonlinconstr
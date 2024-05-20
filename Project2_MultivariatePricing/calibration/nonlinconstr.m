function [c, ceq] = nonlinconstr(x, rho_mkt, eps_th)
% Computation of the non linear constraints for the solver, both the
% equality and the inequality ones
% 
% INPUT:
% x:            [VEC] [K1, T1, S1, K2, T2, S2, Nz, Bz, Gz, a1, a2]
% rho_mkt:      correlation coefficients of the market
% eps:          arbitrary threshold
% 
% OUTPUT:
% c:            inequality constraints
% ceq:          equality constraints
% 
% USES:

    %% Constraints on the equalities

    % Constraint to create the entire equality given on the final
    % parameters
    ceq = [];
%     ceq = [ x(1)*x(2)-x(7)*x(10)*x(8);
%         x(1)*x(3)^2-x(7)*x(10)^2*x(9)^2];

    % ceq = [ x(3)^2/(x(1)*x(2)^2) - x(6)^2/(x(4)*x(5)^2)];
    %% Constraints on the inequalities

    % Constraints to set the correlation cefficients

    rho_model = x(10) * x(11) * (x(9)^2 + x(8)^2 * x(7)) ...
        /(sqrt((x(3)^2 + x(2)^2 * x(1))) * sqrt((x(6)^2 + x(5)^2 * x(4))));
    c = abs(rho_model - rho_mkt) - eps_th;
    % c = [];
end
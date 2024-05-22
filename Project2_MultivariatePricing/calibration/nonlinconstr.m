function [c, ceq] = nonlinconstr(x)
% Computation of the non linear constraints for the solver, both the
% equality and the inequality ones
% 
% INPUT:
% x:            [VEC] [K1, T1, S1, K2, T2, S2, Nz, Bz, Gz, a1, a2]
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
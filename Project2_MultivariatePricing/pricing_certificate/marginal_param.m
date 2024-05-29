function [sol] = marginal_param(params,nu_z)
%
% This function computes the common and idiosyncratic parameters for the
% USA and EU markets
%
% INPUTS
% params: the calibrated parameters
% nu_z: the calibrated nu_z
%
% OUTPUTS
% sol: the solution of the optimization problem
%
% USES: marginal_param()

    kappa = params(1);
    theta = params(2);
    sigma = params(3);

    x = optimvar('x',5);
    eq1 = x(1) * x(2) - kappa * theta / nu_z == 0;
    eq2 = (x(3)*x(1))^2 - kappa * sigma^2 / nu_z == 0;
    eq3 = x(4) +x(1)*x(2) - theta == 0;
    eq4 = x(5)^2 + x(3)^2 * x(3)^2 - sigma^2 == 0;
    eq5 = x(3)^2/x(2)^2 - (sigma^2 * nu_z) / (kappa * theta^2) == 0;

    prob = eqnproblem;
    prob.Equations.eq1 = eq1;
    prob.Equations.eq2 = eq2;
    prob.Equations.eq3 = eq3;
    prob.Equations.eq4 = eq4;
    prob.Equations.eq5 = eq5;


    x0.x = 0.1 * ones(5,1);

    sol = solve(prob,x0);


end % function marginal_param(params,nu_z)
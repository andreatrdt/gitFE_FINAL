function [c, ceq] = nonlinconstr(a_1, a_2, beta_1, nu_1, gamma_1, beta_2, nu_2, gamma_2, ...
    beta_z, nu_z, gamma_z)
% Computation of the non linear constraints for the solver
% 
% INPUT:
% 
% OUTPUT:
% 

    %% General function handles for the writing

    sigma2_i = @(gamma2_i, a_i) gamma2_i + a_i^2 * gamma_z;
    theta_i = @(beta_i, a_i) beta_i + a_i * beta_z;
    k_i = @(nu_i) (nu_i * nu_z)/(nu_i + nu_z)

    %% Constraints on conditions 10
    
    constraint_1 = sigma2_i(gamma_1, a_1)/(k_i(nu_1) * (theta_i(beta_1, a_1)^2)) - gamma_z/(nu_z * beta_z^2);
%     constraint_2 = sigma2_i(gamma_2, a_2)/(k_i(nu_2) * (theta_i(beta_2, a_2)^2)) - gamma_z/(nu_z * beta_z^2);

    %% Creation of the final vectors on inequalities

    c = [];
    ceq = constraint_1;

end
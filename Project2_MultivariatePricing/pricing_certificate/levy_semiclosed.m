function price = levy_semiclosed(s1_0, rate1, rate2, params_USA, params_EU, idiosync_USA, idiosync_EU, syst_Z, TTM)
% Computation of the semi closed formula obtained by proof
% 
% INPUT:
% s1_0:         [SCALAR] initial stock value
% rate1         [SCALAR] rate USA mkt
% rate2         [SCALAR] rate EU mkt
% sigma1        [SCALAR] sigma USA mkt
% sigma2        [SCALAR] sigma EU mkt
% rho           [SCALAR] historical correlation
% TTM           [SCALAR] time to maturity
% OUTPUT:
% price : price od the option obtained via the semiclosed formula
% USES : none

% Authors:
% M.Maspes, A.Tarditi, M.Torba
    
    %% Unpacking of the parameters
    
    kappa_USA = params_USA(1);
    theta_USA = params_USA(2);
    sigma_USA = params_USA(3);

    kappa_EU = params_EU(1);
    theta_EU = params_EU(2);
    sigma_EU = params_EU(3);

    nu_z = syst_Z(1);
    nu_USA = idiosync_USA(1);
    nu_EU = idiosync_EU(1);
    
    beta_z = syst_Z(2);
    beta_USA = idiosync_USA(2);
    beta_EU = idiosync_EU(2);

    gamma_z = syst_Z(3);
    gamma_USA = idiosync_USA(3);  
    gamma_EU = idiosync_EU(3);

    a_USA = idiosync_USA(4);
    a_EU = idiosync_EU(4);
    
    % beta_z = -(1/2+beta_z)*gamma_z^2;
    % beta_USA = -(1/2+beta_USA)*gamma_USA^2;
    % beta_EU = -(1/2+beta_EU)*gamma_EU^2;

    %% Drift compensators:

    drift_compensator_USA = - 1/kappa_USA * (1 - sqrt(1 - 2*kappa_USA*theta_USA - kappa_USA*sigma_USA^2));
    drift_compensator_EU = - 1/kappa_EU * (1 - sqrt(1 - 2*kappa_EU*theta_EU - kappa_EU*sigma_EU^2));

    %% Initialization of the parameters

    discount = exp(-rate1*TTM);

    %% Density of the Normal Inverse Gaussian

    A = @(beta, gamma) beta/(gamma^2);
    B = @(nu, beta, gamma) sqrt(beta^2 + (gamma^2)/nu)/(gamma^2);
    C = @(nu, beta, gamma) TTM/pi * exp(TTM/nu) * sqrt(beta^2/(nu * gamma^2) + 1/nu^2);

    supp = @(nu, gamma, x) sqrt(x.^2 + TTM^2 * gamma^2/nu);

    density_NIG = @(nu, beta, gamma, x)  C(nu, beta, gamma) .* exp(A(beta, gamma).*x) ...
                    .* besselk(1, B(nu, beta, gamma) .* supp(nu, gamma, x))./supp(nu, gamma, x);

    %% Creation of the integrands

    d1 = @(z) arrayfun(@(zval) integral(@(y) exp(y).*density_NIG(nu_USA, beta_USA, gamma_USA, y),-drift_compensator_USA*TTM-a_USA.*zval-rate1*TTM, 100), z);
    d2 = @(z) 1 - arrayfun(@(zval) integral(@(y) density_NIG(nu_USA, beta_USA, gamma_USA, y),-100, -drift_compensator_USA*TTM - a_USA*zval - rate1*TTM), z); 
    d3 = @(z) arrayfun(@(zval) integral(@(y) density_NIG(nu_EU, beta_EU, gamma_EU, y),-60, log(0.95) - drift_compensator_EU*TTM - a_EU*zval - rate2*TTM), z);

    integrand = @(z) ( ( exp(rate1*TTM + drift_compensator_USA*TTM + a_USA.*z) .* d1(z) - d2(z) ) .* d3(z) .*  density_NIG(nu_z, beta_z, gamma_z, z));

    %% Computation of the closed integral
    
    price = s1_0*discount* integral( integrand, -80, 80);

end % function blk_semiclosed
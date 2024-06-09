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

    % USA

    A_USA = beta_USA/(gamma_USA^2);
    B_USA = sqrt(beta_USA^2 + (gamma_USA^2)/nu_USA)/(gamma_USA^2);
    C_USA = TTM/pi * exp(TTM/nu_USA) * sqrt(beta_USA^2/(nu_USA * gamma_USA^2) + 1/nu_USA^2);

    supp_USA = @(x) sqrt(x.^2 + TTM^2 * gamma_USA^2/nu_USA);

    density_NIG_USA = @(x) C_USA .* exp(A_USA.*x) .* besselk(1,B_USA .* supp_USA(x))./supp_USA(x);

    % EU

    A_EU = beta_EU/(gamma_EU^2);
    B_EU = sqrt(beta_EU^2 + (gamma_EU^2)/nu_EU)/(gamma_EU^2);
    C_EU = TTM/pi * exp(TTM/nu_EU) * sqrt(beta_EU^2/(nu_EU * gamma_EU^2) + 1/nu_EU^2);

    supp_EU = @(x) sqrt(x.^2 + TTM^2 * gamma_EU^2/nu_EU);

    density_NIG_EU = @(x) C_EU .* exp(A_EU.*x) .* besselk(1,B_EU .* supp_EU(x))./supp_EU(x);

    % Z process

    A_z = beta_z/(gamma_z^2);
    B_z = sqrt(beta_z^2 + (gamma_z^2)/nu_z)/(gamma_z^2);
    C_z = TTM/pi * exp(TTM/nu_z) * sqrt(beta_z^2/(nu_z * gamma_z^2) + 1/nu_z^2);

    supp_z = @(x) sqrt(x.^2 + TTM^2 * gamma_z^2/nu_z);

    density_NIG_z = @(x) C_z * exp(A_z.*x) .* besselk(1,B_z .* supp_z(x))./supp_z(x);

    %% Creation of the integrands

    % d1 = @(z) arrayfun(@(zval) integral(@(y) exp(y) .* density_NIG(y), -drift_compensator_USA * TTM - a_USA.*zval -rate1*TTM, 50), z);
    % d2 = @(z) 1 - arrayfun(@(zval) cdf('InverseGaussian', -drift_compensator_USA*TTM - a_USA*zval - rate1*TTM, 1, TTM/nu_USA), z);
    % d3 = @(z) arrayfun(@(zval) cdf('InverseGaussian', log(0.95) - drift_compensator_EU*TTM - a_EU*zval - rate2*TTM, 1, TTM/nu_EU), z);
    % 
    % integrand = @(z) ( ( exp(rate1*TTM + drift_compensator_USA*TTM + a_USA.*z) .* d1(z) - d2(z) ) .* d3(z) .* density_NIG_z(z));
    
    d1 = @(z) arrayfun(@(zval) integral(@(y) exp(y).*density_NIG_USA(y),-drift_compensator_USA*TTM-a_USA.*zval-rate1*TTM, 100), z);
    d2 = @(z) 1 - arrayfun(@(zval) integral(@(y) density_NIG_USA(y),-100, -drift_compensator_USA*TTM - a_USA*zval - rate1*TTM), z); 
    d3 = @(z) arrayfun(@(zval) integral(@(y) density_NIG_EU(y),-60, log(0.95) - drift_compensator_EU*TTM - a_EU*zval - rate2*TTM), z);

    integrand = @(z) ( ( exp(rate1*TTM + drift_compensator_USA*TTM + a_USA.*z) .* d1(z) - d2(z) ) .* d3(z) .*  density_NIG_z(z));

    %% Computation of the closed integral
    
    price = s1_0*discount* integral( integrand, -80, 80);

end % function blk_semiclosed
function callPrices = callIntegral(B0, F0, kappa, eta, sigma, t, log_moneyness)
% Compute the price of a call option using the integral of NIG model
%
% INPUT:
%   B0: discount factor at time 0
%   F0: forward price at time 0
%   kappa: vol of vol
%   eta: skewness
%   sigma: variance of the model
%   t: time to maturity
%   log_moneyness: log of the moneyness to compute the price at
%
% OUTPUT:
%   callPrices: price of the call option (same size as log_moneyness)

% 
alpha = 1/2;

% Compute the Laplace exponent as a function of omega and alpha
ln_L = @(omega) t/kappa * (1 - alpha)/alpha * ...
        (1 - (1 + (omega .* kappa * sigma^2)/(1-alpha)).^alpha );

% compute the laPlace exponent at eta
ln_L_eta = ln_L(eta);

% Compute the characteristic function
phi = @(xi) exp(-1i * xi * ln_L_eta) .* exp( ln_L (0.5 * ((xi.^2) + 1i * (1+2*eta) .* xi)));
% @(u) exp(t.*(1./kappa.*(1-sqrt(1-2i*u.*kappa.*eta+u.^2.*kappa.*sigma.^2))));


% Compute the integral
I = integralQuad(phi, log_moneyness);

% apply the lewis formula
callPrices = B0 .* F0 .* (1 - exp(-log_moneyness/2) .* I);
    
end
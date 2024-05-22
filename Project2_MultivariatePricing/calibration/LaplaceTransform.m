function laplace_trasform = LaplaceTransform (deltaT, k, alpha, w, sigma)
% Computation of the Laplace Transform for a general w: L(w)
%
%INPUT:
% deltaT:                    interval of time
% k:                         volvol of the dynamic
% alpha:                     parameter to define a particular model
% w:                         integrand function
% sigma:                     volatility of the dynamic

    laplace_trasform = exp(deltaT/k * (1-alpha)/alpha * (1 - (1 + (w.*k.*sigma.^2)/(1-alpha)).^alpha));

end % function LaplaceTransform
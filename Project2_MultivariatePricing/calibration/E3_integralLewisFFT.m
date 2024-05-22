function [Z, integral] = E3_integralLewisFFT(M, dz, alpha, sigma, k, eta, timeToMaturity)
% Computation of the integral of Lewis formula in case of FFT
%
%INPUT
% M:              parameter of the grid
% dz:             parameter of the moneyness grid
% alpha:          alpha of the model
% sigma:          volatility of the model
% k:              vol of vol of the model
% eta:            simile asimmetry of the model
% timeToMaturity: TTM
%
%OUTPUT
% X:              grid of moneyness
% integral:       integral values 
% 
%USES:
% funciton parametersFFT(M, dz)
% function LaplaceTransform(deltaT, k, alpha, w, sigma)

    % Computation of FFT parameters
    [N, x1, xN, dx, z1, zN] = parametersFFT(M, dz);

    % Computation of the grid
    X = (x1: dx: xN)';
    Z = (z1: dz: zN)';

    % Computation of prefactor
    preFactor = dx * exp(-1i*x1.*Z);

    % Computation of FFT input
    xi = -X-1i/2;
    charFct = exp(-1i*xi*log(LaplaceTransform(timeToMaturity, k, alpha, eta, sigma))) .* ...
        LaplaceTransform(timeToMaturity, k, alpha, (xi.^2+1i*(1+2*eta).*xi)/2, sigma);
    integrand = 1/(2*pi) * 1./(X.^2+1/4) .* charFct; 
    j = (0:1:N-1)';
    inputFFT = integrand .* exp(-1i*z1*dx.*j);
    
    % Computation of the integral
    integral = preFactor .* fft(inputFFT);

end % function E3_integralLewisFFT
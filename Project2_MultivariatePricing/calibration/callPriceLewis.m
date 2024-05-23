function price = callPriceLewis(B0, F0, log_moneyness, sigma, k, theta, TTM, M, dz)
% Call price with Lewis formula
%
%INPUT
% B0:                 [SCALAR] discount factor B(t0, t)
% F0:                 [FORWARD] initial fwd price
% log_moneyness:      [VECTOR] vector of values of moneyness
% sigma:              [SCALAR] volatility of the model
% k:                  [SCALAR] vol of vol of the model
% theta:              [SCALAR] simile asimmetry of the model
% TTM:                [SCALAR} time to maturity
% M:                  parameter for FFT
% dz:                 parameter of the moneyness grid
% 
%USES:
% function E3_integralLewisFFT(M, dz, alpha, sigma, k, eta, timeToMaturity)

    %% Introduction of the FFT parameters
    
    % Initial constraints for FFt parameters
    N = 2^M;

    dx = (2*pi)/(N*dz); xN = ((N-1)*dx)/2; x1 = -xN;
    zN = ((N-1)*dz)/2; z1 = -zN;

    %% Computation of the grid
    X = (x1: dx: xN)';
    Z = (z1: dz: zN)';

    % Computation of prefactor
    preFactor = dx * exp(-1i*x1.*Z);

    % Computation of FFT input
    xi = -X-1i/2;

    %% Computation of the characteristic function
    charFct = exp(TTM .* (1/k * (1 - sqrt(1 - 2i .* xi .* k .* theta + xi.^2 .* k .* sigma.^2))));

    %% Computation of the integral
    integrand = 1/(2*pi) * 1./(X.^2+1/4) .* charFct; 
    j = (0:1:N-1)';
    inputFFT = integrand .* exp(-1i*z1*dx.*j);
    
    % Computation of the integral
    IntLewis = preFactor .* fft(inputFFT);

    IntLewis = real(IntLewis);
    IntLewis = interp1(Z, IntLewis, log_moneyness);

    %% Computation of the Call Price through Lewis
    
    price = B0 * F0 * (1 - exp(-log_moneyness./2) .* IntLewis);
    
end % function callPriceLewis
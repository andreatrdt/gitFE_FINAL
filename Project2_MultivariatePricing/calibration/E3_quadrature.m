function I = E3_quadrature(x_grid, sigma, TTM, k, eta, alpha)
% Computation of the Call Prices through the Quadrature method
% 
%INPUT:
% x_grid:                     logmoneyness grid of the quadrature
% sigma:                      volatility of the contract
% TTM:                        time to maturity
% k:                          volvol of the contract
% eta:                        simmetry of the contract
% alpha:                      parameter to define a particular model
% 
%USES:
% function LaplaceTransform (TTM, k, alpha, eta, sigma)

    %% Creation of the z grid
    
    % The function integral used for the quadrature already takes
    % consideration of the entire space, thus is not required to write a
    % discrete grid
    
    %% Computation of the integral - Riemann

    % Computation of L(eta)
    laplace_transform_eta = LaplaceTransform (TTM, k, alpha, eta, sigma);
    
    % Computation of Characteristic Funciton - Phi
    w = @(xi) (1i .* xi .* (1/2 + eta) + 1/2 .* xi.^2);
    phi = @(xi) exp(-1i .* xi .* log(laplace_transform_eta)) .* LaplaceTransform(TTM, k, alpha, w(xi), sigma);

    % Computation of the Lewis integral
    I = zeros(length(x_grid), 1);
    
    for j = 1:length(x_grid)
        I(j) = 1/(2*pi) * real(integral(@(xi) phi(-xi -1i/2) .* exp(-1i .* xi .* x_grid(j)) ./ (xi.^2 + 0.25), -50, 50));
    end
    
end % function E3_quadrature
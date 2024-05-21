function I = integralQuad(phi, queryPoints)
% Compute the price of the integral using the quadrature method between lower and upper
%
% Inputs:
%   phi: function to integrate
%   queryPoints: points to evaluate the integral
%
% Outputs:
%   I: integral of the function phi
%

% Transform the function to match the Lewis definition
f = @(xi,x) exp(-1i * xi .* x) / (2 * pi) .* phi(-xi-1i/2) .* 1 ./ (xi.^2 + 1/4);

% Define the function I(x) which computes the integral with respect to u
Integral_function = @(x) integral(@(u) f(u, x), -200, 200); 

% Compute the function I(x) for a specific value of x
I = real((arrayfun(@(x) Integral_function(x), queryPoints)));

end
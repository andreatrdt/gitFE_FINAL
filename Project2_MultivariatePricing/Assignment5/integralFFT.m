function I = integralFFT(phi, M, dz, queryPoints)
% I = FFT(integrand, M, xi_1)
%
% This function computes the Fourier Transform of the input integrand
% using the FFT algorithm.
%
% INPUTS:
%   f: function handle to the integrand
%   M: N = 2^M number of points to use in the FFT
%   xi_1: the fourier inferior limit
%
% OUTPUTS:
%   I: The integral of the integrand
%

% compute N
N = 2^M;

% compute the x values
z_1 = -(N-1)/2 * dz;
z = z_1:dz:-z_1;

% compute the dxi value
d_xi = 2 * pi / (N * dz);
xi_1 = -(N-1)/2 * d_xi;
xi = xi_1:d_xi:-xi_1;


% use the lewis formula to compute the function to integrate
f = 1 / (2*pi) *  phi(-xi - 1i/2) ./ (xi.^2 + 1/4);
f_tilde = f .* exp(-1i * z_1 * d_xi .* (0:N-1));

% compute the FFT
FFT = fft(f_tilde);

% compute the prefactor
prefactor = d_xi * exp(-1i * xi_1 * z);

% compute the integral by multiplying by the prefactor
I = prefactor .* FFT;

% get only the real part
I = real(I);

% interpolate the values
I = interp1(z, I, queryPoints);

end
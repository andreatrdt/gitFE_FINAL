function [N, x1, xN, dx, z1, zN] = parametersFFT (M, dz)
% Computation of FFT parameters
%
%INPUT
% M:              parameter of the grid
% dz:             parameter of the moneyness grid
%
%OUTPUT
% N:              number of steps
% x1:             lower bound of grid x
% xN:             upper bound of grid x
% dx:             step grid x
% z1:             lower bound of grid moneyness
% zN:             upper bound of grid moneyness
    
    N = 2^M;
    dx = (2*pi)/(N*dz);
    xN = ((N-1)*dx)/2;
    x1 = -xN;
    zN = ((N-1)*dz)/2;
    z1 = -zN;

end % function parametersFFT
function price = E3_callPriceLewis(discount, F0, alpha, moneyness, sigma, k, eta, timeToMaturity, flag, M, dz)
% Call price with Lewis formula
%
%INPUT
% discount:           discount factor B(t0, t)
% F0:                 initial fwd price
% alpha:              alpha of the model
% moneyness:          vector of values of moneyness
% sigma:              volatility of the model
% k:                  vol of vol of the model
% eta:                simile asimmetry of the model
% timeToMaturity:     TTM
% flag:               1 = FFT, 2 = Quadrature
% M:                  parameter for FFT
% dz:                 parameter of the moneyness grid
% 
%USES:
% function E3_integralLewisFFT(M, dz, alpha, sigma, k, eta, timeToMaturity)
% function gridExtract (oldGrid, oldLewis, newGrid)
% function E3_quadrature(x_grid, sigma, TTM, k, eta, alpha)


    %% Initialization

    if (nargin < 10)
        M = 0;
        dz = 0;
    end

    %% Computation of the Lewis integral
    
    switch (flag)

        case 1  % FFT
            [Z, IntLewis] = E3_integralLewisFFT(M, dz, alpha, sigma, k, eta, timeToMaturity);
            IntLewis = real(IntLewis);
            IntLewis = gridExtract(Z, IntLewis, moneyness);

            % Grid chosen to get 2 values for each 1% of moneyness

        case 2  % Quadrature
            IntLewis = E3_quadrature(moneyness, sigma, timeToMaturity, k, eta, alpha);

        otherwise
            IntLewis = 0;
    end

    %% Computation of the Call Price through Lewis
    
    price = discount .* F0 .* (1 - exp(-moneyness./2) .* IntLewis');
    

end % function E3_callPriceLewis
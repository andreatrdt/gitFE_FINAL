function delta = blkdelta(F0, K, sigma, t0, t, flag)
% Computation of the delta of a European Call or Put Option via the Black
% formula.
% 
% INPUT:
% F0:               Forward price F(0,t)
% K:                Strike price
% sigma:            Black implied volatility of the option
% t0:               settlement date
% t:                expiry date
% flag:             [0: European Put Option; 1: European Call Option]
% 
% OUTPUT:
% delta:            delta of the option
    
    % Convenction Act/365 for yearfrac:
    Act365 = 3;
    % Time to maturity:
    TTM = yearfrac(t0,t,Act365);
    % Computation of d1:
    d1 = log(F0./K)./(sigma.*sqrt(TTM))+sigma.*sqrt(TTM)/2;
    
    % Computation of the Delta:
    if flag
        delta = normcdf(d1);    % delta for a European Call Option
    else
        delta = -normcdf(-d1);  % delta for a European Put Option
    end

end % function blkdelta
function delta = blkdelta(F0, K, B0, sigma, TTM, flag)
% Computation of the delta of a European Call or Put Option via the Black
% formula.
% 
% INPUT:
% F0:                     Forward price F(0,t)
% K:                      Strike price
% B0:                     discount factor
% sigma:                  Black implied volatility of the option
% TTM:                    time to maturity
% flag:                   [0: European Put Option; 1: European Call Option]
% 
% OUTPUT:
% delta:                  delta of the option

    %% Computation of d1:
    d1 = log(F0./K)./(sigma.*sqrt(TTM))+sigma.*sqrt(TTM)/2;
    
    %% Computation of the Delta:
    if flag
        delta = B0 * normcdf(d1);    % delta for a European Call Option
    else
        delta = -B0 *normcdf(-d1);  % delta for a European Put Option
    end

end % function blkdelta
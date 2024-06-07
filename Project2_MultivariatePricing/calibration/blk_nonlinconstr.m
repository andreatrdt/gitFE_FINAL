function [c, ceq] = blk_nonlinconstr(dataset, date_settlement, B0, F0,sigma)
% Computation of the non linear constraints for the solver, both the
% equality and the inequality ones
% 
% INPUT:
% sigma:        [SCALAR] volatility sigma
% 
% OUTPUT:
% c:            inequality constraints
% ceq:          equality constraints
%
% USES:         none

% Authors:
% M.Maspes, A.Tarditi, M.Torba

    %% Constraints on the equalities

    ceq = [];

    %% Constraints on the inequalities

    % Initialization:
    c = zeros(min(length(dataset.datesExpiry),19),1);

    % Computation of the current Put prices
    for ii = 1:min(length(dataset.datesExpiry),19)

        %% Initialization

        % Vector of strikes:
        strikes = dataset.strikes(ii).value;
        % Year frac convenction Act/365:
        conv_ACT365 = 3;           
        % Time to maturity
        TTM = yearfrac(date_settlement, datenum(dataset.datesExpiry(ii)), conv_ACT365);
        % Interest rate
        interest_rate = -log(B0(ii))/TTM;

        %% Pricing 
        
        % Price of Call/Puts through the Black formula
        call_prices = blkprice(F0(ii), strikes, interest_rate, TTM, sigma);
        put_prices = call_prices - B0(ii)*(F0(ii) - strikes);
        

        % Constraint avoid to have negative put prices:
        c(ii) = -min(put_prices);

        
    end


    
    
end % function blk_nonlinconstr
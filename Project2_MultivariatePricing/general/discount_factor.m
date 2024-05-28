function discount =  discount_factor(B_bar , data , date_settlement)
% Computation of the Discount Factor after the USA spread
% 
% INPUT:
% B_bar :               SCALAR [1x1] : discount factor at time 0
% data :                STRUCT : market data
% date_settlement :     SCALAR [1x1] : settlement date
% 
% OUTPUT
% discount : vector [1xN] : discount factor
%
% USES : discount_factor(B_bar , data , date_settlement)

    %% Conventions
    conv_ACT365 = 3;

    %% Computation of the dates
    dates = datenum(data.datesExpiry);
    dt = yearfrac(date_settlement ,   dates(1:end), conv_ACT365);
    
    % mean spread
    spread = 34*1e-4;
    
    %% Computation of the discount factor
    discount = B_bar .* exp( spread .* dt);

end % function discount_factor
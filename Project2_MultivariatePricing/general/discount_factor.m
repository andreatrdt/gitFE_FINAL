function discount =  discount_factor(B_bar , data , date_settlement)
%
% INPUT
% B_bar             [SCALAR]discount factor at time 0
% data :            [STRUCT]market data
% date_settlement   [SCALAR]settlement date
% OUTPUT
% discount
%
% USES : discount_factor()

dates = datenum(data.datesExpiry);

dt = yearfrac(date_settlement ,   dates(1:end));

% mean spread
spread = 34*1e-4;

% Compute the discount factor

discount = B_bar .* exp( spread .* dt);


end % function discount_factor
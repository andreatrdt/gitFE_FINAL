function discount =  discount_factor(B_bar , data , date_settlement)
%
% INPUT
% B_bar : scalar [1x1] : discount factor at time 0
% data : structure : market data
% date_settlement : scalar [1x1] : settlement date
% OUTPUT
% discount : vector [1xN] : discount factor
%
% USES : discount_factor(B_bar , data , date_settlement)

dates = datenum(data.datesExpiry);

dt = yearfrac(date_settlement ,   dates(1:end));

% mean spread
spread = 34*1e-4;

% Compute the discount factor

discount = B_bar .* exp( spread .* dt);


end % function discount_factor
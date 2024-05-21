function prices = prices_calibration(data, date_settlement, date,dates)
% Calibration of the prices for the options
% 
% INPUT:
% 
% OUTPUT:
% 
% USES:
% 

    Act365 = 3;

    idx = find_idx(data, date);

    % compute the forward in 0:
    [F_0, ~ , discount_at_expiry] = forward_prices(data, date, 1);

    % compute the log moneyess from the strikes
    log_moneyness = log(F_0(1,:) ./ data.strikes(idx).value);

    % time to maturity
    t = yearfrac(date_settlement,dates(idx), Act365);

    % create a function that the prices of the call options given the strikes
    prices = @(p) callIntegral(discount_at_expiry, F_0(1,:), p(1), p(2), p(3), t, log_moneyness);

    % compute the implied volatilities:
    volatility = @(p) blkimpv(F_0(1,:), data.strikes(idx).value, -log(discount_at_expiry)/t, t, prices_EU(p));

end
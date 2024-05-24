function prices = cp(data_calib_EU, F0_EU, B_bar_EU, params_marginals, idx, date_settlement)


put_length = length(data_calib_EU.putAsk(idx).prices);

log_moneyness = log(F0_EU(idx) ./ data_calib_EU.strikes(idx).value);
TTM = yearfrac(date_settlement, datenum(data_calib_EU.datesExpiry(idx)), 3);
prices = callPriceLewis(B_bar_EU(idx), F0_EU(idx), log_moneyness, params_marginals(6), params_marginals(4), params_marginals(5), TTM, 16, 0.0025);
call_prices = prices(put_length+1:end);
mean_call_price = (data_calib_EU.callAsk(idx).prices + data_calib_EU.callBid(idx).prices)/2;

put_prices = prices(1:put_length) - F0_EU(idx).* B_bar_EU(idx) + data_calib_EU.strikes(idx).value(1:put_length) .* B_bar_EU(idx);
mean_put_price = (data_calib_EU.putAsk(idx).prices + data_calib_EU.putBid(idx).prices)/2;

perc_var_call = (call_prices - mean_call_price)./mean_call_price * 100;
perc_var_put = (put_prices - mean_put_price)./mean_put_price * 100;

perc = [perc_var_put perc_var_call]

call_prices = prices
end
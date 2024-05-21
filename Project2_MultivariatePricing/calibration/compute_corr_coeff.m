function rho = compute_corr_coeff(data, date_sett, date_maturity)
% Compute the correlation coefficient between two series [EU and USA]
%
% INPUTS
% data:            matrix of the returns
% date_sett:       data di settlement
% date_maturity:   data della maturity
% 
% OUTPUTS
% rho: correlation coefficient between the two series

    %% Fixing of the orders of the returns

    % Computation of the daily returns
    returns_daily = flip(data.Daily);

    %% Creation of the interval

    days = datenum(date_maturity) - datenum(date_sett);

    %% Choice of the interval

    USA = returns_daily(1:days, 1);
    EU = returns_daily(1:days, 2);

    %% Correlation

    C = corrcoef(EU,USA);
    rho = C(1,2);

end % function compute_corr_coeff
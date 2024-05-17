function rho = compute_corr_coeff(data_EU,data_USA,date)
% compute_corr_coeff: compute the correlation coefficient between two series
%
%   INPUTS
%   F_EU: forward prices of the European options
%   F_USA: forward prices of the American options
%
%   OUTPUTS
%   rho: correlation coefficient between the two series

    % load data

    [F_vector, ~] = forward_prices(data_USA, date);
    F_USA = F_vector(1,:);
    [F_vector, ~] = forward_prices(data_EU, date);
    F_EU = F_vector(1,:);

    % find the index of the date

    idx = find_idx(data_EU.datesExpiry,date);

    % interpolate data for the mising strikes
    F_EU = interp1(data_EU.strikes(idx).value,F_EU,data_USA.strikes(idx).value,"linear","extrap");
    
    % compute the correlation coefficient

    C = corrcoef(F_EU,F_USA);

    rho = C(1,2);
end
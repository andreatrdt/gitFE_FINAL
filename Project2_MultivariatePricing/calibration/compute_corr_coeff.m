function rho = compute_corr_coeff(data_EU,data_USA,date)
% Compute the correlation coefficient between two series [EU and USA]
%
% INPUTS
% data_EU: data of the European market
% data_USA: data of the American market
% date: date until to compute the correlation
% 
% OUTPUTS
% rho: correlation coefficient between the two series

    %% 

    % load data


    % find the index of the date

    idx_EU = find_idx(data_EU.datesExpiry,date);
    idx_USA = find_idx(data_USA.datesExpiry,date);

    % interpolate data for the mising strikes
    F_EU = interp1(data_EU.strikes(idx_EU).value,F_EU,data_USA.strikes(idx_USA).value,"linear","extrap");
    
    % compute the correlation coefficient

    C = corrcoef(F_EU,F_USA);

    rho = C(1,2);
end
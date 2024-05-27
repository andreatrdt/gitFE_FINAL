function [rate, interp_F0] = interp_pricing_params(dates, F0, B0, date_settlement)
% Extraction of the interest rate and initial forwards
% 
% INPUT:
% dates:                [DATENUM VEC] dates of the dataset
% F0:                   [VECTOR] initial forward value
% B0:                   [VECTOR] initial discount value
% date_settlement:      [DATENUM] initial date of the certificate
% 
% OUTPUT:
% rate:                 [SCALAR] value of the interest rate
% interp_F0:            [SCALAR] initial forwards F(0, 1y)
% 
% USES:
% function rate_interpolation()

    %% Conventions

    conv_ACT365 = 3;
    
    %% Computation of the TTM
    interp_date = datenum(busdate(datetime(date_settlement, 'ConvertFrom', 'datenum') - caldays(1) + calyears(1), 1, eurCalendar));
    TTM = yearfrac(date_settlement, interp_date, conv_ACT365);

    %% Computation of interest rate

    rate = rate_interpolation(dates, B0, date_settlement);

    %% Computation of the interpolated F(0, 1y)

    interp_F0 = interp1(dates, F0, interp_date);

end % function interp_pricing_params
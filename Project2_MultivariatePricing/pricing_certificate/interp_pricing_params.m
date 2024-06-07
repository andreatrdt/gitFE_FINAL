function [rate, TTM] = interp_pricing_params(dates, B0, date_settlement, year_to_maturity)
% Extraction of the interest rate and time to maturity
% 
% INPUT:
% dates:                [DATENUM VEC] dates of the dataset
% B0:                   [VECTOR] initial discount value
% date_settlement:      [DATENUM] initial date of the certificate
% year_to_maturity:
% 
% OUTPUT:
% rate:                 [SCALAR] value of the interest rate
% TTM:                  [SCALAR] time to maturity
% 
% USES: rate_interpolation()

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Conventions

    conv_ACT365 = 3;
    
    %% Computation of the TTM
    interp_date = datenum(busdate(datetime(date_settlement, 'ConvertFrom', 'datenum') - caldays(1) + calyears(year_to_maturity), 1, eurCalendar));
    TTM = yearfrac(date_settlement, interp_date, conv_ACT365);

    %% Computation of interest rate

    rate = rate_interpolation(dates, B0, date_settlement, interp_date);

end % function interp_pricing_params
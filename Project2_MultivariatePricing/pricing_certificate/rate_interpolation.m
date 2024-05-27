function rate = rate_interpolation(dates, discounts, date_settlement)
% Interpolation of the 1y zero rate
% 
% INPUT:
% dates:                 dates of required points
% discounts:             discounts of the required points
% date_settlement:       initial date of computation
% 
% OUTPUT:
% rate:                  zero rate for the GK formula
% 
% USES:
% function zeroRates()

    %% Creation of the required support vector

    dates = [date_settlement; dates];
    discounts = [1; discounts];

    %% Computation of the Zero Rates

    rates = zeroRates(dates, discounts);

    %% Choice of the date

    interp_date = datenum(busdate(datetime(dates(1), 'ConvertFrom', 'datenum') - caldays(1) + calyears(1), 1, eurCalendar));

    %% Interpolation of the required rate

    rate = interp1(dates, rates, interp_date);

end % function rate_interpolation
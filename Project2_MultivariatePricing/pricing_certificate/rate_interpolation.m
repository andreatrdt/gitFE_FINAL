function rate = rate_interpolation(dates, discounts, date_settlement, interp_date)
% Interpolation of the 1y zero rate
% 
% INPUT:
% dates:                 dates of required points
% discounts:             discounts of the required points
% date_settlement:       initial date of computation
% interp_date:           date of interpolation
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

    %% Interpolation of the required rate

    rate = interp1(dates, rates, interp_date);

end % function rate_interpolation
function prices = stock_simulation(dataset, params, F0, B0, date_settlement)
% Pricing of the underlying process Si(t)
% 
% INPUT:
% 
% OUTPUT:
% 
% USES:
% 

    %% Conventions

    conv_ACT365 = 3;
    
    %% Unpacking of the parameters

    k = params(1); 
    theta = params(2);
    sigma = params(3);

    %% Computation of the support params

    drift_compensator = - 1/k * (1 - sqrt(1 - 2*k*theta - k*sigma^2)),

    % Computation of the TTM
    interp_date = datenum(busdate(datetime(dates(1), 'ConvertFrom', 'datenum') - caldays(1) + calyears(1), 1, eurCalendar));
    TTM = yearfrac(date_settlement, interp_date, conv_ACT365);

    % Computation of interest rate

    dates = datenum(dataset.datesExpiry);
    rate = rate_interpolation(dates, B0, date_settlement);

    % Computation of the initial stock
    


end
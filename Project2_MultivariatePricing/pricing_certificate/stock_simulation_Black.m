function [prices, S0] = stock_simulation_Black(sigmas, F0, rates, rho, date_settlement)
% Pricing of the underlying process Si(t)
% 
% INPUT:
% sigmas:               [VECTOR] volatilities of the Black process
% F0:                   [VECTOR] initial forward value
% B0:                   [VECTOR] initial discount value
% rho:                  [SCALAR] historical correlation
% date_settlement:      [DATENUM] initial date of the certificate
% 
% OUTPUT:
% prices:               [MATRIX] underlying stock to be simulated
% S0:                   [VECTOR] initial value of the stock
% 
% USES:
% function rate_interpolation()

    %% Conventions

    conv_ACT365 = 3;

    %% Unpacking of parameters

    %% Initialization

    nSim = 1e6;

    % Computation of the TTM
    interp_date = datenum(busdate(datetime(date_settlement, 'ConvertFrom', 'datenum') - caldays(1) + calyears(1), 1, eurCalendar));
    TTM = yearfrac(date_settlement, interp_date, conv_ACT365);

    %% Computation of the initial stock price S(0)

    S0 = F0 .* exp(- rates .* TTM);
    
    %% Simulation of the NIG process

    % Stochastic part
    correlationMatrix = [1 rho; rho 1];
    meanVector = [0; 0];
    
    g = mvnrnd(meanVector,correlationMatrix, nSim);
    
    % Creation of Xt dynamic

    Xt = -0.5 .* sigmas'.^2 .* TTM + sigmas' .* sqrt(TTM) .* g;

    %% Computation of the initial stock

    prices = S0' .* exp(rates' + Xt);

end % function stock_simulation_Black
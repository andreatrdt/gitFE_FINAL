function RMSE_total = RMSE_total(params, dataset, F0, B0, date_settlement)
% Computation of the total RMSE for the dataset given
% 
% INPUT:
% params:            [VECTOR] of params [k, theta, sigma]
% dataset:           [STRUCT] dataset preprocessed
% F0:                [VECTOR] forward values
% B0:                [VECTOR] discount values
% date_settlement:   [DATENUM] initial date
% 
% OUTPUT:
% RMSE:              error to minimize
% 
% USES:
% function callIntegral()

    %% Unpack the parameters
    k = params(1);
    theta = params(2);
    sigma = params(3);

    %% Conventions
    conv_ACT365 = 3;

    %% Initial parameters
    RMSE = zeros(length(dataset.datesExpiry), 1);
    N_options = 0;
    
    % FFT parameters
    M = 15;
    dz = 0.001;

    %% Computation 
    
    for ii = 1:length(dataset.datesExpiry)

        %% Initialization
        put_length = length(dataset.putAsk(ii).prices);

        N_options = N_options + length(dataset.strikes(ii).value);

        % Logmoneyness values
        log_moneyness = log(F0(ii) ./ dataset.strikes(ii).value);

        % Time to maturity
        TTM = yearfrac(date_settlement, datenum(dataset.datesExpiry(ii)), conv_ACT365);

        %% Pricing 
        % Price of Call/Puts through the FFT and Lewis formula

        prices = callPriceLewis(B0(ii), F0(ii), log_moneyness, sigma, k, theta, TTM, M, dz);
        
%         call_prices = max(prices(put_length+1:end), 0);
%         put_prices = max(prices(1:put_length) - F0(ii).* B0(ii) + dataset.strikes(ii).value(1:put_length) .* B0(ii), 0); 

        call_prices = prices(put_length+1:end);
        put_prices = prices(1:put_length) - F0(ii).* B0(ii) + dataset.strikes(ii).value(1:put_length) .* B0(ii); 
        
        mean_call_price = (dataset.callAsk(ii).prices + dataset.callBid(ii).prices)/2;
        mean_put_price = (dataset.putAsk(ii).prices + dataset.putBid(ii).prices)/2;

        %% Computation of RMSE

        RMSE(ii) = sum((call_prices - mean_call_price).^2) + sum((put_prices - mean_put_price).^2);
    end

    %% Final adjusting of RMSE
    RMSE_total = sqrt(sum(RMSE))/N_options;

end % function RMSE_total
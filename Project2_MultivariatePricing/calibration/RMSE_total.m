function RMSE = RMSE_total(params, dataset, F0, B0, date_settlement)
% Computation of the total RMSE for the dataset given
% 
% INPUT:
% params:            [k, theta, sigma]
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

    %% Computations
    RMSE = 0;
    N_options = 0;

    for ii = 1:3

        %% Initialization
        put_length = length(dataset.putAsk(ii).prices);
        N_options = N_options + length(dataset.strikes(ii).value);

        % Logmoneyness values
        log_moneyness = log(F0(ii) ./ dataset.strikes(ii).value);

        % Time to maturity
        TTM = yearfrac(date_settlement, dataset.datesExpiry(ii), conv_ACT365);

        %% Pricing 
        % Price of Call/Puts through the FFT and Lewis formula

%         prices = callIntegral(B0(ii), F0(ii), k, theta, sigma, TTM, log_moneyness);
        prices = E3_callPriceLewis(B0(ii), F0(ii), 1/2, log_moneyness, sigma, k, theta, TTM, 1, 19, 0.0048);
        
        call_prices = prices(put_length+1:end);
        put_prices = prices(1:put_length) - F0(ii) + dataset.strikes(ii).value(1:put_length) .* B0(ii); 

        mean_call_price = (dataset.callAsk(ii).prices + dataset.callBid(ii).prices)/2;
        mean_put_price = (dataset.putAsk(ii).prices + dataset.putBid(ii).prices)/2;

        RMSE = RMSE + sum((call_prices - mean_call_price).^2) + sum((put_prices - mean_put_price).^2);
    end

    %% Final rebalancing
    RMSE = sqrt(RMSE/N_options);

end % function RMSE_total
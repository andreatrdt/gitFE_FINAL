function plot_calls_puts(dataset, F0, B0, params, date_settlement)
% Plotting of the calls and put prices after the calibration
% 
% INPUT:
% dataset:            [STRUCT] initial dataset
% F0:                 [VECTOR] initial forward value F(0, T)
% B0:                 [VECTOR] initial discounts B(0, T)
% params:             [VECTOR] [k, theta, sigma]
% date_settlement:    [DATENUM] initial date
% 
% OUTPUT: none
% 
% USES: callPriceLewis_pref(), 
% 
% Authors:
% M.Maspes, A.Tarditi, M.Torba

    for ii = 1:length(dataset.datesExpiry)

        %% Initialization
    
        put_length = length(dataset.putAsk(ii).prices);
    
        conv_ACT365 = 3;

        %% Unpacking of the parameters

        k = params(1);
        theta = params(2);
        sigma = params(3);

        %% Computation of the prices
        
        % Parameters FFT
        M = 15;
        dz = 0.0025;

        % Parameters pricing
        strikes = dataset.strikes(ii).value;

        log_moneyness = log(F0(ii) ./ strikes);
        TTM = yearfrac(date_settlement, datenum(dataset.datesExpiry(ii)), conv_ACT365);
        
        prices = callPriceLewis_pref(B0(ii), F0(ii), log_moneyness, sigma, k, theta, TTM, M, dz);
        
        call_prices = prices(put_length+1:end);
        
        put_prices = prices(1:put_length) - F0(ii).* B0(ii) + strikes(1:put_length) .* B0(ii);
        
        %% Plots

        figure();
        subplot(1, 2, 1);
        plot(strikes(put_length+1:end), call_prices, '*-'); hold on;
        plot(strikes(put_length+1:end), dataset.callAsk(ii).prices); hold on;
        plot(strikes(put_length+1:end), dataset.callBid(ii).prices); grid on;
        title('Call prices'); xlabel('Strikes'); ylabel('Prices');
        legend('Calibrated prices','Location','best');
        
        subplot(1, 2, 2);
        plot(strikes(1:put_length), put_prices, '*-'); hold on;
        plot(strikes(1:put_length), dataset.putAsk(ii).prices); hold on;
        plot(strikes(1:put_length), dataset.putBid(ii).prices); grid on;
        title('Put prices'); xlabel('Strikes'); ylabel('Prices');
        legend('Calibrated prices','Location','best');
    
    end

end % function plot_calls_puts
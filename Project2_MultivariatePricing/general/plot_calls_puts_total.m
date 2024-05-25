function plot_calls_puts_total(dataset, F0, B0, params, date_settlement)
% Plotting of the calls and put prices after the calibration
% 
% INPUT:
% dataset:            [STRUCT] initial dataset
% F0:                 [VECTOR] initial forward value F(0, T)
% B0:                 [VECTOR] initial discounts B(0, T)
% params:             [VECTOR] [k, theta, sigma]
% date_settlement:    [DATENUM] initial date
% 
% USES:
% function callPriceLewis()

    for ii = 1:length(dataset.datesExpiry)

        %% Initialization
    
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
        
        call_prices = callPriceLewis(B0(ii), F0(ii), log_moneyness, sigma, k, theta, TTM, M, dz);
        put_prices = call_prices - B0(ii)*(F0(ii) - strikes);
        
        % Parameters comparison
        mean_call = (dataset.callBid(ii).prices + dataset.callAsk(ii).prices)/2;
        mean_put = (dataset.putBid(ii).prices + dataset.putAsk(ii).prices)/2;

        %% Plots

        figure();

        subplot(1, 2, 1);
        plot(strikes, call_prices, '*-'); hold on;
        plot(strikes, mean_call, 'o-'); grid on;
        xline(F0(ii), '--', 'LineWidth', 2, 'Color', 'r');
        title('Call prices at expiry ', dataset.datesExpiry(ii)); xlabel('Strikes'); ylabel('Prices');
        legend('Calibrated prices', 'Mean prices', 'Strike ATM');

        subplot(1, 2, 2);
        plot(strikes, put_prices, '*-'); hold on;
        plot(strikes, mean_put, 'o-'); grid on;
        xline(F0(ii), '--', 'LineWidth', 2, 'Color', 'r');
        title('Put prices at expiry ', dataset.datesExpiry(ii)); xlabel('Strikes'); ylabel('Prices');
        legend('Calibrated prices', 'Mean prices', 'Strike ATM');
        
    end

end % function plot_calls_puts_total
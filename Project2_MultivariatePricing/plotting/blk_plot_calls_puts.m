function blk_plot_calls_puts(dataset_calib, F0, B0, sigma, date_settlement)
% Plotting of the calls and put prices after the calibration
% 
% INPUT:
% dataset_calib:      [STRUCT] dataset of the options used for the
%                     calibration
% F0:                 [VECTOR] initial forward value F(0, T)
% B0:                 [VECTOR] initial discounts B(0, T)
% sigma:              [SCALAR] calibrated volatility
% date_settlement:    [DATENUM] initial date
% OUTPUT:
% None
% USES:  blkprice()

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    for ii = 1:length(dataset_calib.datesExpiry)

        %% Initialization
        
        put_length = length(dataset_calib.putAsk(ii).prices);

        % Vector of strikes:
        strikes = dataset_calib.strikes(ii).value;

        % Year frac convenction Act/365:
        conv_ACT365 = 3;

        %% Computation of the prices
             
        % Time to maturity
        TTM = yearfrac(date_settlement, datenum(dataset_calib.datesExpiry(ii)), conv_ACT365);
        
        % Interest rate
        interest_rate = -log(B0(ii))/TTM;

        %% Pricing 
        
        % Price of Call/Puts through the Black formula
        prices = blkprice(F0(ii), strikes, interest_rate, TTM, sigma);
        
        call_prices = prices(put_length+1:end);
        put_prices = prices(1:put_length) - F0(ii).* B0(ii) + strikes(1:put_length) .* B0(ii);

        % put_prices = call_prices - B0(ii)*(F0(ii) - strikes);
        
        % Parameters comparison
        mean_call = (dataset_calib.callBid(ii).prices + dataset_calib.callAsk(ii).prices)/2;
        mean_put = (dataset_calib.putBid(ii).prices + dataset_calib.putAsk(ii).prices)/2;

        %% Computation and update of the error:
        % [error_call_prices, error_put_prices] = error_calibration(call_prices, put_prices, ...
        %    dataset.callBid(ii).prices, dataset.callAsk(ii).prices, dataset.putBid(ii).prices, dataset.putAsk(ii).prices);

        %% Plots

        figure();

        subplot(1, 2, 1);
        plot(strikes(put_length+1:end), call_prices, '*-'); hold on;
        plot(strikes(put_length+1:end), mean_call, 'o-'); grid on;
        plot(strikes(put_length+1:end), dataset_calib.callAsk(ii).prices, '--');
        plot(strikes(put_length+1:end), dataset_calib.callBid(ii).prices, '--');
        xline(F0(ii), '--', 'LineWidth', 2, 'Color', 'r');
        title(['Black model Call prices at expiry ', datestr(dataset_calib.datesExpiry(ii))]); xlabel('Strikes'); ylabel('Prices');
        legend('Black calibrated prices', 'Mean prices', 'Call Ask', 'Call Bid', 'Strike ATM', 'Location', 'best');

        subplot(1, 2, 2);
        plot(strikes(1:put_length), put_prices, '*-'); hold on;
        plot(strikes(1:put_length), mean_put, 'o-'); grid on;
        plot(strikes(1:put_length), dataset_calib.putAsk(ii).prices, '--');
        plot(strikes(1:put_length), dataset_calib.putBid(ii).prices, '--');
        xline(F0(ii), '--', 'LineWidth', 2, 'Color', 'r');
        title(['Black Put prices at expiry ', datestr(dataset_calib.datesExpiry(ii))]); xlabel('Strikes'); ylabel('Prices');
        legend('Black calibrated prices', 'Mean prices', 'Put Ask', 'Put Bid', 'Strike ATM','Location', 'best');
        
    end

end % function plot_calls_puts_total
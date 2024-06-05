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
% USES:     none

% Authors:
% M.Maspes, A.Tarditi, M.Torba



%% Initialization of the error vectors:
error_call_prices_vec = zeros(length(dataset.datesExpiry),1);
error_put_prices_vec = zeros(length(dataset.datesExpiry),1);

    for ii = 1:min(length(dataset.datesExpiry),19)

        %% Initialization
    
        conv_ACT365 = 3;

        %% Unpacking of the parameters

        k = params(1);
        theta = params(2);
        sigma = params(3);

        %% Computation of the prices
        
        % Parameters FFT
        M = 15;
        dz = 0.001;

        % Parameters pricing
        strikes = dataset.strikes(ii).value;

        log_moneyness = log(F0(ii) ./ strikes);
        TTM = yearfrac(date_settlement, datenum(dataset.datesExpiry(ii)), conv_ACT365);
        
        call_prices = callPriceLewis_pref(B0(ii), F0(ii), log_moneyness, sigma, k, theta, TTM, M, dz);
        put_prices = call_prices - B0(ii)*(F0(ii) - strikes);
        
        % Parameters comparison
        mean_call = (dataset.callBid(ii).prices + dataset.callAsk(ii).prices)/2;
        mean_put = (dataset.putBid(ii).prices + dataset.putAsk(ii).prices)/2;

        %% Computation and update of the error:

        % Computation of the errors:
        [error_call_prices, error_put_prices] = error_calibration(call_prices, put_prices, ...
           dataset.callBid(ii).prices, dataset.callAsk(ii).prices, dataset.putBid(ii).prices, dataset.putAsk(ii).prices);
        
        % Update of the mean error vector:
        error_call_prices_vec(ii) = mean(error_call_prices);
        error_put_prices_vec(ii) = mean(error_put_prices);

        %% Plots

        figure();

        subplot(1, 2, 1);
        plot(strikes, call_prices, '*-'); hold on;
        plot(strikes, mean_call, 'o-'); grid on;
        plot(strikes, dataset.callAsk(ii).prices, '--');
        plot(strikes, dataset.callBid(ii).prices, '--');
        xline(F0(ii), '--', 'LineWidth', 2, 'Color', 'r');
        title(['Call prices at expiry ', datestr(dataset.datesExpiry(ii))]); xlabel('Strikes'); ylabel('Prices');
        legend('Calibrated prices', 'Mean prices', 'Call Ask', 'Call Bid', 'Strike ATM');

        subplot(1, 2, 2);
        plot(strikes, put_prices, '*-'); hold on;
        plot(strikes, mean_put, 'o-'); grid on;
        plot(strikes, dataset.putAsk(ii).prices, '--');
        plot(strikes, dataset.putBid(ii).prices, '--');
        xline(F0(ii), '--', 'LineWidth', 2, 'Color', 'r');
        title(['Put prices at expiry ', datestr(dataset.datesExpiry(ii))]); xlabel('Strikes'); ylabel('Prices');
        legend('Calibrated prices', 'Mean prices', 'Put Ask', 'Put Bid', 'Strike ATM');
        
    end

    %% Display of the errors in the prices:

    % For American options
    if length(dataset.datesExpiry)>13 
        fprintf('\nMEAN ERROR USA PRICES:\n')
        disp('--------------------------------------------------------------')
        fprintf('Expiry         | Call Prices error | Put prices error\n')
        disp('--------------------------------------------------------------')
        for ii=1:length(dataset.datesExpiry)
            fprintf('%s     |  %f%%       |    %f%%\n', datestr(dataset.datesExpiry(ii)), error_call_prices_vec(ii), error_put_prices_vec(ii))
        end

    % For European options
    else
        fprintf('\nMEAN ERROR EU PRICES:\n')
        disp('--------------------------------------------------------------')
        fprintf('Expiry         | Call Prices error | Put prices error\n')
        disp('--------------------------------------------------------------')
        for ii=1:length(dataset.datesExpiry)
            fprintf('%s     |  %f%%       |    %f%%\n', datestr(dataset.datesExpiry(ii)), error_call_prices_vec(ii), error_put_prices_vec(ii))
        end
    end
end % function plot_calls_puts_total
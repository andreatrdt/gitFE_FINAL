function blk_plot_calls_puts_total(dataset, F0, B0, sigma, date_settlement)
% Plotting of the calls and put prices after the calibration
% 
% INPUT:
% dataset:            [STRUCT] initial dataset
% F0:                 [VECTOR] initial forward value F(0, T)
% B0:                 [VECTOR] initial discounts B(0, T)
% sigma:              [SCALAR] volatility
% date_settlement:    [DATENUM] initial date
% OUTPUT:
% None
% USES:  blkprice()

% Authors:
% M.Maspes, A.Tarditi, M.Torba

%% Initialization of the error vectors:
    error_call_prices_vec = zeros(length(dataset.datesExpiry),1);
    error_put_prices_vec = zeros(length(dataset.datesExpiry),1);

    count_prices_neg = 0;
    count_prices_neg_put = 0;

    for ii = 1:min(length(dataset.datesExpiry),19)

        %% Initialization

        % Vector of strikes:
        strikes = dataset.strikes(ii).value;

        % Year frac convenction Act/365:
        conv_ACT365 = 3;

        %% Computation of the prices
             
        % Time to maturity
        TTM = yearfrac(date_settlement, datenum(dataset.datesExpiry(ii)), conv_ACT365);
        
        % Interest rate
        interest_rate = -log(B0(ii))/TTM;

        %% Pricing 
        
        % Price of Call/Puts through the Black formula
        [call_prices,put_prices] = blkprice(F0(ii), strikes, interest_rate, TTM, sigma);
        
        % Check for negative prices
        count_prices_neg = count_prices_neg + length(find(call_prices < 0));
        count_prices_neg_put = count_prices_neg_put + length(find(put_prices < 0));
        
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
        title(['Black model Call prices at expiry ', datestr(dataset.datesExpiry(ii))]); xlabel('Strikes'); ylabel('Prices');
        legend('Black calibrated prices', 'Mean prices', 'Call Ask', 'Call Bid', 'Strike ATM', 'Location', 'best');

        subplot(1, 2, 2);
        plot(strikes, put_prices, '*-'); hold on;
        plot(strikes, mean_put, 'o-'); grid on;
        plot(strikes, dataset.putAsk(ii).prices, '--');
        plot(strikes, dataset.putBid(ii).prices, '--');
        xline(F0(ii), '--', 'LineWidth', 2, 'Color', 'r');
        title(['Black Put prices at expiry ', datestr(dataset.datesExpiry(ii))]); xlabel('Strikes'); ylabel('Prices');
        legend('Black calibrated prices', 'Mean prices', 'Put Ask', 'Put Bid', 'Strike ATM','Location', 'best');
        
    end

    %% Display of the errors in the prices:

    % For American options
    if length(dataset.datesExpiry)>13 

        fprintf('\nMEAN ERROR USA PRICES WITH BLACK MODEL:\n')
        disp('--------------------------------------------------------------')
        fprintf('The number of negative prices in the USA mkt is: %d\n', count_prices_neg)
        fprintf('The number of negative Put prices in the USA mkt is: %d\n', count_prices_neg_put)
        disp('--------------------------------------------------------------')
        
        fprintf('Expiry         | Black Call Prices error | Black Put prices error\n')
        disp('--------------------------------------------------------------')
        for ii=1:min(length(dataset.datesExpiry),19)
            fprintf('%s     |  %f%%            |    %f%%\n', datestr(dataset.datesExpiry(ii)), error_call_prices_vec(ii), error_put_prices_vec(ii))
        end
        disp('--------------------------------------------------------------')
        fprintf('\n AVERAGE MEAN ERROR USA PRICES WITH BLACK:\n')
        fprintf('Call Prices error | Put prices error\n')
        fprintf(' %f%%       |    %f%%\n', mean(error_call_prices_vec), mean(error_put_prices_vec))

    % For European options
    else

        fprintf('\nMEAN ERROR EU PRICES WITH BLACK MODEL:\n')
        disp('--------------------------------------------------------------')
        fprintf('The number of negative prices in the EU mkt is: %d\n', count_prices_neg)
        fprintf('The number of negative Put prices in the EU mkt is: %d\n', count_prices_neg_put)
        disp('--------------------------------------------------------------')
        fprintf('Expiry         | Black Call Prices error | Black Put prices error\n')
        disp('--------------------------------------------------------------')
        for ii=1:length(dataset.datesExpiry)
            fprintf('%s     |  %f%%            |    %f%%\n', datestr(dataset.datesExpiry(ii)), error_call_prices_vec(ii), error_put_prices_vec(ii))
        end
        disp('--------------------------------------------------------------')
        fprintf('\n AVERAGE MEAN ERROR EU PRICES WITH BLACK:\n')
        fprintf('Call Prices error | Put prices error\n')
        fprintf(' %f%%       |    %f%%\n', mean(error_call_prices_vec), mean(error_put_prices_vec))
    end

end % function plot_calls_puts_total
function blk_plot_volatility_smiles(dataset, F0, B0, sigma, date_settlement)
% Plotting of the calls and put prices after the calibration
% 
% INPUT:
% dataset:            [STRUCT] initial dataset
% F0:                 [VECTOR] initial forward value F(0, T)
% B0:                 [VECTOR] initial discounts B(0, T)
% params:             [VECTOR] [k, theta, sigma]
% date_settlement:    [DATENUM] initial date
% 
% USES:         blkimpv() , blkprice()

% Authors:
% M.Maspes, A.Tarditi, M.Torba

    for ii = 1:length(dataset.datesExpiry)

        %% Initialization
    
        conv_ACT365 = 3;

        %% Computation of the prices
        
        % Parameters pricing
        strikes = dataset.strikes(ii).value;
        
        put_length = length(dataset.putAsk(ii).prices);

        log_moneyness = log(F0(ii) ./ strikes);

        TTM = yearfrac(date_settlement, datenum(dataset.datesExpiry(ii)), conv_ACT365);

        interest_rate = -log(B0(ii))/TTM;

        % Pricing the calls through FFT
        call_prices = blkprice(F0(ii), strikes, interest_rate, TTM, sigma);
        
        % Compute the implied volatilities OLD
        mid_price_put = (dataset.putAsk(ii).prices + dataset.putBid(ii).prices)/2;
        mid_price_put = mid_price_put + B0(ii) .* (F0(ii) - strikes(1:put_length));

        mid_price_call = (dataset.callAsk(ii).prices + dataset.callBid(ii).prices)/2;
        impvol_OLD= blkimpv(F0(ii), strikes, interest_rate, TTM, [mid_price_put mid_price_call], 'Class', {'Call'});

        % Compute the implied volatilities NEW
        impvol_NEW= blkimpv(F0(ii), strikes, interest_rate, TTM, call_prices, 'Class', {'Call'});
        
        
        %% Plot of the volatilities

        figure();
        plot(strikes, impvol_OLD, 'o-'); hold on;
        plot(strikes, impvol_NEW, '*-'); grid on;
%         plot(log_moneyness, impvol_OLD, 'o-'); hold on;
%         plot(log_moneyness, impvol_NEW, '*-'); grid on;

        title(['Implied volatilities at expiry ', datestr(dataset.datesExpiry(ii))]);
        xlabel('Log-moneyness'); ylabel('Volatilities');
        legend('Imp Vol dataset', 'Imp Vol calibrated','Location','best');

    end

end % function plot_volatility_smiles
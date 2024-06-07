function [dataset] = dataset_preprocessing(dataset, F0, B0, date_settlement, flag)
% Preprocessing of the dataset for each maturity in order to remove all
% those options without delta_Black between [10%, 90%]
% 
% INPUT:
% dataset:            [STRUCT] initial dataset
% F0:                 [VECTOR] initial forward F(0, T)
% B0:                 [VECTOR] initial discount B(0, T)
% date_settlement:    [SCALAR] settlement date
% flag:               [0: no plots, 1: with plots]
% 
% OUTPUT:
% dataset:            modified version of the dataset
%
% USES:     blkimpv() , blsdelta()

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Conventions
    conv_ACT365 = 3;

    %% Computation
    
    for ii=1:length(dataset.datesExpiry)

        %% Quantities of interest
        
        % Initial values
        spot_ATM = dataset.spot;

        strike_ATM = F0(ii);
        strikes = dataset.strikes(ii).value;

        TTM = yearfrac(date_settlement, datenum(dataset.datesExpiry(ii)), conv_ACT365);
        interest_rate = -log(B0(ii))/TTM;

        % Indexes OTM option
        idx_call_OTM = find(strikes > strike_ATM);
        idx_put_OTM = find(strikes <= strike_ATM);

        mid_price_call = (dataset.callAsk(ii).prices + dataset.callBid(ii).prices)/2;
        mid_price_call = mid_price_call(idx_call_OTM);

        mid_price_put = (dataset.putAsk(ii).prices + dataset.putBid(ii).prices)/2;
        mid_price_put = mid_price_put(idx_put_OTM);

        %% Comparison for the last implied volatility
        % We modified the strikes of the last American maturity in order to
        % find a better shape of the implied volatilities (the code is
        % brute force since we need to work directly on this specific
        % dataset)
        
        if length(dataset.datesExpiry) == 20 && ii == length(dataset.datesExpiry)

            % Strikes to remove the jumps
            strikes_removal_put = [2000 2200 2400 2600 2700 4900 5000 5100 5200];
            strikes_removal_call = [5300 7200];

            % Creation of the indexes
            [~, idx_removed_strikes_put] = ismember(strikes(idx_put_OTM), strikes_removal_put);
            [~, idx_removed_strikes_call] = ismember(strikes(idx_call_OTM), strikes_removal_call);

            idx_put = find(idx_removed_strikes_put == 0);
            idx_call = find(idx_removed_strikes_call == 0) + length(idx_put_OTM);

            % Upgrade of the quantities
            strikes = strikes([idx_put idx_call]);

            mid_price_put = mid_price_put(idx_put);
            mid_price_call = mid_price_call(idx_call - length(idx_put_OTM));

            idx_call_OTM = find(strikes > strike_ATM);
            idx_put_OTM = find(strikes <= strike_ATM);

            mid_price_call = mid_price_call(idx_call_OTM - length(idx_put_OTM));
            mid_price_put = mid_price_put(idx_put_OTM);
        end

        %% Compute the Black implied volatilities
        
        impvol_call_i= blkimpv(F0(ii), strikes(idx_call_OTM), interest_rate, TTM, mid_price_call, 'Class', {'Call'});
        impvol_put_i = blkimpv(F0(ii), strikes(idx_put_OTM), interest_rate, TTM, mid_price_put, 'Class', {'Put'});
        
        %% Plot of the matched volatilities
        
        if flag
            figure();
            plot(strikes([idx_call_OTM(1)-1 idx_call_OTM]), [impvol_put_i(end) impvol_call_i], 'o-'); hold on;
            plot(strikes(idx_put_OTM ), impvol_put_i, '*-'); grid on;

            title('Implied volatilities');
            xlabel('Strikes'); ylabel('Volatilities');
            legend('Implied vol Call', 'Implied vol Put');
        end

        %% Computation of the delta
        
        [delta_call, ~] = blsdelta(spot_ATM, strikes(idx_call_OTM), interest_rate, TTM, impvol_call_i);
        [~, delta_put] = blsdelta(spot_ATM, strikes(idx_put_OTM), interest_rate, TTM, impvol_put_i);

        %% Restructuring of the dataset

        % Find the indexes to cut for the put
        indicator_10 = delta_put <= -0.1;
        indicator_90 = delta_put >= -0.9;
        indicator = indicator_10 .* indicator_90;

        idx_put = find(indicator > 0);

        % Find the indexes to cut for the call
        indicator_10 = delta_call >= 0.1;
        indicator_90 = delta_call <= 0.9;
        indicator = indicator_10 .* indicator_90;

        idx_call = find(indicator > 0) + length(delta_put);

        % Cut of the structs
        dataset.callBid(ii).prices = dataset.callBid(ii).prices(idx_call);
        dataset.callAsk(ii).prices = dataset.callAsk(ii).prices(idx_call);
        dataset.callAsk(ii).impvol = dataset.callAsk(ii).impvol(idx_call);
        dataset.callAsk(ii).impvol = impvol_call_i(idx_call - length(delta_put));

        dataset.putAsk(ii).prices = dataset.putAsk(ii).prices(idx_put);
        dataset.putBid(ii).prices = dataset.putBid(ii).prices(idx_put);
        dataset.putBid(ii).impvol = dataset.putBid(ii).impvol(idx_put);
        dataset.putBid(ii).impvol = impvol_put_i(idx_put);

        idx_combined = unique([idx_put idx_call]);
        dataset.strikes(ii).value = dataset.strikes(ii).value(idx_combined);

        dataset.Volume_call(ii).volume = dataset.Volume_call(ii).volume(idx_call);
        dataset.Volume_put(ii).volume = dataset.Volume_put(ii).volume(idx_put);
        
    end

end % function dataset_preprocessing
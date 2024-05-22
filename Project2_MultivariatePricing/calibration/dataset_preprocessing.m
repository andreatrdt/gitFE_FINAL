function dataset_preprocessing(dataset, F0, B0, date_settlement, flag)
% Preprocessing of the dataset for each maturity in order to remove all
% those options without delta_Black between [10%, 90%]
% 
% INPUT:
% 
% OUTPUT:
% 
% USES:
% 
    
    for ii=1:length(dataset.datesExpiry)

        %% Quantities of interest

        TTM = yearfrac(date_settlement, dataset.datesExpiry(ii));
        interest_rate = -log(B0(ii))/TTM;
        
        mid_price_call = (dataset.callAsk(ii).prices + dataset.callBid(ii).prices)/2;
        mid_price_put = (dataset.putAsk(ii).prices + dataset.putBid(ii).prices)/2;

        strikes = dataset.strikes(ii).value;

        %% Compute the Black implied volatilities
        
        impvol_call_i_vector = blkimpv(F0(ii), strikes, interest_rate, TTM, mid_price_call, 'Class', {'Call'});
        impvol_put_i_vector = blkimpv(F0(ii), strikes, interest_rate, TTM, mid_price_put, 'Class', {'Put'});
        
        % Plot of the matched volatilities

        if flag
            figure();
            plot(strikes, impvol_call_i_vector, 'o-'); hold on;
            plot(strikes, impvol_put_i_vector, '*-');
        end

        % Keep only Delta which fall in the interval [10% , 90%]
        [delta_call, delta_put] = blsdelta(F0(ii)/B0(ii), strikes, interest_rate, TTM, impvol_call_i_vector)
%         delta_call = blkdelta(F0(ii), strikes, B0(ii), impvol_call_i_vector, TTM, 1);
%         delta_put = blkdelta(F0(ii), strikes, B0(ii), impvol_put_i_vector, TTM, 0);

        % Consider only the OTM options
        a
        % 
    
    end
end
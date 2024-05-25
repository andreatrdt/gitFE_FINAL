function dataset = OTM_preprocessing(dataset, B0, F0)
% Computation of the OTM options only
% 
% INPUT:
% dataset:          initial dataset
% B0:               initial discount 
% F0:               forward price at time 0
% 
% OUTPUT:
% dataset:          dataset modified for the OTM options 
    
    for ii=1:length(dataset.datesExpiry)

        %% Quantities of interest
        
%         strike_ATM = dataset.spot/B0(ii);
        strike_ATM = F0(ii);

        strikes = dataset.strikes(ii).value;
        
        %% Compute the related index for the strike ATM
        
        indicator_call_OTM = strikes > strike_ATM;
        indicator_put_OTM = strikes <= strike_ATM;
        
        %% Restructuring of the dataset

        dataset.callBid(ii).prices = dataset.callBid(ii).prices .* indicator_call_OTM;
        dataset.callAsk(ii).prices = dataset.callAsk(ii).prices .* indicator_call_OTM;

        dataset.putAsk(ii).prices = dataset.putAsk(ii).prices .* indicator_put_OTM;
        dataset.putBid(ii).prices = dataset.putBid(ii).prices .* indicator_put_OTM;
    
    end
    
end % function OTM_preprocessing
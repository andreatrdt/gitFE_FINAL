function dataset_renovated = removal_expiry(dataset, indexes)
% Removal of a specific expiry after preprocessing decisions
% 
% INPUT:
% dataset:                 [STRUCT] dataset from the mkt
% indexes:                 [VECTOR] expiry to remove
% 
% OUTPUT:
% dataset_renovated:       [STRUCT] dataset modified

    %% Introductive values
    length_dataset = length(dataset.datesExpiry);

    indexes = sort(unique(indexes), 'descend');

    for ii = 1:length(indexes)
        if indexes(ii) > length_dataset
            error("Idx required doesn't exist");
        end
    end

    %% Removal of the components and creation of the new dataset

    for ii = 1:length(indexes)

        % Choose the index
        idx = indexes(ii);

        % Initializa the new dataset
        dataset_renovated = struct();
        
        % Expiries
        new_dates = [dataset.datesExpiry(1:idx-1); dataset.datesExpiry(idx+1:end)];
        dataset_renovated.datesExpiry = new_dates;
    
        % Call Ask
        new_callAsk = [dataset.callAsk(1:idx-1) dataset.callAsk(idx+1:end)];
        dataset_renovated.callAsk = new_callAsk;
    
        % Call Bid
        new_callBid = [dataset.callBid(1:idx-1) dataset.callBid(idx+1:end)];
        dataset_renovated.callBid = new_callBid;
    
        % Put Ask
        new_putAsk = [dataset.putAsk(1:idx-1) dataset.putAsk(idx+1:end)];
        dataset_renovated.putAsk = new_putAsk;
    
        % Put Bid
        new_putBid = [dataset.putBid(1:idx-1) dataset.putBid(idx+1:end)];
        dataset_renovated.putBid = new_putBid;
    
        % Volume Call
        new_Volume_call = [dataset.Volume_call(1:idx-1) dataset.Volume_call(idx+1:end)];
        dataset_renovated.Volume_call = new_Volume_call;
    
        % Volume Put
        new_Volume_put = [dataset.Volume_put(1:idx-1) dataset.Volume_put(idx+1:end)];
        dataset_renovated.Volume_put = new_Volume_put;
    
        % Strikes
        new_strikes = [dataset.strikes(1:idx-1) dataset.strikes(idx+1:end)];
        dataset_renovated.strikes = new_strikes;
    
        % Spot
        dataset_renovated.spot = dataset.spot;

        % Final change for new for computations
        dataset = dataset_renovated;
    end

end % function removal_expiry
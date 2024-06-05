function [dataset_renovated, F0_renovated, B0_renovated] = removal_expiry(dataset, F0, B0, indexes)
% Removal of a specific expiry after preprocessing decisions
% 
% INPUT:
% dataset:                 [STRUCT] dataset from the mkt
% F0:                      [VECTOR] initial fwd 
% B0:                      [VECTOR] initial discount 
% indexes:                 [VECTOR] expiry to remove
% 
% OUTPUT:
% dataset_renovated:       [STRUCT] dataset modified
% F0_renovated:            [VECTOR] initial fwd modified
% B0_renovated:            [VECTOR] initial discount modified
%
% USES: none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    %% Introductive values
    length_dataset = length(dataset.datesExpiry);

    indexes = sort(unique(indexes), 'descend');

    for ii = 1:length(indexes)
        if indexes(ii) > length_dataset
            error("Idx required doesn't exist");
        end
    end

    %% Controls
    % Control that F0 and B0 are columns
    
    if (size(F0, 2) > size(F0, 1))
        F0 = F0';
    end

    if (size(B0, 2) > size(B0, 1))
        B0 = B0';
    end

    %% Removal of the components and creation of the new dataset

    for ii = 1:length(indexes)

        %% Choose the index
        idx = indexes(ii);

        %% New dataset
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

        %% New Forwards and discounts

        F0_renovated = [F0(1:idx-1); F0(idx+1:end)];
        B0_renovated = [B0(1:idx-1); B0(idx+1:end)];

        %% Final changes for new computations

        dataset = dataset_renovated;
        F0 = F0_renovated;
        B0 = B0_renovated;

    end

end % function removal_expiry
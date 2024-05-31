function dataset_renovated = removal_expiry(dataset, idx)
% Removal of a specific expiry after preprocessing decisions
% 
% INPUT:
% dataset:            [STRUCT] dataset from the mkt
% idx:                [SCALAR] expiry to remove
% 
% OUTPUT:
% dataset_renovated:            [STRUCT] dataset modified

    
    %% Introductive values
    length_dataset = length(dataset.datesExpiry);

    if idx < length_dataset
        error("Idx required doesn't exist");
    end

    %% Removal of the components

    dataset_renovated = struct;
    
    % Expiries

    for ii = 1: idx
        dataset_renovated.datesExpiry = [dataset_renovated.datesExpiry; dataset.datesExpiry(ii)];
    end

    for ii = idx+1: length_dataset
        dataset_renovated.datesExpiry = [dataset_renovated.datesExpiry; dataset.datesExpiry(ii)];
    end




end % function removal_expiry
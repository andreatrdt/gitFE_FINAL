function idx = find_idx(data, date)
% Computation that finds the index in a struct given a certain date
%
% INPUT:
% data:         [STRUCT]struct containing the data
% date:         [DATENUM]date to be found
% 
% OUTPUT:
% idx: index of the date in the struct
%
% USES: find_idx()

% Authors:
% M.Maspes, A.Tarditi, M.Torba


    idx = find(data == date);

end % function find_idx